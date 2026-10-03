"""Conservative alignment-span intervals. No neighbour filling or component labels."""
from collections import defaultdict

def regions(alignments, minimum):
    events = defaultdict(lambda: defaultdict(int))
    for lo, hi, chrom, *_ in alignments:
        if hi <= lo:
            continue
        events[lo][chrom] += 1
        events[hi][chrom] -= 1
    active, pieces, previous = {}, [], None
    for position in sorted(events):
        labels = [ch for ch, count in active.items() if count > 0]
        if previous is not None and position > previous and len(labels) == 1:
            ch = labels[0]
            if pieces and pieces[-1][1] == previous and pieces[-1][2] == ch:
                pieces[-1] = (pieces[-1][0], position, ch)
            else:
                pieces.append((previous, position, ch))
        for ch, delta in events[position].items():
            active[ch] = active.get(ch, 0) + delta
        previous = position
    # Overlapping alignments to the same chromosome count once.
    return [p for p in pieces if p[1] - p[0] >= minimum]

def intervals(alignments, minimum):
    anchors = regions(alignments, minimum)
    result = []
    previous = None
    for lo, hi, ch in anchors:
        if previous and previous[2] != ch:
            result.append(dict(lo=previous[1], hi=lo, left=previous[2], right=ch,
                               left_start=previous[0], right_end=hi))
        previous = (lo, hi, ch)
    return result

def compatible_gaps(interval, gaps):
    # gaps are (midpoint, source record). No nearest-gap snapping.
    return [(pos, row) for pos, row in gaps if interval["lo"] <= pos <= interval["hi"]]


def arm_context(alignments, interval, minimum):
    """Descriptive arm support; no absolute chromosome-size gate or cut authorization.

    Consecutive same-label anchors form an arm even across unaligned sequence.
    Support is union aligned bases, not the arm's bounding span. All other arms
    remain visible: an A-B-A excursion is not silently converted into an A-A join.
    """
    arms = []
    for lo, hi, chrom in regions(alignments, minimum):
        if arms and arms[-1]["chrom"] == chrom:
            arms[-1]["end"] = hi
            arms[-1]["bp"] += hi-lo
        else:
            arms.append(dict(start=lo, end=hi, chrom=chrom, bp=hi-lo))
    total = sum(arm["bp"] for arm in arms)
    for i, (left, right) in enumerate(zip(arms, arms[1:])):
        if left["end"] != interval["lo"] or right["start"] != interval["hi"]:
            continue
        excursion = ((i > 0 and arms[i-1]["chrom"] == right["chrom"]) or
                     (i+2 < len(arms) and arms[i+2]["chrom"] == left["chrom"]))
        return dict(left_arm_aligned_bp=left["bp"], right_arm_aligned_bp=right["bp"],
                    left_arm_fraction=left["bp"]/total, right_arm_fraction=right["bp"]/total,
                    assigned_union_bp=total, arm_count=len(arms),
                    arm_pattern="out_and_back" if excursion else "chromosome_transition")
    raise ValueError("Transition has no matching unambiguous arms")

def verdict_for_interval(candidate, interval, number):
    # Scaffold-wide votes cannot authorize one of several unadjudicated junctions.
    # Preserve only a single two-member candidate's conservative decision.
    members = candidate.get("members", "").split("+")
    named = {n.rsplit("_", 1)[0] for n in candidate.get("name", "").split("+")}
    return (candidate.get("verdict", "REVIEW")
            if number == 1 and len(members) == 2 and
            named == {interval["left"], interval["right"]} else "REVIEW")

