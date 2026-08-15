#!/usr/bin/env python3
"""
harmonize_names.py -- give homologous chromosomes consistent names across a set
of same-species assemblies, using one in-batch assembly as the naming reference.

RENAME-AND-REORIENT ONLY. This never re-scaffolds: scaffold count, boundaries and
internal gaps are invariant. The reference's chromosome-scale scaffolds define the
chromosome slots (chrN, N = descending-size rank); every other assembly inherits
those names by minimap2 homology.

Name grammar (always suffixed, so the output is uniform for later graph partitioning):

    chr{K}_{part}                  one scaffold covering chromosome K, piece <part>
    chr{A}_{i}+chr{B}_{j}[+...]    one scaffold spanning >1 chromosome (fusion / chimera)
    unplaced_{n}                   no confident chromosome placement

  * Split case (2 scaffolds cover chr K): both get primary K, ordered by their
    reference start -> chrK_1, chrK_2. Each is reverse-complemented to match K.
  * Merge / fusion (1 scaffold covers chr A and chr B): composite name carrying a
    per-chromosome part token for each chromosome it touches; NOT reoriented (its two
    segments can disagree on strand), so left native. Flagged, with batch concordance.

Placement is decided by alignment, not by size: a small scaffold that aligns
confidently to chr1 is chr1; a large scaffold that aligns nowhere is unplaced.

Inputs via a manifest TSV (header: id  role  fai  paf):
    - one 'ref' row  (paf column = '-')
    - one 'query' row per other assembly (paf = minimap2 <reference> <query> PAF; .gz ok)
Each 'fai' is a samtools .fai (name<TAB>length<TAB>...); only the first two cols are read.

Outputs (into --outdir):
    {id}.harmonized_name_map.tsv   per assembly (incl. the reference)
    {species}.harmonization_report.tsv   one row per (assembly, scaffold) + flags
"""

import argparse
import csv
import gzip
import os
import re
import statistics
import sys
from collections import defaultdict
from itertools import combinations

HAP_SUFFIX_RE = re.compile(r"_hap[0-9]+$", re.IGNORECASE)


def individual_of(assembly_id, hap_re=HAP_SUFFIX_RE):
    """Collapse an assembly id to its individual by stripping a haplotype suffix.

    Presence signals must be counted per individual, not per haplotype: a sex-limited
    chromosome (W / Y) is carried by exactly one haplotype of a heterogametic individual,
    so counting haplotypes halves its apparent presence.
    """
    return hap_re.sub("", assembly_id)


def scaffold_n50(lengths):
    """N50 of a list of scaffold lengths (0 for an empty list)."""
    ls = sorted(lengths, reverse=True)
    tot = sum(ls)
    if tot == 0:
        return 0
    half, cum = tot / 2.0, 0
    for x in ls:
        cum += x
        if cum >= half:
            return x
    return ls[-1]


def assign_roles(metrics, min_cut_ratio, min_genome_frac, min_n50_ratio,
                 require_dropoff, min_voters, reference_id):
    """Split assemblies into consensus 'voter's and naming-only 'passenger's.

    A passenger is still aligned, placed, named and emitted -- it simply does not vote in
    the batch-consensus coverage count or the fusion concordance counts. This keeps a
    failed assembly's name map as a deliverable while stopping it from corroborating
    reference fragments or inflating concordance denominators.

    Every criterion is contiguity / self-consistency based. None encodes a chromosome
    number, so this is safe for species with no known karyotype.

    Returns (roles, reasons, flags).
    """
    n50s = [m["n50"] for m in metrics.values() if m["n50"] > 0]
    med_n50 = statistics.median(n50s) if n50s else 0
    roles, reasons, flags = {}, {}, []
    for rid, m in metrics.items():
        why = []
        if require_dropoff and m["method"] != "dropoff":
            why.append("chrom_set_method=%s" % m["method"])
        if min_cut_ratio > 0 and (m["cut_ratio"] is None or m["cut_ratio"] < min_cut_ratio):
            why.append("cut_ratio=%s<%.2f" % (
                "NA" if m["cut_ratio"] is None else "%.2f" % m["cut_ratio"], min_cut_ratio))
        if m["genome_fraction"] < min_genome_frac:
            why.append("genome_frac=%.3f<%.2f" % (m["genome_fraction"], min_genome_frac))
        if med_n50 > 0 and m["n50"] < min_n50_ratio * med_n50:
            why.append("n50_ratio=%.3f<%.2f" % (m["n50"] / med_n50, min_n50_ratio))
        roles[rid] = "passenger" if why else "voter"
        reasons[rid] = ";".join(why) if why else "-"
    voters = [r for r, v in roles.items() if v == "voter"]
    if len(voters) < min_voters:
        flags.append("voter_filter_disabled(only_%d_of_%d_passed)" % (len(voters), len(roles)))
        for rid in roles:
            roles[rid] = "voter"
            reasons[rid] = (reasons[rid] + ";overridden_min_voters").lstrip("-;")
    elif roles.get(reference_id) == "passenger":
        # The PAFs are already computed against this reference, so we cannot re-select
        # here. Keep it as a voter but make the contradiction loud.
        flags.append("reference_failed_voter_criteria")
        roles[reference_id] = "voter"
        reasons[reference_id] = reasons[reference_id] + ";forced_voter_is_reference"
    return roles, reasons, flags



def openf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fai(path):
    """Return [(name, length), ...] in file order."""
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out.append((f[0], int(f[1])))
    return out


def parse_paf(path, chrom_of):
    """Accumulate per (query_scaffold, chromosome) alignment evidence.

    chrom_of maps a *reference scaffold name* -> chromosome number, and only
    contains chromosome-scale reference scaffolds. Alignments to sub-threshold
    reference sequences are ignored (they confer no chromosome name).

    'bp' is the sum of PAF block lengths, which double-counts overlapping alignments and
    can exceed the true footprint (observed at 120-130% of query length on repeat-rich
    sequence). 'ivals' keeps the raw target intervals so a merged, non-double-counted
    coverage can be computed where correctness matters.

    Returns:
        ev:   qname -> chr_num -> {bp, tstart, tend, plus, minus, ivals}
    """
    ev = defaultdict(lambda: defaultdict(
        lambda: {"bp": 0, "tstart": None, "tend": 0, "plus": 0, "minus": 0,
                 "ivals": []}))
    if not path or path == "-" or not os.path.exists(path):
        return ev
    with openf(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            q = c[0]
            strand = c[4]
            t = c[5]
            ts = int(c[7])
            te = int(c[8])
            blk = int(c[10])          # alignment block length -> aligned bp
            if t not in chrom_of:
                continue
            k = chrom_of[t]
            e = ev[q][k]
            e["bp"] += blk
            e["tstart"] = ts if e["tstart"] is None else min(e["tstart"], ts)
            e["tend"] = max(e["tend"], te)
            e["ivals"].append((ts, te))
            if strand == "-":
                e["minus"] += blk
            else:
                e["plus"] += blk
    return ev



GRAPH_RULES = ("permissive", "majority", "plurality")


def edge_keep(rule, n_f, n_s, min_fused=2):
    """Should a fused pair become a join-graph edge under `rule`?

    min_fused is a hard floor applied BEFORE any rule: a junction seen fused in fewer than
    this many voters is a single observation, not evidence. Without it, `plurality`
    (n_f >= n_s) accepts a 1f/1s pair -- one assembly's mis-join against one assembly's
    split -- and union-find then chains it into an unrelated component. On the Sde run that
    turned chr2 + chr4 + chr20 into a single 185 Mb "chromosome", larger than the largest
    real chromosome (95 Mb). `permissive`'s `n_s == 0` branch has the same hole: 1f/0s would
    pass on one observation.

    permissive -- >1 voter fused, or nobody splits it. Passes 2f/10s, so it over-joins.
    majority   -- strictly more voters fused than split. Conservative; leaves ties open.
    plurality  -- fused >= split. Treats a genuine tie as a join.

    Note the rules are NOT nested: permissive's own `n_f > 1` clause makes it stricter than
    plurality at n_f == 1. All three are reported so tie cases stay visible instead of being
    silently resolved by whichever rule happens to be wired into naming.
    """
    if n_f < max(1, min_fused):
        return False
    if rule == "permissive":
        return n_f > 1 or n_s == 0
    if rule == "majority":
        return n_f > n_s
    if rule == "plurality":
        return n_f >= n_s
    raise ValueError("unknown graph rule %r" % rule)


def components(nodes, edges, node_len=None):
    """Union-find over reference-frame piece numbers -> consensus chromosomes.

    Consensus ids are 1..N ordered by DESCENDING component total length (ties broken by
    smallest member), so chr1 is the largest chromosome -- the universal convention.

    Ordering by smallest member instead, as this did originally, is wrong once any merge
    happens: the id inherits the rank of the component's smallest member while the total
    length comes from all of them, so a merged component lands in a low-numbered slot with
    a large total while mid-rank singletons keep their original positions. On the Sde run
    that put chr9 at 49.3 Mb between chr8 at 78.0 Mb and chr10 at 86.3 Mb, and left 13 of
    15 ids out of size order. Indefensible in a published assembly.

    node_len maps piece number -> length; omit it to fall back to smallest-member ordering.

    NB consensus ids are not comparable across graph rules -- a different edge set gives a
    different component count and different totals, so the numbering changes. The rule must
    therefore be settled before any name derived from it ships.
    """
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for x, y in edges:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[max(rx, ry)] = min(rx, ry)
    groups = defaultdict(list)
    for n in nodes:
        groups[find(n)].append(n)
    if node_len:
        ordered = sorted(groups.values(),
                         key=lambda g: (-sum(node_len.get(n, 0) for n in g), min(g)))
    else:
        ordered = sorted(groups.values(), key=lambda g: min(g))
    comp_of, comp_members = {}, {}
    for cid, g in enumerate(ordered, start=1):
        comp_members[cid] = sorted(g)
        for n in g:
            comp_of[n] = cid
    return comp_of, comp_members


def merged_length(ivals):
    """Total length of the union of [start, end) intervals."""
    if not ivals:
        return 0
    iv = sorted(ivals)
    tot = 0
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s > ce:
            tot += ce - cs
            cs, ce = s, e
        elif e > ce:
            ce = e
    return tot + (ce - cs)


def classify(chr_ev, scaf_len, min_frac, sec_frac, chrom_len=None,
             member_cover_frac=0.0):
    """Decide placement of one query scaffold from its per-chromosome evidence.

    Only called for scaffolds that passed the per-query chromosome-scale eligibility set
    (select_chromosome_set applied to the query's own lengths); placement here is purely
    alignment-decided.

    A chromosome K becomes a composite MEMBER of this scaffold when either:

      share  -- K contributes >= sec_frac of the scaffold's total aligned bp, or
      cover  -- >= member_cover_frac of reference chromosome K's own length is inside
                this scaffold (merged, non-double-counted coverage).

    The share test alone has a size-asymmetry blind spot: a small chromosome joined to a
    large one is invisible because its share of the combined alignment is small by
    construction. An 8.7 Mb chromosome inside an 84.6 Mb chromosome contributes
    8.7/93.3 = 0.093 of the aligned total -- under a 0.2 threshold -- so the scaffold was
    classified as a plain chromosome and the smaller chromosome silently vanished from the
    frame with no flag. The smaller the chromosome, the more invisible it became.

    The cover test is size-invariant: it asks whether this scaffold holds most of K,
    regardless of how large the scaffold is. It uses merged interval coverage rather than
    summed block length, so repeat-driven double counting cannot satisfy it.
    """
    total_bp = sum(e["bp"] for e in chr_ev.values())
    if scaf_len == 0 or total_bp == 0 or (total_bp / scaf_len) < min_frac:
        return {"class": "unplaced",
                "aligned_frac": (total_bp / scaf_len) if scaf_len else 0.0}
    cov = {k: merged_length(e.get("ivals") or []) for k, e in chr_ev.items()}
    by_bp = sorted(chr_ev.items(), key=lambda kv: kv[1]["bp"], reverse=True)
    primary = by_bp[0][0]
    partners, reasons = [], {}
    for k, e in by_bp[1:]:
        by_share = e["bp"] >= sec_frac * total_bp
        klen = (chrom_len or {}).get(k, 0)
        by_cover = (member_cover_frac > 0 and klen > 0
                    and cov.get(k, 0) >= member_cover_frac * klen)
        if by_share or by_cover:
            partners.append(k)
            if by_cover and not by_share:
                # only reachable via the size-invariant test -- worth surfacing, because
                # every one of these was silently dropped before
                reasons[k] = "size_asymmetric(ref%d:share=%.3f,cover=%.2f)" % (
                    k, e["bp"] / total_bp, (cov.get(k, 0) / klen) if klen else 0.0)
    members = sorted([primary] + partners)
    orient = {k: ("-" if chr_ev[k]["minus"] > chr_ev[k]["plus"] else "+")
              for k in members}
    info = {
        "class": "composite" if partners else "chromosome",
        "primary": primary,
        "members": members,
        "tstart": {k: (chr_ev[k]["tstart"] or 0) for k in members},
        "tend": {k: chr_ev[k]["tend"] for k in members},
        "orient": orient,
        "aligned_frac": total_bp / scaf_len,
        "cover_frac": sum(cov.get(k, 0) for k in chr_ev) / scaf_len if scaf_len else 0.0,
        "member_reasons": reasons,
    }
    return info


def select_chromosome_set(ref_fai, min_scaffold_bp, method, dropoff_ratio, dropoff_min_frac,
                          min_chrom_frac=0.0):
    """Pick the reference chromosome set (descending by length).

    method='threshold' -> every reference scaffold >= min_scaffold_bp.
    method='dropoff'   -> cut at the sharpest adjacent length ratio among boundaries where
                          the last chromosome is >= min_scaffold_bp (floor) and the
                          cumulative genome fraction has reached dropoff_min_frac (so we
                          don't cut after a single dominant chromosome). The floor also
                          excludes deep-tail cuts. If no boundary reaches dropoff_ratio,
                          fall back to the plain threshold and flag it.

    min_chrom_frac     -> relative-size guardrail: after the cut, drop any member shorter
                          than this fraction of the median selected-chromosome length
                          (0 disables). Chromosomes in a karyotype are within an order of
                          magnitude of each other, so a member far below its peers is a
                          fragment/haplotig pulled in when two adjacent size cliffs are
                          near-equal, not a real chromosome.

    Returns (chrom_list, meta) where meta records method / n / cut ratio / genome fraction
    / flags so the choice is auditable in the report.
    """
    ordered = sorted(ref_fai, key=lambda x: (-x[1], x[0]))
    n = len(ordered)
    total = sum(L for _, L in ordered) or 1
    meta = {"method": method, "n_chrom": 0, "cut_ratio": None,
            "genome_fraction": 0.0, "flags": []}

    def threshold_set():
        return [(nm, L) for nm, L in ordered if L >= min_scaffold_bp]

    if method == "threshold" or n <= 1:
        cs = threshold_set()
        meta["method"] = "threshold"
    else:
        cum, cf = 0.0, []
        for _, L in ordered:
            cum += L
            cf.append(cum / total)
        best_k, best_ratio = None, -1.0
        for k in range(1, n):                       # ordered[k] must exist for the ratio
            last_len, nxt_len = ordered[k-1][1], ordered[k][1]
            if last_len < min_scaffold_bp or cf[k-1] < dropoff_min_frac:
                continue
            ratio = last_len / nxt_len if nxt_len > 0 else float("inf")
            if ratio > best_ratio:
                best_ratio, best_k = ratio, k
        if best_k is not None and best_ratio >= dropoff_ratio:
            cs = ordered[:best_k]
            meta.update(method="dropoff", cut_ratio=best_ratio)
        else:
            cs = threshold_set()
            meta.update(method="threshold_fallback",
                        cut_ratio=(best_ratio if best_k is not None else None))
            meta["flags"].append("no_sharp_dropoff")

    # relative-size guardrail: cs is sorted descending, so any member below the floor is a
    # trailing tail -> trim it (it falls through to unplaced). Guards the near-equal-cliff
    # case where a mid-sized scaffold gets pulled in as a spurious extra chromosome.
    if min_chrom_frac and min_chrom_frac > 0 and len(cs) > 1:
        med = statistics.median([L for _, L in cs])
        floor = min_chrom_frac * med
        kept = [(nm, L) for nm, L in cs if L >= floor]
        n_trim = len(cs) - len(kept)
        if n_trim > 0:
            meta["flags"].append(f"size_floor_trimmed={n_trim}")
            cs = kept

    meta["n_chrom"] = len(cs)
    meta["genome_fraction"] = sum(L for _, L in cs) / total
    if cs and meta["genome_fraction"] < 0.7:
        meta["flags"].append("low_genome_fraction")
    return cs, meta


def interval_coverage(b0, b1, ivals):
    """Length of [b0, b1) covered by the union of ivals (each clipped to [b0, b1))."""
    clipped = sorted((max(s, b0), min(e, b1)) for s, e in ivals if min(e, b1) > max(s, b0))
    if not clipped:
        return 0
    cov, cs, ce = 0, clipped[0][0], clipped[0][1]
    for s, e in clipped[1:]:
        if s > ce:
            cov += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return cov + (ce - cs)


def apply_containment(place, contained_frac, demote):
    """Within one assembly, find chromosome-class pieces whose reference territory on their
    chromosome is >= contained_frac covered by LARGER pieces on the same chromosome -- i.e.
    redundant haplotigs sitting inside a bigger scaffold. This is the size-independent analog
    of the reference chromosome-set cut, applied per query: it catches a large haplotig that a
    size floor cannot (a haplotig is the same size as a real fragment) while keeping genuine
    splits (whose pieces cover disjoint reference territory and so are not contained).

    demote=True  -> reclassify contained pieces as unplaced (they lose the chromosome name).
    demote=False -> keep them as chromosomes but tag 'contained(chrK)'.
    """
    by_chr = defaultdict(list)
    for nm, info in place.items():
        if info["class"] == "chromosome":
            by_chr[info["primary"]].append(nm)
    demote_to, flag_to = {}, {}                # nm -> chr  (decisions, applied after the scan)
    for k, names in by_chr.items():
        names.sort(key=lambda nm: -place[nm]["length"])     # largest keeper first
        ivals = {nm: (place[nm]["tstart"][k], place[nm]["tend"][k]) for nm in names}
        for i, nm in enumerate(names):
            b0, b1 = ivals[nm]
            span = b1 - b0
            if span <= 0:
                continue
            larger = [ivals[o] for o in names[:i]]           # every larger piece on this chr
            if interval_coverage(b0, b1, larger) / span >= contained_frac:
                (demote_to if demote else flag_to)[nm] = k
    # mutate only after all decisions are made (reading a demoted dict mid-scan would KeyError)
    for nm, k in demote_to.items():
        place[nm] = {"class": "unplaced", "length": place[nm]["length"],
                     "aligned_frac": 0.0, "contained_in": k}
    for nm, k in flag_to.items():
        place[nm].setdefault("extra_flags", []).append(f"contained(chr{k})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="TSV: id  role  fai  paf")
    ap.add_argument("--reference-id", required=True)
    ap.add_argument("--species", default="species", help="label for the report file")
    ap.add_argument("--min-scaffold-bp", type=int, default=1_000_000,
                    help="floor for a chromosome-scale scaffold; also the threshold used "
                         "when the drop-off method finds no clear boundary")
    ap.add_argument("--chromosome-set-method", choices=("dropoff", "threshold"),
                    default="dropoff",
                    help="pick the reference chromosome set by data-driven length "
                         "drop-off (default) or a plain --min-scaffold-bp threshold")
    ap.add_argument("--dropoff-ratio", type=float, default=2.0,
                    help="min adjacent length ratio to accept a drop-off boundary")
    ap.add_argument("--dropoff-min-frac", type=float, default=0.5,
                    help="cumulative genome fraction reached before a drop-off cut is "
                         "allowed (guards against cutting after a dominant chromosome)")
    ap.add_argument("--min-chrom-frac", type=float, default=0.1,
                    help="relative-size guardrail: drop a chromosome-set member shorter "
                         "than this fraction of the median chromosome length (0 disables). "
                         "Lower it for genomes with genuine small microchromosomes.")
    ap.add_argument("--batch-consensus", dest="batch_consensus", action="store_true",
                    help="when >=3 assemblies are present and a strict majority agree on a "
                         "chromosome count below the reference's, demote the reference's "
                         "smallest extra chromosomes to unplaced so refContigs match the "
                         "batch (default on)")
    ap.add_argument("--no-batch-consensus", dest="batch_consensus", action="store_false")
    ap.set_defaults(batch_consensus=True)
    ap.add_argument("--min-aligned-frac", type=float, default=0.5,
                    help="min (aligned-to-chromosomes bp / scaffold length) to place")
    ap.add_argument("--contained-frac", type=float, default=0.9,
                    help="a chromosome piece whose reference territory is >= this fraction "
                         "covered by larger pieces on the same chromosome is treated as a "
                         "contained (redundant) haplotig")
    ap.add_argument("--demote-contained", dest="demote_contained", action="store_true",
                    default=True, help="move contained pieces to unplaced (default)")
    ap.add_argument("--no-demote-contained", dest="demote_contained", action="store_false",
                    help="keep contained pieces as chromosomes but tag 'contained(chrK)'")
    ap.add_argument("--secondary-frac", type=float, default=0.2,
                    help="a 2nd chromosome contributing >= this share -> composite")
    ap.add_argument("--overlap-tol", type=int, default=0,
                    help="bp of reference overlap between two pieces before flagging")
    ap.add_argument("--graph-min-fused", type=int, default=2,
                    help="join-graph support floor: a junction must be seen fused in at "
                         "least this many voters before any rule can keep it. Guards "
                         "against a single assembly's mis-join chaining unrelated "
                         "chromosomes together through union-find")
    ap.add_argument("--graph-rule", choices=GRAPH_RULES, default="majority",
                    help="which join-graph rule the PROVISIONAL consensus chromosome map "
                         "is written under. Reporting only in this version -- scaffold "
                         "names are still derived from the reference frame; all three "
                         "rules are reported side by side regardless")
    ap.add_argument("--member-cover-frac", type=float, default=0.5,
                    help="size-invariant composite-member test: a reference chromosome "
                         "with at least this fraction of its own length inside a scaffold "
                         "is a member of that scaffold, regardless of how small its share "
                         "of the scaffold's total alignment is (0 disables, restoring the "
                         "share-only behaviour that hides small chromosomes joined to "
                         "large ones)")
    ap.add_argument("--inflated-aln-tol", type=float, default=1.05,
                    help="flag a placed scaffold when summed block length / scaffold "
                         "length exceeds this (repeat-driven double counting); diagnostic "
                         "only, does not change placement")
    # ---- voter / passenger role filter -------------------------------------------
    ap.add_argument("--voter-min-cut-ratio", type=float, default=0.0,
                    help="voter criterion: min chromosome-set drop-off cut ratio "
                         "(0 disables; --voter-require-dropoff already implies "
                         "cut_ratio >= --dropoff-ratio)")
    ap.add_argument("--voter-min-genome-frac", type=float, default=0.8,
                    help="voter criterion: min fraction of the assembly held in its own "
                         "chromosome-scale set")
    ap.add_argument("--voter-min-n50-ratio", type=float, default=0.2,
                    help="voter criterion: min scaffold N50 as a fraction of the batch "
                         "median scaffold N50")
    ap.add_argument("--voter-require-dropoff", dest="voter_require_dropoff",
                    action="store_true", default=True,
                    help="voter criterion: the assembly's own chromosome set must come "
                         "from a real length drop-off, not the threshold fallback "
                         "(default on)")
    ap.add_argument("--no-voter-require-dropoff", dest="voter_require_dropoff",
                    action="store_false")
    ap.add_argument("--min-voters", type=int, default=3,
                    help="if fewer assemblies than this pass the voter criteria, disable "
                         "the filter (everything votes) and flag it")
    ap.add_argument("--batch-consensus-action", choices=("flag", "demote"), default="flag",
                    help="what to do with a reference chromosome no voter covers: 'flag' "
                         "keeps the name and tags it (default; safe for sex-limited "
                         "chromosomes), 'demote' moves it to unplaced (pre-patch "
                         "behaviour)")
    ap.add_argument("--restricted-presence-min-individuals", type=int, default=2,
                    help="a chromosome covered by at least this many individuals but not "
                         "all of them is tagged restricted_presence (a sex-limited "
                         "chromosome candidate) rather than reference_specific")
    ap.add_argument("--concordance-exclude-reference", action="store_true", default=False,
                    help="exclude the reference from the fusion split count; a reference "
                         "that is split at a junction otherwise votes against the very "
                         "fusion that reveals it is split")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    with open(a.manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    by_id = {r["id"]: r for r in rows}
    if a.reference_id not in by_id:
        sys.exit(f"reference id {a.reference_id!r} not in manifest")

    # ---- reference chromosome set (chrN = descending-size rank) ----
    ref_fai = read_fai(by_id[a.reference_id]["fai"])
    chrom, chrom_meta = select_chromosome_set(
        ref_fai, a.min_scaffold_bp, a.chromosome_set_method,
        a.dropoff_ratio, a.dropoff_min_frac, a.min_chrom_frac)

    # ---- per-query chromosome-scale sets (computed once; reused for eligibility below
    # and for the batch-consensus check) ----
    qset_by_id = {}
    set_meta = {a.reference_id: dict(chrom_meta, flags=list(chrom_meta["flags"]))}
    fai_n50 = {a.reference_id: scaffold_n50([L for _, L in ref_fai])}
    fai_nscaf = {a.reference_id: len(ref_fai)}
    for rid, row in by_id.items():
        if rid == a.reference_id:
            continue
        qfai = read_fai(row["fai"])
        qcs, qmeta = select_chromosome_set(
            qfai, a.min_scaffold_bp, a.chromosome_set_method,
            a.dropoff_ratio, a.dropoff_min_frac, a.min_chrom_frac)
        qset_by_id[rid] = {nm for nm, _ in qcs}
        set_meta[rid] = qmeta
        fai_n50[rid] = scaffold_n50([L for _, L in qfai])
        fai_nscaf[rid] = len(qfai)

    # ---- voter / passenger roles. Each assembly's own chromosome-set metrics are already
    # computed above; this just surfaces them and applies contiguity-based thresholds. A
    # passenger is named like everyone else but does not corroborate reference chromosomes
    # or contribute to fusion concordance counts.
    metrics = {rid: {"n_chrom": set_meta[rid]["n_chrom"],
                     "method": set_meta[rid]["method"],
                     "cut_ratio": set_meta[rid]["cut_ratio"],
                     "genome_fraction": set_meta[rid]["genome_fraction"],
                     "set_flags": list(set_meta[rid]["flags"]),
                     "n_scaffolds": fai_nscaf[rid],
                     "n50": fai_n50[rid]} for rid in by_id}
    roles, role_reasons, role_flags = assign_roles(
        metrics, a.voter_min_cut_ratio, a.voter_min_genome_frac, a.voter_min_n50_ratio,
        a.voter_require_dropoff, a.min_voters, a.reference_id)
    voters = {rid for rid, v in roles.items() if v == "voter"}
    chrom_meta["flags"].extend(role_flags)
    for rid in sorted(by_id):
        if roles[rid] == "passenger":
            sys.stderr.write(
                "[harmonize %s] passenger (named, but excluded from consensus voting): "
                "%s -- %s\n" % (a.species, rid, role_reasons[rid]))
    for f in role_flags:
        sys.stderr.write("[harmonize %s] WARNING: %s\n" % (a.species, f))
    sys.stderr.write("[harmonize %s] voters=%d/%d individuals=%d\n" % (
        a.species, len(voters), len(by_id),
        len({individual_of(r) for r in voters})))

    # per-assembly chromosome-scale scaffold counts (post size-floor) -> report audit.
    # NB these are scaffold counts, so a chromosome split across two scaffolds counts twice;
    # they are recorded for provenance, not used to decide the chromosome number (see below).
    assembly_chrom_counts = {a.reference_id: len(chrom)}
    for rid in by_id:
        if rid != a.reference_id:
            assembly_chrom_counts[rid] = len(qset_by_id[rid])

    chrom_of = {n: i for i, (n, _) in enumerate(chrom, start=1)}  # ref scaffold -> chr num
    n_chrom = len(chrom)
    if n_chrom == 0:
        sys.stderr.write(
            f"[harmonize {a.species}] WARNING: reference {a.reference_id} yields no "
            f"chromosome-scale scaffolds; nothing will be placed. Lower --min-scaffold-bp "
            f"or use a chromosome-scale reference.\n")

    # ---- resolve every assembly into placements (reusable so we can re-place if the
    # batch-consensus check trims a reference chromosome) ----
    ref_len = dict(ref_fai)

    def place_all(chrom_of):
        place, present_chr, fused_pairs, covered_chr = {}, {}, {}, {}
        # chromosome number -> reference chromosome length, for the size-invariant
        # member-cover test in classify(). Rebuilt per call because chrom_of is
        # renumbered if batch-consensus demotes.
        chrom_len = {i: ref_len[n] for n, i in chrom_of.items() if n in ref_len}
        for rid, row in by_id.items():
            fai = read_fai(row["fai"])
            if rid == a.reference_id:
                p = {}
                for n, L in fai:
                    if n in chrom_of:
                        k = chrom_of[n]
                        p[n] = {"class": "chromosome", "primary": k, "members": [k],
                                "tstart": {k: 0}, "tend": {k: L}, "orient": {k: "+"},
                                "aligned_frac": 1.0, "length": L}
                    else:
                        p[n] = {"class": "unplaced", "length": L, "aligned_frac": 0.0}
            else:
                ev = parse_paf(row["paf"], chrom_of)
                # per-query chromosome-scale eligibility: only scaffolds in the query's own
                # chromosome-scale set may claim a chromosome.
                eligible = qset_by_id[rid]
                p = {}
                for n, L in fai:
                    if n in eligible:
                        info = classify(ev.get(n, {}), L, a.min_aligned_frac,
                                        a.secondary_frac, chrom_len,
                                        a.member_cover_frac)
                    else:
                        info = {"class": "unplaced", "aligned_frac": 0.0}
                    info["length"] = L
                    p[n] = info
                # demote/flag redundant haplotigs contained within a larger piece
                apply_containment(p, a.contained_frac, a.demote_contained)
            place[rid] = p
            pres, fused, cov = set(), set(), set()
            for info in p.values():
                if info["class"] == "chromosome":
                    pres.add(info["primary"])
                    cov.add(info["primary"])
                elif info["class"] == "composite":
                    # composite membership counts as COVERAGE but not as PRESENCE: the
                    # assembly carries the sequence (coverage) while not holding it as a
                    # standalone chromosome (presence), and the split/fused distinction in
                    # concordance() depends on presence staying standalone-only.
                    cov.update(info["members"])
                    for x, y in combinations(info["members"], 2):
                        fused.add(frozenset((x, y)))
            present_chr[rid] = pres
            fused_pairs[rid] = fused
            covered_chr[rid] = cov
        return place, present_chr, fused_pairs, covered_chr

    place, present_chr, fused_pairs, covered_chr = place_all(chrom_of)

    # ---- batch consensus. A reference chromosome is corroborated when a VOTER covers it,
    # either as a standalone chromosome or as a member of a composite. Composite membership
    # is the strongest corroboration available: it means the voter carries the sequence AND
    # that the reference is split exactly where the voter is joined. Counting standalone
    # chromosomes only (the pre-patch rule) scored every chromosome the reference is split
    # at as unsupported, so such chromosomes survived only when some low-quality assembly
    # happened to drop a fragment on them -- and were destroyed the moment that assembly was
    # removed from the cohort.
    #
    # Zero voter coverage no longer demotes by default. A sex-limited chromosome (W / Y) is
    # carried by one haplotype of the heterogametic sex only, so it is legitimately
    # uncorroborated in a cohort where that sex is a minority; demoting it silently deletes
    # the most interesting sequence in the assembly. Coverage is therefore counted per
    # INDIVIDUAL as well, and a clean partial partition is tagged restricted_presence -- a
    # nomination for read-depth / telomere follow-up, not a verdict.
    consensus_notes = {}
    if a.batch_consensus and len(by_id) >= 3 and n_chrom > 0:
        vote_pool = [r for r in sorted(voters) if r != a.reference_id]
        support = {k: 0 for k in chrom_of.values()}
        ind_cov = defaultdict(set)
        for rid in sorted(voters):
            for k in covered_chr[rid]:
                if k in support:
                    ind_cov[k].add(individual_of(rid))
        for rid in vote_pool:
            for k in covered_chr[rid]:
                if k in support:
                    support[k] += 1
        n_ind = len({individual_of(r) for r in voters}) or 1
        unsupported = sorted(k for k, s in support.items() if s == 0)
        for k in sorted(set(chrom_of.values())):
            n_i = len(ind_cov[k])
            if support[k] == 0 and n_i < a.restricted_presence_min_individuals:
                consensus_notes[k] = "reference_specific(ref%d:%d/%d_individuals)" % (
                    k, n_i, n_ind)
            elif n_i < n_ind:
                consensus_notes[k] = "restricted_presence(ref%d:%d/%d_individuals)" % (
                    k, n_i, n_ind)
        if unsupported:
            tags = ",".join("chr%d" % k for k in unsupported)
            chrom_meta["flags"].append(
                "batch_consensus_%s=%d(%s)" % (
                    "demoted" if a.batch_consensus_action == "demote" else "flagged",
                    len(unsupported), tags))
            sys.stderr.write(
                "[harmonize %s] batch-consensus: reference chromosome(s) %s are covered by "
                "no other voter (of %d); action=%s. A sex-limited chromosome is EXPECTED to "
                "look uncorroborated -- check the restricted_presence tags and per-"
                "chromosome read depth before treating these as fragments.\n" % (
                    a.species, tags, len(vote_pool), a.batch_consensus_action))
            if a.batch_consensus_action == "demote":
                num_to_ref = {i: n for n, i in chrom_of.items()}
                drop = {num_to_ref[k] for k in unsupported}
                chrom = [(n, L) for (n, L) in chrom if n not in drop]
                chrom_of = {n: i for i, (n, _) in enumerate(chrom, start=1)}
                n_chrom = len(chrom)
                place, present_chr, fused_pairs, covered_chr = place_all(chrom_of)
                consensus_notes = {}

    # ---- fusion concordance over VOTERS only. A passenger assembly fragments into hundreds
    # of pieces and manufactures spurious fusion pairs, which both invent phantom junctions
    # and inflate the split denominator for real ones.
    conc_pool = [r for r in sorted(voters)
                 if not (a.concordance_exclude_reference and r == a.reference_id)]

    def concordance(x, y):
        pair = frozenset((x, y))
        n_f = sum(1 for r in conc_pool if pair in fused_pairs[r])
        n_s = sum(1 for r in conc_pool if x in present_chr[r] and y in present_chr[r])
        return n_f, n_s

    # ---- assign per-chromosome part indices, names, order ----
    # ---- CONSENSUS NAMING FRAME ---------------------------------------------------------
    # The reference is now an ALIGNMENT COORDINATE SYSTEM ONLY. Scaffold names come from the
    # join graph's connected components, not from whichever assembly won reference
    # selection. A reference that is split at a junction therefore accumulates _2 parts in
    # ITSELF, instead of exporting its fragmentation into every other assembly's names.
    chrom_nums = sorted(set(chrom_of.values()))
    cand_pairs = set()
    for r in conc_pool:
        cand_pairs |= set(fused_pairs[r])
    pair_stats = {}
    for pr in cand_pairs:
        x, y = sorted(pr)
        pair_stats[(x, y)] = concordance(x, y)

    # refN = reference-frame piece (size rank in the reference); chrN = consensus
    # chromosome. Two namespaces, never the same token for both.
    rs = {i: n for n, i in chrom_of.items()}
    rl = {i: ref_len.get(n, 0) for n, i in chrom_of.items()}

    graph = {}
    for rule in GRAPH_RULES:
        edges = [(x, y) for (x, y), (n_f, n_s) in pair_stats.items()
                 if edge_keep(rule, n_f, n_s, a.graph_min_fused)]
        # rl orders consensus ids by descending total length -> chr1 is the largest
        comp_of, comp_members = components(chrom_nums, edges, rl)
        graph[rule] = {"edges": edges, "comp_of": comp_of, "comp_members": comp_members}

    frame = graph[a.graph_rule]
    comp_of = frame["comp_of"]
    comp_members = frame["comp_members"]
    # Order of reference pieces WITHIN a consensus chromosome: descending length, then
    # piece number. Arbitrary but deterministic -- the reference cannot tell us the true
    # order of pieces it failed to join, and part indices only have to be stable.
    piece_rank = {}
    for _cid, _mem in comp_members.items():
        for _r, _m in enumerate(sorted(_mem, key=lambda m: (-rl.get(m, 0), m)), start=1):
            piece_rank[_m] = _r
    n_consensus = len(comp_members)
    sys.stderr.write(
        "[harmonize %s] naming frame: rule=%s -> %d consensus chromosomes from %d "
        "reference pieces\n" % (a.species, a.graph_rule, n_consensus, len(chrom_nums)))

    results = {}
    for rid in by_id:
        p = place[rid]
        # --- reference-piece level: overlap detection only. Comparing coordinates across
        #     two different reference pieces is meaningless, so overlap stays in the piece
        #     frame with exactly the pre-Step-4 semantics.
        rclaim = defaultdict(list)
        for n, info in p.items():
            if info["class"] in ("chromosome", "composite"):
                for k in info["members"]:
                    rclaim[k].append((info["tstart"][k], n))
        overlap = defaultdict(set)
        for k, lst in rclaim.items():
            lst.sort(key=lambda x: (x[0], x[1]))
            prev_end, prev = None, None
            for ts, n in lst:
                if prev_end is not None and ts < prev_end - a.overlap_tol:
                    overlap[n].add(k)
                    overlap[prev].add(k)
                prev_end, prev = p[n]["tend"][k], n

        # --- consensus level: part indices. A scaffold claims a consensus chromosome ONCE,
        #     positioned by its earliest-ranked reference piece within that chromosome --
        #     otherwise a scaffold carrying two pieces of the same chromosome would claim it
        #     twice and take two part numbers.
        cclaim = defaultdict(list)
        for n, info in p.items():
            if info["class"] not in ("chromosome", "composite"):
                continue
            best = {}
            for k in info["members"]:
                c = comp_of.get(k)
                if c is None:
                    continue
                key = (piece_rank.get(k, 999), info["tstart"][k])
                if c not in best or key < best[c]:
                    best[c] = key
            for c, key in best.items():
                cclaim[c].append((key[0], key[1], n))
        part = defaultdict(dict)               # scaffold -> {consensus chr: index}
        for c, lst in cclaim.items():
            lst.sort()
            for idx, (_r, _t, n) in enumerate(lst, start=1):
                part[n][c] = idx

        unplaced = sorted([(n, info["length"]) for n, info in p.items()
                           if info["class"] == "unplaced"],
                          key=lambda x: (-x[1], x[0]))
        unplaced_idx = {n: i for i, (n, _) in enumerate(unplaced, start=1)}

        out = []
        for n, info in p.items():
            flags = []
            if info["class"] in ("chromosome", "composite"):
                mem = list(info["members"])
                cons = sorted({comp_of[k] for k in mem if k in comp_of})
                span = ";".join(f"ref{k}:{info['tstart'][k]}-{info['tend'][k]}"
                                for k in mem)
                if len(cons) == 1:
                    # ONE consensus chromosome. This is where the frame decoupling pays:
                    # a scaffold spanning ref8+ref17 is a whole chromosome, not a
                    # composite, and gets a plain chrN_p name that downstream
                    # /^chr[0-9]+_[0-9]+$/ matchers -- cactus_pangenome.nf refContigs,
                    # pangenome.nf, dotplot_paf.R, riparian_paf.R -- actually match.
                    c = cons[0]
                    new = f"chr{c}_{part[n][c]}"
                    rcls = "chromosome"
                    key = (0, c, part[n][c])
                    # orient by the piece contributing most alignment, and say so when the
                    # pieces disagree. Pre-Step-4 a composite was NEVER reoriented, leaving
                    # 22-26% of some assemblies unnormalised and inflating apparent SV in
                    # the pangenome.
                    prim = max(mem, key=lambda k: info["tend"][k] - info["tstart"][k])
                    orient = "rev" if info["orient"][prim] == "-" else "fwd"
                    if len(mem) > 1:
                        flags.append("spans(" + "+".join("ref%d" % k for k in mem) + ")")
                        if len({info["orient"][k] for k in mem}) > 1:
                            flags.append("orient_conflict(" + ",".join(
                                "ref%d:%s" % (k, info["orient"][k]) for k in mem) + ")")
                elif len(cons) > 1:
                    # spans several consensus chromosomes -- a genuine chimera in this
                    # frame, and the only case that still earns a '+' name
                    toks = [f"chr{c}_{part[n][c]}" for c in cons]
                    new = "+".join(toks)
                    rcls = "composite"
                    orient = "fwd"
                    key = (0, cons[0], part[n][cons[0]])
                    flags.append("fusion")
                else:
                    new, rcls, orient, key = f"unframed_{n}", "unplaced", "fwd", (2, 0, 0)
                    flags.append("no_consensus_chromosome")
                for x, y in combinations(mem, 2):
                    n_f, n_s = concordance(x, y)
                    if n_f == 0:
                        # no voter carries this junction -- only reachable for a passenger,
                        # whose fusions are unverifiable by construction
                        tag = "unsupported"
                    elif n_f > 1 or n_s == 0:
                        tag = "concordant"
                    else:
                        tag = "chimera_suspect"
                    flags.append(f"{tag}(ref{min(x, y)}+ref{max(x, y)}:{n_f}f/{n_s}s)")
                for k in mem:
                    if k in consensus_notes:
                        flags.append(consensus_notes[k])
                for k in sorted(info.get("member_reasons") or {}):
                    flags.append(info["member_reasons"][k])
                if info.get("aligned_frac", 0) > a.inflated_aln_tol:
                    flags.append("inflated_aln(bp=%.2f,cov=%.2f)" % (
                        info["aligned_frac"], info.get("cover_frac", 0.0)))
                if any(k in overlap[n] for k in mem):
                    flags.append("overlap")
                flags += info.get("extra_flags", [])
            else:
                u = unplaced_idx[n]
                new = f"unplaced_{u}"
                rcls, orient, span = "unplaced", "fwd", "-"
                key = (1, u, 0)
                if info.get("contained_in") is not None:
                    flags.append(f"contained(ref{info['contained_in']})")
                elif info.get("aligned_frac", 0) > 0:
                    flags.append(f"low_cov({info['aligned_frac']:.2f})")
            out.append({"old": n, "new": new, "orient": orient, "length": info["length"],
                        "class": rcls, "span": span,
                        "flags": ";".join(flags) if flags else "-", "_key": key})
        out.sort(key=lambda d: d["_key"])
        for i, d in enumerate(out, start=1):
            d["order"] = i
        results[rid] = out

    # ---- write outputs ----
    os.makedirs(a.outdir, exist_ok=True)
    for rid, out in results.items():
        with open(os.path.join(a.outdir, f"{rid}.harmonized_name_map.tsv"), "w") as fh:
            fh.write("old_name\tnew_name\torient\torder\tlength\tclass\tref_span\tflags\n")
            for d in out:
                fh.write("\t".join(str(x) for x in (
                    d["old"], d["new"], d["orient"], d["order"],
                    d["length"], d["class"], d["span"], d["flags"])) + "\n")

    # ---- per-assembly chromosome-set audit table. select_chromosome_set() runs on every
    # assembly but its metrics were previously discarded for everything except the
    # reference, which made the voter decision and the frame quality invisible.
    csets = os.path.join(a.outdir, f"{a.species}.chromosome_sets.tsv")
    with open(csets, "w") as fh:
        fh.write("id\trole\tis_reference\tn_chrom_set\tchrom_set_method\tcut_ratio\t"
                 "genome_fraction\tscaffold_n50\tn_scaffolds\tchrom_set_flags\t"
                 "role_reasons\n")
        for rid in sorted(by_id):
            m = metrics[rid]
            fh.write("\t".join(str(x) for x in (
                rid, roles[rid], 1 if rid == a.reference_id else 0, m["n_chrom"],
                m["method"],
                "NA" if m["cut_ratio"] is None else "%.4f" % m["cut_ratio"],
                "%.4f" % m["genome_fraction"], m["n50"], m["n_scaffolds"],
                ";".join(m["set_flags"]) if m["set_flags"] else "-",
                role_reasons[rid])) + "\n")

    # graph, rs and rl are computed above -- the naming frame is derived from them, so
    # they cannot be rebuilt here without risking the report describing a different frame
    # from the one the names came out of.
    gpath = os.path.join(a.outdir, f"{a.species}.chromosome_graph.tsv")
    with open(gpath, "w") as fh:
        fh.write("# refN = reference-frame piece (size rank in %s); chrN = consensus "
                 "chromosome. min_fused=%d\n" % (a.reference_id, a.graph_min_fused))
        fh.write("ref_a\tref_b\tref_scaffold_a\tref_scaffold_b\tlen_a\tlen_b\t"
                 "n_fused\tn_split\tn_voters\t"
                 + "\t".join("edge_%s" % r for r in GRAPH_RULES) + "\n")
        for (x, y) in sorted(pair_stats, key=lambda p: (-pair_stats[p][0], p)):
            n_f, n_s = pair_stats[(x, y)]
            fh.write("\t".join(str(v) for v in (
                "ref%d" % x, "ref%d" % y, rs.get(x, "?"), rs.get(y, "?"),
                rl.get(x, 0), rl.get(y, 0), n_f, n_s, len(conc_pool),
                *(1 if edge_keep(r, n_f, n_s, a.graph_min_fused) else 0
                  for r in GRAPH_RULES))) + "\n")

    cpath = os.path.join(a.outdir, f"{a.species}.chromosome_components.tsv")
    with open(cpath, "w") as fh:
        fh.write("# chrN = consensus chromosome, numbered by descending total length; "
                 "ties by smallest member. members are reference-frame pieces (refN, size "
                 "rank in %s). Ids are NOT comparable across rules.\n" % a.reference_id)
        fh.write("rule\tn_chromosomes\tn_edges\tchromosome\tn_members\t"
                 "members\tmember_scaffolds\ttotal_bp\n")
        for rule in GRAPH_RULES:
            g = graph[rule]
            for cid, mem in g["comp_members"].items():
                fh.write("\t".join(str(v) for v in (
                    rule, len(g["comp_members"]), len(g["edges"]),
                    "chr%d" % cid, len(mem),
                    "+".join("ref%d" % m for m in mem),
                    "+".join(rs.get(m, "?") for m in mem),
                    sum(rl.get(m, 0) for m in mem))) + "\n")

    # presence matrix: how each assembly holds each reference chromosome. Makes genuine
    # between-haplotype differences explicit rather than absorbing them silently.
    ppath = os.path.join(a.outdir, f"{a.species}.presence_matrix.tsv")
    ids_sorted = sorted(by_id)
    with open(ppath, "w") as fh:
        fh.write("# rows are reference-frame pieces (refN, size rank in %s), not consensus "
                 "chromosomes -- see chromosome_components.tsv for the mapping\n"
                 % a.reference_id)
        fh.write("ref\tref_scaffold\tlength\tn_chromosome\tn_composite\tn_absent\t"
                 "n_individuals_present\tn_individuals\t"
                 + "\t".join(ids_sorted) + "\n")
        n_ind_all = len({individual_of(r) for r in by_id})
        for k in chrom_nums:
            cells, nc, nx, na, inds = [], 0, 0, 0, set()
            for rid in ids_sorted:
                if k in present_chr[rid]:
                    parts = sum(1 for info in place[rid].values()
                                if info["class"] == "chromosome"
                                and info.get("primary") == k)
                    cells.append("chromosome" + (":%dparts" % parts if parts > 1 else ""))
                    nc += 1
                    inds.add(individual_of(rid))
                elif k in covered_chr[rid]:
                    cells.append("composite")
                    nx += 1
                    inds.add(individual_of(rid))
                else:
                    cells.append("absent")
                    na += 1
            fh.write("\t".join(str(v) for v in (
                "ref%d" % k, rs.get(k, "?"), rl.get(k, 0), nc, nx, na,
                len(inds), n_ind_all, *cells)) + "\n")

    # provisional consensus chromosome map -- the schema a per-chromosome pangenome build
    # would consume. Written under --graph-rule and labelled provisional because names are
    # NOT yet derived from it.
    gr = graph[a.graph_rule]
    mpath = os.path.join(a.outdir, f"{a.species}.consensus_chromosome_map.tsv")
    with open(mpath, "w") as fh:
        fh.write("# Consensus frame under graph rule '%s' (min_fused=%d). This IS the "
                 "naming frame: new_name is derived from consensus_chrom. refN identifies "
                 "a reference-frame piece, chrN a consensus chromosome -- one namespace "
                 "each, never the same token for both.\n"
                 % (a.graph_rule, a.graph_min_fused))
        fh.write("assembly\trole\told_name\tcurrent_new_name\tclass\tlength\torient\t"
                 "consensus_chrom\tconsensus_members\tflags\n")
        for rid in ids_sorted:
            for d in results[rid]:
                cc, cm = "-", "-"
                info = place[rid].get(d["old"], {})
                mem = info.get("members") or []
                cids = sorted({gr["comp_of"][m] for m in mem if m in gr["comp_of"]})
                # consensus ids are prefixed 'cons' so they can never be confused with
                # reference chromosome numbers -- they are a different numbering system
                if len(cids) == 1:
                    cc = "chr%d" % cids[0]
                    cm = "+".join("ref%d" % m
                                  for m in gr["comp_members"][cids[0]])
                elif len(cids) > 1:
                    # spans >1 consensus chromosome: under this rule the scaffold reads as
                    # a chimera, i.e. the rule disagrees with this assembly's topology
                    cc = "MULTI:" + "+".join("chr%d" % c for c in cids)
                    cm = "-"
                fh.write("\t".join(str(v) for v in (
                    rid, roles[rid], d["old"], d["new"], d["class"], d["length"],
                    d["orient"], cc, cm, d["flags"])) + "\n")

    for rule in GRAPH_RULES:
        g = graph[rule]
        sys.stderr.write(
            "[harmonize %s] join graph rule=%-10s edges=%2d consensus_chromosomes=%2d%s\n"
            % (a.species, rule, len(g["edges"]), len(g["comp_members"]),
               "   <- --graph-rule" if rule == a.graph_rule else ""))
    multi = sorted(m for m in graph[a.graph_rule]["comp_members"].values() if len(m) > 1)
    if multi:
        sys.stderr.write(
            "[harmonize %s] reference is split across %d consensus chromosome(s) under "
            "rule '%s': %s. The consensus frame now names it -- those pieces become _1/_2 "
            "parts of one chrN in the reference's own name map, and assemblies that made "
            "the join get a plain chrN_1.\n" % (
                a.species, len(multi), a.graph_rule,
                ", ".join("+".join("ref%d" % m for m in g) for g in multi)))

    rep = os.path.join(a.outdir, f"{a.species}.harmonization_report.tsv")
    with open(rep, "w") as fh:
        fh.write(f"# reference_id\t{a.reference_id}\n")
        fh.write(f"# n_reference_chromosomes\t{n_chrom}\n")
        fh.write(f"# chromosome_set_method\t{chrom_meta['method']}\n")
        fh.write(f"# chromosome_set_cut_ratio\t"
                 f"{('%.2f' % chrom_meta['cut_ratio']) if chrom_meta['cut_ratio'] is not None else 'NA'}\n")
        fh.write(f"# chromosome_set_genome_fraction\t{chrom_meta['genome_fraction']:.3f}\n")
        fh.write(f"# chromosome_set_flags\t{';'.join(chrom_meta['flags']) if chrom_meta['flags'] else '-'}\n")
        fh.write(f"# params\tmin_scaffold_bp={a.min_scaffold_bp}\t"
                 f"dropoff_ratio={a.dropoff_ratio}\tdropoff_min_frac={a.dropoff_min_frac}\t"
                 f"min_chrom_frac={a.min_chrom_frac}\tbatch_consensus={a.batch_consensus}\t"
                 f"contained_frac={a.contained_frac}\tdemote_contained={a.demote_contained}\t"
                 f"min_aligned_frac={a.min_aligned_frac}\t"
                 f"secondary_frac={a.secondary_frac}\t"
                 f"batch_consensus_action={a.batch_consensus_action}\t"
                 f"voter_min_genome_frac={a.voter_min_genome_frac}\t"
                 f"voter_min_n50_ratio={a.voter_min_n50_ratio}\t"
                 f"voter_min_cut_ratio={a.voter_min_cut_ratio}\t"
                 f"voter_require_dropoff={a.voter_require_dropoff}\t"
                 f"min_voters={a.min_voters}\t"
                 f"restricted_presence_min_individuals="
                 f"{a.restricted_presence_min_individuals}\t"
                 f"concordance_exclude_reference={a.concordance_exclude_reference}\t"
                 f"member_cover_frac={a.member_cover_frac}\t"
                 f"inflated_aln_tol={a.inflated_aln_tol}\t"
                 f"graph_rule={a.graph_rule}\t"
                 f"graph_min_fused={a.graph_min_fused}\n")
        fh.write("# assembly_chromosome_counts\t"
                 + ";".join(f"{k}={v}" for k, v in assembly_chrom_counts.items()) + "\n")
        fh.write(f"# voters\t{len(voters)}/{len(by_id)}\t"
                 + ";".join(f"{r}={roles[r]}" for r in sorted(by_id)) + "\n")
        for rule in GRAPH_RULES:
            g = graph[rule]
            fh.write("# chromosome_graph\t%s\tedges=%d\tconsensus_chromosomes=%d\t%s\n"
                     % (rule, len(g["edges"]), len(g["comp_members"]),
                        ";".join("+".join("ref%d" % m for m in mem)
                                 for mem in g["comp_members"].values()
                                 if len(mem) > 1) or "-"))
        fh.write("# chromosome_graph_rule_reported\t%s\n" % a.graph_rule)
        fh.write("# n_consensus_chromosomes\t%d\tfrom\t%d\treference_pieces\trule\t%s\n"
                 % (n_consensus, len(chrom_nums), a.graph_rule))
        # how many '+' names each rule would leave. A '+' name now means a scaffold spans
        # several CONSENSUS chromosomes, i.e. a genuine chimera -- and it is the thing that
        # breaks /^chr[0-9]+_[0-9]+$/ downstream, so the count is worth stating per rule.
        for rule in GRAPH_RULES:
            co = graph[rule]["comp_of"]
            n_plus = sum(1 for rid in by_id for _n, inf in place[rid].items()
                         if inf["class"] in ("chromosome", "composite")
                         and len({co[k] for k in inf["members"] if k in co}) > 1)
            fh.write("# composite_names_under\t%s\tchromosomes\t%d\tplus_names\t%d\n"
                     % (rule, len(graph[rule]["comp_members"]), n_plus))
        fh.write("# frame_vocabulary\tnew_name/ref_span/flags=reference-frame chrN_p "
                 "(interim)\tchromosome_graph+components=refN pieces -> chrN consensus\n")
        fh.write("# chromosome_consensus_notes\t"
                 + (";".join(consensus_notes[k] for k in sorted(consensus_notes))
                    if consensus_notes else "-") + "\n")
        fh.write("# assembly_chromosome_set_header\tid\trole\tn_chrom_set\tmethod\t"
                 "cut_ratio\tgenome_fraction\tscaffold_n50\tchrom_set_flags\t"
                 "role_reasons\n")
        for rid in sorted(by_id):
            m = metrics[rid]
            fh.write("# assembly_chromosome_set\t" + "\t".join(str(x) for x in (
                rid, roles[rid], m["n_chrom"], m["method"],
                "NA" if m["cut_ratio"] is None else "%.2f" % m["cut_ratio"],
                "%.3f" % m["genome_fraction"], m["n50"],
                ";".join(m["set_flags"]) if m["set_flags"] else "-",
                role_reasons[rid])) + "\n")
        fh.write("assembly\told_name\tnew_name\tclass\tlength\torient\tref_span\tflags\n")
        for rid, out in results.items():
            for d in out:
                fh.write("\t".join(str(x) for x in (
                    rid, d["old"], d["new"], d["class"], d["length"],
                    d["orient"], d["span"], d["flags"])) + "\n")

    n_flagged = sum(1 for out in results.values() for d in out if d["flags"] != "-")
    sys.stderr.write(
        f"[harmonize {a.species}] reference={a.reference_id} "
        f"consensus_chromosomes={n_consensus} reference_pieces={n_chrom} "
        f"({chrom_meta['method']}, "
        f"{chrom_meta['genome_fraction']:.2f} of genome) "
        f"assemblies={len(results)} "
        f"voters={len(voters)} passengers={len(by_id) - len(voters)} "
        f"flagged_scaffolds={n_flagged}\n")


if __name__ == "__main__":
    main()
