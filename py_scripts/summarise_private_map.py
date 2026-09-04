#!/usr/bin/env python3
# ======================================================================================
# summarise_private_map.py
#
# Turns a PAF of one haplotype's private segments aligned against all assemblies into a
# per-segment table answering: does this "private" sequence exist elsewhere?
#
# WHY THIS EXISTS
# ---------------
# "Private" is defined by GRAPH COVERAGE -- a node walked by exactly one haplotype. That is
# not the same as sequence absent from other assemblies. A segment can be private in the
# graph because the aligner failed to merge it, or because minigraph-cactus collapsed a
# repeat, while the sequence itself sits plainly in three other assemblies. Nothing else in
# the pipeline can tell those apart, and every private-sequence figure depends on the
# distinction.
#
# The live question it is aimed at: chr8 and chr9 carry 23.6% and 22.9% private sequence
# against an ~11% floor across other chromosomes, consistently across all ten haplotypes and
# both graph flavours. Real divergence, or collapsed repeat? If those segments align cleanly
# to other assemblies, it is collapse.
#
# TWO THINGS DONE CAREFULLY
# -------------------------
# 1. SELF-HITS ARE EXCLUDED BY ASSEMBLY TAG, NOT BY IDENTITY. A private segment is present in
#    its own assembly by construction and will align there at ~100%. Filtering on identity
#    would also discard genuine near-identical copies elsewhere -- exactly the collapsed-repeat
#    case we are trying to detect. Target names are <assembly_id>::<contig> and the query's
#    own assembly id is passed in explicitly, because it cannot be recovered from the PanSN
#    haplotype name (sampleOf() maps dots to underscores and is not invertible).
#
# 2. ALIGNED FRACTION USES MERGED QUERY INTERVALS. Summing PAF block lengths double-counts
#    overlapping alignments, which is the same error as measuring an alignment envelope
#    instead of a merged footprint -- it has appeared three times in this project already, so
#    it is merged here per (segment, target assembly) and again across assemblies.
#
# Dependencies: python3 stdlib.
#
# USAGE
#   summarise_private_map.py --paf x.paf.gz --self-assembly Sde-CBau_104_hap1 \
#       --haplotype 'Sde-CBau_104#1' --label 373251.full --outdir .
# ======================================================================================

import argparse
import gzip
import os
import sys
from collections import defaultdict


def popen(p):
    with open(p, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return (gzip.open if gz else open)(p, "rt", encoding="utf-8", errors="replace")


def merged_len(iv):
    """Union length of [start, end) intervals. NOT the sum of block lengths."""
    if not iv:
        return 0
    iv = sorted(iv)
    total, s, e = 0, iv[0][0], iv[0][1]
    for a, b in iv[1:]:
        if a > e:
            total += e - s
            s, e = a, b
        elif b > e:
            e = b
    return total + (e - s)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True)
    p.add_argument("--fasta", required=True,
                   help="the query FASTA. Read for its headers ONLY, so that every segment "
                        "gets a row -- minimap2 omits a query entirely when it aligns "
                        "nowhere, including to its own assembly, and 41 segments on this "
                        "cohort did exactly that. Emitting only PAF-observed queries left "
                        "them absent from the table and broke the downstream join.")
    p.add_argument("--self-assembly", required=True,
                   help="assembly id of the query haplotype; its hits are EXCLUDED")
    p.add_argument("--haplotype", required=True, help="PanSN haplotype key, for the output")
    p.add_argument("--set", dest="segset", default="private", choices=["private", "control"],
                   help="which sequence set this FASTA is. Carried into the table so the "
                        "join can contrast private against its matched control -- the whole "
                        "point of the control is that it is measured identically.")
    p.add_argument("--label", required=True)
    p.add_argument("--outdir", required=True)
    # 0 = DERIVE FROM THE CONTROL. Same convention as --single-copy in
    # summarise_private_kmer.py, so the two behave alike.
    p.add_argument("--min-identity", type=float, default=0.0,
                   help="gap-compressed identity floor for a hit to count. 0 = derive from "
                        "this run's CONTROL set (only valid with --set control). Control "
                        "windows are non-private by construction, so the identity "
                        "distribution of their hits IS the empirical cross-haplotype "
                        "homology identity for this cohort -- 0.80 fits this clupeid but a "
                        "plant or copepod cohort would differ, and a hardcoded value would "
                        "be wrong invisibly.")
    p.add_argument("--min-frac", type=float, default=0.0,
                   help="merged aligned fraction at or above which a segment is NOT_PRIVATE. "
                        "0 = derive from the CONTROL. Gives the cutoff a defensible meaning: "
                        "'covered at least as well as most known-shared sequence'.")
    p.add_argument("--derive-percentile", type=float, default=5.0,
                   help="percentile of the control distribution to use as each floor "
                        "(default 5.0 = retain 95%% of known-true homology)")
    p.add_argument("--max-control-no-hit-frac", type=float, default=0.5,
                   help="if more than this fraction of CONTROL segments have no non-self hit "
                        "at all, the thresholds are underivable and the run FAILS rather "
                        "than falling back to a guess (default 0.5). This is the check that "
                        "would have caught both the split minimap2 index and the 0.90 "
                        "identity floor immediately.")
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    stem = "%s.%s.%s" % (a.label, a.haplotype.replace("#", "_").replace("/", "_"), a.segset)

    # every segment in the FASTA, in file order, so the table is COMPLETE
    all_q = []
    fa_len = {}
    with open(a.fasta, encoding="utf-8", errors="replace") as fh:
        name = None
        n = 0
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    fa_len[name] = n
                name = line[1:].strip().split()[0] if len(line) > 1 else ""
                all_q.append(name)
                n = 0
            else:
                n += len(line.strip())
        if name is not None:
            fa_len[name] = n

    # ---- derive thresholds from the CONTROL, if asked -------------------------------
    # Two passes over the control PAF: the first with NO identity filter, to measure what
    # cross-haplotype homology actually looks like in this cohort; the second applies the
    # derived floors. Private runs receive those same values so the two sets stay comparable.
    derive = (a.min_identity <= 0.0) or (a.min_frac <= 0.0)
    derived_note = "supplied"
    if derive and a.segset != "control":
        sys.exit("ERROR: --min-identity/--min-frac of 0 means DERIVE, which is only valid "
                 "with --set control. For the private set, pass the values the control run "
                 "wrote into its audit (derived_min_identity / derived_min_frac).")
    if derive:
        best_id_c = {}
        merged_c = defaultdict(list)
        clen = {}
        n_seen = set()
        with popen(a.paf) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
                if f[5].split("::", 1)[0] == a.self_assembly:
                    continue
                nm, bl = int(f[9]), int(f[10])
                clen[q] = ql
                n_seen.add(q)
                idt = (nm / bl) if bl else 0.0
                if idt > best_id_c.get(q, 0.0):
                    best_id_c[q] = idt
                merged_c[q].append((qs, qe))

        # every control segment, including those with no non-self hit at all
        n_ctrl = 0
        with open(a.fasta, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    n_ctrl += 1
        n_nohit = n_ctrl - len(n_seen)
        frac_nohit = (n_nohit / n_ctrl) if n_ctrl else 1.0
        sys.stderr.write("[private_map] derive: %d control segments, %d with no non-self hit "
                         "(%.3f)\n" % (n_ctrl, n_nohit, frac_nohit))
        if frac_nohit > a.max_control_no_hit_frac:
            sys.stderr.write(
                "ERROR: %.1f%% of control segments have NO non-self alignment. Control "
                "sequence is non-private by construction, so this cannot be right -- the "
                "index is probably wrong (e.g. split into parts, or missing assemblies) or "
                "the control windows are not what they claim. Thresholds are not derivable "
                "from this; refusing to guess.\n" % (100.0 * frac_nohit))
            sys.exit(1)

        def pct(vals, q):
            if not vals:
                return 0.0
            vs = sorted(vals)
            i = max(0, min(len(vs) - 1, int(round((q / 100.0) * (len(vs) - 1)))))
            return vs[i]

        if a.min_identity <= 0.0:
            a.min_identity = round(pct(list(best_id_c.values()), a.derive_percentile), 4)
        if a.min_frac <= 0.0:
            fr = [min(1.0, merged_len(merged_c[q]) / clen[q]) for q in merged_c if clen.get(q)]
            a.min_frac = round(pct(fr, a.derive_percentile), 4)
        derived_note = "derived_p%.1f" % a.derive_percentile
        sys.stderr.write("[private_map] derived min_identity=%.4f min_frac=%.4f "
                         "(p%.1f of %d control segments)\n"
                         % (a.min_identity, a.min_frac, a.derive_percentile, len(best_id_c)))

    qlen = {}
    # (query, target assembly) -> merged query intervals
    per_asm = defaultdict(list)
    best_id = defaultdict(float)
    n_rows = n_self = n_lowid = n_kept = 0
    # identity histogram over NON-SELF rows, kept or not. Without it a badly set
    # --min-identity is invisible: the post-filter counts look like "these segments have no
    # homologue elsewhere" when the real answer is "every homologue was filtered out". Observed
    # cross-haplotype identity for shared sequence in this species is 0.75-0.90, so a 0.90
    # cutoff discarded essentially all of it and 70% of CONTROL windows -- non-private by
    # construction -- came back with n_other_assemblies = 0.
    id_hist = defaultdict(int)
    seen_q = set()

    with popen(a.paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            n_rows += 1
            q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            tname = f[5]
            nmatch, blen = int(f[9]), int(f[10])
            qlen[q] = ql
            seen_q.add(q)

            asm = tname.split("::", 1)[0]
            if asm == a.self_assembly:
                n_self += 1
                continue
            ident = (nmatch / blen) if blen else 0.0
            id_hist["%.2f" % (int(ident * 20) / 20.0)] += 1
            if ident < a.min_identity:
                n_lowid += 1
                continue
            n_kept += 1
            per_asm[(q, asm)].append((qs, qe))
            if ident > best_id[q]:
                best_id[q] = ident

    # every segment in the FASTA should appear; those with no PAF row aligned nowhere
    fout = os.path.join(a.outdir, "%s.private_map.tsv" % stem)
    with open(fout, "w") as out:
        out.write("# per private segment: does this sequence exist in ANY OTHER assembly?\n")
        out.write("# Self-hits excluded by assembly tag (%s), not by identity -- a genuine\n"
                  % a.self_assembly)
        out.write("#   near-identical copy elsewhere is precisely the collapsed-repeat case\n")
        out.write("#   this is meant to detect, so an identity filter would hide it.\n")
        out.write("# aligned_frac is MERGED query coverage across all other assemblies;\n")
        out.write("#   summing PAF block lengths would double-count overlaps.\n")
        out.write("# verdict PRIVATE_CONFIRMED = no other assembly covers >= %.2f of it.\n"
                  % a.min_frac)
        out.write("haplotype\tset\tsegment\tcontig\tstart\tend\tsegment_bp\t"
                  "n_other_assemblies\tbest_identity\taligned_frac_merged\t"
                  "max_frac_single_assembly\tverdict\n")

        rows = 0
        n_no_paf = 0
        for q in all_q:
            # length from the PAF when available, else from the FASTA -- a segment with no
            # PAF row still has a real length and must not report 0
            ql = qlen.get(q, fa_len.get(q, 0))
            if q not in seen_q:
                n_no_paf += 1
            # header is <haplotype>|<contig>|<start>|<end>|seg<N>
            parts = q.split("|")
            contig = parts[1] if len(parts) > 3 else "."
            start = parts[2] if len(parts) > 3 else "0"
            end = parts[3] if len(parts) > 3 else "0"
            seg = parts[4] if len(parts) > 4 else q

            asms = [k[1] for k in per_asm if k[0] == q]
            all_iv = []
            per_max = 0.0
            for asm in asms:
                iv = per_asm[(q, asm)]
                all_iv.extend(iv)
                fr = merged_len(iv) / ql if ql else 0.0
                per_max = max(per_max, fr)
            frac = (merged_len(all_iv) / ql) if ql else 0.0
            # no PAF row at all means it aligned nowhere, not even to itself. That is not
            # the same as "aligned nowhere else", so it gets its own verdict rather than
            # being silently counted as confirmed-private.
            if q not in seen_q:
                verdict = "NO_ALIGNMENT"
            elif frac >= a.min_frac:
                verdict = "NOT_PRIVATE"
            else:
                verdict = "PRIVATE_CONFIRMED"
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%s\n"
                      % (a.haplotype, a.segset, seg, contig, start, end, ql,
                         len(asms), best_id.get(q, 0.0), min(1.0, frac),
                         min(1.0, per_max), verdict))
            rows += 1

    with open(os.path.join(a.outdir, "%s.private_map_audit.tsv" % stem), "w") as out:
        out.write("metric\tvalue\n")
        out.write("haplotype\t%s\nset\t%s\nself_assembly\t%s\n"
                  % (a.haplotype, a.segset, a.self_assembly))
        out.write("paf_rows\t%d\n" % n_rows)
        out.write("rows_self_excluded\t%d\n" % n_self)
        out.write("rows_below_identity_%.2f\t%d\n" % (a.min_identity, n_lowid))
        out.write("rows_kept\t%d\n" % n_kept)
        out.write("derived_min_identity\t%.4f\n" % a.min_identity)
        out.write("derived_min_frac\t%.4f\n" % a.min_frac)
        out.write("threshold_source\t%s\n" % derived_note)
        out.write("min_identity_threshold\t%.4f\n" % a.min_identity)
        out.write("# derived_* are what the PRIVATE run must be given, so both sets are\n")
        out.write("# filtered identically. threshold_source says whether they came from the\n")
        out.write("# control distribution or were supplied by hand.\n")
        for b in sorted(id_hist):
            out.write("identity_ge_%s\t%d\n" % (b, id_hist[b]))
        out.write("# identity_ge_* is the gap-compressed identity of every NON-SELF row, in\n")
        out.write("# 0.05 bins, BEFORE the min_identity filter. If the mass sits below the\n")
        out.write("# threshold the filter is the reason segments look private, not the data.\n")
        out.write("segments_in_fasta\t%d\n" % len(all_q))
        out.write("segments_with_any_paf_row\t%d\n" % len(seen_q))
        out.write("segments_with_no_alignment\t%d\n" % n_no_paf)
        out.write("# segments_with_no_alignment aligned NOWHERE, not even to their own\n")
        out.write("# assembly -- low-complexity or masked sequence, usually. They are\n")
        out.write("# reported as NO_ALIGNMENT rather than omitted, because omitting them\n")
        out.write("# left the downstream join with rows present in one stream and not the\n")
        out.write("# other, which looks identical to a failed task.\n")
        out.write("# segments_with_any_paf_row counts only segments minimap2 reported at all.\n")
        out.write("# Segments absent from the PAF aligned NOWHERE, including to their own\n")
        out.write("# assembly, which should not happen -- rows_self_excluded near zero is a\n")
        out.write("# warning that the --self-assembly tag does not match the index tags.\n")

    sys.stderr.write("[private_map] %s: %d PAF rows, %d self-excluded, %d kept, %d segments\n"
                     % (a.haplotype, n_rows, n_self, n_kept, len(seen_q)))
    if n_rows and n_self == 0:
        sys.stderr.write("WARNING: no self-hits excluded -- --self-assembly '%s' may not "
                         "match the index's assembly tags\n" % a.self_assembly)


if __name__ == "__main__":
    main()
