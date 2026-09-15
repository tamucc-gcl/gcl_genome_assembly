#!/usr/bin/env python3
# ======================================================================================
# pangenome_input_coverage.py
#
# How much of each input haplotype actually made it into the pangenome graph, per
# chromosome -- and, where a chromosome is missing, why.
#
# WHY THIS EXISTS
# ---------------
# Both halves of this answer were already being computed and neither was surfaced:
#
#   cactus writes  <outName>.stats.tgz  containing sample-stats.tsv: bp of each haplotype
#                  present in each reference chromosome's graph.
#   harmonization  writes the report, flagging composite scaffolds with a cross-haplotype
#                  concordance vote (`chimera_suspect(ref5+ref9:1f/7s)` = one other
#                  haplotype carries this junction, seven keep those chromosomes apart).
#
# Nothing connected them, so this went unnoticed for every pangenome run built so far:
#
#   chr9_1  Sde-CTlk_104#1  1,281,066 bp        -- 1.9% of the 68 Mb chromosome median
#
# `Sde-CTlk_104_hap1` has no chr9 scaffold at all. Its chromosome 9 is fused into a
# 111.6 Mb scaffold named `chr5_1+chr9_1`, and `cactus-graphmap-split` assigns each contig
# to a SINGLE chromosome -- so the whole scaffold went to chr5's graph, 45.2 Mb of it
# aligned there, and 73.6 Mb of genuine chromosome 9 is absent from every graph.
#
# Cactus does not break chimeric contigs. Inspector cannot see this either: it works from
# read-to-contig alignment and a Hi-C scaffold join is an N-gap with no reads spanning it,
# so there is nothing to disagree with. It found two structural errors on this haplotype,
# neither on the fused scaffold.
#
# THE JOIN IS THE POINT
# ---------------------
# A depressed chromosome count on its own says something is wrong. A composite flag on its
# own says a scaffold spans two reference chromosomes. Together they say: THIS haplotype's
# chrB is missing from the graph BECAUSE it is fused to chrA, and N other haplotypes keep
# them separate. That is a statement someone can act on.
#
# THE CHROMOSOME-SCALE FLOOR IS DERIVED, NOT SET
# ----------------------------------------------
# Deciding which scaffolds are large enough to matter reuses harmonization's own dropoff
# rule rather than a second hardcoded threshold that could drift from it:
#
#   sort descending; cut at the sharpest adjacent length ratio (>= dropoff_ratio) among
#   boundaries where the last member is >= min_scaffold_bp AND cumulative genome fraction
#   has reached dropoff_min_frac; then drop members below min_chrom_frac x the median of
#   the kept set. Fall back to a plain threshold if no boundary is sharp enough.
#
# Applied to the REFERENCE, giving ONE floor for the cohort. Applied per assembly it would
# pick a lower cut in a fragmented assembly -- `Sde-CPla_115_hap1` has 76 composites of
# 1-14 Mb, and some would qualify as "chromosome-scale for this assembly" while being
# nothing of the kind. The reference frame is also the frame harmonization uses to name
# composites in the first place, so the two agree by construction.
#
# Dependencies: python3 stdlib.
#
# USAGE
#   pangenome_input_coverage.py --sample-stats sample-stats.tsv \
#       --harmonization <taxid>.harmonization_report.tsv --ref-fai <reference>.fasta.fai \
#       --label <taxid> --outdir . [--low-frac 0.5] [--min-scaffold-bp 1000000]
# ======================================================================================

import argparse
import os
import re
import statistics
import sys

FLAG_VOTE = re.compile(r"(chimera_suspect|concordant|unsupported)\(ref(\d+)\+ref(\d+):(\d+)f/(\d+)s\)")


def read_rows(path, sep="\t"):
    """'#'-commented table -> list of dict rows.

    Comments are stripped by hand and comment.char stays OFF: '#' is a legitimate
    character inside a PanSN haplotype name (Sde-CBau_104#1), and treating it as a comment
    marker truncates every key at the separator and shifts every column left.
    """
    out = []
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return out
    hdr = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split(sep)
            if hdr is None:
                hdr = f
                continue
            if len(f) != len(hdr):
                continue
            out.append(dict(zip(hdr, f)))
    return out


def select_chromosome_set(ref_fai, min_scaffold_bp, method, dropoff_ratio, dropoff_min_frac,
                          min_chrom_frac=0.0):
    """Verbatim port of harmonize_names.py's rule. Kept identical on purpose.

    Duplicating it is deliberate: this script must run BEFORE the pangenome is built and
    cannot import from a module that lives beside it in the pipeline, but the two must not
    disagree about what a chromosome is. If harmonize_names.py's version changes, change
    this one in the same commit -- the audit records the derived floor so a divergence shows
    up as a changed number rather than as silence.
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
        for k in range(1, n):
            last_len, nxt_len = ordered[k - 1][1], ordered[k][1]
            if last_len < min_scaffold_bp or cf[k - 1] < dropoff_min_frac:
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

    if min_chrom_frac and min_chrom_frac > 0 and len(cs) > 1:
        med = statistics.median([L for _, L in cs])
        floor = min_chrom_frac * med
        kept = [(nm, L) for nm, L in cs if L >= floor]
        n_trim = len(cs) - len(kept)
        if n_trim > 0:
            meta["flags"].append("size_floor_trimmed=%d" % n_trim)
        cs = kept

    meta["n_chrom"] = len(cs)
    meta["genome_fraction"] = sum(L for _, L in cs) / total
    return cs, meta


def chrom_of(contig):
    """chr10_1 -> chr10; chr10_17+chr11_12 -> chr10 (the first member)."""
    c = str(contig).split("#")[-1]
    return c.split("+")[0].split("_")[0]


def indiv_of(hap):
    """PanSN haplotype -> individual. Sde-CBau_104#1 -> Sde-CBau_104;
    Sde-CMat_203_hap1#0 -> Sde-CMat_203. Grouping only, so the irreversible dot
    substitution in the PanSN naming does not matter here."""
    base = str(hap).split("#")[0]
    return re.sub(r"_hap[0-9]+$", "", base)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample-stats", required=True,
                   help="sample-stats.tsv from cactus's stats bundle: "
                        "<ref_chrom> <haplotype> <bp_in_graph>")
    p.add_argument("--harmonization", default="",
                   help="harmonization report. Optional: without it the coverage table is "
                        "still emitted, but a depressed chromosome has no explanation.")
    p.add_argument("--ref-fai", required=True,
                   help="reference .fai, for deriving the chromosome-scale floor")
    p.add_argument("--label", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--low-frac", type=float, default=0.5,
                   help="flag a (haplotype, chromosome) whose bp is below this fraction of "
                        "that chromosome's median across haplotypes (default 0.5). The "
                        "chr9 case is at 0.019, so this is not a tight threshold -- it is "
                        "meant to catch whole missing chromosomes, not size variation.")
    # harmonization's own defaults, so the derived floor matches
    p.add_argument("--min-scaffold-bp", type=int, default=1000000)
    p.add_argument("--chromosome-set-method", choices=("dropoff", "threshold"),
                   default="dropoff")
    p.add_argument("--dropoff-ratio", type=float, default=2.0)
    p.add_argument("--dropoff-min-frac", type=float, default=0.5)
    p.add_argument("--min-chrom-frac", type=float, default=0.0)
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    op = lambda s: os.path.join(a.outdir, a.label + s)

    # ---- the chromosome-scale floor, from the reference -----------------------------
    fai = []
    with open(a.ref_fai, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].isdigit():
                fai.append((f[0], int(f[1])))
    if not fai:
        sys.exit("ERROR: no usable rows in %s" % a.ref_fai)
    cs, meta = select_chromosome_set(fai, a.min_scaffold_bp, a.chromosome_set_method,
                                     a.dropoff_ratio, a.dropoff_min_frac, a.min_chrom_frac)
    floor_bp = min(L for _, L in cs) if cs else a.min_scaffold_bp
    ref_len = {chrom_of(nm): L for nm, L in cs}
    sys.stderr.write("[input_coverage] chromosome set: %d members by %s, floor %d bp, "
                     "genome fraction %.3f\n"
                     % (meta["n_chrom"], meta["method"], floor_bp, meta["genome_fraction"]))

    # ---- coverage: bp of each haplotype in each chromosome's graph -------------------
    rows = []
    with open(a.sample_stats, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.lstrip().startswith("#"):
                continue
            f = line.split()
            if len(f) >= 3 and f[2].isdigit():
                rows.append((f[0], f[1], int(f[2])))
    if not rows:
        sys.exit("ERROR: no usable rows in %s" % a.sample_stats)

    by_chrom, by_hap = {}, {}
    for c, h, bp in rows:
        by_chrom.setdefault(chrom_of(c), {})[h] = by_chrom.setdefault(chrom_of(c), {}).get(h, 0) + bp
        by_hap[h] = by_hap.get(h, 0) + bp
    haps = sorted(by_hap)

    # ---- harmonization: composites and their concordance votes ----------------------
    # key: (assembly, chromosome absorbed INTO a composite) -> the evidence
    absorbed = {}
    comp_rows = []
    for r in read_rows(a.harmonization):
        vals = list(r.values())
        asm = vals[0] if vals else ""
        name = ""
        for v in vals:
            if "+chr" in str(v):
                name = str(v)
                break
        if not name:
            continue
        flags = vals[-1] if vals else ""
        members = [m for m in name.split("+") if m.startswith("chr")]
        try:
            span = max(int(v) for v in vals if str(v).isdigit())
        except ValueError:
            span = 0
        votes = {}
        for m in FLAG_VOTE.finditer(str(flags)):
            tag, x, y, nf, ns = m.groups()
            votes["ref%s+ref%s" % (x, y)] = (tag, int(nf), int(ns))
        comp_rows.append((asm, name, span, members, votes, str(flags)))
        # every member EXCEPT the first is at risk of being lost: cactus assigns the whole
        # contig to one chromosome, and the name's first member is what it is sorted under
        for m in members[1:]:
            absorbed.setdefault((asm, chrom_of(m)), []).append((name, span, votes, str(flags)))

    # PanSN haplotype <-> assembly id. The PanSN key is <individual>#<hap> for most and
    # <assembly>#0 for the reference individual, so match on the individual and the hap
    # number rather than trying to invert the naming.
    def asm_for(hap):
        ind, _, h = str(hap).partition("#")
        ind_only = re.sub(r"_hap[0-9]+$", "", ind)
        cands = [asm for asm, *_ in comp_rows if asm.startswith(ind_only)]
        if not cands:
            return ""
        if ind.endswith(tuple("_hap%d" % i for i in range(1, 10))):
            return ind                          # reference individual: PanSN base IS the id
        want = "_hap%s" % h
        for asm in cands:
            if asm.endswith(want):
                return asm
        # NO FALLBACK. Returning cands[0] here made every haplotype of an individual
        # inherit its sibling's composites -- the chr9 row for Sde-CTlk_104#2 showed the
        # chr5+chr9 fusion, which is hap1's alone. An unmatched haplotype has no composite.
        return ""

    # ---- the joined table ------------------------------------------------------------
    n_low = n_explained = 0
    with open(op(".input_coverage.tsv"), "w") as out:
        out.write("# bp of each input haplotype present in each reference chromosome's "
                  "graph, from cactus sample-stats.tsv, joined to the harmonization\n")
        out.write("# composite flags. A chromosome far below its cohort median with a\n")
        out.write("# composite explanation is a FUSED scaffold: cactus assigns a contig to\n")
        out.write("# a single chromosome, so the other member's sequence is lost.\n")
        out.write("# vote n_f/n_s: other haplotypes carrying this junction / keeping the\n")
        out.write("#   chromosomes separate. Low n_f with high n_s means a scaffolding error.\n")
        out.write("# chromosome_scale: does the reference chromosome clear the derived floor\n")
        out.write("#   (%d bp, %s over %d members)? Only these are candidates for breaking.\n"
                  % (floor_bp, meta["method"], meta["n_chrom"]))
        out.write("chromosome\thaplotype\tindividual\tbp_in_graph\tchrom_median_bp\t"
                  "frac_of_median\tref_chrom_bp\tchromosome_scale\tstatus\t"
                  "composite\tcomposite_bp\tvote\tflags\n")
        for c in sorted(by_chrom, key=lambda x: (len(x), x)):
            vals = sorted(by_chrom[c].values())
            med = statistics.median(vals) if vals else 0
            for h in haps:
                bp = by_chrom[c].get(h, 0)
                frac = (bp / med) if med else 0.0
                is_scale = c in ref_len and ref_len[c] >= floor_bp
                asm = asm_for(h)
                exp = absorbed.get((asm, c), [])
                low = frac < a.low_frac
                if low:
                    n_low += 1
                if low and exp:
                    n_explained += 1
                status = ("LOW_EXPLAINED_BY_COMPOSITE" if (low and exp)
                          else "LOW_UNEXPLAINED" if low else "ok")
                comp = exp[0] if exp else None
                vote = ""
                if comp:
                    vote = ";".join("%s:%s(%df/%ds)" % (k, v[0], v[1], v[2])
                                    for k, v in sorted(comp[2].items()))
                out.write("%s\t%s\t%s\t%d\t%d\t%.4f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                          % (c, h, indiv_of(h), bp, med, frac,
                             ref_len.get(c, "NA"), "yes" if is_scale else "no", status,
                             comp[0] if comp else ".", comp[1] if comp else ".",
                             vote or ".", (comp[3] if comp else ".")))

    # ---- break candidates: composites on chromosome-scale members only ---------------
    # The size floor does two jobs. It keeps the 111.6 Mb chr5+chr9 fusion in and the
    # 1-14 Mb composites of a fragmented assembly out -- and that is the right outcome,
    # because breaking a 2 Mb chimeric fragment accomplishes nothing: both halves stay
    # unplaced. `Sde-CPla_115` has 135 such composites across its two haplotypes.
    with open(op(".break_candidates.tsv"), "w") as out:
        out.write("# composites whose members are chromosome-scale in the REFERENCE frame.\n")
        out.write("# Candidates only -- nothing is broken by this script. Review, then write\n")
        out.write("# the breakpoints file for the second round.\n")
        out.write("# vote n_f/n_s is the primary evidence: it is the only signal independent\n")
        out.write("#   of how this scaffold was built. Hi-C cannot be used to validate a\n")
        out.write("#   break of a Hi-C-made join -- measured on the chr5+chr9 fusion,\n")
        out.write("#   cross-junction contact was 1.209x matched distance, i.e. ELEVATED,\n")
        out.write("#   because the junction sits in subtelomeric repeat. That is how the\n")
        out.write("#   error was made, not evidence against breaking.\n")
        out.write("# A candidate must clear the floor TWICE: the scaffold's own span, and\n")
        out.write("#   each member naming a chromosome-scale reference chromosome.\n")
        out.write("assembly\tcomposite\tspan_bp\tn_members\tmembers_chromosome_scale\t"
                  "vote\tflags\n")
        n_cand = 0
        for asm, name, span, members, votes, flags in sorted(comp_rows):
            # TWO conditions, and the first is the one that matters. Testing only that the
            # MEMBERS name chromosome-scale reference chromosomes let a 2.6 Mb scaffold
            # fusing fragments of chr1 and chr5 through, because chr1 is 92 Mb. The
            # SCAFFOLD must itself be chromosome-scale: a 2.6 Mb chimeric fragment is not a
            # mis-joined chromosome, and breaking it achieves nothing because both halves
            # stay unplaced.
            if span < floor_bp:
                continue
            scale = [m for m in members if chrom_of(m) in ref_len
                     and ref_len[chrom_of(m)] >= floor_bp]
            if len(scale) < 2:
                continue
            n_cand += 1
            out.write("%s\t%s\t%d\t%d\t%s\t%s\t%s\n"
                      % (asm, name, span, len(members), ",".join(scale),
                         ";".join("%s:%s(%df/%ds)" % (k, v[0], v[1], v[2])
                                  for k, v in sorted(votes.items())) or ".", flags))

    # ---- audit ------------------------------------------------------------------------
    with open(op(".input_coverage_audit.tsv"), "w") as out:
        out.write("metric\tvalue\n")
        out.write("label\t%s\n" % a.label)
        out.write("haplotypes\t%d\nchromosomes\t%d\n" % (len(haps), len(by_chrom)))
        out.write("chromosome_set_method\t%s\n" % meta["method"])
        out.write("chromosome_set_n\t%d\n" % meta["n_chrom"])
        out.write("chromosome_scale_floor_bp\t%d\n" % floor_bp)
        out.write("chromosome_set_genome_fraction\t%.4f\n" % meta["genome_fraction"])
        if meta["cut_ratio"] is not None:
            out.write("chromosome_set_cut_ratio\t%.3f\n" % meta["cut_ratio"])
        for f in meta["flags"]:
            out.write("chromosome_set_flag\t%s\n" % f)
        out.write("low_frac_threshold\t%.3f\n" % a.low_frac)
        out.write("cells_low\t%d\n" % n_low)
        out.write("cells_low_explained_by_composite\t%d\n" % n_explained)
        out.write("cells_low_unexplained\t%d\n" % (n_low - n_explained))
        out.write("composites_total\t%d\n" % len(comp_rows))
        out.write("break_candidates\t%d\n" % n_cand)
        for h in haps:
            out.write("hap_bp_in_graph.%s\t%d\n" % (h, by_hap[h]))
        out.write("# cells_low_unexplained is the number worth reading: a chromosome far\n")
        out.write("# below its cohort median with NO composite to explain it is a different\n")
        out.write("# problem -- a genuinely absent chromosome, or a chromosome-assignment\n")
        out.write("# failure that left no composite name behind.\n")

    sys.stderr.write("[input_coverage] %d low cells, %d explained by a composite, "
                     "%d unexplained; %d break candidates of %d composites\n"
                     % (n_low, n_explained, n_low - n_explained, n_cand, len(comp_rows)))


if __name__ == "__main__":
    main()
