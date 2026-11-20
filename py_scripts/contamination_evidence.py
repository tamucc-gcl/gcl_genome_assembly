#!/usr/bin/env python3
"""Translate search coordinates and assemble reviewable contamination evidence."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
from shlex import quote
import shutil
from types import SimpleNamespace

from contamination_review import FIELDS, digest, records, write_record, write_table


HIT_FIELDS = ["query", "sequence", "start0", "end0", "strand", "subject", "subject_title", "taxids",
              "evalue", "bitscore", "identity_pct", "lineage_at_host_rank", "host_relation"]
WINDOW_FIELDS = ["query", "sequence", "start0", "end0", "protein_evidence", "megablast_evidence",
                 "nucleotide_target"]
TARGET_FIELDS = ["query", "sequence", "start0", "end0", "hit_start0", "hit_end0", "evidence"]
CONTIG_REVIEW_FIELDS = ["non_host_hit_bp", "distinct_non_host_bp", "distinct_non_host_pct",
                        "competing_hit_bp", "host_aligned_bp", "host_aligned_pct",
                        "dominant_non_host_group", "dominant_non_host_bp", "whole_contig_candidate"]
SEARCH_LABELS = {"protein": "BLASTX", "megablast": "MEGABLAST", "nucleotide": "BLASTN"}
LOCUS_REVIEW_FIELDS = ["non_host_match", "host_match", "alternative_non_host_taxa",
                       "search_agreement", "blastn_followup", "review_reason"]
PLOT_FIELDS = ["action_id", "sequence", "interval_start0", "interval_end0", "contig_bp",
               "method", "role", "group", "subject", "start0", "end0", "identity_pct", "bitscore",
               "blastn_start0", "blastn_end0"]


def table(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def protein_identity(path):
    identity = json.loads(path.read_text())
    resolution = identity.get('taxonomy_resolution')
    if resolution:
        if resolution['audit_file'] != 'taxonomy_resolution.tsv.gz':
            raise ValueError('Unexpected protein taxonomy audit filename')
        audit = path.resolve().parent / resolution['audit_file']
        if digest(audit) != resolution['audit_sha256']:
            raise ValueError('Protein taxonomy audit does not match its manifest')
    return identity


def protein_taxonomy_note(identity):
    resolution = identity.get('taxonomy_resolution')
    if not resolution:
        return []
    return [f"Protein reference taxonomy: {resolution['mapped_protein_records']:,} proteins mapped; "
            f"**{resolution['unresolved_protein_records']:,} proteins retained with unresolved taxonomy**. "
            "Unresolved hits remain searchable and are recorded as unresolved evidence. "
            "[Database taxonomy-resolution audit (gzipped TSV)](protein_taxonomy_resolution.tsv.gz).", ""]


def taxonomy(taxdump, rank, host, ids):
    with (taxdump / "nodes.dmp").open() as handle:
        nodes = {}
        for line in handle:
            tid, parent, level, *_ = line.split("|")
            nodes[int(tid)] = (int(parent), level.strip() == rank)
    merged = {}
    if (taxdump / "merged.dmp").is_file():
        with (taxdump / "merged.dmp").open() as handle:
            for line in handle:
                old, new, *_ = line.split("|")
                merged[int(old)] = int(new)

    def ancestor(tid):
        tid = merged.get(tid, tid)
        seen = set()
        while tid in nodes and tid not in seen:
            parent, match = nodes[tid]
            if match:
                return tid
            seen.add(tid)
            tid = parent
        return None

    host_group = ancestor(host)
    if host_group is None:
        raise ValueError(f"Host taxid {host} has no lineage at rank {rank}")
    groups = {tid: ancestor(tid) for tid in ids}
    needed = {value for value in groups.values() if value is not None}
    names = {}
    with (taxdump / "names.dmp").open() as handle:
        for line in handle:
            tid, name, _, kind, *_ = line.split("|")
            if kind.strip() == "scientific name" and int(tid) in needed:
                names[int(tid)] = name.strip()
    return host_group, groups, names


def blob_taxonomy(args):
    """Write the lineage files BlobToolKit needs for the host and observed hits."""
    with (args.taxdump / "nodes.dmp").open() as handle:
        parents = {}
        for line in handle:
            tid, parent, *_ = line.split("|")
            parents[int(tid)] = int(parent)
    aliases = {}
    if (args.taxdump / "merged.dmp").is_file():
        with (args.taxdump / "merged.dmp").open() as handle:
            for line in handle:
                old, new, *_ = line.split("|")
                aliases[int(old)] = int(new)

    args.out.mkdir()
    needed = {aliases.get(args.taxid, args.taxid)}
    with args.hits.open() as hits, (args.out / "hits.tsv").open("w") as output:
        for line in hits:
            columns = line.rstrip("\n").split("\t")
            ids = {aliases.get(int(tid), int(tid)) for tid in columns[1].split(";") if tid.isdigit()}
            ids = {tid for tid in ids if tid > 1 and tid in parents}
            if ids:
                needed.update(ids)
                columns[1] = ";".join(str(tid) for tid in sorted(ids))
                output.write("\t".join(columns) + "\n")
    lineages = {}
    for tid in needed:
        lineage, ancestor = [], tid
        while parents[ancestor] != ancestor:
            ancestor = parents[ancestor]
            if ancestor in lineage:
                raise ValueError(f"Taxonomy parent cycle for {tid}")
            lineage.append(ancestor)
        lineages[tid] = lineage[::-1]
    needed.update(ancestor for lineage in lineages.values() for ancestor in lineage)
    with (args.out / "taxidlineage.dmp").open("w") as output:
        for tid, lineage in sorted(lineages.items()):
            if lineage:
                output.write(f"{tid}\t|\t{' '.join(map(str, lineage))}\t|\n")
    for filename in ("nodes.dmp", "names.dmp"):
        with (args.taxdump / filename).open() as source, (args.out / filename).open("w") as output:
            for line in source:
                if int(line.split("|", 1)[0]) in needed:
                    output.write(line)


def search_hits(paths, windows, args):
    hits = []
    taxa = set()
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                query, subject, taxids, qstart, qend, score, evalue, identity, title = line.rstrip("\n").split("\t", 8)
                window = windows[query]
                start, end = sorted((int(qstart), int(qend)))
                if start < 1 or end > int(window["end0"]) - int(window["start0"]):
                    raise ValueError(f"Search hit outside its query window: {query}")
                ids = {int(value) for value in taxids.split(";") if value.isdigit() and int(value) > 0}
                taxa.update(ids)
                start0 = int(window["start0"])
                hits.append((ids, dict(zip(HIT_FIELDS, [
                    query, window["sequence"], start0 + start - 1, start0 + end,
                    "+" if int(qstart) <= int(qend) else "-", subject, title, taxids,
                    evalue, float(score), identity, "", "",
                ]))))
    host, groups, names = taxonomy(args.taxdump, args.host_rank, args.taxid, taxa)
    for ids, hit in hits:
        lineages = {groups[tid] for tid in ids}
        hit["lineage_at_host_rank"] = "; ".join(
            sorted(names.get(tid, str(tid)) for tid in lineages if tid is not None))
        if None in lineages or not lineages:
            relation = "unresolved"
        elif lineages == {host}:
            relation = "host-compatible"
        elif host in lineages:
            relation = "mixed"
        else:
            relation = "non-host"
        hit["host_relation"] = relation
    return sorted((hit for _, hit in hits),
                  key=lambda hit: (hit["query"], -hit["bitscore"], hit["start0"], hit["subject"]))


def informative_hits(hits):
    """Drop clearly weaker alternatives over the same bases, preserving separate loci."""
    keep = []
    for hit in hits:
        span = hit["end0"] - hit["start0"]
        if not any(other["bitscore"] > hit["bitscore"] / .9 and
                   min(hit["end0"], other["end0"]) - max(hit["start0"], other["start0"]) >= .8 * span
                   for other in keep):
            keep.append(hit)
    return keep


def relationship(hits):
    states = {hit["host_relation"] for hit in hits}
    if "mixed" in states or {"host-compatible", "non-host"} <= states:
        return "mixed"
    if "unresolved" in states and "host-compatible" in states:
        return "unresolved"
    if "non-host" in states:
        return "non-host"
    if "host-compatible" in states:
        return "host-compatible"
    return "unresolved" if hits else "no-hit"


def group_by_query(hits):
    result = defaultdict(list)
    for hit in hits:
        result[hit["query"]].append(hit)
    return result


def grouped_hits(hits):
    return {query: informative_hits(values) for query, values in group_by_query(hits).items()}


def overlap(left, right):
    return min(left["end0"], right["end0"]) - max(left["start0"], right["start0"])


def comparable_hit(reference, other):
    return (other["bitscore"] >= .9 * reference["bitscore"]
            and overlap(reference, other) >= .8 * (reference["end0"] - reference["start0"]))


def candidate_hits(hits):
    for hit in hits:
        if hit["host_relation"] not in {"non-host", "mixed"}:
            continue
        competing = hit["host_relation"] == "mixed" or any(
            other["host_relation"] == "host-compatible"
            and comparable_hit(hit, other)
            for other in hits
        )
        yield hit, competing


def span_union(spans):
    total = end = 0
    for lo, hi in sorted(spans):
        total += max(0, hi - max(lo, end))
        end = max(end, hi)
    return total


def ranked_loci(protein, mega, nucleotide=()):
    """Group overlapping evidence and rank loci without treating it as a removal decision."""
    evidence = sorted(
        ((hit, method, competing) for method, hits in
         (("protein", protein), ("megablast", mega), ("nucleotide", nucleotide))
         for hit, competing in candidate_hits(hits)),
        key=lambda item: (item[0]["start0"], item[0]["end0"], item[1]),
    )
    groups = []
    for item in evidence:
        hit = item[0]
        if groups and hit["start0"] <= groups[-1]["end0"]:
            groups[-1]["end0"] = max(groups[-1]["end0"], hit["end0"])
            groups[-1]["hits"].append(item)
        else:
            groups.append(dict(start0=hit["start0"], end0=hit["end0"], hits=[item]))
    for locus in groups:
        hits = locus["hits"]
        distinct = [(hit, method) for hit, method, competing in hits if not competing]
        locus["distinct_spans"] = [(hit["start0"], hit["end0"]) for hit, _ in distinct]
        locus["competing_spans"] = [(hit["start0"], hit["end0"]) for hit, _, competing in hits if competing]
        locus["distinct_bp"] = span_union(locus["distinct_spans"])
        agreement = {
            frozenset((method, other_method))
            for i, (hit, method) in enumerate(distinct)
            for other, other_method in distinct[i + 1:]
            if other_method != method and hit["lineage_at_host_rank"]
            and hit["lineage_at_host_rank"] == other["lineage_at_host_rank"]
            and ";" not in hit["lineage_at_host_rank"] and overlap(hit, other) >= 250
        }
        locus["concordant"] = any("protein" in pair for pair in agreement)
        locus["nucleotide_agreement"] = frozenset(("megablast", "nucleotide")) in agreement
        locus["methods"] = {method for _, method, _ in hits}
        locus["lineages"] = {hit["lineage_at_host_rank"] for hit, _, _ in hits} - {""}
        locus["context"] = (
            ({"distinct non-host hit"} if distinct else set())
            | ({"same-locus host competitor"} if locus["competing_spans"] else set())
        )
        if locus["concordant"]:
            locus["priority"] = "concordant non-host"
        elif locus["distinct_bp"] >= 1000:
            locus["priority"] = "extended non-host"
        elif distinct:
            locus["priority"] = "short non-host signal"
        else:
            locus["priority"] = "competing hits only"
    return groups


def locus_comparison(locus, hosts):
    non_host_rows, host_rows, alternative_rows = [], [], []
    plot_hits = []
    any_host = False
    for method, label in SEARCH_LABELS.items():
        candidates = [hit for hit, source, _ in locus["hits"] if source == method]
        best = max(candidates, key=lambda hit: (hit["bitscore"], overlap(locus, hit)), default=None)
        overlapping = [hit for hit in hosts.get(method, ()) if overlap(best or locus, hit) > 0]
        host = max(overlapping, key=lambda hit: (
            bool(best and comparable_hit(best, hit)), hit["bitscore"], overlap(best or locus, hit)
        ), default=None)
        if best:
            plot_hits.append((label, "Mixed taxids" if best["host_relation"] == "mixed" else "Non-host", best))
            non_host_rows.append(f"{label} {describe_hit(best, overlap(locus, best))}")
            alternatives = {}
            for hit in sorted(candidates, key=lambda hit: -hit["bitscore"]):
                group = hit["lineage_at_host_rank"]
                if group and group != best["lineage_at_host_rank"] and comparable_hit(best, hit):
                    alternatives.setdefault(group, hit)
            if alternatives:
                alternative_rows.append(f"{label} " + ", ".join(sorted(alternatives)))
                plot_hits.extend((label, "Alternative taxon", hit) for hit in list(alternatives.values())[:2])
        if host:
            plot_hits.append((label, "Host-compatible", host))
            any_host = True
            span = "bp overlapping non-host hit" if best else "bp in interval"
            qualifier = " (comparable score and span)" if best and comparable_hit(best, host) else ""
            host_rows.append(f"{label} {describe_hit(host, overlap(best or locus, host), span)}{qualifier}")
        elif best:
            host_rows.append(f"{label} no overlapping host hit found")

    if not locus["distinct_bp"]:
        reason = "Lower priority — host evidence competes at the same locus"
    elif locus["competing_spans"]:
        reason = "Compare alignments — competing host evidence at this locus"
    elif any_host:
        reason = "Compare alignments — an overlapping host hit was also found"
    elif locus["priority"] == "short non-host signal":
        reason = "Lower priority — short non-host alignment"
    else:
        reason = "Inspect interval — no overlapping host hit found"
    agreement = []
    if locus["concordant"]:
        agreement.append("Protein and nucleotide agree")
    if locus["nucleotide_agreement"]:
        agreement.append("Both nucleotide searches agree (same nt database)")
    return dict(non_host_match=" · ".join(non_host_rows), host_match=" · ".join(host_rows),
                alternative_non_host_taxa=" · ".join(alternative_rows),
                search_agreement="; ".join(agreement) or "No matching cross-search evidence",
                review_reason=reason, plot_hits=plot_hits)


def describe_hit(hit, bp, span="bp in interval"):
    return (f"{hit['lineage_at_host_rank'] or 'Unresolved taxon'} ({hit['subject']}), "
            f"{bp:,} {span}, {float(hit['identity_pct']):.1f}% identity, {hit['bitscore']:g} bits")


def markdown_cell(value):
    return value.replace("|", "&#124;").replace(" · ", "<br>")


INTERVAL_RANKING = (
    "The table ranks concordant protein and nucleotide evidence first, then extended non-host evidence. "
    "Within each group, more distinct non-host bases rank higher. Ties favor a larger fraction of "
    "the contig covered by that interval, then a longer interval. The remaining concordant or extended "
    "intervals also qualify for review. Short hits and competing-only hits rank lower in the full sheet."
)


def interval_review_guide(host_rank):
    return [
        f"**Where host matches come from.** Each search returns subject taxids. The taxonomy lookup traces them "
        f"to the configured host rank ({host_rank}). A hit is host-compatible when its usable taxids all resolve "
        f"to the same group as the assembly taxid at that rank. A match at phylum can come from another species "
        "in that phylum. The host column shows overlapping host-compatible hits from the searches named in the row.",
        "",
        "**How to decide.**",
        "",
        "1. **Compare the same bases.** A host hit covering most of the non-host hit at a similar score weakens "
        "the case for removal. Compare bit scores within the same search. Identity is measured in amino acids "
        "for BLASTX and nucleotides for MEGABLAST and BLASTN.",
        "2. **Check whether the non-host assignment is consistent.** Look for a coherent stretch of aligned "
        "sequence and a consistent group across matches. Agreement between protein and nucleotide evidence adds "
        "support. Agreement between MEGABLAST and BLASTN tests sensitivity against the same database. "
        "Several similarly scored taxa or sparse coverage favor leaving the interval unselected for further review.",
        "3. **Choose the scope of any removal.** A strong local signal calls for review of that interval. Use "
        "the whole-contig section for evidence spanning most of a contig. "
        "[Confirm the coordinates and record a decision](#confirm-coordinates-and-record-a-decision).",
        "",
        "The table shows a representative non-host hit and the best overlapping host hit per search. "
        "An absent host hit means that search returned no overlapping host-compatible alignment.",
        "",
    ]


REVIEW_ACTION_GUIDE = [
    "## Confirm coordinates and record a decision", "",
    "1. Open the [review sheet (TSV)](review.tsv) and find the **Review ID** from this report. "
    "Check the contig name and the `start0` and `end0` columns against the evidence.",
    "2. Use the original **post-FCS FASTA** for a post-FCS report. Coordinates are 0-based and half-open. "
    "The interval `[100, 200)` corresponds to bases 101 through 200 in a 1-based viewer. "
    "Confirm the interval on the source contig and edit the coordinates if the evidence supports a different boundary.",
    "3. For a local removal use `REMOVE_INTERVAL`. For an entire contig use `REMOVE_SEQUENCE` with "
    "`start0=0` and `end0` equal to the contig length. Set `selected` to `YES` only on the chosen rows "
    "and fill in `reviewer` and `reason`. Save a reviewed copy as TSV. Leave unresolved intervals at `NO`.",
    "",
    "To add a manual removal, use a new action ID and copy the assembly, stage and "
    "source FASTA checksum from the matching source assembly. Selected intervals on the same contig must be disjoint.",
    "",
]


def add_read_support_guide(review_dir):
    readme = review_dir / "README.md"
    text = readme.read_text(encoding="utf-8")
    readme.write_text(text.replace(
        "### Intervals to review first\n\n",
        "### Intervals to review first\n\n"
        "Locus plots include [local read depth](locus_read_depth.tsv) and "
        "[reads spanning each interval edge](locus_read_boundaries.tsv). "
        "Counts use primary alignments with mapping quality at least 20, with at least 30 aligned bases "
        "in each 50 bp flank. Edges within 50 bp of the contig end have insufficient flanking sequence "
        "for a spanning test. Read support measures sequence continuity. Taxonomic assignments come from the search tracks.\n\n",
        1), encoding="utf-8")


def review_apply_guide(source):
    return [
        "### Apply the reviewed decisions", "",
        "Run this from the assembly project directory with the exact FASTA used to create the report. "
        "Replace the FASTA and reviewed TSV paths and choose a new output directory.", "",
        "```bash",
        "python3 gcl_genome_assembly/py_scripts/contamination_review.py apply \\",
        '  --fasta "/path/to/original.fasta" \\',
        f"  --assembly {quote(source['assembly'])} --stage {quote(source['stage'])} \\",
        '  --review "/path/to/reviewed.tsv" \\',
        '  --out "/path/to/cleaned-output"',
        "```", "",
        "For a post-FCS report, add `--contamination_review /path/to/reviewed.tsv` and `-resume` "
        "to the original assembly command. Set `--stop_after_contig_review false` to continue into scaffolding. "
        "Keep the original assembly inputs and FCS settings so the reviewed coordinates still match.", "",
        "The standalone command above checks the assembly, stage, and FASTA checksum before writing "
        "`retained.fasta` and the removal audit.",
        "",
    ]


def whole_contig_evidence(spans_by_group, distinct_bp, host_bp, length):
    group, bp = max(((group, span_union(spans)) for group, spans in spans_by_group.items()),
                    key=lambda item: (item[1], item[0]), default=("", 0))
    candidate = bool(length and 2 * bp >= length and 5 * bp >= 4 * distinct_bp
                     and 10 * host_bp <= length)
    return dict(dominant_non_host_group=group, dominant_non_host_bp=bp,
                whole_contig_candidate="YES" if candidate else "NO")


def blastn_result(locus, target, hits, searched):
    if not searched:
        return "", False
    if not target or overlap(locus, target) <= 0:
        return "BLASTN: outside follow-up", False
    covered = overlap(locus, target)
    if any(method == "nucleotide" for _, method, _ in locus["hits"]):
        status = "non-host or mixed hit listed"
    elif any(overlap(locus, hit) > 0 and hit["host_relation"] == "host-compatible" for hit in hits):
        status = "overlapping host hit found"
    elif any(overlap(locus, hit) > 0 for hit in hits):
        status = "overlapping hits did not meet review criteria"
    else:
        status = "no overlapping hit found"
    nominated = dict(start0=target["hit_start0"], end0=target["hit_end0"])
    return f"BLASTN: {status} ({covered:,} bp of interval searched)", overlap(locus, nominated) > 0


PRIORITY = {name: i for i, name in enumerate(
    ("concordant non-host", "extended non-host", "short non-host signal", "competing hits only"))}


def followup_plan(windows, protein, mega, max_bp):
    eligible = []
    for query, window in windows.items():
        for locus in ranked_loci(protein.get(query, []), mega.get(query, [])):
            if locus["priority"] in {"concordant non-host", "extended non-host"}:
                eligible.append((query, window, locus))
    eligible.sort(key=lambda item: (PRIORITY[item[2]["priority"]], -item[2]["distinct_bp"],
                                    item[0], item[2]["start0"]))
    targets = {}
    chunk, chunk_bp = 1, 0
    for query, window, locus in eligible:
        if query in targets:
            continue
        # Retain the original window ID but record the exact cropped query coordinates.
        lo = max(int(window["start0"]), min(start for start, _ in locus["distinct_spans"]) - 500)
        hi = min(int(window["end0"]), max(end for _, end in locus["distinct_spans"]) + 500)
        if hi - lo > max_bp:
            continue
        max_bp -= hi - lo
        if chunk_bp and chunk_bp + hi - lo > 25_000:
            chunk += 1
            chunk_bp = 0
        chunk_bp += hi - lo
        targets[query] = dict(query=query, sequence=window["sequence"], start0=lo, end0=hi,
                              hit_start0=min(start for start, _ in locus["distinct_spans"]),
                              hit_end0=max(end for _, end in locus["distinct_spans"]),
                              evidence=locus["priority"], chunk=f"{chunk:04d}.fasta")
    return targets


def select(args):
    windows = {row["query"]: row for row in table(args.prepared / "windows.tsv")}
    hits = search_hits(args.hits, windows, args)
    grouped = grouped_hits(hits)
    mega_hits = search_hits(args.megablast_hits, windows, args) if args.megablast_hits else []
    mega_grouped = grouped_hits(mega_hits)
    plan = followup_plan(windows, grouped, mega_grouped, args.followup_max_bp)
    targets = defaultdict(list)
    for row in plan.values():
        targets[row["sequence"]].append(row)
    assignments = []
    for query, window in windows.items():
        status = relationship(grouped.get(query, [])) if window["status"] == "search" else "insufficient-acgt"
        mega = "not-run"
        if args.megablast:
            mega = (relationship(mega_grouped.get(query, [])) if window["status"] == "search"
                    else "insufficient-acgt")
        assignments.append(dict(window, protein_evidence=status, megablast_evidence=mega,
                                nucleotide_target="YES" if query in plan else "NO"))
    args.out.mkdir()
    write_table(args.out / "protein_hits.tsv", hits, HIT_FIELDS)
    write_table(args.out / "megablast_hits.tsv", mega_hits, HIT_FIELDS)
    write_table(args.out / "windows.tsv",
                ({field: row[field] for field in WINDOW_FIELDS} for row in assignments), WINDOW_FIELDS)
    write_table(args.out / "nucleotide_targets.tsv",
                ({field: row[field] for field in TARGET_FIELDS} for row in plan.values()), TARGET_FIELDS)
    with (args.out / "targets.fasta").open("w") as output:
        for header, sequence in records(args.fasta):
            for row in targets[header.split()[0]]:
                fragment = sequence[int(row["start0"]):int(row["end0"])]
                write_record(output, row["query"], fragment)
                with (args.out / f"target_{row['chunk']}").open("a") as chunk:
                    write_record(chunk, row["query"], fragment)


def retained_hits(path):
    hits = table(path)
    for hit in hits:
        for key, convert in (("bitscore", float), ("start0", int), ("end0", int)):
            hit[key] = convert(hit[key])
    return hits


def lineage_group(hits, no_hits):
    by_lineage = defaultdict(list)
    for hit in hits:
        lineage = hit["lineage_at_host_rank"]
        if lineage and ";" not in lineage:
            by_lineage[lineage].append((hit["start0"], hit["end0"]))
    weights = {}
    for lineage, spans in by_lineage.items():
        covered, end = 0, 0
        for lo, hi in sorted(spans):
            covered += max(0, hi - max(lo, end))
            end = max(hi, end)
        weights[lineage] = covered
    if not weights:
        return "no-hit" if no_hits else "unresolved"
    best = max(weights, key=weights.get)
    return best if weights[best] >= .8 * sum(weights.values()) else "conflicting"


def taxon_summary(metrics):
    groups = defaultdict(list)
    for row in metrics:
        groups[row["plot_taxon"]].append(int(row["length_bp"]))
    for group, sizes in sorted(groups.items(), key=lambda item: (-sum(item[1]), item[0])):
        total, cumulative = sum(sizes), 0
        for length in sorted(sizes, reverse=True):
            cumulative += length
            if 2 * cumulative >= total:
                yield dict(group=group, contigs=len(sizes), span_bp=total, n50_bp=length)
                break


def summarize(args, translated_nucleotide_hits=None, previous=None):
    assigned = table(args.selected / "windows.tsv")
    targets = {row["query"]: row for row in table(args.selected / "nucleotide_targets.tsv")}
    for target in targets.values():
        for field in ("start0", "end0", "hit_start0", "hit_end0"):
            target[field] = int(target[field])
    if translated_nucleotide_hits is None:
        nucleotide = search_hits(args.nucleotide_hits, targets, args) if args.nucleotide_hits else []
    else:
        nucleotide = translated_nucleotide_hits
        for hit in nucleotide:
            target = targets[hit["query"]]
            if (hit["sequence"] != target["sequence"] or
                    not target["start0"] <= hit["start0"] < hit["end0"] <= target["end0"]):
                raise ValueError(f"Translated BLASTN hit outside its planned region: {hit['query']}")
    raw_nt_by_query = group_by_query(nucleotide)
    nt_by_query = grouped_hits(nucleotide)
    protein = retained_hits(args.selected / "protein_hits.tsv")
    mega = retained_hits(args.selected / "megablast_hits.tsv")
    protein_by_query = grouped_hits(protein)
    mega_by_query = grouped_hits(mega)
    raw_hosts = {method: group_by_query(hit for hit in hits if hit["host_relation"] == "host-compatible")
                 for method, hits in (("protein", protein), ("megablast", mega), ("nucleotide", nucleotide))}
    plot_hits = defaultdict(list)
    mega_sequences = {hit["sequence"] for hits in mega_by_query.values() for hit in hits
                      if hit["lineage_at_host_rank"]}
    by_sequence = defaultdict(set)
    loci_by_sequence = defaultdict(list)
    host_spans = defaultdict(list)
    non_host_group_spans = defaultdict(lambda: defaultdict(list))
    for row in assigned:
        query = row["query"]
        nt = relationship(nt_by_query.get(query, []))
        if row["nucleotide_target"] != "YES":
            nt = "not-targeted"
        elif not args.targeted:
            nt = "not-run"
        row["nucleotide_evidence"] = nt
        name = row["sequence"]
        by_sequence[name].update([row["protein_evidence"], row["megablast_evidence"], nt])
        plotting = mega_by_query if name in mega_sequences else protein_by_query
        plot_hits[name].extend(plotting.get(query, []))
        protein_hits = protein_by_query.get(query, [])
        mega_hits = mega_by_query.get(query, [])
        nt_hits = nt_by_query.get(query, [])
        host_spans[name].extend((hit["start0"], hit["end0"]) for hits in (protein_hits, mega_hits, nt_hits)
                                for hit in hits if hit["host_relation"] == "host-compatible")
        hosts = {method: hits.get(query, []) for method, hits in raw_hosts.items()}
        for locus in ranked_loci(protein_hits, mega_hits, nt_hits):
            locus.update(locus_comparison(locus, hosts))
            locus["blastn_followup"], locus["planned_hit"] = blastn_result(
                locus, targets.get(query), raw_nt_by_query.get(query, ()), args.targeted)
            locus["blastn_target"] = targets.get(query) if locus["blastn_followup"] and query in targets else None
            for hit, _, competing in locus["hits"]:
                group = hit["lineage_at_host_rank"]
                if not competing and group and ";" not in group:
                    non_host_group_spans[name][group].append((hit["start0"], hit["end0"]))
            del locus["hits"]
            loci_by_sequence[name].append(locus)

    if previous:
        source = {key: previous[key] for key in (
            "schema", "assembly", "stage", "source_sha256", "coordinates", "window_bp", "chunk_bp",
            "minimum_acgt_bp", "source_bp", "windows", "search_queries", "search_query_bp", "chunks",
            "boundary_policy")}
        metrics = table(args.selected / "sequence_evidence.tsv")
        review = [{key: row[key] for key in FIELDS} for row in table(args.selected / "review.tsv")
                  if row["action"] == "REMOVE_SEQUENCE"]
        if any(row["selected"] != "NO" for row in review):
            raise ValueError("Use a generated report directory to rebuild the report")
    else:
        source = json.loads((args.prepared / "source.json").read_text())
        metrics = table(args.metrics / "sequence_metrics.tsv")
        review = table(args.metrics / "review.tsv")
    depths = {}
    if args.coverage:
        with args.coverage.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                name = row["#rname"]
                if name in depths or int(row["startpos"]) != 1:
                    raise ValueError(f"Invalid coverage row: {name}")
                length, covered = int(row["endpos"]), row.get("covbases")
                if covered is not None and not 0 <= int(covered) <= length:
                    raise ValueError(f"Invalid covered bases for {name}: {covered}")
                depths[name] = (length, row["meandepth"], covered)
        if {name: values[0] for name, values in depths.items()} != {
            row["sequence"]: int(row["length_bp"]) for row in metrics
        }:
            raise ValueError("Coverage BAM sequence names and lengths differ from the evidence FASTA")

    whole_sequence = {row["sequence"]: row for row in review}
    planned_hit_ids = set()
    plot_rows = {}
    for row in metrics:
        name = row["sequence"]
        states = by_sequence[name]
        loci = loci_by_sequence[name]
        length = int(row["length_bp"])
        if loci:
            status = min((locus["priority"] for locus in loci), key=PRIORITY.get)
        elif "host-compatible" in states:
            status = "host-compatible"
        elif "unresolved" not in states and "no-hit" in states:
            status = "no-hit"
        else:
            status = "unresolved"
        if previous:
            depth, covered = row["mean_depth"], row["read_covered_bp"]
        else:
            _, depth, covered = depths.get(name, (0, "NA", None))
        distinct_bp = sum(locus["distinct_bp"] for locus in loci)
        host_bp = span_union(host_spans[name])
        row.update(mean_depth=depth, read_covered_bp=covered if covered is not None else "NA", evidence=status,
                   plot_taxon=lineage_group(plot_hits[name], "no-hit" in states and not plot_hits[name]),
                   non_host_hit_bp=sum(locus["end0"] - locus["start0"] for locus in loci),
                   distinct_non_host_bp=distinct_bp,
                   distinct_non_host_pct=100 * distinct_bp / length if length else 0,
                   competing_hit_bp=sum(span_union(locus["competing_spans"]) for locus in loci),
                   host_aligned_bp=host_bp, host_aligned_pct=100 * host_bp / length if length else 0)
        row.update(whole_contig_evidence(non_host_group_spans[name], distinct_bp, host_bp, length))
        for locus in loci:
            lo, hi = locus["start0"], locus["end0"]
            if lo == 0 and hi == length:
                action = whole_sequence[name]
            else:
                action = dict(zip(FIELDS, [
                    f"A{len(review)+1:05d}", source["assembly"], source["stage"], source["source_sha256"],
                    name, "REMOVE_INTERVAL", lo, hi, "NO", "", "",
                ]))
                review.append(action)
            action.update(support=", ".join(sorted(locus["methods"])),
                          hit_lineages="; ".join(sorted(locus["lineages"])),
                          hit_context="; ".join(sorted(locus["context"])),
                          priority=locus["priority"], distinct_hit_bp=locus["distinct_bp"],
                          distinct_contig_pct=100 * locus["distinct_bp"] / length)
            action.update({key: locus[key] for key in LOCUS_REVIEW_FIELDS})
            target = locus["blastn_target"] or {}
            plot_rows[action["action_id"]] = [
                dict(action_id=action["action_id"], sequence=name, interval_start0=lo, interval_end0=hi,
                     contig_bp=length, method=method, role=role, group=hit["lineage_at_host_rank"],
                     blastn_start0=target.get("start0", ""), blastn_end0=target.get("end0", ""),
                     **{key: hit[key] for key in ("subject", "start0", "end0", "identity_pct", "bitscore")})
                for method, role, hit in locus["plot_hits"]
            ]
            if locus["planned_hit"]:
                planned_hit_ids.add(action["action_id"])

    args.out.mkdir()
    write_table(args.out / "windows.tsv", assigned, WINDOW_FIELDS + ["nucleotide_evidence"])
    write_table(args.out / "protein_hits.tsv", protein, HIT_FIELDS)
    write_table(args.out / "megablast_hits.tsv", mega, HIT_FIELDS)
    write_table(args.out / "nucleotide_hits.tsv", nucleotide, HIT_FIELDS)
    write_table(args.out / "sequence_evidence.tsv", metrics, list(metrics[0]))
    summaries = list(taxon_summary(metrics))
    write_table(args.out / "taxon_summary.tsv", summaries, ["group", "contigs", "span_bp", "n50_bp"])
    write_table(args.out / "blob_sequence_evidence.tsv", (row for row in metrics if int(row["acgt_bp"])),
                list(metrics[0]))
    write_table(args.out / "blobdir_exclusions.tsv", (
        dict(sequence=row["sequence"], length_bp=row["length_bp"], reason="no-ACGT")
        for row in metrics if not int(row["acgt_bp"])
    ), ["sequence", "length_bp", "reason"])
    contigs = {row["sequence"]: row for row in metrics}
    for row in review:
        contig = contigs[row["sequence"]]
        row["sequence_evidence"] = contig["evidence"]
        row.update({key: contig[key] for key in CONTIG_REVIEW_FIELDS})
    supported = sorted((row for row in review if row.get("support")), key=lambda row: (
        PRIORITY[row["priority"]], -row["distinct_hit_bp"], -row["distinct_contig_pct"],
        -(int(row["end0"]) - int(row["start0"])), row["sequence"], int(row["start0"])))
    write_table(args.out / "review.tsv", supported + [row for row in review if not row.get("support")],
                FIELDS + ["priority", "distinct_hit_bp", "distinct_contig_pct", "sequence_evidence",
                          *CONTIG_REVIEW_FIELDS, *LOCUS_REVIEW_FIELDS, "support", "hit_lineages", "hit_context"])
    target_fasta = "nucleotide_targets.fasta" if previous else "targets.fasta"
    shutil.copyfile(args.selected / target_fasta, args.out / "nucleotide_targets.fasta")
    shutil.copyfile(args.selected / "nucleotide_targets.tsv", args.out / "nucleotide_targets.tsv")
    if args.nucleotide_manifest:
        shutil.copyfile(args.nucleotide_manifest, args.out / "nucleotide_database.json")
    with (args.out / "blob_hits.tsv").open("w") as handle:
        for hit in (hit for hits in plot_hits.values() for hit in hits):
            if any(value.isdigit() and int(value) > 0 for value in hit["taxids"].split(";")):
                handle.write("\t".join(str(hit[key]) for key in
                    ["sequence", "taxids", "bitscore", "subject", "evalue"]) +
                    f"\t{hit['start0']+1}\t{hit['end0']}\n")
    counts = Counter(row["protein_evidence"] for row in assigned)
    priorities = Counter(row["priority"] for row in supported)
    shortlist = [row for row in supported if PRIORITY[row["priority"]] <= PRIORITY["extended non-host"]]
    interval_preview = shortlist[:30]
    shown_ids = {row["action_id"] for row in interval_preview}
    interval_preview.extend(row for row in supported
                            if row["action_id"] in planned_hit_ids and row["action_id"] not in shown_ids)
    write_table(args.out / "locus_plot_hits.tsv",
                (hit for row in interval_preview for hit in plot_rows[row["action_id"]]), PLOT_FIELDS)
    contig_shortlist = sorted((row for row in metrics if row["whole_contig_candidate"] == "YES"),
                              key=lambda row: (-row["dominant_non_host_bp"] / int(row["length_bp"]),
                                               row["sequence"]))
    target_bp = sum(int(row["end0"]) - int(row["start0"]) for row in targets.values())
    flagged_bp = sum(int(row["end0"]) - int(row["start0"]) for row in supported)
    covered_bp = (sum(int(row["read_covered_bp"]) for row in metrics)
                  if all(row["read_covered_bp"] != "NA" for row in metrics) else None)
    mega_counts = Counter(row["megablast_evidence"] for row in assigned)
    metadata = dict(source, host_taxid=args.taxid, host_rank=args.host_rank,
                    protein_windows=dict(counts), nucleotide_target_windows=len(targets),
                    nucleotide_target_bp=target_bp,
                    megablast_windows=dict(mega_counts), flagged_interval_bp=flagged_bp,
                    interval_priorities=dict(priorities),
                    whole_contig_candidates=len(contig_shortlist),
                    read_covered_bp=covered_bp,
                    nucleotide_status=("searched" if targets else "no-targets") if args.targeted else "not-run",
                    coverage=previous["coverage"] if previous else
                    (args.coverage_source if args.coverage else "unavailable"),
                    protein_database=previous["protein_database"] if previous else protein_identity(args.database),
                    nucleotide_database=args.nucleotide_db or None)
    resolution = metadata['protein_database'].get('taxonomy_resolution')
    if resolution:
        audit = (args.selected / 'protein_taxonomy_resolution.tsv.gz' if previous else
                 args.database.resolve().parent / resolution['audit_file'])
        if digest(audit) != resolution['audit_sha256']:
            raise ValueError('Protein taxonomy audit checksum mismatch')
        shutil.copyfile(audit, args.out / 'protein_taxonomy_resolution.tsv.gz')
    (args.out / "evidence.json").write_text(json.dumps(metadata, indent=2) + "\n")
    lines = [
        f"# Contamination evidence — {source['assembly']}", "",
        f"Coordinate stage: `{source['stage']}`. Coordinates are 0-based, half-open.", "",
        f"Host comparison: NCBI taxid {args.taxid}, {args.host_rank}. "
        f"Nucleotide follow-up: {metadata['nucleotide_status']} "
        f"({len(targets):,} planned regions, {target_bp:,} of {source['source_bp']:,} bp "
        f"or {100 * target_bp / source['source_bp']:.3f}%). "
        f"Read coverage: {metadata['coverage']}.", "",
        *protein_taxonomy_note(metadata['protein_database']),
        *(["[Managed nucleotide database provenance](nucleotide_database.json).", ""]
          if args.nucleotide_manifest else []),
        *([f"**Bases covered by reads: {covered_bp:,} of {source['source_bp']:,} bp "
           f"({100 * covered_bp / source['source_bp']:.1f}%).**", ""] if covered_bp is not None else []),
        "| Window evidence | Protein | MEGABLAST |", "| --- | ---: | ---: |",
        *(f"| {status} | {counts[status]:,} | {mega_counts[status]:,} |"
          for status in sorted(counts.keys() | mega_counts.keys())), "",
        "## Review priorities", "",
        f"The searches flagged {flagged_bp:,} bp across {len(supported):,} evidence intervals. "
        "Every removal option starts unselected.", "",
        "| Evidence | Intervals |", "| --- | ---: |",
        *(f"| {name} | {priorities[name]:,} |" for name in PRIORITY), "",
        "Concordant means distinct protein and nucleotide hits agree on the taxonomic group "
        "and overlap by at least 250 bp. MEGABLAST and BLASTN agreement is shown separately because "
        "both search nt. Extended means at least 1 kb of distinct hit coverage in a locus.", "",
        "The review reasons guide inspection. Host and non-host spans can overlap. "
        "The host-aligned total counts the union of spans from retained host-compatible hits. "
        "Bases outside those spans may lack informative alignments.", "",
        "### Whole-contig removal candidates", "",
        (f"{len(contig_shortlist):,} contigs met the whole-contig screen. "
         f"Showing {min(len(contig_shortlist), 20):,}.")
        if contig_shortlist else "**No contigs met the whole-contig screen.**", "",
        "A candidate needs distinct hits to one non-host group across at least half the contig. "
        "That group must cover at least 80% of the distinct non-host span, and host-compatible hits "
        "must cover no more than 10% of the contig. A local hit can still qualify for interval review below.", "",
        *(
            ["| Review ID | Contig | Length (bp) | Dominant non-host group | Group-aligned bp (% of contig) | Host-aligned bp (% of contig) | Mean depth |",
             "| --- | --- | ---: | --- | ---: | ---: | ---: |"]
            + [f"| {whole_sequence[row['sequence']]['action_id']} | {markdown_cell(row['sequence'])} | {int(row['length_bp']):,} | "
               f"{markdown_cell(row['dominant_non_host_group'])} | {row['dominant_non_host_bp']:,} "
               f"({100 * row['dominant_non_host_bp'] / int(row['length_bp']):.2f}%) | "
               f"{row['host_aligned_bp']:,} ({row['host_aligned_pct']:.2f}%) | {row['mean_depth']} |"
               for row in contig_shortlist[:20]]
            + [""]
            if contig_shortlist else []
        ),
        "### Intervals to review first", "",
        (f"Showing the top {min(len(shortlist), 30):,} of {len(shortlist):,} concordant or extended intervals.")
        if shortlist else "No intervals met the concordant or extended evidence criteria.", "",
        "[**Open the full interval review sheet (TSV)**](review.tsv) for every evidence interval "
        "and manual whole-contig removal options.", "",
        *([INTERVAL_RANKING, ""] if shortlist else []),
        *([f"Included {len(interval_preview) - min(len(shortlist), 30):,} additional intervals "
           "with targeted BLASTN results.", ""] if len(interval_preview) > min(len(shortlist), 30) else []),
        *interval_review_guide(args.host_rank),
        *(
            ["| Review ID | Contig interval (bp) | Non-host match | Host match at the same locus | Agreement across searches | Review reason |",
             "| --- | --- | --- | --- | --- | --- |"]
            + [f"| [{row['action_id']}](../diagnostics/loci/{row['action_id']}.png) | {markdown_cell(row['sequence'])} "
               f"[{int(row['start0']):,}, {int(row['end0']):,})<br>{row['distinct_hit_bp']:,} distinct non-host bp | "
               f"{markdown_cell(row['non_host_match'])}"
               + (f"<br>Similarly scored taxa: {markdown_cell(row['alternative_non_host_taxa'])}"
                  if row["alternative_non_host_taxa"] else "")
               + f" | {markdown_cell(row['host_match'])} | {row['search_agreement']}"
               + (f"<br>{row['blastn_followup']}" if row["blastn_followup"] else "")
               + f" | {row['review_reason']} |"
               for row in interval_preview]
            + [""]
            if interval_preview else []
        ),
        "See the [window assignments](windows.tsv), [protein hits](protein_hits.tsv) and "
        "[MEGABLAST hits](megablast_hits.tsv) and [targeted nucleotide hits](nucleotide_hits.tsv) "
        "for intervals and supporting matches. [All contig summaries](sequence_evidence.tsv).", "",
        "[Sensitive BLASTN plan](nucleotide_targets.tsv) · [Target FASTA](nucleotide_targets.fasta).", "",
        *REVIEW_ACTION_GUIDE,
        *review_apply_guide(source),
        "## Diagnostics", "",
        "Taxon groups use the union of aligned intervals per lineage. A lineage must account for "
        "at least 80% of the summed lineage spans to color a contig. MEGABLAST is used where it "
        "provides a ranked hit, with protein as fallback. Candidate span above counts only flagged intervals.", "",
        "| Taxon group | Contigs | Contig span (bp) | N50 (bp) |", "| --- | ---: | ---: | ---: |",
        *(f"| {row['group'].replace('|', '&#124;')} | {row['contigs']:,} | "
          f"{row['span_bp']:,} | {row['n50_bp']:,} |" for row in summaries[:12]), "",
        "![Blob diagnostic](../diagnostics/blob.png)", "",
        "![Cumulative sequence by evidence](../diagnostics/cumulative.png)", "",
        "[Taxon counts, span and N50](taxon_summary.tsv). "
        "[BlobToolKit dataset](../blobdir/) with MEGABLAST taxonomy where available, "
        "protein evidence as fallback, and the separate review categories. "
        "[Sequences without measurable GC](blobdir_exclusions.tsv) remain in the evidence tables.", "",
    ]
    (args.out / "README.md").write_text("\n".join(lines), encoding="utf-8")


def rebuild_report(args):
    previous = json.loads((args.previous / "evidence.json").read_text())
    if digest(args.fasta) != previous["source_sha256"]:
        raise ValueError("Original FASTA checksum differs from the previous evidence report")
    nt_manifest = args.previous / "nucleotide_database.json"
    summarize(SimpleNamespace(
        selected=args.previous, out=args.out, coverage=None,
        targeted=previous["nucleotide_status"] in {"searched", "no-targets"},
        nucleotide_manifest=nt_manifest if nt_manifest.is_file() else None,
        nucleotide_db=previous["nucleotide_database"], taxid=previous["host_taxid"],
        host_rank=previous["host_rank"]),
        translated_nucleotide_hits=retained_hits(args.previous / "nucleotide_hits.tsv"),
        previous=previous)
    read_tables = ["locus_read_depth.tsv", "locus_read_boundaries.tsv"]
    if all((args.previous / name).is_file() for name in read_tables):
        # Carry forward read evidence only when the review IDs still identify the same loci.
        def coordinates(directory):
            return {tuple(row[key] for key in ("action_id", "sequence", "interval_start0", "interval_end0"))
                    for row in table(directory / "locus_plot_hits.tsv")}
        if (args.previous / "locus_plot_hits.tsv").is_file() and coordinates(args.previous) == coordinates(args.out):
            for name in read_tables:
                shutil.copyfile(args.previous / name, args.out / name)
            if table(args.out / "locus_read_boundaries.tsv"):
                add_read_support_guide(args.out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    verify = commands.add_parser("verify-database")
    for name in ("db", "database", "taxdump"):
        verify.add_argument("--" + name, type=Path, required=True)
    rebuild = commands.add_parser("rebuild-report")
    for name in ("previous", "fasta", "out"):
        rebuild.add_argument("--" + name, type=Path, required=True)
    blob = commands.add_parser("blob-input")
    for name in ("fasta", "exclusions", "out"):
        blob.add_argument("--" + name, type=Path, required=True)
    blob_tax = commands.add_parser("blob-taxonomy")
    for name in ("taxdump", "hits", "out"):
        blob_tax.add_argument("--" + name, type=Path, required=True)
    blob_tax.add_argument("--taxid", type=int, required=True)
    for command in ("select", "summarize"):
        sub = commands.add_parser(command)
        sub.add_argument("--prepared", type=Path, required=True)
        sub.add_argument("--taxdump", type=Path, required=True)
        sub.add_argument("--taxid", type=int, required=True)
        sub.add_argument("--host-rank", default="phylum")
        sub.add_argument("--out", type=Path, required=True)
        if command == "select":
            sub.add_argument("--fasta", type=Path, required=True)
            sub.add_argument("--hits", type=Path, nargs="*", default=[])
            sub.add_argument("--megablast-hits", type=Path, nargs="*", default=[])
            sub.add_argument("--megablast", action="store_true")
            sub.add_argument("--followup-max-bp", type=int, default=100_000)
        else:
            for name in ("selected", "metrics", "database"):
                sub.add_argument("--" + name, type=Path, required=True)
            sub.add_argument("--coverage", type=Path)
            sub.add_argument("--coverage-source", default="provided-bam")
            sub.add_argument("--nucleotide-hits", type=Path, nargs="*", default=[])
            sub.add_argument("--targeted", action="store_true")
            sub.add_argument("--nucleotide-db")
            sub.add_argument("--nucleotide-manifest", type=Path)
    args = parser.parse_args()
    if args.command == "verify-database":
        identity = protein_identity(args.database)
        if digest(args.db) != identity["database_sha256"]:
            raise ValueError("DIAMOND index does not match its database manifest")
        for name, checksum in identity["taxonomy_sha256"].items():
            if digest(args.taxdump / name) != checksum:
                raise ValueError(f"Search taxonomy differs from the database manifest: {name}")
    elif args.command == "rebuild-report":
        rebuild_report(args)
    elif args.command == "blob-input":
        # BlobToolKit's FASTA parser cannot compute GC for all-ambiguous sequences.
        excluded = {row["sequence"] for row in table(args.exclusions)}
        if not excluded:
            args.out.symlink_to(args.fasta.resolve())
        else:
            with args.out.open("w") as output:
                for header, sequence in records(args.fasta):
                    if header.split()[0] not in excluded:
                        write_record(output, header, sequence)
    elif args.command == "blob-taxonomy":
        blob_taxonomy(args)
    elif args.command == "select":
        if args.followup_max_bp < 1:
            parser.error("--followup-max-bp must be a positive integer")
        select(args)
    else:
        summarize(args)


if __name__ == "__main__":
    main()
