#!/usr/bin/env python3
"""Assembly identities and species cohorts; FAI-only, no alignments or graph work."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from harmonize_names import select_chromosome_set, scaffold_n50, assign_roles


def representation(sample):
    if sample["ploidy"] == 1 and sample["n_hap"] == 1:
        return "haploid", True
    if sample["ploidy"] == 2 and sample["n_hap"] == 2 and sample["assembler"] == "hifiasm":
        return ("hic_phased" if sample["hic_phasing"] else "partially_phased"), True
    if sample["n_hap"] == 1:
        return "collapsed_or_primary", False
    return "unsupported_ploidy_representation", False


def evaluate(data):
    cfg = data["settings"]
    if cfg["min_individuals"] < 2:
        raise ValueError("pangenome_min_individuals must be at least 2")
    expected, observed = data["expected"], data["observed"]
    by_id = {r["meta"]["id"]: r for r in observed}
    if len(by_id) != len(observed):
        raise ValueError("Duplicate finalized assembly identity")
    records, by_species = [], defaultdict(list)
    seen = set()
    for s in sorted(expected, key=lambda r: r["sample"]):
        rep, eligible = representation(s)
        suffixes = ["primary"] if s["n_hap"] == 1 else ["hap1", "hap2"]
        for hap in suffixes:
            aid = s["sample"] + "_" + hap
            if aid in seen:
                raise ValueError("Duplicate expected assembly identity: " + aid)
            seen.add(aid)
            obs = by_id.get(aid)
            row = dict(id=aid, sample=s["sample"], taxid=str(s["taxid"]),
                       species=s["species"], ploidy=s["ploidy"], haplotype=hap,
                       representation=rep, representation_eligible=eligible,
                       completion="finalized" if obs else "missing",
                       qc="unassessed" if cfg["qc_mode"] == "none" else "see_QC_report",
                       chromosome_scale=False, voter=False, voter_reason="not_assessed",
                       reference_candidate=False, graph_member=False)
            if obs:
                if str(obs["meta"]["taxid"]) != row["taxid"] or obs["meta"]["sample"] != row["sample"]:
                    raise ValueError("Finalized identity disagrees with inputs: " + aid)
                fai = []
                with open(obs["fai_name"]) as handle:
                    for line in handle:
                        f = line.rstrip("\n").split("\t")
                        if len(f) < 2 or int(f[1]) <= 0:
                            raise ValueError("Invalid FAI for " + aid)
                        fai.append((f[0], int(f[1])))
                if not fai or len({n for n, _ in fai}) != len(fai):
                    raise ValueError("Empty/duplicate FAI sequence identities: " + aid)
                chromosomes, metrics = select_chromosome_set(
                    fai, cfg["min_scaffold_bp"], cfg["chromosome_method"],
                    cfg["dropoff_ratio"], cfg["dropoff_min_frac"], cfg["min_chrom_frac"])
                metrics["n50"] = scaffold_n50([n for _, n in fai])
                row.update(metrics=metrics, reference_contigs=[n for n, _ in chromosomes],
                           fasta=obs["fasta"], fai=obs["fai"])
            records.append(row)
            by_species[row["taxid"]].append(row)
    if set(by_id) - seen:
        raise ValueError("Unexpected finalized assembly IDs: " + ",".join(sorted(set(by_id)-seen)))
    refs = {}
    for taxid, aid in data["harmonization_references"]:
        key = str(taxid)
        if key in refs:
            raise ValueError('Duplicate harmonization reference for taxid ' + key)
        if aid not in seen or next(r for r in records if r['id'] == aid)['taxid'] != key:
            raise ValueError('Harmonization reference is absent or belongs to another species: ' + aid)
        refs[key] = aid
    cohorts = []
    for taxid, members in sorted(by_species.items()):
        metrics = {m["id"]: m["metrics"] for m in members if "metrics" in m}
        roles, reasons, flags = assign_roles(
            metrics, cfg["min_cut_ratio"], cfg["min_genome_frac"], cfg["min_n50_ratio"],
            cfg["require_dropoff"], cfg["min_voters"], None)
        for m in members:
            m["voter"] = roles.get(m["id"]) == "voter" and bool(m.get("reference_contigs"))
            m["chromosome_scale"] = m["voter"]
            m["voter_reason"] = reasons.get(m["id"], "missing")
            # Composite names describe the selected reference frame, not proven errors.
            m["reference_candidate"] = m["voter"] and m["representation_eligible"]
            m["reference_warning"] = (
                "composite_chromosome_assignment;review_structure"
                if any("+" in n for n in m.get("reference_contigs", [])) else "")
        kept = [m for m in members if m["representation_eligible"] and m["completion"] == "finalized"]
        candidates = [m for m in kept if m["reference_candidate"]]
        preferred = refs.get(taxid)
        candidates.sort(key=lambda m: (m["id"] != preferred, -m["metrics"]["genome_fraction"],
                                       -m["metrics"]["n50"], m["id"]))
        missing = [m["id"] for m in members if m["completion"] != "finalized"]
        individuals = sorted({m["sample"] for m in kept})
        reasons_skip = []
        if missing:
            reasons_skip.append("incomplete_species_assembly_set")
        if len(individuals) < cfg["min_individuals"]:
            reasons_skip.append("fewer_than_minimum_eligible_individuals")
        if not candidates:
            reasons_skip.append("no_suitable_chromosome_scale_reference")
        ready = not reasons_skip
        cohort = dict(taxid=taxid, species=members[0]["species"], ready=ready,
                      status=("ready" if cfg["run_pangenome"] else "disabled_but_eligible") if ready else "withheld",
                      reasons=reasons_skip, missing=missing, individuals=individuals,
                      flags=flags, members=[])
        cohort['harmonization_reference_id'] = preferred
        cohort['reference_selection'] = 'not_selected_cohort_withheld'
        if ready:
            ref = candidates[0]
            cohort['reference_selection'] = (
                'retained_harmonization_reference' if preferred == ref['id'] else
                'replacement_harmonization_reference_not_eligible' if preferred else
                'ranked_final_assembly_no_harmonization_reference')
            refname = "REF_" + ref["id"].encode().hex()
            cohort.update(reference_id=ref["id"], reference_name=refname,
                          reference_contigs=ref["reference_contigs"])
            for m in sorted(kept, key=lambda m: (m["id"] != ref["id"], m["id"])):
                m["graph_member"] = True
                # Explicit reversible names avoid punctuation collisions and reference-prefix restrictions.
                name = "IND_" + m["sample"].encode().hex()
                if m["haplotype"] in ("hap1", "hap2"):
                    name += "." + m["haplotype"][-1]
                if m["id"] == ref["id"]:
                    name = refname
                m["graph_name"] = name
                cohort["members"].append(dict(id=m["id"], sample=m["sample"],
                    haplotype=m["haplotype"], graph_name=name, fasta=m["fasta"]))
        cohorts.append(cohort)
    return dict(schema_version=2, settings=cfg, assemblies=records, cohorts=cohorts)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    a = p.parse_args()
    result = evaluate(json.loads(Path(a.input).read_text()))
    Path("assembly_eligibility.json").write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    columns = ["taxid", "sample", "id", "ploidy", "haplotype", "representation", "completion",
               "qc", "chromosome_scale", "voter", "voter_reason", "reference_candidate", "reference_warning",
               "graph_member", "graph_name"]
    with open("assembly_eligibility.tsv", "w") as out:
        writer = csv.DictWriter(out, fieldnames=columns, extrasaction="ignore", delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result["assemblies"])
    with open("assembly_eligibility.md", "w") as out:
        out.write("<details>\n<summary>Pangenome eligibility (final assemblies)</summary>\n\n")
        out.write("| Species/taxid | Status | Eligible individuals | Reason |\n|---|---|---:|---|\n")
        for c in result["cohorts"]:
            reason = "; ".join(c["reasons"]) or "Cohort complete; suitable reference available"
            out.write(f'| {c["taxid"]} | {c["status"]} | {len(c["individuals"])} | {reason} |\n')
        out.write("\nEligibility uses declared ploidy and actual assembly mode; it never infers haploidy "
                  "from the other samples. Fragmented phased assemblies may be graph members without "
                  "being voters or reference candidates. Chromosome-scale status is an assembly-size "
                  "heuristic, not a karyotype or completeness guarantee. QC does not silently exclude inputs.\n\n")
        out.write('| Taxid | Harmonization reference | Planned graph reference | Selection reason |\n|---|---|---|---|\n')
        for c in result['cohorts']:
            out.write('| %s | %s | %s | %s |\n' % (c['taxid'], c.get('harmonization_reference_id') or 'none',
                      c.get('reference_id', 'none'), c['reference_selection']))
        out.write('\n')
        out.write("Graph membership in this audit means planned eligibility, even when graph construction is disabled. Graph names map explicitly to biological individuals/haplotypes in "
                  "[the identity table](assembly_eligibility.tsv). A reference path's special graph name "
                  "does not create a new biological individual. Synthetic GREF paths are not biological samples.\n\n")
        out.write("</details>\n")


if __name__ == "__main__":
    main()
