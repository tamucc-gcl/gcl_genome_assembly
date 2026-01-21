process FCS_BLOB_EVIDENCE_REPORT {
  tag "fcs_blob_evidence_report"

  input:
    path action_report
    path taxonomy_report optional true
    path blobtools_outputs_dir   // expects blobtable.tsv inside if available
    path decontaminated_fasta
    path contaminants_fasta

  output:
    path "fcs_contig_actions.tsv", emit: actions_tsv
    path "fcs_blob_annotated.tsv", emit: annotated_tsv
    path "decontam_evidence.md",   emit: report_md

  script:
  """
  set -euo pipefail

  python - << 'PY'
  import csv
  import re
  from pathlib import Path

  action_report = Path("${action_report}")
  taxonomy_report = Path("${taxonomy_report}") if "${taxonomy_report}" != "null" else None
  blobdir = Path("${blobtools_outputs_dir}")
  blobtable = blobdir / "blobtable.tsv"

  # Parse FASTA headers to get contig names (robust, no dependencies)
  def fasta_ids(fp: Path):
    ids = []
    if not fp.exists():
      return ids
    with fp.open("rt", errors="ignore") as f:
      for line in f:
        if line.startswith(">"):
          ids.append(line[1:].strip().split()[0])
    return ids

  kept_ids = set(fasta_ids(Path("${decontaminated_fasta}")))
  removed_ids = set(fasta_ids(Path("${contaminants_fasta}")))

  # Parse FCS action report into per-contig actions.
  # FCS formats vary; we handle common patterns:
  # - tabular with a header including seq_id/sequence/contig
  # - lines containing contig id + action keyword
  actions = {}
  header = None
  rows = []

  txt = action_report.read_text(errors="ignore").splitlines()
  # find a header row if present
  for i, line in enumerate(txt):
    if not line.strip() or line.lstrip().startswith("#"):
      continue
    # candidate header if it includes seq_id/contig/sequence and action
    low = line.lower()
    if ("seq" in low or "contig" in low or "sequence" in low) and ("action" in low or "result" in low):
      header = re.split(r"\\t+", line.strip())
      # consume remaining lines as tsv-like
      for j in range(i+1, len(txt)):
        ln = txt[j].strip()
        if not ln or ln.startswith("#"):
          continue
        parts = re.split(r"\\t+", ln)
        if len(parts) < 2:
          continue
        rows.append(parts)
      break

  if header:
    # map columns
    hlow = [h.strip().lower() for h in header]
    def col(*names):
      for n in names:
        if n in hlow:
          return hlow.index(n)
      return None
    i_id = col("seq_id", "sequence", "contig", "sequence_id", "seqid", "name", "id")
    i_action = col("action", "result", "decision")
    i_reason = col("reason", "comment", "note", "details")
    for parts in rows:
      if i_id is None or i_action is None: 
        continue
      cid = parts[i_id].strip()
      act = parts[i_action].strip()
      reason = parts[i_reason].strip() if (i_reason is not None and i_reason < len(parts)) else ""
      actions[cid] = (act, reason)
  else:
    # fallback heuristic parse
    # Look for lines like: <contig> <...> <ACTION>
    # We'll record only if we can identify an ID-like token and an action-like token.
    action_words = {"keep","trim","remove","exclude","split","drop","contam","contaminant"}
    for line in txt:
      if not line.strip() or line.lstrip().startswith("#"):
        continue
      parts = re.split(r"[\\t ]+", line.strip())
      if len(parts) < 2:
        continue
      cid = parts[0]
      # action candidate: first token in remainder that looks like an action
      act = ""
      for tok in parts[1:]:
        low = tok.lower()
        if low in action_words:
          act = tok
          break
      if act:
        actions[cid] = (act, "")

  # Write actions table
  # Determine status by FASTA membership (ground truth of what got removed)
  out_actions = Path("fcs_contig_actions.tsv")
  with out_actions.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\\t")
    w.writerow(["contig", "status", "fcs_action", "fcs_reason"])
    all_ids = sorted(kept_ids | removed_ids | set(actions.keys()))
    for cid in all_ids:
      status = "kept" if cid in kept_ids else ("removed" if cid in removed_ids else "unknown")
      act, reason = actions.get(cid, ("", ""))
      w.writerow([cid, status, act, reason])

  # Merge with blobtools table if present
  annotated = Path("fcs_blob_annotated.tsv")
  if blobtable.exists():
    # Read blobtable as dict by contig id (common columns: 'name' or 'contig_id')
    with blobtable.open("r", newline="") as f:
      r = csv.DictReader(f, delimiter="\\t")
      blobrows = list(r)
      # find id column
      id_col = None
      for c in r.fieldnames or []:
        if c.lower() in ("contig_id", "name", "seqid", "id"):
          id_col = c
          break
      blob_by_id = {}
      if id_col:
        for row in blobrows:
          blob_by_id[row.get(id_col, "")] = row

    # Write merged table
    with out_actions.open("r", newline="") as f_in, annotated.open("w", newline="") as f_out:
      r = csv.DictReader(f_in, delimiter="\\t")
      # choose a few useful blob columns if they exist
      wanted = []
      candidate_cols = ["length","gc","gc_content","coverage","cov","cov_mean","taxid","besttaxid","taxname","phylum","class","order","family","genus","species"]
      # discover from first blob row
      sample_blob = next(iter(blob_by_id.values()), {})
      for c in candidate_cols:
        for k in sample_blob.keys():
          if k.lower() == c:
            wanted.append(k)
            break
      fieldnames = r.fieldnames + [f"blob_{c}" for c in wanted] + ["blob_present"]
      w = csv.DictWriter(f_out, delimiter="\\t", fieldnames=fieldnames)
      w.writeheader()
      for row in r:
        cid = row["contig"]
        b = blob_by_id.get(cid)
        row2 = dict(row)
        row2["blob_present"] = "yes" if b else "no"
        if b:
          for c in wanted:
            row2[f"blob_{c}"] = b.get(c, "")
        else:
          for c in wanted:
            row2[f"blob_{c}"] = ""
        w.writerow(row2)
  else:
    # If no blobtable, just copy actions table as annotated output
    annotated.write_text(out_actions.read_text())

  # Write a short markdown report
  kept_n = len(kept_ids)
  rem_n = len(removed_ids)
  unknown_n = sum(1 for _ in (set(actions.keys()) - kept_ids - removed_ids))

  md = Path("decontam_evidence.md")
  lines = []
  lines.append("# Decontamination evidence report (FCS-GX + blobtools2)")
  lines.append("")
  lines.append("## Outputs")
  lines.append("- `fcs_contig_actions.tsv`: per-contig kept/removed + FCS action/reason")
  lines.append("- `fcs_blob_annotated.tsv`: above + blobtools evidence columns when available")
  lines.append("")
  lines.append("## Counts (from FASTA outputs)")
  lines.append(f"- Kept contigs: **{kept_n}**")
  lines.append(f"- Removed contigs: **{rem_n}**")
  if unknown_n:
    lines.append(f"- Contigs present in action report but not in FASTAs: **{unknown_n}** (check parsing / naming)")
  lines.append("")
  if blobtable.exists():
    lines.append("## Blobtools evidence")
    lines.append("- `blobtools_outputs/blobtable.tsv` was found and merged into `fcs_blob_annotated.tsv`.")
  else:
    lines.append("## Blobtools evidence")
    lines.append("- `blobtools_outputs/blobtable.tsv` was **not** found (your blobtools2 build may not support `blobtools table`).")
    lines.append("- You still have `blobtools view/plot` outputs for visual evidence.")
  md.write_text("\\n".join(lines) + "\\n")

  PY
  """
}
