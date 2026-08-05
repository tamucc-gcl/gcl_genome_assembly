process DIAMOND_DB_GET {
  tag "diamond_db_get"
  label 'diamond' 

  input:
    val profile           // "uniprot_ref_proteomes" | "custom"
    path db_dir
    val db_name
    val fasta_url
    val taxonmap_url
    path taxdump_dir
    val force_download

  output:
    path "${db_dir}/${db_name}.dmnd", emit: dmnd
    path "${taxdump_dir}",            emit: taxdump_out

  container "quay.io/biocontainers/diamond:2.1.10--h43eeafb_0"

  script:
  """
  set -euo pipefail

  mkdir -p "${db_dir}"
  mkdir -p "${taxdump_dir}"

  DMND="${db_dir}/${db_name}.dmnd"
  SENTINEL="${db_dir}/.${db_name}.ready"

  # --- get taxdump if missing ---
  if [ "${force_download}" = "true" ] || [ ! -s "${taxdump_dir}/names.dmp" ] || [ ! -s "${taxdump_dir}/nodes.dmp" ]; then
    echo "[DIAMOND_DB_GET] Downloading NCBI taxdump into ${taxdump_dir}"
    # NCBI taxonomy FTP referenced by NCBI docs; taxdump is standard. :contentReference[oaicite:6]{index=6}
    curl -L -o "${taxdump_dir}/taxdump.tar.gz" "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    tar -xzf "${taxdump_dir}/taxdump.tar.gz" -C "${taxdump_dir}" names.dmp nodes.dmp
    rm -f "${taxdump_dir}/taxdump.tar.gz"
  else
    echo "[DIAMOND_DB_GET] Taxdump present; skipping."
  fi

  # --- determine download URLs ---
  if [ "${profile}" = "uniprot_ref_proteomes" ]; then
    # NOTE: UniProt reference proteome packaging can change by release.
    # For fully reproducible builds, set fasta_url and taxonmap_url explicitly in config.
    : "${fasta_url:?Set params.diamond_fasta_url for uniprot_ref_proteomes (release-specific URL)}"
    : "${taxonmap_url:?Set params.diamond_taxonmap_url for uniprot_ref_proteomes (mapping file URL)}"
  else
    : "${fasta_url:?Set params.diamond_fasta_url for profile=custom}"
    : "${taxonmap_url:?Set params.diamond_taxonmap_url for profile=custom}"
  fi

  # --- build db if missing ---
  if [ "${force_download}" = "true" ] || [ ! -s "${DMND}" ] || [ ! -f "${SENTINEL}" ]; then
    echo "[DIAMOND_DB_GET] Downloading FASTA + taxonmap and building ${DMND}"

    curl -L -o "${db_dir}/${db_name}.fasta.gz" "${fasta_url}"
    curl -L -o "${db_dir}/${db_name}.taxonmap.tsv.gz" "${taxonmap_url}"

    gunzip -f "${db_dir}/${db_name}.fasta.gz"
    gunzip -f "${db_dir}/${db_name}.taxonmap.tsv.gz"

    diamond makedb \
      --in "${db_dir}/${db_name}.fasta" \
      --db "${db_dir}/${db_name}" \
      --taxonmap "${db_dir}/${db_name}.taxonmap.tsv" \
      --taxonnodes "${taxdump_dir}/nodes.dmp" \
      --taxonnames "${taxdump_dir}/names.dmp" \
      --threads ${task.cpus}

    date -Is > "${SENTINEL}"
  else
    echo "[DIAMOND_DB_GET] DIAMOND db present at ${DMND}; skipping build."
  fi
  """
}
