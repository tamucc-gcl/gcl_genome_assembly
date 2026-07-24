/*
========================================================================================
    FIND MITO REFERENCE MODULE
========================================================================================
    Uses MitoHiFi's findMitoReference.py to automatically download the closest
    related mitochondrial genome from NCBI for a given species.

    Runs at the very start of the pipeline (no dependencies on reads or assembly)
    so the reference is ready by the time MITOHIFI needs it.

    Input:
    - species_name: Scientific name (e.g., "Spratelloides delicatulus")

    Output:
    - ref_fasta: Reference mitogenome FASTA
    - ref_gb: Reference mitogenome GenBank annotation
    - ref_info: TSV with reference species, accession, and length
========================================================================================
*/

process FIND_MITO_REFERENCE {
    tag "${species_name}"
    label 'mitohifi'

    publishDir "${params.outdir}/mitogenome", mode: params.publish_dir_mode, saveAs: { filename ->
        filename == 'mito_reference_info.tsv' ? "${taxid}_mito_reference_info.tsv" : null
    }

    input:
    tuple val(taxid), val(species_name)

    output:
    tuple val(taxid), path("*.fasta"),                 emit: ref_fasta
    tuple val(taxid), path("*.gb"),                    emit: ref_gb
    tuple val(taxid), path("mito_reference_info.tsv"), emit: ref_info

    script:
    def min_length = params.mitohifi_ref_min_bp ?: 14000
    """
    set -euo pipefail

    echo "[FIND_MITO_REFERENCE] Searching NCBI for closest mitogenome to: ${species_name}"
    echo "[FIND_MITO_REFERENCE] Minimum reference length: ${min_length} bp"
    echo "[FIND_MITO_REFERENCE] Started: \$(date)"

    findMitoReference.py \\
        --species "${species_name}" \\
        --outfolder . \\
        --min_length ${min_length}

    # Identify downloaded files
    REF_FASTA=\$(ls *.fasta 2>/dev/null | head -n 1)
    REF_GB=\$(ls *.gb 2>/dev/null | head -n 1)

    if [ -z "\${REF_FASTA}" ] || [ -z "\${REF_GB}" ]; then
        echo "[FIND_MITO_REFERENCE] ERROR: findMitoReference.py did not produce expected outputs"
        echo "[FIND_MITO_REFERENCE] Files in output directory:"
        ls -la .
        exit 1
    fi

    # Extract accession and species from the downloaded filename
    # findMitoReference.py names files as <accession>.fasta and <accession>.gb
    ACCESSION=\$(basename "\${REF_FASTA}" .fasta)

    # Parse the GenBank file for organism name
    REF_SPECIES=\$(grep -m1 'ORGANISM' "\${REF_GB}" | sed 's/.*ORGANISM  *//' | xargs)
    REF_LENGTH=\$(grep -v '^>' "\${REF_FASTA}" | tr -d '\\n' | wc -c)

    # Write reference info TSV
    cat > mito_reference_info.tsv <<EOF
accession\tspecies\tlength\tquery_species
\${ACCESSION}\t\${REF_SPECIES}\t\${REF_LENGTH}\t${species_name}
EOF

    echo "[FIND_MITO_REFERENCE] Reference found:"
    echo "  Accession: \${ACCESSION}"
    echo "  Species:   \${REF_SPECIES}"
    echo "  Length:    \${REF_LENGTH} bp"
    echo "[FIND_MITO_REFERENCE] Complete: \$(date)"
    """

    stub:
    """
    echo ">NC_000000.0 stub mitogenome" > NC_000000.0.fasta
    echo "ATCG" >> NC_000000.0.fasta
    touch NC_000000.0.gb
    printf 'accession\\tspecies\\tlength\\tquery_species\\nNC_000000.0\\tStub species\\t16000\\t${species_name}\\n' > mito_reference_info.tsv
    """
}