process ASSEMBLY_ELIGIBILITY {
    label 'assembly_eligibility'
    publishDir "${params.outdir}", mode: params.publish_dir_mode
    input:
    val expected
    tuple val(observed), path(fais)
    val harmonization_references
    path auditor
    path harmonizer
    output:
    path 'assembly_eligibility.json', emit: manifest
    path 'assembly_eligibility.tsv', emit: table
    path 'assembly_eligibility.md', emit: report
    script:
    def nz = { key, fallback -> params[key] != null ? params[key] : fallback }
    def settings = [
        min_individuals: nz('pangenome_min_individuals', 2) as int,
        run_pangenome: params.run_pangenome as boolean, qc_mode: params.qc_mode,
        min_scaffold_bp: nz('harmonize_min_scaffold_bp', nz('finalize_min_scaffold_bp', 1000000)) as long,
        chromosome_method: nz('harmonize_chromosome_method', 'dropoff'),
        dropoff_ratio: nz('harmonize_dropoff_ratio', 2.0),
        dropoff_min_frac: nz('harmonize_dropoff_min_frac', 0.5),
        min_chrom_frac: nz('harmonize_min_chrom_frac', 0.1),
        min_cut_ratio: nz('harmonize_voter_min_cut_ratio', 0.0),
        min_genome_frac: nz('harmonize_voter_min_genome_frac', 0.8),
        min_n50_ratio: nz('harmonize_voter_min_n50_ratio', 0.2),
        require_dropoff: params.harmonize_voter_require_dropoff != false,
        min_voters: nz('harmonize_min_voters', 3) as int
    ]
    def payload = groovy.json.JsonOutput.toJson([expected: expected, observed: observed,
        harmonization_references: harmonization_references, settings: settings])
    """
    set -euo pipefail
    cat > eligibility_input.json <<'ELIGIBILITY_INPUT'
${payload}
ELIGIBILITY_INPUT
    python3 ${auditor} --input eligibility_input.json
    """
}