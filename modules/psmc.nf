/* Optional PSMC: stable expensive inputs, separate masks/fits/bootstrap/plots. */
process PSMC_REFERENCE {
    tag "${selection.sample}"
    cpus 1
    memory '2 GB'
    time '2h'
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/post_assembly/psmc/${selection.taxid}/${selection.sample}/reference", mode: params.publish_dir_mode
    input:
    tuple val(selection), path(fasta, stageAs:'reference.fa'), path(fai, stageAs:'reference.fa.fai')
    path verifier
    path fasta_inspector
    output:
    tuple val(selection), path('reference.json'), emit: reference
    path 'reference.md', emit: report
    script:
    """
    set -euo pipefail
    cat > selection.json <<'PSMC_REFERENCE_SELECTION'
${groovy.json.JsonOutput.toJson(selection)}
PSMC_REFERENCE_SELECTION
    python3 ${verifier} --selection selection.json --fasta ${fasta} --fai ${fai}
    """
}

process PSMC_PLAN {
    cpus 1
    memory '2 GB'
    time '1h'
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/post_assembly/psmc", mode: params.publish_dir_mode
    input:
    path eligibility
    val references
    val reads
    path manifest
    path planner
    output:
    path 'psmc_plan.json', emit: plan
    path 'psmc_plan.md', emit: report
    script:
    def options = params.psmc_mutation_rate == null ? '' : "--mutation-rate ${params.psmc_mutation_rate}"
    options += params.psmc_generation_time == null ? '' : " --generation-time ${params.psmc_generation_time}"
    """
    set -euo pipefail
    cat > references.json <<'PSMC_REFERENCES'
${groovy.json.JsonOutput.toJson(references)}
PSMC_REFERENCES
    cat > reads.json <<'PSMC_READS'
${groovy.json.JsonOutput.toJson(reads)}
PSMC_READS
    python3 ${planner} --eligibility ${eligibility} --references references.json --reads reads.json \\
        --manifest '${manifest}' --mode '${params.psmc_mode}' --min-contig-bp ${params.psmc_min_contig_bp} ${options} --out psmc_plan.json
    """
}

process PSMC_ASSEMBLY_CALL {
    tag "${meta.id}"
    cpus { Math.min(16, params.max_cpus as int) }
    memory '32 GB'
    time '48h'
    conda "${projectDir}/environments/psmc.yaml"
    publishDir "${params.outdir}/post_assembly/psmc/${meta.taxid}/${meta.id}/genotypes", mode: params.publish_dir_mode,
        pattern: 'calls/*', saveAs: { it.replaceFirst('^calls/','') }
    input:
    tuple val(meta), path(reference, stageAs:'reference.fa'), path(haplotypes, stageAs:'haplotypes/*')
    val alignment_preset
    path caller
    path helper
    output:
    tuple val(meta), path('reference.fa'), path('calls/diploid.dip.vcf.gz'), path('calls/diploid.dip.vcf.gz.tbi'),
        path('calls/diploid.dip.bed'), path('calls/NO_DEPTH'), emit: calls
    path 'calls/dipcall_settings.json', emit: settings
    path 'calls/dipcall.log', emit: log
    path 'versions.tsv', emit: versions
    script:
    """
    set -euo pipefail
    samtools faidx reference.fa
    python3 ${caller} --preset '${alignment_preset}' --threads-per-hap ${Math.max(1, (task.cpus as int).intdiv(2))} --haplotypes ${haplotypes[0]} ${haplotypes[1]}
    bcftools index --tbi diploid.dip.vcf.gz
    mkdir calls
    cp diploid.dip.vcf.gz diploid.dip.vcf.gz.tbi diploid.dip.bed dipcall.log dipcall_settings.json calls/
    touch calls/NO_DEPTH
    python3 ${helper} versions --tools dipcall minimap2 samtools bcftools
    """
}

process PSMC_MAP_HIFI {
    tag "${meta.id}"
    cpus { Math.min(16, params.max_cpus as int) }
    memory '32 GB'
    time '48h'
    conda "${projectDir}/environments/psmc.yaml"
    input:
    tuple val(meta), path(reference, stageAs:'reference.fa'), path(reads, stageAs:'reads.fastq.gz')
    path helper
    output:
    tuple val(meta), path('reference.fa'), path('hifi.bam'), path('hifi.bam.bai'), emit: bam
    path 'versions.tsv', emit: versions
    script:
    """
    set -euo pipefail
    samtools faidx reference.fa
    minimap2 -ax map-hifi -t ${Math.max(1,(task.cpus as int)-3)} -R '@RG\\tID:${meta.id}\\tSM:${meta.id}' reference.fa reads.fastq.gz \\
        | samtools sort -@ ${Math.min(2,Math.max(0,(task.cpus as int)-2))} -m 1G -o hifi.bam -
    samtools index hifi.bam
    python3 ${helper} versions --tools minimap2 samtools
    """
}

process PSMC_READ_CALL {
    tag "${meta.id}"
    cpus { Math.min(2*(params.psmc_read_call_workers as int), params.max_cpus as int) }
    memory { "${Math.max(16, (task.cpus as int)*2)} GB" }
    time '72h'
    conda "${projectDir}/environments/psmc.yaml"
    publishDir "${params.outdir}/post_assembly/psmc/${meta.taxid}/${meta.id}/genotypes", mode: params.publish_dir_mode,
        pattern: 'calls/*', saveAs: { it.replaceFirst('^calls/','') }
    input:
    tuple val(meta), path(reference, stageAs:'reference.fa'), path(bam), path(bai)
    val filters
    val workers
    path caller
    path helper
    output:
    tuple val(meta), path('reference.fa'), path('calls/diploid.vcf.gz'), path('calls/diploid.vcf.gz.tbi'),
        path('calls/NO_CONFIDENT_BED'), path('calls/depth.tsv'), emit: calls
    path 'calls/calling_manifest.json', emit: manifest
    path 'calls/calling.log', emit: log
    path 'versions.tsv', emit: versions
    script:
    """
    set -euo pipefail
    samtools faidx reference.fa
    python3 ${caller} --reference reference.fa --bam ${bam} --out calls \\
        --min-contig-bp ${filters.min_contig_bp} --mapq ${filters.mapq} --baseq ${filters.baseq} \\
        --max-pileup-depth ${filters.max_pileup_depth} --cpus ${task.cpus} --workers ${workers}
    python3 ${helper} versions --tools samtools bcftools
    """
}

process PSMC_INPUT {
    tag "${meta.id}"
    cpus 1
    memory '24 GB'
    time '48h'
    conda "${projectDir}/environments/psmc.yaml"
    publishDir "${params.outdir}/post_assembly/psmc/${meta.taxid}/${meta.id}", mode: params.publish_dir_mode,
        pattern: 'input'
    input:
    tuple val(meta), path(reference, stageAs:'reference.fa'), path(vcf), path(tbi), path(confident), path(depth),
        val(filters), path(include_bed, stageAs:'scope/*'), path(exclude_bed, stageAs:'exclude/*')
    path adapter
    path helper
    output:
    tuple val(meta), path('input'), emit: input
    path 'versions.tsv', emit: versions
    script:
    def provenance=[identity:meta,filters:filters,vcf:vcf.toString(),confident_bed:confident.toString(),
        caller:meta.mode=='assemblies' ? 'Dipcall 0.3' : 'bcftools 1.21 multiallelic diploid caller; complete HiFi read set']
    """
    set -euo pipefail
    samtools faidx reference.fa
    python3 ${adapter} prepare --reference reference.fa --reference-sha256 '${meta.reference_sha256}' \\
        --vcf ${vcf} --mode ${meta.mode} --confident ${confident} --depth ${depth} \\
        --include ${include_bed} --exclude ${exclude_bed} --min-contig-bp ${filters.min_contig_bp} \\
        --min-depth ${filters.min_depth} --max-depth ${filters.max_depth} --min-gq ${filters.min_gq} \\
        --indel-padding ${filters.indel_padding} --min-allele-balance ${filters.min_allele_balance} --out input
    fq2psmcfa -s ${filters.bin_size} -g ${filters.min_good_bases} input/diploid.fa.gz > input/input.psmcfa
    python3 ${adapter} finalize --out input --bin-size ${filters.bin_size} \\
        --min-callable-bp ${filters.min_callable_bp} --min-callable-fraction ${filters.min_callable_fraction} --min-hets ${filters.min_hets}
    cat > input/provenance.json <<'PSMC_INPUT_PROVENANCE'
${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(provenance))}
PSMC_INPUT_PROVENANCE
    python3 ${helper} versions --tools psmc samtools python
    """
}

process PSMC_FIT {
    tag "${meta.id}"
    cpus 1
    memory '8 GB'
    time '48h'
    conda "${projectDir}/environments/psmc.yaml"
    publishDir "${params.outdir}/post_assembly/psmc/${meta.taxid}/${meta.id}/fits", mode: params.publish_dir_mode,
        pattern:'observed.psmc'
    input:
    tuple val(meta), path(input)
    val settings
    path helper
    path version_helper
    output:
    tuple val(meta), path('observed.psmc'), emit: model
    path 'versions.tsv', emit: versions
    script:
    """
    set -euo pipefail
    psmc -I0 -N ${settings.iterations} -t ${settings.max_time} -r ${settings.ratio} -p '${settings.pattern}' -o observed.psmc ${input}
    python3 ${helper} validate --input observed.psmc --last-iteration ${settings.iterations}
    python3 ${version_helper} versions --tools psmc
    """
}

process PSMC_SPLIT {
    tag "${meta.id}"
    cpus 1
    memory '2 GB'
    time '1h'
    conda "${projectDir}/environments/psmc.yaml"
    input:
    tuple val(meta), path(input)
    val block_bins
    output:
    tuple val(meta), path('split.psmcfa'), emit: split
    script:
    """
    set -euo pipefail
    splitfa ${input} ${block_bins} > split.psmcfa
    test -s split.psmcfa
    """
}

process PSMC_BOOTSTRAP {
    tag "${meta.id}:bootstrap_${replicate}"
    cpus 1
    memory '8 GB'
    time '48h'
    maxForks (params.psmc_bootstrap_max_forks as int)
    conda "${projectDir}/environments/psmc.yaml"
    publishDir "${params.outdir}/post_assembly/psmc/${meta.taxid}/${meta.id}/fits", mode: params.publish_dir_mode,
        pattern:'bootstrap_*'
    input:
    tuple val(meta), path(split), val(replicate)
    val settings
    val seed
    path helper
    path version_helper
    output:
    tuple val(meta), path('bootstrap_*.psmc'), emit: model
    tuple val(meta), path('bootstrap_*.resampling.json'), emit: resampling
    path 'versions.tsv', emit: versions
    script:
    def name=String.format('bootstrap_%04d',replicate as int)
    """
    set -euo pipefail
    python3 ${helper} resample --input ${split} --sample '${meta.id}' --replicate ${replicate} --seed ${seed} --out resampled.psmcfa
    mv resampling.json ${name}.resampling.json
    psmc -I0 -N ${settings.iterations} -t ${settings.max_time} -r ${settings.ratio} -p '${settings.pattern}' -o ${name}.psmc resampled.psmcfa
    python3 ${helper} validate --input ${name}.psmc --last-iteration ${settings.iterations}
    python3 ${version_helper} versions --tools psmc python
    """
}

process PSMC_EXPORT {
    tag "${meta.id}"
    cpus 1
    memory '2 GB'
    time '1h'
    conda 'conda-forge::python=3.11'
    input:
    tuple val(meta), path(models), val(scaling)
    path helper
    output:
    tuple val(meta), path('curves'), emit: curves
    script:
    def options=scaling.mutation_rate==null ? '' : "--mutation-rate ${scaling.mutation_rate} --generation-time ${scaling.generation_time}"
    options += params.psmc_show_last_interval ? ' --show-last' : ''
    """
    set -euo pipefail
    python3 ${helper} export --files ${models} --out curves --iteration ${params.psmc_plot_iteration} \\
        --bin-size ${params.psmc_bin_size} --expected-bootstraps ${params.psmc_bootstraps} ${options}
    """
}

process PSMC_RENDER {
    tag "${meta.id}"
    cpus 1
    memory '4 GB'
    time '1h'
    conda "${projectDir}/environments/psmc_plot.yaml"
    publishDir "${params.outdir}/post_assembly/psmc/${meta.taxid}/${meta.id}", mode: params.publish_dir_mode,
        pattern:'report'
    input:
    tuple val(meta), path(input), path(curves), val(context)
    path renderer
    output:
    tuple val(meta), path('report'), emit: report
    path 'versions.tsv', emit: versions
    script:
    """
    set -euo pipefail
    cat > context.json <<'PSMC_CONTEXT'
${groovy.json.JsonOutput.toJson(context)}
PSMC_CONTEXT
    Rscript ${renderer} context.json ${input} ${curves} report
    printf 'R\\t%s\\n' "\$(Rscript --version 2>&1)" > versions.tsv
    """
}

process PSMC_INDEX {
    cpus 1
    memory '2 GB'
    time '1h'
    conda "${projectDir}/environments/psmc_plot.yaml"
    publishDir "${params.outdir}/post_assembly/psmc", mode: params.publish_dir_mode,
        pattern:'psmc/*', saveAs:{it.replaceFirst('^psmc/','')}
    input:
    path plan
    path directories, stageAs:'reports/report??'
    path helper
    path renderer
    output:
    path 'psmc/README.md', emit: report
    path 'psmc/summary.tsv', emit: summary
    path 'psmc/*.png', optional:true, emit: figures
    path 'psmc/*.pdf', optional:true, emit: pdfs
    path 'psmc/combined_curves.tsv', optional:true, emit: curves
    path 'psmc/psmc_plan.json', emit: plan_copy
    path 'versions.tsv', emit: versions
    script:
    """
    set -euo pipefail
    python3 ${helper} index --plan ${plan} --directories ${directories} --out psmc
    Rscript ${renderer} psmc
    printf 'R\\t%s\\n' "\$(Rscript --version 2>&1)" > versions.tsv
    """
}
