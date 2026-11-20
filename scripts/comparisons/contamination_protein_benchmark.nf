nextflow.enable.dsl = 2

params.benchmark_project = null
params.benchmark_trace = null
params.benchmark_protein_db_dir = null
params.benchmark_samples = 12
params.benchmark_jobs = 8
params.benchmark_threads = 8
params.benchmark_scratch = '/scratch'
params.benchmark_node = null

def benchQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

process SELECT_PROTEIN_BENCHMARK {
    conda 'conda-forge::python=3.11'
    cpus 1
    memory '2 GB'
    time '30m'
    scratch false
    publishDir "${params.outdir}", mode: 'copy', pattern: 'queries/selection.json'
    input:
    val project
    val trace
    path helper
    output:
    path 'queries'
    script:
    """
    python3 ${helper} select --project ${benchQuote(project)} --trace ${benchQuote(trace)} \
      --samples ${params.benchmark_samples} --out queries
    """
}

process RUN_PROTEIN_BENCHMARK {
    conda 'bioconda::diamond=2.1.14 conda-forge::python=3.11 conda-forge::time=1.9'
    cpus { (params.benchmark_jobs as int) * (params.benchmark_threads as int) }
    memory '0 B'
    clusterOptions { '--exclusive --mem=0' + (params.benchmark_node ? " --nodelist=${params.benchmark_node}" : '') }
    time '48h'
    scratch false
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path queries
    val bundle
    path helper
    output:
    path 'benchmark'
    script:
    """
    python3 ${helper} run --queries ${benchQuote(queries)} --bundle ${benchQuote(bundle)} \
      --scratch ${benchQuote(params.benchmark_scratch)} --jobs ${params.benchmark_jobs} \
      --threads ${params.benchmark_threads} --out benchmark
    """
}

workflow {
    if (!params.benchmark_project || !params.benchmark_trace || !params.benchmark_protein_db_dir)
        error 'Supply --benchmark_project, --benchmark_trace and --benchmark_protein_db_dir'
    ['benchmark_samples', 'benchmark_jobs', 'benchmark_threads'].each { name ->
        if (!(params[name].toString() ==~ /[1-9][0-9]*/)) error "Supply a positive integer for ${name}"
    }
    if (params.benchmark_node && !(params.benchmark_node.toString() ==~ /[A-Za-z0-9._-]+/))
        error 'Supply one Slurm node name for --benchmark_node'
    def helper = file("${projectDir}/contamination_protein_benchmark.py", checkIfExists: true)
    SELECT_PROTEIN_BENCHMARK(file(params.benchmark_project).toAbsolutePath().toString(),
        file(params.benchmark_trace, checkIfExists: true).toAbsolutePath().toString(), helper)
    RUN_PROTEIN_BENCHMARK(SELECT_PROTEIN_BENCHMARK.out,
        file(params.benchmark_protein_db_dir, checkIfExists: true).toAbsolutePath().toString(), helper)
}
