include { PSMC_PLAN; PSMC_ASSEMBLY_CALL; PSMC_MAP_HIFI; PSMC_READ_CALL; PSMC_INPUT;
          PSMC_FIT; PSMC_SPLIT; PSMC_BOOTSTRAP; PSMC_EXPORT; PSMC_RENDER; PSMC_INDEX } from '../modules/psmc.nf'

workflow PSMC_ANALYSIS {
    take:
    eligibility
    references
    hifi
    main:
    reports=Channel.empty();versions=Channel.empty()
    if(params.run_psmc) {
        if(!(params.psmc_mode in ['auto','assemblies','hifi'])) error 'psmc_mode must be auto, assemblies or hifi'
        ['psmc_min_contig_bp','psmc_bin_size','psmc_min_good_bases','psmc_min_callable_bp','psmc_min_hets',
         'psmc_iterations','psmc_max_time','psmc_theta_rho_ratio','psmc_bootstrap_block_bp','psmc_bootstrap_max_forks',
         'psmc_max_pileup_depth'].each { key ->
            if((params[key] as double)<=0) error "${key} must be positive"
        }
        ['psmc_min_depth','psmc_max_depth','psmc_min_mapq','psmc_min_baseq','psmc_min_gq','psmc_indel_padding',
         'psmc_bootstraps','psmc_plot_iteration','psmc_seed'].each { key ->
            if((params[key] as double)<0) error "${key} must be nonnegative"
        }
        ['psmc_min_contig_bp','psmc_bin_size','psmc_min_good_bases','psmc_min_callable_bp','psmc_min_hets',
         'psmc_iterations','psmc_bootstraps','psmc_plot_iteration','psmc_bootstrap_block_bp','psmc_bootstrap_max_forks',
         'psmc_seed','psmc_min_depth','psmc_max_depth','psmc_min_mapq','psmc_min_baseq','psmc_min_gq',
         'psmc_indel_padding','psmc_max_pileup_depth'].each { key ->
            def value=params[key] as double
            if(!Double.isFinite(value) || value!=Math.rint(value)) error "${key} must be a finite integer"
        }
        ['psmc_max_time','psmc_theta_rho_ratio','psmc_min_callable_fraction','psmc_min_allele_balance'].each { key ->
            if(!Double.isFinite(params[key] as double)) error "${key} must be finite"
        }
        if((params.psmc_plot_iteration as int)>(params.psmc_iterations as int)) error 'PSMC plot iteration exceeds fitting iterations'
        if((params.psmc_min_callable_fraction as double)<0 || (params.psmc_min_callable_fraction as double)>1)
            error 'psmc_min_callable_fraction must be between 0 and 1'
        if((params.psmc_min_allele_balance as double)<0 || (params.psmc_min_allele_balance as double)>.5)
            error 'psmc_min_allele_balance must be between 0 and 0.5'
        if(!(params.psmc_pattern.toString()==~ /\d+(?:\*\d+)?(?:\+\d+(?:\*\d+)?)*\z/)) error 'Invalid PSMC interval pattern'
        def groups=params.psmc_pattern.toString().split(/\+/).collectMany { token ->
            def fields=token.split(/\*/).collect { it as int }
            if(fields.any { it<=0 }) error 'PSMC interval pattern values must be positive'
            Collections.nCopies(fields.size()==2 ? fields[0] : 1,fields[-1])
        }
        if(groups.sum()>255) error 'PSMC interval pattern exceeds 255 atomic intervals'
        if(groups.size()<2 && !params.psmc_show_last_interval) error 'A single PSMC parameter group requires psmc_show_last_interval=true'
        if((params.psmc_bootstrap_block_bp as int)%(params.psmc_bin_size as int)!=0)
            error 'psmc_bootstrap_block_bp must be divisible by psmc_bin_size'
        if((params.psmc_mutation_rate==null)!=(params.psmc_generation_time==null))
            error 'Set both PSMC scaling parameters or leave both null'
        def helper=file("${projectDir}/py_scripts/psmc_results.py",checkIfExists:true)
        def version_helper=file("${projectDir}/py_scripts/psmc_versions.py",checkIfExists:true)
        def index_helper=file("${projectDir}/py_scripts/psmc_report_index.py",checkIfExists:true)
        def absent=file("${projectDir}/assets/NO_PSMC_INPUT",checkIfExists:true)
        def mask_settings=[min_contig_bp:params.psmc_min_contig_bp,bin_size:params.psmc_bin_size,
            min_good_bases:params.psmc_min_good_bases,min_callable_bp:params.psmc_min_callable_bp,
            min_callable_fraction:params.psmc_min_callable_fraction,min_hets:params.psmc_min_hets,
            min_depth:params.psmc_min_depth,max_depth:params.psmc_max_depth,min_gq:params.psmc_min_gq,
            min_allele_balance:params.psmc_min_allele_balance,indel_padding:params.psmc_indel_padding,
            min_mapq:params.psmc_min_mapq,min_baseq:params.psmc_min_baseq,max_pileup_depth:params.psmc_max_pileup_depth]
        def fit_settings=[iterations:params.psmc_iterations,max_time:params.psmc_max_time,
                          ratio:params.psmc_theta_rho_ratio,pattern:params.psmc_pattern]
        PSMC_PLAN(eligibility,
            references.map { taxid,path -> new groovy.json.JsonSlurperClassic().parseText(path.text) }.toList(),
            hifi.map { meta,fq -> [sample:meta.sample.toString(),fastq:fq.toString()] }.toList(),
            params.psmc_manifest ? file(params.psmc_manifest,checkIfExists:true) : absent,
            file("${projectDir}/py_scripts/psmc_plan.py",checkIfExists:true))
        jobs=PSMC_PLAN.out.plan.flatMap { path -> new groovy.json.JsonSlurperClassic().parseText(path.text).jobs }
        def meta_for={ j -> [id:j.id,sample:j.sample,taxid:j.taxid,species:j.species,mode:j.mode,
            reference_id:j.reference.id,reference_sha256:j.reference.sha256,representation:j.representation] }
        PSMC_ASSEMBLY_CALL(jobs.filter { it.mode=='assemblies' }.map { j ->
            tuple(meta_for(j),file(j.reference.fasta,checkIfExists:true),j.assemblies.collect { file(it,checkIfExists:true) })
        },version_helper)
        PSMC_MAP_HIFI(jobs.filter { it.mode=='hifi' }.map { j ->
            tuple(meta_for(j),file(j.reference.fasta,checkIfExists:true),file(j.reads,checkIfExists:true))
        },version_helper)
        PSMC_READ_CALL(PSMC_MAP_HIFI.out.bam,
            [mapq:params.psmc_min_mapq,baseq:params.psmc_min_baseq,max_pileup_depth:params.psmc_max_pileup_depth,
             min_contig_bp:params.psmc_min_contig_bp],version_helper)
        call_inputs=PSMC_ASSEMBLY_CALL.out.calls.mix(PSMC_READ_CALL.out.calls)
            .map { m,fa,vcf,tbi,bed,depth -> tuple(m.id,m,fa,vcf,tbi,bed,depth) }
            .join(jobs.map { j -> tuple(j.id,
                j.include_bed ? file(j.include_bed,checkIfExists:true) : absent,
                j.exclude_bed ? file(j.exclude_bed,checkIfExists:true) : absent) })
            .map { id,m,fa,vcf,tbi,bed,depth,inc,exc -> tuple(m,fa,vcf,tbi,bed,depth,mask_settings,inc,exc) }
        PSMC_INPUT(call_inputs,file("${projectDir}/py_scripts/psmc_input.py",checkIfExists:true),version_helper)
        assessed=PSMC_INPUT.out.input.branch { m,dir ->
            ready: new groovy.json.JsonSlurperClassic().parseText(dir.resolve('input_qc.json').text).status=='ready'
            not_fitted: true
        }
        fit_inputs=assessed.ready.map { m,dir -> tuple(m,dir.resolve('input.psmcfa')) }
        PSMC_FIT(fit_inputs,fit_settings,helper,version_helper)
        models=PSMC_FIT.out.model
        versions=PSMC_ASSEMBLY_CALL.out.versions.mix(PSMC_MAP_HIFI.out.versions,PSMC_READ_CALL.out.versions,
            PSMC_INPUT.out.versions,PSMC_FIT.out.versions)
        if((params.psmc_bootstraps as int)>0) {
            PSMC_SPLIT(fit_inputs,(params.psmc_bootstrap_block_bp as int).intdiv(params.psmc_bin_size as int))
            PSMC_BOOTSTRAP(PSMC_SPLIT.out.split.combine(Channel.fromList((1..(params.psmc_bootstraps as int)).toList())),
                fit_settings,params.psmc_seed,helper,version_helper)
            models=models.mix(PSMC_BOOTSTRAP.out.model)
            versions=versions.mix(PSMC_BOOTSTRAP.out.versions)
        }
        scaling=jobs.map { j -> tuple(j.id,[mutation_rate:j.mutation_rate,generation_time:j.generation_time]) }
        PSMC_EXPORT(models.map { m,path -> tuple(groupKey(m.id,(params.psmc_bootstraps as int)+1),m,path) }
            .groupTuple().map { key,metas,paths -> tuple(key.getGroupTarget(),metas[0],paths) }
            .join(scaling).map { id,m,paths,s -> tuple(m,paths,s) },helper)
        plot_context=jobs.map { j -> tuple(j.id,meta_for(j)+[mutation_rate:j.mutation_rate,generation_time:j.generation_time,
            iteration:params.psmc_plot_iteration,bootstrap_count:params.psmc_bootstraps]) }
        ready_reports=PSMC_EXPORT.out.curves.map { m,dir -> tuple(m.id,m,dir) }
            .join(assessed.ready.map { m,dir -> tuple(m.id,dir) })
            .map { id,m,curves,input -> tuple(id,m,input,curves) }
        skipped_reports=assessed.not_fitted.map { m,input -> tuple(m.id,m,input,absent) }
        PSMC_RENDER(ready_reports.mix(skipped_reports).join(plot_context)
            .map { id,m,input,curves,context -> tuple(m,input,curves,context) },
            file("${projectDir}/r_scripts/plot_psmc.R",checkIfExists:true))
        versions=versions.mix(PSMC_RENDER.out.versions)
        PSMC_INDEX(PSMC_PLAN.out.plan,PSMC_RENDER.out.report.map { m,dir -> dir }.toList(),index_helper)
        reports=PSMC_INDEX.out.report.mix(PSMC_INDEX.out.summary)
    }
    emit:
    report=reports
    versions=versions
}
