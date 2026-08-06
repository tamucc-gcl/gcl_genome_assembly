/*
========================================================================================
    HARMONIZE SPECIES MODULE
========================================================================================
    Repo location: modules/harmonize_species.nf

    One task per species (meta.taxid) over ALL its long-read assemblies (>= 2).
    Picks an in-batch reference -- the assembly whose chromosome-scale scaffold count is
    nearest the batch median (N50 tie-break), which avoids N50's bias toward spuriously
    fused assemblies, or a pinned id -- aligns every other assembly to it with minimap2,
    and runs py_scripts/harmonize_names.py to emit a per-assembly homology name map.
    Alignments are currently serial within the task
    (fine at a handful of same-species individuals; parallelise later if needed).

    Inputs preserve positional correspondence between `ids` and `assemblies` (both come
    from the same groupTuple), so ids[i] <-> assemblies[i].

    Input : tuple(taxid, [ids], [assemblies]) ; path(resolver)
    Output: name_maps (per assembly) / report (per species) / versions / logs
========================================================================================
*/

process HARMONIZE_SPECIES {
    tag "taxid_${taxid}"
    label 'harmonize_scaffolds'

    publishDir "${params.outdir}/assembly/harmonization", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.tsv' ? null : fn }

    input:
    tuple val(taxid), val(ids), path(assemblies)
    path(resolver)

    output:
    tuple val(taxid), path("*.harmonized_name_map.tsv"), emit: name_maps
    tuple val(taxid), path("${taxid}.reference_id.txt"), emit: reference_id
    path("*.harmonization_report.tsv"),                  emit: report
    path("*.reference_selection.tsv"),                   emit: reference_selection
    path("versions.tsv"),                                emit: versions
    path("*.minimap2.log"),                              emit: logs, optional: true

    script:
    def preset       = params.harmonize_minimap2_preset ?: 'asm5'
    def min_scaffold = params.harmonize_min_scaffold_bp ?: (params.finalize_min_scaffold_bp ?: 1000000)
    def chr_method   = params.harmonize_chromosome_method ?: 'dropoff'
    def drop_ratio   = params.harmonize_dropoff_ratio ?: 2.0
    def drop_minfrac = params.harmonize_dropoff_min_frac ?: 0.5
    def min_frac     = params.harmonize_min_aligned_frac ?: 0.5
    def contained_frac = params.harmonize_contained_frac ?: 0.9
    def demote       = (params.harmonize_demote_contained == false) ? '--no-demote-contained' : '--demote-contained'
    def sec_frac     = params.harmonize_secondary_frac ?: 0.2
    def ref_pref     = params.harmonize_reference_ids ?: ''
    """
    set -euo pipefail

    ids=(${ids.join(' ')})
    fastas=(${assemblies})
    PREF="${ref_pref}"

    # ------------------------------------------------------------------
    # id-based symlinks + index (avoids filename collisions / decouples
    # output names from upstream staging names)
    # ------------------------------------------------------------------
    for i in "\${!ids[@]}"; do
        ln -sf "\${fastas[\$i]}" "\${ids[\$i]}.input.fasta"
        samtools faidx "\${ids[\$i]}.input.fasta"
    done

    # ------------------------------------------------------------------
    # choose reference:
    #   1) a pinned id (params.harmonize_reference_ids) if present in this group; else
    #   2) the assembly whose chromosome-scale scaffold count is nearest the batch
    #      median -- a fusion lowers the count and inflates N50, so N50-selection is
    #      biased toward fused assemblies; count-nearest-median removes that bias.
    #      Ties (incl. the k=2 case, where median == mean) broken by higher N50.
    # ------------------------------------------------------------------
    n50 () {
        awk '{print \$2}' "\$1" | sort -rn \\
          | awk '{a[NR]=\$1; s+=\$1} END{h=s/2; c=0; for(i=1;i<=NR;i++){c+=a[i]; if(c>=h){print a[i]; exit}}}'
    }
    nchrom () { awk -v m=${min_scaffold} '\$2>=m{c++} END{print c+0}' "\$1"; }

    # per-assembly selection metrics (also published as provenance)
    printf 'id\\tn_chrom_scaffolds\\tscaffold_n50\\n' > ${taxid}.reference_selection.tsv
    for id in "\${ids[@]}"; do
        printf '%s\\t%s\\t%s\\n' "\$id" \\
            "\$(nchrom "\${id}.input.fasta.fai")" "\$(n50 "\${id}.input.fasta.fai")" \\
            >> ${taxid}.reference_selection.tsv
    done

    REF=""
    if [ -n "\${PREF}" ]; then
        IFS=',' read -ra PREFA <<< "\${PREF}"
        for p in "\${PREFA[@]}"; do
            for id in "\${ids[@]}"; do
                if [ "\$id" = "\$p" ]; then REF="\$id"; break 2; fi
            done
        done
    fi
    AMB=0
    if [ -z "\${REF}" ]; then
        SEL=\$(tail -n +2 ${taxid}.reference_selection.tsv | awk 'BEGIN{FS=OFS="\\t"}
            { id[NR]=\$1; cnt[NR]=\$2; nfifty[NR]=\$3; n=NR }
            END{
                for(i=1;i<=n;i++) s[i]=cnt[i]
                for(i=2;i<=n;i++){ v=s[i]; j=i-1; while(j>=1 && s[j]>v){ s[j+1]=s[j]; j-- } s[j+1]=v }
                med = (n%2==1) ? s[(n+1)/2] : (s[n/2]+s[n/2+1])/2.0
                bi=0; bd=-1; bn=-1
                for(i=1;i<=n;i++){
                    d=cnt[i]-med; if(d<0) d=-d
                    if(bd<0 || d<bd || (d==bd && nfifty[i]>bn)){ bd=d; bn=nfifty[i]; bi=i }
                }
                # ambiguous if (k>=3 and) another assembly ties at min distance with a
                # different chromosome count -> no majority, N50 tie-break decided it
                amb=0
                for(i=1;i<=n;i++){ d=cnt[i]-med; if(d<0) d=-d; if(n>=3 && d==bd && cnt[i]!=cnt[bi]) amb=1 }
                printf "%s\\n%s\\n", id[bi], amb
            }')
        REF=\$(printf '%s\\n' "\${SEL}" | sed -n '1p')
        AMB=\$(printf '%s\\n' "\${SEL}" | sed -n '2p')
    fi
    echo "[HARMONIZE taxid ${taxid}] reference = \${REF} (see ${taxid}.reference_selection.tsv)"
    if [ "\${AMB}" = "1" ]; then
        echo "[HARMONIZE taxid ${taxid}] WARNING: no clear chromosome-count consensus (even split at the median); reference chosen by N50 tie-break and may sit on a fused or fragmented frame. Inspect ${taxid}.reference_selection.tsv or pin a reference via harmonize_reference_ids." >&2
    fi
    printf '%s\\n' "\${REF}" > ${taxid}.reference_id.txt

    # ------------------------------------------------------------------
    # align each non-reference assembly to the reference (query = 2nd arg,
    # so PAF col1 = query scaffold, col6 = reference chromosome)
    # ------------------------------------------------------------------
    printf 'id\\trole\\tfai\\tpaf\\n' > manifest.tsv
    printf '%s\\tref\\t%s\\t-\\n' "\${REF}" "\${REF}.input.fasta.fai" >> manifest.tsv
    for id in "\${ids[@]}"; do
        [ "\$id" = "\${REF}" ] && continue
        minimap2 -t ${task.cpus} -x ${preset} --eqx --secondary=no \\
            "\${REF}.input.fasta" "\${id}.input.fasta" \\
            2> "\${id}.minimap2.log" > "\${id}.paf"
        printf '%s\\tquery\\t%s\\t%s\\n' "\$id" "\${id}.input.fasta.fai" "\${id}.paf" >> manifest.tsv
    done

    # ------------------------------------------------------------------
    # resolve harmonized names
    # ------------------------------------------------------------------
    python3 ${resolver} \\
        --manifest manifest.tsv \\
        --reference-id "\${REF}" \\
        --species "${taxid}" \\
        --min-scaffold-bp ${min_scaffold} \\
        --chromosome-set-method ${chr_method} \\
        --dropoff-ratio ${drop_ratio} \\
        --dropoff-min-frac ${drop_minfrac} \\
        --min-aligned-frac ${min_frac} \\
        --contained-frac ${contained_frac} ${demote} \\
        --secondary-frac ${sec_frac} \\
        --outdir .

    # ------------------------------------------------------------------
    # versions (awk drains input -> no SIGPIPE under pipefail; no `head`)
    # ------------------------------------------------------------------
    MM=\$(minimap2 --version 2>&1 | awk 'NR==1{print}')
    ST=\$(samtools --version 2>&1 | awk 'NR==1{print \$2}')
    PY=\$(python3 --version 2>&1 | awk '{print \$2}')
    {
        printf 'process\\ttool\\tversion\\n'
        printf '%s\\tminimap2\\t%s\\n' "${task.process}" "\${MM}"
        printf '%s\\tsamtools\\t%s\\n' "${task.process}" "\${ST}"
        printf '%s\\tpython\\t%s\\n'   "${task.process}" "\${PY}"
    } > versions.tsv
    """

    stub:
    """
    ids=(${ids.join(' ')})
    for id in "\${ids[@]}"; do
        printf 'old_name\\tnew_name\\torient\\torder\\tlength\\tclass\\tref_span\\tflags\\n' \\
            > "\${id}.harmonized_name_map.tsv"
    done
    printf 'assembly\\told_name\\tnew_name\\tclass\\tlength\\torient\\tref_span\\tflags\\n' \\
        > "${taxid}.harmonization_report.tsv"
    printf 'id\\tn_chrom_scaffolds\\tscaffold_n50\\n' > "${taxid}.reference_selection.tsv"
    printf 'id\\tn_chrom_scaffolds\\tscaffold_n50\\n' > "${taxid}.reference_selection.tsv"
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
