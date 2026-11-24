include { FCS_ADAPTOR as FINAL_ADAPTOR } from '../modules/fcs_adaptor.nf'
include { FCS_GX_SCREEN as FINAL_GX } from '../modules/fcs_gx_screen.nf'

// Audit exact finalized records. Never invalidate naming/coordinate evidence by editing here.
workflow FCS_FINAL_SCREEN {
    take:
    assemblies
    database
    main:
    FINAL_ADAPTOR(assemblies.map { m,fa -> tuple(m,fa,params.decon_fcsadaptor_mode,'nextflow','final') },
        file("${moduleDir}/../py_scripts/run_fcs_adaptor.py",checkIfExists:true))
    adaptor = FINAL_ADAPTOR.out.action_report
    versions = FINAL_ADAPTOR.out.versions
    FINAL_GX(assemblies.map { m,fa -> tuple(m,fa,m.taxid ?: params.taxid) }
        .combine(database).map { m,fa,tax,db -> tuple(m,fa,tax,db,'final') })
    // Release each assembly as soon as both of its screens finish, without a cohort barrier.
    ch_audited = assemblies.map { m,fa -> tuple(m.id.toString(),m,fa) }
        .join(FINAL_GX.out.action_report.map { m,p -> tuple(m.id.toString(),true) },failOnDuplicate:true,failOnMismatch:true)
        .join(adaptor.map { m,p -> tuple(m.id.toString(),true) },failOnDuplicate:true,failOnMismatch:true)
        .map { id,m,fa,gx,ad -> tuple(m,fa) }
    ch_all_versions = versions.mix(FINAL_GX.out.versions)
    emit:
    action_report = FINAL_GX.out.action_report
    taxonomy_report = FINAL_GX.out.taxonomy_report
    adaptor_report = adaptor
    assembly = ch_audited
    versions = ch_all_versions
}
