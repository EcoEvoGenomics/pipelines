include { GET_GENOMICS_GENERAL; GENOMICS_GENERAL_VCF_TO_GENO; GENOMICS_GENERAL_POPGEN_WINDOWS } from "../modules/genomicsgeneral.nf"

nextflow.preview.output = true

workflow {
    main:
    metadata = file(params.metadata)
    vcfs = Channel.fromPath("${params.pw_vcfdir}/**.vcf.gz")
    gg_repo = GET_GENOMICS_GENERAL()
    geno_file = GENOMICS_GENERAL_VCF_TO_GENO(gg_repo, vcfs)
    GENOMICS_GENERAL_POPGEN_WINDOWS(
        gg_repo,
        geno_file,
        metadata,
        params.pw_sliding_window,
        params.pw_step_size,
        params.pw_min_sites,
        params.pw_input_format,
        params.pw_compare_species
    )

    publish:
    popgen = GENOMICS_GENERAL_POPGEN_WINDOWS.out
}

output {
    popgen { path "popgen_windows" }
}
