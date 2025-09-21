include { GET_GENOMICS_GENERAL; GENOMICS_GENERAL_VCF_TO_GENO; GENOMICS_GENERAL_POPGEN_WINDOWS } from "../modules/genomicsgeneral.nf"

nextflow.preview.output = true

workflow {
    main:
    metadata = file(params.metadata)
    vcfs = Channel.fromPath("${params.pw_vcfdir}/**.vcf.gz")
    GG = GET_GENOMICS_GENERAL()
    geno = GENOMICS_GENERAL_VCF_TO_GENO(GG, vcfs)
    GENOMICS_GENERAL_POPGEN_WINDOWS(GG, geno, metadata, params.gg_sliding_window, params.gg_min_sites, params.gg_format, params.gg_compare_species)

    publish:
    popgen = GENOMICS_GENERAL_POPGEN_WINDOWS.out
}

output {
    popgen { path "popgen_windows" }
}
