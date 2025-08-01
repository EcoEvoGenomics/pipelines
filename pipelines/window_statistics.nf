include { PLINK_INIT_BEDFILES; PLINK_PAIRWISE_LD; PARSE_PLINK_LD_DECAY; PLOT_PLINK_LD_DECAY } from "../modules/plink.nf"
include { GET_GENOMICS_GENERAL; GENOMICS_GENERAL_VCF_TO_GENO; GENOMICS_GENERAL_POPGEN_WINDOWS } from "../modules/genomicsgeneral.nf"

nextflow.preview.output = true

workflow {
    main:
    vcfs = Channel.fromPath("${params.vcfdir}/**.vcf.gz")
    PLINK_INIT_BEDFILES(vcfs, params.n_chroms)
    PLINK_PAIRWISE_LD(PLINK_INIT_BEDFILES.out, params.ld_thin, params.ld_window, params.ld_window_kb)
    PARSE_PLINK_LD_DECAY(PLINK_PAIRWISE_LD.out, params.scaffold_name, params.ld_bin_size)
    PLOT_PLINK_LD_DECAY(PARSE_PLINK_LD_DECAY.out, params.ld_window_kb)

    GG = GET_GENOMICS_GENERAL()
    geno = GENOMICS_GENERAL_VCF_TO_GENO(GG, vcfs)
    GENOMICS_GENERAL_POPGEN_WINDOWS(GG, geno, params.metadata, params.gg_sliding_window, params.gg_min_sites, params.gg_format, params.gg_compare_species)

    publish:
    ld_stats = PLINK_PAIRWISE_LD.out
    ld_parsed = PARSE_PLINK_LD_DECAY.out
    ld_plot = PLOT_PLINK_LD_DECAY.out
    popgen = GENOMICS_GENERAL_POPGEN_WINDOWS.out
}

output {
    ld_stats { path "window_statistics/ld_decay" }
    ld_parsed { path "window_statistics/ld_decay" }
    ld_plot { path "window_statistics/ld_decay" }
    popgen { path "window_statistics/popgen_windows" }
}
