include { PLINK_INIT_BEDFILES; PLINK_PAIRWISE_LD; PARSE_PLINK_LD_DECAY; PLOT_PLINK_LD_DECAY } from "../modules/plink.nf"

nextflow.preview.output = true

workflow {
    main:
    vcf = file(params.diag_vcf)

    PLINK_INIT_BEDFILES(vcf, params.n_chroms)
    PLINK_PAIRWISE_LD(PLINK_INIT_BEDFILES.out, params.ld_thin, params.ld_window, params.ld_window_kb)
    PARSE_PLINK_LD_DECAY(PLINK_PAIRWISE_LD.out, params.scaffold_name, params.ld_bin_size)
    PLOT_PLINK_LD_DECAY(PARSE_PLINK_LD_DECAY.out, params.ld_window_kb)

    publish:
    ld_stats = PLINK_PAIRWISE_LD.out
    ld_parsed = PARSE_PLINK_LD_DECAY.out
    ld_plot = PLOT_PLINK_LD_DECAY.out
}

output {
    ld_stats { path "vcf_diagnostics/ld_decay" }
    ld_parsed { path "vcf_diagnostics/ld_decay" }
    ld_plot { path "vcf_diagnostics/ld_decay" }
}
