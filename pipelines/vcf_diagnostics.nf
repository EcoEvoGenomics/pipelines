include { PLINK_INIT_BEDFILES; PLINK_PAIRWISE_LD; PARSE_PLINK_LD_DECAY; PLOT_PLINK_LD_DECAY } from "../modules/plink.nf"
include { BCFTOOLS_INDEX; BCFTOOLS_SAMPLE_VCF } from "../modules/bcftools.nf"
include { VCFTOOLS_SNP_DENSITY; PLOT_VCFTOOLS_SNP_DENSITY } from "../modules/vcftools.nf"
include { VCFTOOLS_VCF_STATS; PLOT_VCFTOOLS_VCF_STATS } from "../modules/vcftools.nf"

nextflow.preview.output = true

workflow {
    main:
    vcf = file(params.diag_vcf)
    
    VCFTOOLS_SNP_DENSITY(vcf, params.snpden_binsize, params.scaffold_name)
    PLOT_VCFTOOLS_SNP_DENSITY(VCFTOOLS_SNP_DENSITY.out)

    PLINK_INIT_BEDFILES(vcf, params.n_chroms)
    PLINK_PAIRWISE_LD(PLINK_INIT_BEDFILES.out, params.ld_thin, params.ld_window, params.ld_window_kb)
    PARSE_PLINK_LD_DECAY(PLINK_PAIRWISE_LD.out, params.scaffold_name, params.ld_bin_size)
    PLOT_PLINK_LD_DECAY(PARSE_PLINK_LD_DECAY.out, params.ld_window_kb)

    vcf_indexed = BCFTOOLS_INDEX(vcf)
    vcf_sampled = BCFTOOLS_SAMPLE_VCF(vcf_indexed, params.n_sampled_sites)

    VCFTOOLS_VCF_STATS(vcf_sampled)
    frq = VCFTOOLS_VCF_STATS.out.frq
    idepth = VCFTOOLS_VCF_STATS.out.idepth
    imiss = VCFTOOLS_VCF_STATS.out.imiss
    ldepth_mean = VCFTOOLS_VCF_STATS.out.ldepth_mean
    lqual = VCFTOOLS_VCF_STATS.out.lqual
    lmiss = VCFTOOLS_VCF_STATS.out.lmiss
    het = VCFTOOLS_VCF_STATS.out.het
    PLOT_VCFTOOLS_VCF_STATS(vcf_sampled, frq, idepth, imiss, ldepth_mean, lqual, lmiss, het)

    publish:
    ld_stats = PLINK_PAIRWISE_LD.out
    ld_decay = PARSE_PLINK_LD_DECAY.out
    ld_decay_plot = PLOT_PLINK_LD_DECAY.out
    snp_density = VCFTOOLS_SNP_DENSITY.out
    snp_density_plot = PLOT_VCFTOOLS_SNP_DENSITY.out
    frq = frq
    idepth = idepth
    imiss = imiss
    ldepth_mean = ldepth_mean
    lqual = lqual
    lmiss = lmiss
    het = het
    vcf_stats_plot = PLOT_VCFTOOLS_VCF_STATS.out
}

output {
    ld_stats { path "vcf_diagnostics/ld_decay" }
    ld_decay { path "vcf_diagnostics/ld_decay" }
    ld_decay_plot { path "vcf_diagnostics/ld_decay" }
    snp_density { path "vcf_diagnostics/snp_density" }
    snp_density_plot { path "vcf_diagnostics/snp_density" }
    frq { path "vcf_diagnostics/vcf_stats" }
    idepth { path "vcf_diagnostics/vcf_stats" }
    imiss { path "vcf_diagnostics/vcf_stats" }
    ldepth_mean { path "vcf_diagnostics/vcf_stats" }
    lqual { path "vcf_diagnostics/vcf_stats" }
    lmiss { path "vcf_diagnostics/vcf_stats" }
    het { path "vcf_diagnostics/vcf_stats" }
    vcf_stats_plot { path "vcf_diagnostics/vcf_stats" }
}
