include { BCFTOOLS_INDEX; BCFTOOLS_SAMPLE_VCF } from "../modules/bcftools.nf"
include { VCFTOOLS_SNP_DENSITY; PLOT_VCFTOOLS_SNP_DENSITY } from "../modules/vcftools.nf"
include { VCFTOOLS_CALCULATE_RELATEDNESS; PLOT_VCFTOOLS_RELATEDNESS } from "../modules/vcftools.nf"
include { VCFTOOLS_VCF_STATS; PLOT_VCFTOOLS_VCF_STATS } from "../modules/vcftools.nf"

nextflow.preview.output = true

workflow {
    main:
    vcf = file(params.vd_vcf)
    
    // Get full, non-downsampled SNP density
    VCFTOOLS_SNP_DENSITY(vcf, params.vd_snpden_binsize, params.ref_scaffold_name)
    PLOT_VCFTOOLS_SNP_DENSITY(VCFTOOLS_SNP_DENSITY.out)

    // For efficiency the remaining stats are run on a downsampled version of the VCF
    vcf_indexed = BCFTOOLS_INDEX(vcf)
    vcf_sampled = BCFTOOLS_SAMPLE_VCF(vcf_indexed, params.vd_n_sampled_sites)

    // Relatedness between individuals
    relatedness = VCFTOOLS_CALCULATE_RELATEDNESS(vcf_sampled)
    PLOT_VCFTOOLS_RELATEDNESS(relatedness)

    // General stats for variants and individuals
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
    snp_density = VCFTOOLS_SNP_DENSITY.out
    snp_density_plot = PLOT_VCFTOOLS_SNP_DENSITY.out
    relatedness = relatedness
    relatedness_plot = PLOT_VCFTOOLS_RELATEDNESS.out
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
    snp_density { path "vcf_diagnostics/snp_density" }
    snp_density_plot { path "vcf_diagnostics/snp_density" }
    relatedness { path "vcf_diagnostics/relatedness" }
    relatedness_plot { path "vcf_diagnostics/relatedness" }
    frq { path "vcf_diagnostics/vcf_stats" }
    idepth { path "vcf_diagnostics/vcf_stats" }
    imiss { path "vcf_diagnostics/vcf_stats" }
    ldepth_mean { path "vcf_diagnostics/vcf_stats" }
    lqual { path "vcf_diagnostics/vcf_stats" }
    lmiss { path "vcf_diagnostics/vcf_stats" }
    het { path "vcf_diagnostics/vcf_stats" }
    vcf_stats_plot { path "vcf_diagnostics/vcf_stats" }
}
