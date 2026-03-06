include { BCFTOOLS_INDEX; BCFTOOLS_SAMPLE_VCF } from "../modules/bcftools.nf"
include { VCFTOOLS_SNP_DENSITY; PLOT_VCFTOOLS_SNP_DENSITY } from "../modules/vcftools.nf"
include { VCFTOOLS_CALCULATE_RELATEDNESS; PLOT_VCFTOOLS_RELATEDNESS } from "../modules/vcftools.nf"
include { VCFTOOLS_VCF_STATS; PLOT_VCFTOOLS_VCF_STATS } from "../modules/vcftools.nf"

nextflow.preview.output = true

workflow {
    main:
    vcf = file(params.vd_vcf)
    
    snp_density_whole_vcf = VCFTOOLS_SNP_DENSITY(vcf, params.vd_snpden_binsize, params.ref_scaffold_name)
    PLOT_VCFTOOLS_SNP_DENSITY(snp_density_whole_vcf)

    vcf_indexed = BCFTOOLS_INDEX(vcf)
    vcf_downsampled = BCFTOOLS_SAMPLE_VCF(vcf_indexed, params.vd_n_sampled_sites)

    relatedness = VCFTOOLS_CALCULATE_RELATEDNESS(vcf_downsampled)
    PLOT_VCFTOOLS_RELATEDNESS(relatedness)

    VCFTOOLS_VCF_STATS(vcf_downsampled)
    frq = VCFTOOLS_VCF_STATS.out.frq
    idepth = VCFTOOLS_VCF_STATS.out.idepth
    imiss = VCFTOOLS_VCF_STATS.out.imiss
    ldepth_mean = VCFTOOLS_VCF_STATS.out.ldepth_mean
    lqual = VCFTOOLS_VCF_STATS.out.lqual
    lmiss = VCFTOOLS_VCF_STATS.out.lmiss
    het = VCFTOOLS_VCF_STATS.out.het
    PLOT_VCFTOOLS_VCF_STATS(vcf_downsampled, frq, idepth, imiss, ldepth_mean, lqual, lmiss, het)

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
