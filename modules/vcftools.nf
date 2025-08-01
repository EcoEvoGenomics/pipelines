process VCFTOOLS_VCF_TO_HAP {

    label "VCFTOOLS"

    input:
    path(vcf)

    output:
    path("${vcf.simpleName}.hap")

    script:
    """
    vcftools --vcf ${vcf} --out ${vcf.simpleName} --IMPUTE
    """
}
