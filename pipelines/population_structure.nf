include { PLINK_INIT_BEDFILES; PLINK_EXCLUDE_CHROMS; PLINK_LD_PRUNE; PLINK_PCA; PLINK_MISSINGNESS; PLINK_TO_VCF } from "../modules/plink.nf"
include { PLINK_EXTRACT_SITES as PLINK_EXTRACT_PRUNED; PLINK_EXTRACT_SITES as PLINK_EXTRACT_AIMS } from "../modules/plink.nf"
include { PLINK_FILTER; PLINK_FILTER as PLINK_REFILTER } from "../modules/plink.nf"
include { ADMIXTURE; ADMIXTURE_AIMS } from "../modules/admixture.nf"

// Output syntax current as of version 25.04.6
nextflow.preview.output = true

workflow {
    main:
    vcf = file(params.vcf)
    PLINK_INIT_BEDFILES(vcf, params.n_chroms)
    PLINK_EXCLUDE_CHROMS(PLINK_INIT_BEDFILES.out, params.exclude_chroms)
    PLINK_FILTER(PLINK_EXCLUDE_CHROMS.out, params.filter_mind, params.filter_geno, params.filter_maf)
    PLINK_LD_PRUNE(PLINK_FILTER.out, params.prune_window, params.prune_step, params.prune_threshold)
    PLINK_EXTRACT_PRUNED(PLINK_FILTER.out, PLINK_LD_PRUNE.out.prune_in)
    plinkpruned = PLINK_REFILTER(PLINK_EXTRACT_PRUNED.out, params.filter_mind, params.filter_geno, params.filter_maf)

    PLINK_MISSINGNESS(plinkpruned)
    PLINK_PCA(plinkpruned)

    k_values = Channel.of(params.admixture_kmin..params.admixture_kmax)
    ADMIXTURE(plinkpruned, k_values)
    ADMIXTURE_AIMS(ADMIXTURE.out.pfile, PLINK_LD_PRUNE.out.prune_in, params.aim_variance_threshold)
    plinkaims = PLINK_EXTRACT_AIMS(plinkpruned, ADMIXTURE_AIMS.out.flatten())

    PLINK_TO_VCF(plinkaims)

    publish:
    missingness = PLINK_MISSINGNESS.out
    pca = PLINK_PCA.out
    admixture = ADMIXTURE.out.concat()
    aims = PLINK_TO_VCF.out
}

output {
    missingness { path "population_structure/missingness" }
    pca { path "population_structure/pca" }
    admixture { path "population_structure/admixture" }
    aims { path "population_structure/aims" }
}
