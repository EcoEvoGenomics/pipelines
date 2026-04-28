include { WRITE_POPULATION_CENSUS_LIST } from "../modules/system.nf"
include { VCFTOOLS_CALCULATE_PAIRWISE_FST } from "../modules/vcftools.nf"
include { GET_WINPCA; WINPCA_CHROM; WINPCA_GENOMEPLOT } from "../modules/winpca.nf"
include { PLINK_INIT_BEDFILES; PLINK_EXCLUDE_CHROMS; PLINK_LD_PRUNE; PLINK_PCA; PLINK_MISSINGNESS; PLINK_TO_VCF } from "../modules/plink.nf"
include { PLINK_EXTRACT_SITES as PLINK_EXTRACT_PRUNED; PLINK_EXTRACT_SITES as PLINK_EXTRACT_AIMS } from "../modules/plink.nf"
include { ADMIXTURE; ADMIXTURE_AIMS } from "../modules/admixture.nf"

nextflow.preview.output = true

workflow {
    main:
    vcf = file(params.ps_vcf)
    pops = Channel.from(params.ps_fst_populations)

    // GENOME-WIDE FST
    pop_lists = WRITE_POPULATION_CENSUS_LIST(pops, params.metadata)

    // Pairwise channel self-comparison without item self-comparison by David Mas-Ponte
    // https://github.com/nextflow-io/nextflow/discussions/2109

    pairwise_pop_lists = pop_lists \
    | combine(pop_lists) \
    | filter { scan -> scan[0] != scan[1] } \
    | map { scan -> scan.sort() } \
    | unique

    VCFTOOLS_CALCULATE_PAIRWISE_FST(vcf, pairwise_pop_lists)

    // CHROMOSOMAL AND GENOME-WIDE WINDOWED PCA
    wpca = GET_WINPCA()

    genome_index = Channel.fromPath(params.ref_genome_index) \
    | splitCsv( sep:"\t") \
    | map { cols -> ref_index: [cols[0], cols[1]] } \
    | filter { index_entry -> ! index_entry[0].toString().contains(params.ref_scaffold_name) }

    comma_separated_chrom_list = genome_index \
    | map { index_entry -> chrom_name: index_entry[0] } \
    | collectFile(sort:true) { chrom_name -> chrom_name + "," } \
    | splitText \
    | map { list -> list.toString().substring(0, list.toString().length() - 2) } // Remove trailing comma

    WINPCA_CHROM(wpca, vcf, params.metadata, genome_index)
    WINPCA_GENOMEPLOT(wpca, vcf, params.metadata, WINPCA_CHROM.out.data.collect(), comma_separated_chrom_list)

    // LD-PRUNING AND WHOLE-GENOME PCA
    PLINK_INIT_BEDFILES(vcf, params.ref_n_chroms)
    PLINK_EXCLUDE_CHROMS(PLINK_INIT_BEDFILES.out, params.ps_exclude_chroms)
    PLINK_LD_PRUNE(PLINK_EXCLUDE_CHROMS.out, params.ps_prune_window_kb, params.ps_prune_step_snps, params.ps_prune_threshold)
    plinkpruned = PLINK_EXTRACT_PRUNED(PLINK_EXCLUDE_CHROMS.out, PLINK_LD_PRUNE.out.prune_in)
    PLINK_PCA(plinkpruned)

    // PERFORM ADMIXTURE ON LD-PRUNED PLINK FILES, FIND AIMs
    k_values = Channel.of(params.ps_admixture_kmin..params.ps_admixture_kmax)
    ADMIXTURE(plinkpruned, k_values)
    ADMIXTURE_AIMS(ADMIXTURE.out.pfile, PLINK_LD_PRUNE.out.prune_in, params.ps_aim_variance_threshold)
    plinkaims = PLINK_EXTRACT_AIMS(plinkpruned, ADMIXTURE_AIMS.out.flatten())
    PLINK_TO_VCF(plinkaims)

    // REPORT MISSINGNESS OF LD-PRUNED DATA
    PLINK_MISSINGNESS(plinkpruned)

    publish:
    pairwise_fst_full = VCFTOOLS_CALCULATE_PAIRWISE_FST.out.full
    pairwise_fst_means = VCFTOOLS_CALCULATE_PAIRWISE_FST.out.means
    winpca_data = WINPCA_CHROM.out.data
    winpca_plot = WINPCA_CHROM.out.plot
    winpca_genomeplot = WINPCA_GENOMEPLOT.out
    pca = PLINK_PCA.out
    admixture = ADMIXTURE.out.concat()
    aims = PLINK_TO_VCF.out
    missingness = PLINK_MISSINGNESS.out
}

output {
    pairwise_fst_full { path "population_structure/pairwise_fst" }
    pairwise_fst_means { path "population_structure/pairwise_fst" }
    winpca_data { path "population_structure/windowed_pca"}
    winpca_plot { path "population_structure/windowed_pca"}
    winpca_genomeplot { path "population_structure/windowed_pca"}
    pca { path "population_structure/ld_pruned/pca" }
    admixture { path "population_structure/ld_pruned/admixture" }
    aims { path "population_structure/ld_pruned/aims" }
    missingness { path "population_structure/ld_pruned/missingness" }
}
