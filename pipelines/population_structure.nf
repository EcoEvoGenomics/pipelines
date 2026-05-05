include { WRITE_POPULATION_CENSUS_LIST } from "../modules/system.nf"
include { VCFTOOLS_CALCULATE_PAIRWISE_FST } from "../modules/vcftools.nf"
include { GET_WINPCA; WINPCA_CHROM; WINPCA_GENOMEPLOT } from "../modules/winpca.nf"
include { PLINK_INIT_BEDFILES; PLINK_EXCLUDE_CHROMS; PLINK_LD_PRUNE; PLINK_PCA; PLINK_MISSINGNESS; PLINK_TO_VCF; PLINK_TO_VCF as AIMS_TO_VCF } from "../modules/plink.nf"
include { PLINK_EXTRACT_SITES as PLINK_EXTRACT_PRUNED; PLINK_EXTRACT_SITES as PLINK_EXTRACT_AIMS } from "../modules/plink.nf"
include { ADMIXTURE; ADMIXTURE_AIMS } from "../modules/admixture.nf"

nextflow.preview.output = true

workflow {
    main:
    // WRITE STRING FOR CHROMS TO EXCLUDE WITH PLINK
    comma_separated_chrom_exclude_list = channel.from(params.ps_exclude_chroms) \
    | collectFile(sort:true) { chrom_name -> chrom_name + "," } \
    | splitText \
    | map { list -> list.toString().substring(0, list.toString().length() - 2) } // Remove trailing comma
    
    // PARSE GENOME INDEX WITHOUT EXCLUDED CHROMS AND SCAFFOLDS
    genome_index = Channel.fromPath(params.ref_genome_index) \
    | splitCsv( sep:"\t") \
    | map { cols -> ref_index: [cols[0], cols[1]] } \
    | filter { index_entry -> ! index_entry[0].toString().contains(params.ref_scaffold_name) } \
    | filter { index_entry -> index_entry[0] !in params.ps_exclude_chroms }

    // DETERMINE NUMBER OF CHROMS FOR PLINK
    n_chroms = genome_index.count()

    // GET FILES
    PLINK_INIT_BEDFILES(file(params.ps_vcf), n_chroms)
    bed = PLINK_EXCLUDE_CHROMS(PLINK_INIT_BEDFILES.out, comma_separated_chrom_exclude_list)
    vcf = PLINK_TO_VCF(bed)

    // GENOME-WIDE FST
    pops = Channel.from(params.ps_fst_populations)
    pop_lists = WRITE_POPULATION_CENSUS_LIST(pops, params.metadata)

    // Pairwise channel self-comparison without item self-comparison by David Mas-Ponte
    // https://github.com/nextflow-io/nextflow/discussions/2109

    pairwise_pop_lists = pop_lists \
    | combine(pop_lists) \
    | filter { scan -> scan[0] != scan[1] } \
    | map { scan -> scan.sort() } \
    | unique

    VCFTOOLS_CALCULATE_PAIRWISE_FST(vcf.combine(pairwise_pop_lists))

    // CHROMOSOMAL AND GENOME-WIDE WINDOWED PCA
    wpca = GET_WINPCA()

    comma_separated_chrom_list = genome_index \
    | map { index_entry -> chrom_name: index_entry[0] } \
    | collectFile(sort:true) { chrom_name -> chrom_name + "," } \
    | splitText \
    | map { list -> list.toString().substring(0, list.toString().length() - 2) } // Remove trailing comma

    WINPCA_CHROM(wpca, params.metadata, vcf.combine(genome_index))
    WINPCA_GENOMEPLOT(wpca, vcf, params.metadata, WINPCA_CHROM.out.data.collect(), comma_separated_chrom_list)

    // LD-PRUNING
    PLINK_LD_PRUNE(bed, params.ps_prune_window_kb, params.ps_prune_step_snps, params.ps_prune_threshold)
    plinkpruned = PLINK_EXTRACT_PRUNED(bed, PLINK_LD_PRUNE.out.prune_in)
    
    // WHOLE-GENOME PCA, ADMIXTURE, AIMS
    k_values = Channel.of(params.ps_admixture_kmin..params.ps_admixture_kmax)
    
    ADMIXTURE(plinkpruned, k_values)
    ADMIXTURE_AIMS(ADMIXTURE.out.pfile, PLINK_LD_PRUNE.out.prune_in, params.ps_aim_variance_threshold)
    plinkaims = PLINK_EXTRACT_AIMS(plinkpruned, ADMIXTURE_AIMS.out.flatten())
    AIMS_TO_VCF(plinkaims)

    PLINK_PCA(plinkpruned)

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
    aims = AIMS_TO_VCF.out
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
