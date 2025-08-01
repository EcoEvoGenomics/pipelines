include { BCFTOOLS_CALL_REGION_VARIANTS; BCFTOOLS_NORMALISE; BCFTOOLS_FILTER } from "../modules/bcftools.nf"
include { BCFTOOLS_INDEX; BCFTOOLS_INDEX as BCFTOOLS_INDEX_NORMALISED; BCFTOOLS_INDEX as BCFTOOLS_INDEX_FILTERED } from "../modules/bcftools.nf"
include { BCFTOOLS_MAKE_CONSENSUS_FASTA } from "../modules/bcftools.nf"
include { SAMTOOLS_EXTRACT_REGION; SAMTOOLS_INDEX_FASTA } from "../modules/samtools.nf"
include { CONCATENATE_FILES as CONCATENATE_FASTAS } from "../modules/system.nf"
include { MAFFT_ALIGN } from "../modules/mafft.nf"
include { IQTREE_BUILD_TREE; IQTREE_TO_PLAIN_NEWICK } from "../modules/iqtree.nf"
include { METADATA_TO_SPART } from "../modules/system.nf"

nextflow.preview.output = true

workflow {
    main:
    crams = Channel.fromPath("${params.cramdir}/**.cram")
    BCFTOOLS_CALL_REGION_VARIANTS(crams, params.metadata, params.ref_genome, params.ref_ploidy, params.genome_region)
    BCFTOOLS_INDEX(BCFTOOLS_CALL_REGION_VARIANTS.out.vcf)
    BCFTOOLS_NORMALISE(BCFTOOLS_INDEX.out.indexed_vcf, params.ref_genome)
    BCFTOOLS_INDEX_NORMALISED(BCFTOOLS_NORMALISE.out.normalised_vcf)
    BCFTOOLS_FILTER(BCFTOOLS_INDEX_NORMALISED.out.indexed_vcf, params.filt_indelgap, params.filt_inclusions)
    BCFTOOLS_INDEX_FILTERED(BCFTOOLS_FILTER.out.filtered_vcf)
    BCFTOOLS_MAKE_CONSENSUS_FASTA(BCFTOOLS_INDEX_FILTERED.out.indexed_vcf, params.ref_genome)
    SAMTOOLS_INDEX_FASTA(BCFTOOLS_MAKE_CONSENSUS_FASTA.out.fasta)
    SAMTOOLS_EXTRACT_REGION(SAMTOOLS_INDEX_FASTA.out.indexed_fasta, params.genome_region)
    haplotype_fastas = SAMTOOLS_EXTRACT_REGION.out.extracted.collect()

    CONCATENATE_FASTAS(haplotype_fastas, "${params.genome_region}.fasta")
    MAFFT_ALIGN(CONCATENATE_FASTAS.out.concat, params.mafft_iterations)
    IQTREE_BUILD_TREE(MAFFT_ALIGN.out.aligned, params.iqtree_bootstraps)
    IQTREE_TO_PLAIN_NEWICK(IQTREE_BUILD_TREE.out.contree)

    METADATA_TO_SPART(params.metadata)

    publish:
    fasta = SAMTOOLS_EXTRACT_REGION.out
    iqtree = IQTREE_BUILD_TREE.out.all_treefiles
    alignment = MAFFT_ALIGN.out.aligned
    plain_newick = IQTREE_TO_PLAIN_NEWICK.out
    spart = METADATA_TO_SPART.out
}

output {
    fasta { path "consensus_haplotree/fasta" }
    iqtree { path "consensus_haplotree/iqtree" }
    alignment { path "consensus_haplotree/hapsolutely" }
    plain_newick { path "consensus_haplotree/hapsolutely" }
    spart { path "consensus_haplotree/hapsolutely" }
}
