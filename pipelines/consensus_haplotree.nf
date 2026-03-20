include { BCFTOOLS_CALL_REGION_VARIANTS; BCFTOOLS_NORMALISE } from "../modules/bcftools.nf"
include { BCFTOOLS_INDEX; BCFTOOLS_INDEX as BCFTOOLS_INDEX_NORMALISED} from "../modules/bcftools.nf"
include { BCFTOOLS_MAKE_CONSENSUS_FASTA } from "../modules/bcftools.nf"
include { SAMTOOLS_EXTRACT_REGION; SAMTOOLS_INDEX_FASTA } from "../modules/samtools.nf"
include { CONCATENATE_FILES as CONCATENATE_FASTAS } from "../modules/system.nf"
include { MAFFT_ALIGN } from "../modules/mafft.nf"
include { IQTREE_BUILD_TREE; IQTREE_TO_PLAIN_NEWICK } from "../modules/iqtree.nf"
include { METADATA_TO_SPART } from "../modules/system.nf"

nextflow.preview.output = true

workflow {
    main:
    crams = Channel.fromPath("${params.ch_cramdir}/**.cram")
    BCFTOOLS_CALL_REGION_VARIANTS(crams, params.metadata, params.ref_genome, params.ref_ploidy, params.ch_genome_region)
    BCFTOOLS_INDEX(BCFTOOLS_CALL_REGION_VARIANTS.out.vcf)
    BCFTOOLS_NORMALISE(BCFTOOLS_INDEX.out.indexed_vcf, params.ref_genome)
    BCFTOOLS_INDEX_NORMALISED(BCFTOOLS_NORMALISE.out.normalised_vcf)
    BCFTOOLS_MAKE_CONSENSUS_FASTA(BCFTOOLS_INDEX_NORMALISED.out.indexed_vcf, params.ch_filt_indelgap, params.ch_filt_inclusions, params.ref_genome)
    SAMTOOLS_INDEX_FASTA(BCFTOOLS_MAKE_CONSENSUS_FASTA.out.fasta)
    SAMTOOLS_EXTRACT_REGION(SAMTOOLS_INDEX_FASTA.out.indexed_fasta, params.ch_genome_region)
    haplotype_fastas = SAMTOOLS_EXTRACT_REGION.out.extracted.collect()

    CONCATENATE_FASTAS(haplotype_fastas, "${params.ch_genome_region}.fasta")
    MAFFT_ALIGN(CONCATENATE_FASTAS.out.concat, params.ch_mafft_iterations)
    IQTREE_BUILD_TREE(MAFFT_ALIGN.out.aligned, params.ch_iqtree_bootstraps)
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
