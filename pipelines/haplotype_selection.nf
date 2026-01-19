include { WRITE_POPULATION_CENSUS_LIST; JOIN_GROUPED_CSVS } from "../modules/system.nf"
include { BCFTOOLS_PICK_SAMPLES } from "../modules/bcftools.nf"
include { REHH_LOAD_VCF; REHH_SCAN_HAPLOTYPE_HOMOZYGOSITY } from "../modules/rehh.nf"
include { REHH_CALCULATE_IHS; REHH_CALCULATE_XPEHH } from "../modules/rehh.nf"
include { REHH_PARSE_PLOT_SCAN } from "../modules/rehh.nf"

nextflow.preview.output = true

workflow {
    main:
    vcfs = Channel.fromPath("${params.hs_vcfdir}/**.vcf.gz")
    pops = Channel.from(params.populations)
    
    pop_vcfs = WRITE_POPULATION_CENSUS_LIST(pops, params.metadata) \
    | combine(vcfs) \
    | BCFTOOLS_PICK_SAMPLES

    pop_scans = REHH_LOAD_VCF(pop_vcfs, params.is_polarised) \
    | REHH_SCAN_HAPLOTYPE_HOMOZYGOSITY \
    | flatten \
    | map { scan -> tuple(scan.name.tokenize("_").get(0), scan) } \
    | groupTuple(by: 0) \
    | JOIN_GROUPED_CSVS
    REHH_CALCULATE_IHS(
        pop_scans,
        params.is_polarised,
        params.ihs_freqbin,
        params.cand_pval,
        params.cand_window,
        params.cand_overlap,
        params.cand_min_n_mrk,
        params.cand_min_n_extr_mrk,
        params.cand_min_perc_extr_mrk
    )

    // Pairwise channel self-comparison without item self-comparison by David Mas-Ponte
    // https://github.com/nextflow-io/nextflow/discussions/2109

    pairwise_pop_scans = pop_scans \
    | combine(pop_scans) \
    | filter { scan -> scan[0] != scan[1] } \
    | map { scan -> scan.sort() } \
    | unique
    REHH_CALCULATE_XPEHH(
        pairwise_pop_scans,
        params.cand_pval,
        params.cand_window,
        params.cand_overlap,
        params.cand_min_n_mrk,
        params.cand_min_n_extr_mrk,
        params.cand_min_perc_extr_mrk
    )

    xpehh_resultfile_list = REHH_CALCULATE_XPEHH.out.csv \
    | map { xpehh -> "${xpehh.toString()}" } \
    | collectFile(name: "xpehh.list", newLine: true)

    xpehh_candfile_list = REHH_CALCULATE_XPEHH.out.candidates \
    | map { cand -> "${cand.toString()}" } \
    | collectFile(name: "cand.list", newLine: true)

    REHH_PARSE_PLOT_SCAN(
        file("${launchDir}/utils/parse_plot_rehh.R"),
        xpehh_resultfile_list,
        xpehh_candfile_list,
        params.ref_gff,
        params.cand_pval,
        170,
        25
    )

    publish:
    haplohh = REHH_LOAD_VCF.out.rds
    pop_scans = JOIN_GROUPED_CSVS.out.joined
    ihs_csv = REHH_CALCULATE_IHS.out.csv
    ihs_rds = REHH_CALCULATE_IHS.out.rds
    ihs_candidate_regions = REHH_CALCULATE_IHS.out.candidates
    xpehh = REHH_CALCULATE_XPEHH.out.csv
    xpehh_candidate_regions = REHH_CALCULATE_XPEHH.out.candidates
    xpehh_parsed_main = REHH_PARSE_PLOT_SCAN.output.mainplot
    xpehh_parsed_candplots = REHH_PARSE_PLOT_SCAN.output.candplots
    xpehh_parsed_candgenes = REHH_PARSE_PLOT_SCAN.output.candgenes
}

output {
    haplohh { path "haplotype_selection/chr_haplohh"}
    pop_scans { path "haplotype_selection/gw_ihh" }
    ihs_csv { path "haplotype_selection/gw_ihs" }
    ihs_rds { path "haplotype_selection/gw_ihs" }
    ihs_candidate_regions { path "haplotype_selection/gw_ihs" }
    xpehh { path "haplotype_selection/gw_xpehh" }
    xpehh_candidate_regions { path "haplotype_selection/gw_xpehh" }
    xpehh_parsed_main { path "haplotype_selection/gw_xpehh" }
    xpehh_parsed_candplots { path "haplotype_selection/gw_xpehh" }
    xpehh_parsed_candgenes { path "haplotype_selection/gw_xpehh" }
}
