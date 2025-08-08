process REHH_LOAD_VCF {

    label "REHH"

    input:
    path(vcf)
    val(is_polarised)

    output:
    path("${vcf.simpleName}.haplohh.rds"), emit: rds

    script:
    """
    #!/usr/bin/env Rscript
    print(getwd())
    library("rehh")

    hh <- rehh::data2haplohh(
        hap_file = "${vcf.toString()}",
        polarize_vcf = ifelse(${is_polarised ? 1 : 0}, TRUE, FALSE)
    )

    saveRDS(hh, file = "${vcf.simpleName}.haplohh.rds")
    """
}

process REHH_SCAN_HAPLOTYPE_HOMOZYGOSITY {

    label "REHH"

    input:
    path(haplohh)

    output:
    path("${haplohh.simpleName}.hh.csv"), emit: csv

    script:
    """
    #!/usr/bin/env Rscript
    print(getwd())
    library("rehh")

    scan <- rehh::scan_hh(
        haplohh = readRDS("${haplohh.toString()}")
    )

    write.csv(scan, row.names = FALSE, file = "${haplohh.simpleName}.hh.csv")
    """
}

process REHH_CALCULATE_IHS {

    label "REHH"

    input:
    path(csv)
    val(is_polarised)
    val(ihs_freqbin)
    val(cand_pval)
    val(cand_window)
    val(cand_overlap)

    output:
    path("${csv.simpleName}.ihs.csv"), emit: csv
    path("${csv.simpleName}.ihs.rds"), emit: rds
    path("${csv.simpleName}.cand.csv"), emit: candidates

    script:
    """
    #!/usr/bin/env Rscript
    print(getwd())
    library("rehh")
    
    ihs <- rehh::ihh2ihs(
        scan = read.csv("${csv.toString()}"),
        freqbin = ${is_polarised ? ihs_freqbin : 0}
    )

    cr <- rehh::calc_candidate_regions(
        scan = ihs,
        pval = TRUE,
        threshold = ${cand_pval},
        window_size = ${cand_window},
        overlap = ${cand_overlap}
    )

    write.csv(ihs\$ihs, row.names = FALSE, file = "${csv.simpleName}.ihs.csv")
    write.csv(cr, row.names = FALSE, file = "${csv.simpleName}.cand.csv")
    saveRDS(ihs, file = "${csv.simpleName}.ihs.rds")
    """
}

process REHH_CALCULATE_XPEHH {

    label "REHH"

    input:
    tuple path(csv_a), path(csv_b)

    output:
    path("${csv_a.simpleName}_${csv_b.simpleName}.xpehh.csv"), emit: csv

    script:
    """
    #!/usr/bin/env Rscript
    print(getwd())
    library("rehh")
    
    xpehh <- rehh::ies2xpehh(
        scan_pop1 = read.csv("${csv_a.toString()}"),
        scan_pop2 = read.csv("${csv_b.toString()}"),
        popname1 = "${csv_a.simpleName}",
        popname2 = "${csv_b.simpleName}",
        include_freq = TRUE
    )

    write.csv(xpehh, row.names = FALSE, file = "${csv_a.simpleName}_${csv_b.simpleName}.xpehh.csv")
    """
}
