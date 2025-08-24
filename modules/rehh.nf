process REHH_LOAD_VCF {

    label "REHH"

    cpus 1
    time { 2.h * task.attempt }
    memory { 256.MB * Math.ceil(vcf.size() / 1024 ** 2) * task.attempt }
    errorStrategy "retry"
    maxRetries 2

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

    cpus {
        def haplohh_size_mb = Math.ceil(haplohh.size() / 1024 ** 2)
        haplohh_size_mb < 15
        ? 8
        : haplohh_size_mb < 30
          ? 16
          : 32
    }
    time {
        def haplohh_size_mb = Math.ceil(haplohh.size() / 1024 ** 2)
        task.exitStatus == 140
        ? task.previousTrace.time * 2
        : haplohh_size_mb < 5
          ? 4.h
          : haplohh_size_mb < 20
            ? 8.h
            : 16.h
    }
    memory {
        def haplohh_size_mb = Math.ceil(haplohh.size() / 1024 ** 2)
        task.exitStatus == 137
        ? 256.MB * haplohh_size_mb * task.attempt
        : 256.MB * haplohh_size_mb
    }
    errorStrategy "retry"
    maxRetries 3

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
        haplohh = readRDS("${haplohh.toString()}"),
        threads = ${task.cpus}
    )

    write.csv(scan, row.names = FALSE, file = "${haplohh.simpleName}.hh.csv")
    """
}

process REHH_CALCULATE_IHS {

    label "REHH"

    cpus 1
    time { 2.h * task.attempt }
    memory { 16.MB * Math.ceil(csv.size() / 1024 ** 2) * task.attempt }
    errorStrategy "retry"
    maxRetries 2

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

    cpus 1
    time { 2.h * task.attempt }
    memory { 16.MB * Math.ceil((csv_a.size() + csv_b.size()) / 1024 ** 2) * task.attempt }
    errorStrategy "retry"
    maxRetries 2

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
