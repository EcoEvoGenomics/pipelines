process ADMIXTURE {

    label "ADMIXTURE"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)
    each(K)

    output:
    path("${bed.simpleName}.k${K}.out"), emit: outfile
    path("${bed.simpleName}.k${K}.Q"), emit: qfile
    path("${bed.simpleName}.k${K}.P"), emit: pfile

    script:
    """
    awk '{\$1="0";print \$0}' ${bed.simpleName}.bim > ${bed.simpleName}.bim.tmp
    mv ${bed.simpleName}.bim.tmp ${bed.simpleName}.bim
    admixture --cv -j${task.cpus} ${bed.simpleName}.bed ${K} > ${bed.simpleName}.k${K}.out
    mv ${bed.simpleName}.${K}.Q ${bed.simpleName}.k${K}.Q
    mv ${bed.simpleName}.${K}.P ${bed.simpleName}.k${K}.P
    """
}

process ADMIXTURE_AIMS {

    // Using .P-file between-column variances to find AIMs,
    // inspired by https://doi.org/10.3389/fgene.2019.00043
    
    label "RBASE"

    input:
    path(pfile)
    path(snplist)
    val(variance_threshold)

    output:
    path("*.list"), optional: true

    script:
    """
    #!/usr/bin/env Rscript
    snp_ids <- read.table("${snplist.toString()}")\$V1
    p_table <- read.table("${pfile.toString()}")
    k <- ncol(p_table)
    for (i in seq_len(k)) {
        if (i == k) break
        for (j in seq(i + 1, k)) {
            ij_vars <- apply(p_table[c(i, j)], 1, \\(x) var(x))
            ij_aims <- snp_ids[which(ij_vars > ${variance_threshold})]
            writeLines(ij_aims, paste("aims_k", k, "_p", i, "p", j, ".list", sep = ""))
        }
    }
    """

}
