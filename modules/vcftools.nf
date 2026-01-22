process VCFTOOLS_SNP_DENSITY {

    label "VCFTOOLS"

    input:
    path(vcf)
    val(binsize)
    val(scaffold_name)

    output:
    path("${vcf.simpleName}.snpden")

    script:
    """
    vcftools --gzvcf ${vcf} --SNPdensity ${binsize} --out ${vcf.simpleName}
    grep -v ${scaffold_name} ${vcf.simpleName}.snpden > ${vcf.simpleName}.snpden.tmp
    mv ${vcf.simpleName}.snpden.tmp ${vcf.simpleName}.snpden
    """
}

process VCFTOOLS_VCF_STATS {

    label "VCFTOOLS"

    input:
    path(vcf)

    output:
    path("${vcf.simpleName}.frq"), emit: frq
    path("${vcf.simpleName}.idepth"), emit: idepth
    path("${vcf.simpleName}.imiss"), emit: imiss
    path("${vcf.simpleName}.ldepth.mean"), emit: ldepth_mean
    path("${vcf.simpleName}.lqual"), emit: lqual
    path("${vcf.simpleName}.lmiss"), emit: lmiss
    path("${vcf.simpleName}.het"), emit: het

    script:
    """
    vcftools --gzvcf ${vcf} --freq2 --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --depth --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --missing-indv --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --site-mean-depth --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --site-quality --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --missing-site --out ${vcf.simpleName}
    vcftools --gzvcf ${vcf} --het --out ${vcf.simpleName}
    """
}

process PLOT_VCFTOOLS_SNP_DENSITY {

    label "RBASE"

    input:
    path(snpden)

    output:
    path("${snpden}.pdf")

    script:
    """
    #!/usr/bin/env Rscript
    tbl <- read.table("${snpden.toString()}", header = TRUE)
    chroms <- unique(tbl["CHROM"]) |> as.vector() |> unlist()
    ymax <- tbl["SNP_COUNT"] |> max()
    plotgrid_ncol <- 6
    plotgrid_nrow <- ceiling(length(chroms) / plotgrid_ncol)
    pdf("${snpden}.pdf", width = 20, height = 20)
    par(mfrow = c(plotgrid_nrow, plotgrid_ncol))
    for (chrom in chroms) {
        plot(data = tbl[which(tbl\$CHROM == chrom), ], SNP_COUNT ~ BIN_START, type = "l", col = "red", ylim = c(0, ymax))
        title(main = chrom)
    }
    dev.off()
    """
}

process PLOT_VCFTOOLS_VCF_STATS {

    label "RBASE"

    input:
    path(vcf)
    path(frq)
    path(idepth)
    path(imiss)
    path(ldepth_mean)
    path(lqual)
    path(lmiss)
    path(het)

    output:
    path("${vcf.simpleName}.stats.pdf")

    script:
    """
    #!/usr/bin/env Rscript

    # Count greatest number of alleles in VCF to format frq table for MAF
    frq <- "${frq.toString()}"
    skip_header_line <- 1L

    field_counts <- count.fields(frq)
    field_counts <- field_counts[-seq_len(skip_header_line)]
    n_nonallele_cols <- 4
    n_alleles <- max(field_counts, na.rm = TRUE) - n_nonallele_cols
    col_names <- c("CHROM","POS","N_ALLELES","N_CHR", paste0("A", seq_len(n_alleles)))

    frq <- read.table(frq, skip = skip_header_line, header = FALSE, fill = TRUE, col.names = col_names, na.strings = c("", "NA"), stringsAsFactors = FALSE)
    
    frq\$MAF <- frq[grep("^A", names(frq), value = TRUE)] |> apply(1, \\(x) min(x))
    idepth <- read.table("${idepth.toString()}", header = TRUE)
    imiss <- read.table("${imiss.toString()}", header = TRUE)
    ldepth_mean <- read.table("${ldepth_mean.toString()}", header = TRUE)
    lqual <- read.table("${lqual.toString()}", header = TRUE)
    lmiss <- read.table("${lmiss.toString()}", header = TRUE)
    het <- read.table("${het.toString()}", header = TRUE)

    pdf("${vcf.simpleName}.stats.pdf", width = 20, height = 20)
    par(mfrow = c(3, 3))
    hist(het\$F, main = "INDV INBREEDING COEFFICIENT (F)")
    hist(idepth\$MEAN_DEPTH, main = "INDV MEAN DEPTH")
    hist(imiss\$F_MISS, main = "INDV MISSINGNESS")
    hist(ldepth_mean\$MEAN_DEPTH[ldepth_mean\$MEAN_DEPTH <= quantile(ldepth_mean\$MEAN_DEPTH, 0.99)], main = "SITE MEAN DEPTH (Cutoff at 99th percentile)")
    hist(lqual\$QUAL[lqual\$QUAL <= quantile(lqual\$QUAL, 0.99)], main = "SITE QUALITY (Cutoff at 99th percentile)")
    hist(lmiss\$F_MISS, main = "SITE MISSINGNESS")
    hist(frq\$MAF, main = "SITE MINOR ALLELE FREQUENCY")
    dev.off()
    """
}
