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
    frq <- read.table("${frq.toString()}", skip = 1, col.names = c("CHROM", "POS", "N_ALLELES", "N_CHR", "A1", "A2"))
    frq\$MAF <- frq[c("A1", "A2")] |> apply(1, \\(x) min(x))
    idepth <- read.table("${idepth.toString()}", header = TRUE)
    imiss <- read.table("${imiss.toString()}", header = TRUE)
    ldepth_mean <- read.table("${ldepth_mean.toString()}", header = TRUE)
    lqual <- read.table("${lqual.toString()}", header = TRUE)
    lmiss <- read.table("${lmiss.toString()}", header = TRUE)
    het <- read.table("${het.toString()}", header = TRUE)

    maf_density <- data.frame(MAF = density(frq\$MAF)\$x, DENSITY = density(frq\$MAF)\$y)
    ldepth_mean_density <- data.frame(MEAN_DEPTH = density(ldepth_mean\$MEAN_DEPTH)\$x, DENSITY = density(ldepth_mean\$MEAN_DEPTH)\$y)
    lqual_density <- data.frame(QUAL = density(lqual\$QUAL)\$x, DENSITY = density(lqual\$QUAL)\$y)
    lmiss_density <- data.frame(F_MISS = density(lmiss\$F_MISS)\$x, DENSITY = density(lmiss\$F_MISS)\$y)

    pdf("${vcf.simpleName}.stats.pdf", width = 20, height = 20)
    par(mfrow = c(3, 3))
    hist(het\$F, main = "INDV INBREEDING COEFFICIENT (F)")
    hist(idepth\$MEAN_DEPTH, main = "INDV MEAN DEPTH")
    hist(imiss\$F_MISS, main = "INDV MISSINGNESS")
    plot(data = ldepth_mean_density, DENSITY ~ MEAN_DEPTH, main = "SITE MEAN DEPTH", type = "l", col = "red")
    plot(data = lqual_density, DENSITY ~ QUAL, main = "SITE QUALITY", type = "l", col = "red")
    plot(data = lmiss_density, DENSITY ~ F_MISS, main = "SITE MISSINGNESS", type = "l", col = "red")
    plot(data = maf_density, DENSITY ~ MAF, main = "SITE MINOR ALLELE FREQUENCY", type = "l", col = "red")
    dev.off()
    """
}
