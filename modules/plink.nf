process PLINK_INIT_BEDFILES {

    label "PLINK"

    input:
    path(vcf)
    val(n_chroms)

    output:
    tuple(path("${vcf.simpleName}.bed"), path("${vcf.simpleName}.bim"), path("${vcf.simpleName}.fam"), val(n_chroms))

    script:
    """
    plink \
    --vcf ${vcf} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --double-id \
    --set-missing-var-ids @:# \
    --make-bed --out ${vcf.simpleName}
    """
}

process PLINK_FILTER {

    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)
    val(mind)
    val(geno)
    val(maf)

    output:
    tuple path(bed), path(bim), path(fam), val(n_chroms)

    script:
    """
    plink \
    --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --mind ${mind} \
    --geno ${geno} \
    --maf ${maf} \
    --make-bed --out ${bed.simpleName}
    """
}

process PLINK_LD_PRUNE {

    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)
    val(window_size)
    val(step_size)
    val(r2_threshold)

    output:
    path("prune.in"), emit: prune_in
    path("prune.out"), emit: prune_out

    script:
    """
    plink \
    --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --indep-pairwise ${window_size} ${step_size} ${r2_threshold} \
    --make-bed --out ${bed.simpleName}
    mv ${bed.simpleName}.prune.in prune.in
    mv ${bed.simpleName}.prune.out prune.out   
    """
}

process PLINK_EXCLUDE_CHROMS {
    
    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)
    val(exclude_chroms)

    output:
    tuple path(bed), path(bim), path(fam), val(n_chroms)

    script:
    """
    plink \
    --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --not-chr ${exclude_chroms} \
    --make-bed --out ${bed.simpleName}
    """
}

process PLINK_EXTRACT_SITES {
    
    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)
    path(sitelist)

    output:
    tuple \
    path("${bed.simpleName}_${sitelist.simpleName}.bed"), \
    path("${bed.simpleName}_${sitelist.simpleName}.bim"), \
    path("${bed.simpleName}_${sitelist.simpleName}.fam"), \
    val(n_chroms)

    script:
    """
    plink \
    --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --extract ${sitelist} \
    --make-bed --out ${bed.simpleName}_${sitelist.simpleName}
    """
}

process PLINK_MISSINGNESS {

    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)

    output:
    tuple path("${bed.simpleName}.imiss"), path("${bed.simpleName}.lmiss")

    script:
    """
    plink \
    --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --missing

    mv plink.imiss ${bed.simpleName}.imiss
    mv plink.lmiss ${bed.simpleName}.lmiss
    """
}

process PLINK_PCA {

    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)

    output:
    tuple path("${bed.simpleName}.eigenvec"), path("${bed.simpleName}.eigenval")

    script:
    """
    plink \
    --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --pca

    mv plink.eigenvec ${bed.simpleName}.eigenvec
    mv plink.eigenval ${bed.simpleName}.eigenval
    """
}

process PLINK_TO_VCF {

    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)

    output:
    path("${bed.simpleName}.vcf.gz")

    script:
    """
    plink --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --recode vcf-iid bgz --out ${bed.simpleName}
    """
}

process PLINK_PAIRWISE_LD {

    label "PLINK"

    input:
    tuple path(bed), path(bim), path(fam), val(n_chroms)
    val(thin)
    val(ld_window)
    val(ld_window_kb)

    output:
    path("${bed.simpleName}.ld.gz"), emit: ld_stats

    script:
    """
    plink --bfile ${bed.simpleName} \
    --allow-extra-chr --chr-set ${n_chroms} \
    --thin ${thin} \
    --ld-window ${ld_window} \
    --ld-window-kb ${ld_window_kb} \
    --ld-window-r2 0 \
    --r2 gz \
    --out ${bed.simpleName}
    """
}

process PARSE_PLINK_LD_DECAY {

    // Modified from https://github.com/speciationgenomics/scripts/blob/master/ld_decay_calc.py
    // by Mark Ravinet. See https://speciationgenomics.github.io/ld_decay/

    label "NUMPY"

    input:
    path(ld_stats)
    val(scaffold_name)
    val(bin_size)

    output:
    path("${ld_stats.simpleName}.ld*"), emit: ld_decay_bins

    script:
    """
    #!/usr/bin/env python
    import gzip, sys
    import numpy as np

    file = gzip.open("${ld_stats.toString()}", "r")
    file.readline()

    chroms = dict()

    lines_read = 0
    flush_every = 1000000
    for line in file:

        line = line.strip().split()
        chr = line[0]

        if "${scaffold_name}" in str(chr):
            continue

        pos1 = int(line[1])
        pos2 = int(line[4])
        dist = abs(pos1 - pos2)
        ld = float(line[6])

        if not chr in chroms:
            chroms[chr] = dict()
        if not dist in chroms[chr]:
            chroms[chr][dist] = list()
        chroms[chr][dist].append([ld])

        lines_read += 1
        if (lines_read % flush_every) == 0:
            print("Processed %d lines ..." % lines_read)
            sys.stdout.flush()

    sys.stdout.flush()
    file.close()

    outfile = open("${ld_stats.simpleName.toString()}.ldd", "w")
    outfile.write("CHR\\tDIST\\tAVG_R2\\tSTD\\tN_SNP\\n")

    bin_size = int("${bin_size}")
    for chr in chroms:

        dist_bins = dict()

        for dist in sorted(chroms[chr]):
            bin = int(round(dist - (bin_size / 2), -len(str(bin_size)) + 1)) + (bin_size / 2)
            lds = np.array(chroms[chr][dist])
            if not bin in dist_bins:
                dist_bins[bin] = lds
                continue
            dist_bins[bin] = np.concatenate([dist_bins[bin], lds])
            
        for bin in sorted(dist_bins):
            mean = np.mean(dist_bins[bin])
            std = np.std(dist_bins[bin])
            nsnp = np.shape(dist_bins[bin])[0]
            outfile.write("%s\\t%d\\t%g\\t%g\\t%s\\n" % (chr.decode(), bin, mean, std, nsnp))

    outfile.close()
    """
}

process PLOT_PLINK_LD_DECAY {

    label "RBASE"

    input:
    path(ld_decay_bins)
    val(ld_window_kb)

    output:
    path("${ld_decay_bins}.pdf")

    script:
    """
    #!/usr/bin/env Rscript
    tbl <- read.table("${ld_decay_bins.toString()}", header = TRUE)
    pdf("${ld_decay_bins}.pdf")
    plot(data = tbl, AVG_R2 ~ DIST, col = rgb(0, 0, 0, 0.05), pch = 20)
    dev.off()
    """
}
