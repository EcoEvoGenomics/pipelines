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
