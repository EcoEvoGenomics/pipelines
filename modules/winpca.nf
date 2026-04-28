process GET_WINPCA {

    // Obtains Moritz Blumer's WinPCA software
    // See https://github.com/MoritzBlumer/winpca

    label "SYSTEM"

    output:
    path("winpca-1.2.1/*")

    script:
    """
    curl -L https://github.com/MoritzBlumer/winpca/archive/refs/tags/v1.2.1.tar.gz > winpca.tar.gz
    tar -zxvf winpca.tar.gz
    """
}

process WINPCA_CHROM {

    label "WINPCA"

    // If there is no data for a chromosome, we ignore the resultant error
    errorStrategy "ignore"

    input:
    path(winpca)
    path(vcf)
    path(metadata_csv)
    tuple val(chrom), val(chrom_length)

    output:
    path("${vcf.simpleName}_${chrom}*"), emit: data
    path("*.html"), emit: plot

    script:
    """
    echo -e "sample_id\\tspecies\\tpopulation" > metadata.tsv
    cat ${metadata_csv} | awk -F, '{print \$1"\\t"\$2"\\t"\$3}' >> metadata.tsv
    python winpca pca ${vcf.simpleName}_${chrom} ${vcf} ${chrom}:1-${chrom_length} --np
    python winpca chromplot ${vcf.simpleName}_${chrom} ${chrom}:1-${chrom_length} -m metadata.tsv -g species
    """
}

process WINPCA_GENOMEPLOT {

    label "WINPCA"

    input:
    path(winpca)
    path(vcf)
    path(metadata_csv)
    path(inputs)
    val(chrom_list)

    output:
    path("*.html")

    script:
    """
    echo -e "sample_id\\tspecies\\tpopulation" > metadata.tsv
    cat ${metadata_csv} | awk -F, '{print \$1"\\t"\$2"\\t"\$3}' >> metadata.tsv
    python winpca genomeplot ${vcf.simpleName}_ '${chrom_list}' -m metadata.tsv -g species
    mv ${vcf.simpleName}_.genomeplot.pc_1.html ${vcf.simpleName}.genomeplot.pc_1.html 
    """
}
