process GET_GENOMICS_GENERAL {

    // Obtains Simon Martin's genomics_general repository
    // See https://simonmartinlab.org/software/
    // Also replace "NaN" with "nan" for more recent versions of Numpy

    label "SYSTEM"

    output:
    path("genomics_general-0.5/*")

    script:
    """
    curl -L https://github.com/simonhmartin/genomics_general/archive/refs/tags/v0.5.tar.gz > genomics_general.tar.gz
    tar -zxvf genomics_general.tar.gz
    sed -i -e 's/np.NaN/np.nan/g' genomics_general-0.5/genomics.py
    """
}

process GENOMICS_GENERAL_VCF_TO_GENO {

    label "NUMPY"

    input:
    path(genomics_general)
    path(vcf)

    output:
    path("${vcf.simpleName}.geno.gz")

    script:
    """
    python VCF_processing/parseVCF.py -i ${vcf} -o ${vcf.simpleName}.geno.gz
    """
}

process GENOMICS_GENERAL_POPGEN_WINDOWS {

    label "NUMPY"

    // The script may fail if a single sample is missing for a (chrom) VCF
    errorStrategy "ignore"

    input:
    path(genomics_general)
    path(geno)
    path(metadata)
    val(window)
    val(min_sites)
    val(format)
    val(compare_species)

    output:
    path("${geno.simpleName}.csv")

    script:
    """
    cat ${metadata} | awk -F, '{print \$1 " " \$${compare_species ? 2 : 3}}' > sample.pops

    echo '-w ${window}' >> popgenWindows.args
    echo '-m ${min_sites}' >> popgenWindows.args
    echo '-g ${geno}' >> popgenWindows.args
    echo '-o ${geno.simpleName}.csv' >> popgenWindows.args
    echo '-f ${format}' >> popgenWindows.args
    echo '-T ${task.cpus}' >> popgenWindows.args
    echo '--popsFile sample.pops' >> popgenWindows.args
    cat sample.pops | awk '{print \$2}' | sed 's/^/-p /' | uniq >> popgenWindows.args

    cat popgenWindows.args | xargs python popgenWindows.py
    """
}
