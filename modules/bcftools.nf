process BCFTOOLS_CALL_REGION_VARIANTS {

    label "BCFTOOLS"

    input:
    path(cram)
    path(metadata)
    path(ref_genome)
    path(ref_ploidy)
    val(region)

    output:
    path("${cram.simpleName}.vcf.gz"), emit: vcf

    script:
    """
    bcftools mpileup \
        --threads ${task.cpus} \
        --targets ${region} \
        --annotate DP \
        --output-type b --output ${cram.simpleName}_tmp.bcf \
        --fasta-ref ${ref_genome} \
        ${cram}

    # MAKE SAMPLES FILE TO MATCH PLOIDY FILE IN BCFTOOLS CALL
    bcftools query --list-samples ${cram.simpleName}_tmp.bcf > samples.txt
    while read -r sample; do
        echo "^\$sample" >> samples_greps.txt
    done < "samples.txt"
    grep -f samples_greps.txt ${metadata} | awk -F, '{print \$1, \$4}' > call.samples

    bcftools call \
        --threads ${task.cpus} \
        --ploidy-file ${ref_ploidy} \
        --samples-file call.samples \
        --targets ${region} \
        --variants-only \
        --multiallelic-caller \
        --output-type z --output ${cram.simpleName}.vcf.gz \
        ${cram.simpleName}_tmp.bcf
    """
}

process BCFTOOLS_INDEX {

    label "BCFTOOLS"

    input:
    path(vcf)

    output:
    tuple path(vcf, includeInputs: true), path("${vcf.name}.csi"), emit: indexed_vcf

    script:
    """
    bcftools index --threads ${task.cpus} ${vcf}
    """
}

process BCFTOOLS_NORMALISE {

    label "BCFTOOLS"

    input:
    tuple path(vcf), path(csi)
    path(ref_genome)

    output:
    path("${vcf.simpleName}_norm.vcf.gz"), emit: normalised_vcf

    script:
    """
    bcftools norm \
        --threads ${task.cpus} \
        --fasta-ref ${ref_genome} \
        --output-type z --output ${vcf.simpleName}_norm.vcf.gz \
        ${vcf}
    """
}

process BCFTOOLS_FILTER {

    label "BCFTOOLS"

    input:
    tuple path(vcf), path(csi)
    val(indelgap)
    val(inclusions)

    output:
    path("${vcf.simpleName}_filt.vcf.gz"), emit: filtered_vcf

    script:
    """
    bcftools filter \
        --threads ${task.cpus} \
        --IndelGap ${indelgap} \
        --include "${inclusions}" \
        --output-type z --output ${vcf.simpleName}_filt.vcf.gz \
        ${vcf}
    """
}

process BCFTOOLS_MAKE_CONSENSUS_FASTA {

    label "BCFTOOLS"

    input:
    tuple path(vcf), path(csi)
    path(ref_genome)

    output:
    tuple env("samples"), path("${vcf.simpleName}.fasta.gz"), emit: fasta

    script:
    """
    samples=\$(bcftools query --list-samples ${vcf})
    cat ${ref_genome} \
    | bcftools consensus --samples \$samples ${vcf} \
    | bgzip -c > ${vcf.simpleName}.fasta.gz
    """
}

process BCFTOOLS_PICK_SAMPLES {

    label "BCFTOOLS"

    input:
    tuple path(sample_list), path(vcf)

    output:
    path("${sample_list.simpleName}_${vcf.simpleName}.vcf.gz"), emit: samples_vcf

    script:
    """
    bcftools view \
        --samples-file ${sample_list} \
        --force-samples \
        --output-type z --output ${sample_list.simpleName}_${vcf.simpleName}.vcf.gz \
        ${vcf}
    """
}
