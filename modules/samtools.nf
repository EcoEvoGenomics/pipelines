process SAMTOOLS_INDEX_FASTA {

    label "SAMTOOLS"

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path(fasta, includeInputs: true), path("${fasta.name}.fai"), emit: indexed_fasta

    script:
    """
    samtools faidx ${fasta}
    """
}

process SAMTOOLS_EXTRACT_REGION {

    label "SAMTOOLS"

    input:
    tuple val(sample), path(fasta), path(fai)
    val(region)

    output:
    path("${sample}_${region}.fasta"), emit: extracted

    script:
    """
    samtools faidx ${fasta} ${region} > ${sample}_${region}.fasta
    sed -i -e 's/>${region}/>${sample}/g' ${sample}_${region}.fasta
    """
}
