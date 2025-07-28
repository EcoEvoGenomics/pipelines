process MAFFT_ALIGN {

    label "MAFFT"

    input:
    path(fasta)
    val(max_iterations)

    output:
    path("${fasta.simpleName}_aligned.fasta"), emit: aligned

    script:
    """
    mafft --maxiterate ${max_iterations} --thread ${task.cpus} ${fasta} > ${fasta.simpleName}_aligned.fasta
    """
}
