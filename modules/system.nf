process CONCATENATE_FILES {

    label "SYSTEM"

    input:
    path(files, stageAs: "inputs/*")
    val(catfile)

    output:
    path("${catfile}"), emit: concat

    script:
    """
    find inputs/ -type f,l | xargs cat > ${catfile}
    """
}
