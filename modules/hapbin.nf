process HAPBIN_HAP_TO_HAPBIN {

    label "HAPBIN"

    input:
    path(hap)
    
    output:
    path("${hap.simpleName}.hapbin")

    script:
    """
    hapbinconv --in ${hap} --out ${hap.simpleName}.hapbin
    """
}

process HAPBIN_RUN_EHHBIN {

    label "HAPBIN"

    input:
    path(hapbin)
    path(map)
    val(locus)
    
    output:
    path("${hapbin.simpleName}.csv")

    script:
    """
    ehhbin --hap ${hapbin} --map ${map} --locus ${locus} --out ${map.simpleName}
    """
}

process HAPBIN_RUN_IHSBIN {

    label "HAPBIN"

    input:
    tuple path(hapbin), path(map)
    
    output:
    path("${hapbin.simpleName}.csv")

    script:
    """
    ihsbin --hap ${hapbin} --map ${map} --out ${map.simpleName}
    """
}

process HAPBIN_RUN_XPEHHBIN {

    label "HAPBIN"

    input:
    path(hapbin_a)
    path(hapbin_b)
    path(map)
    
    output:
    path("${hapbin.simpleName}.csv")

    script:
    """
    xpehhbin --hapA ${hapbin_a} --hapB ${hapbin_b} --map ${map} --out ${map.simpleName}
    """
}
