process IQTREE_BUILD_TREE {

    label "IQTREE"

    input:
    path(alignment)
    val(bootstraps)

    output:
    path("${alignment.name}.*"), emit: all_treefiles
    path("${alignment.name}.contree"), emit: contree

    script:
    """
    iqtree -s ${alignment} -bb ${bootstraps} -m MFP -T AUTO --out-format FASTA
    """
}

process IQTREE_TO_PLAIN_NEWICK {

    label "RBASE"

    input:
    path(contree)

    output:
    path("${contree.simpleName}.nwk")

    script:
    """
    #!/usr/bin/env Rscript
    contree <- readLines("${contree.toString()}")
    rm_bootstrap <- gsub(":([0-9.])+", "", contree)
    rm_branchlen <- gsub("\\\\)([0-9.])+", ")", rm_bootstrap)
    writeLines(rm_branchlen, "${contree.simpleName.toString()}.nwk")
    """
}
