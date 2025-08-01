process REHH_TEST {

    label "REHH"

    script:
    """
    #!/usr/bin/env Rscript
    library("rehh")
    library("vcfR")
    """
}
