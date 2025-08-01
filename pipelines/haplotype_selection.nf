include { REHH_TEST } from "../modules/rehh.nf"

nextflow.preview.output = true

workflow {
    main:
    REHH_TEST()
}

output {

}
