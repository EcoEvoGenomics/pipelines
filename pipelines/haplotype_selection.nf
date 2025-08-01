include { HAPBIN_HAP_TO_HAPBIN; HAPBIN_RUN_IHSBIN } from "../modules/hapbin.nf"
include { VCFTOOLS_VCF_TO_HAP } from "../modules/vcftools.nf"

nextflow.preview.output = true

workflow {
    main:
    vcfs = Channel.fromPath("${params.hs_vcfdir}/**.vcf.gz")
    maps = Channel.fromPath("${params.hs_mapdir}/**.map")
    haps = VCFTOOLS_VCF_TO_HAP(vcfs) | HAPBIN_HAP_TO_HAPBIN
    maps.map { map -> return tuple(map.simpleName(), map) }
    haps.map { map -> return tuple(map.simpleName(), map) }
    haps_with_maps = haps.join(maps)

    HAPBIN_RUN_IHSBIN(haps_with_maps)

    publish:
    hapbins = HAPBIN_HAP_TO_HAPBIN.out
    ihsbin = HAPBIN_RUN_IHSBIN.out
}

output {
    hapbins { path "haplotype_selection/hapbins"}
    ihsbin { path "haplotype_selection/ihsbin"}
}
