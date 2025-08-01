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

process METADATA_TO_SPART {

    // For a description of SPART, see Miralles et al. (2021): https://doi.org/10.1111/1755-0998.13470
    // This format works with Hapsolutely (Vences et al. 2024): https://doi.org/10.1093/bioadv/vbae083

    label "SYSTEM"

    input:
    path(metadata)

    output:
    path("${metadata.simpleName}.spart")

    script:
    """
    # PARSE .SPART HEADER METADATA
    n_samples=\$(wc -l < ${metadata} | bc)
    n_species=\$(awk -F, '{print \$2}' ${metadata} | uniq | wc -l | bc)
    n_populations=\$(awk -F, '{print \$3}' ${metadata} | uniq | wc -l | bc)
    
    # CONSTRUCT .SPART HEADER
    echo "begin spart;" >> ${metadata.simpleName}.spart
    echo "Project_name = ${metadata.simpleName};" >> ${metadata.simpleName}.spart
    echo "Date = \$(date +'%Y-%m-%dT%H:%M:%S');" >> ${metadata.simpleName}.spart
    echo "N_spartitions = 2: Species / Population;" >> ${metadata.simpleName}.spart
    echo "N_individuals = \${n_samples} / \${n_samples};" >> ${metadata.simpleName}.spart
    echo "N_subsets = \${n_species} / \${n_populations};" >> ${metadata.simpleName}.spart

    # CONSTRUCT LOOKUP TABLES FOR INDIVIDUAL ASSIGNMENTS
    cat ${metadata} | sort -u -t, -k2,2 | awk -F, '{print NR "," \$2}' >> spp_lookup.csv.tmp
    cat ${metadata} | sort -u -t, -k3,3 | awk -F, '{print NR "," \$3}' >> pop_lookup.csv.tmp

    # PARSE INDIVIDUAL ASSIGNMENTS
    echo -n "Individual_assignment =" >> ${metadata.simpleName}.spart
    while read -r metadata_line; do
        sample=\$(echo \$metadata_line | awk -F, '{print \$1}')
        spp=\$(cat ${metadata} | grep -e "\${sample}" | awk -F, '{print \$2}')
        pop=\$(cat ${metadata} | grep -e "\${sample}" | awk -F, '{print \$3}')
        spp_idx=\$(cat spp_lookup.csv.tmp | grep -e \${spp} | awk -F, '{print \$1}')
        pop_idx=\$(cat pop_lookup.csv.tmp | grep -e \${pop} | awk -F, '{print \$1}')
        echo -n "\n\${sample}: \${spp_idx} / \${pop_idx}" >> ${metadata.simpleName}.spart
    done < ${metadata}
    echo ";\n" >> ${metadata.simpleName}.spart

    echo "end;" >> ${metadata.simpleName}.spart
    """
}
