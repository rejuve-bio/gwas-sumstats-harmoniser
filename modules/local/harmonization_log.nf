process harmonization_log {
    tag "$GCST"
    
    conda (params.enable_conda ? "${task.ext.conda}" : null)

    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

    input:
    val chr
    tuple val(GCST), val(mode), path(all_hm), path(qc_result), path(delete_sites), path(count), path(raw_yaml), path(input), path(unmapped_files)

    output:
    tuple val(chr), val(GCST), path(raw_yaml), path("${GCST}.h.tsv.gz"), path("${GCST}.h.tsv.gz.tbi"), path ("${GCST}.running.log"), env(result)

    shell:
    """
    # Merge per-chromosome unmapped files into one (keep header from first file)
    files="!{unmapped_files instanceof List ? unmapped_files.join(' ') : unmapped_files}"
    first=\$(echo \$files | tr ' ' '\\n' | head -1)
    head -1 "\$first" > combined_unmapped.tsv
    for f in \$files; do tail -n+2 "\$f"; done >> combined_unmapped.tsv

    # Generating running log
    log_script.sh \
    -r "${params.ref}/homo_sapiens-${chr}.vcf.gz" \
    -i $input \
    -c $count \
    -d $delete_sites \
    -h $all_hm \
    -u combined_unmapped.tsv \
    -o ${GCST}.running.log \
    -p ${params.version}

    N=\$(awk -v RS='\t' '/hm_code/{print NR; exit}' $qc_result)
    sed 1d $qc_result| awk -F "\t" '{print \$'"\$N"'}' | creat_log.py >> ${GCST}.running.log
    
    # extract harmonise result
    result=\$(grep Result ${GCST}.running.log | cut -f2)

    # Prepare the gzip data
    chr=\$(awk -v RS='\t' '/chromosome/{print NR; exit}' $qc_result)
    pos=\$(awk -v RS='\t' '/base_pair_location/{print NR; exit}' $qc_result)

    cat $qc_result | bgzip -c > ${GCST}.h.tsv.gz
    tabix -c N -S 1 -f -s \$chr -b \$pos -e \$pos ${GCST}.h.tsv.gz
    """
}