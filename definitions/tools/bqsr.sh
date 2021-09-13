#!/bin/bash
# mutect 0.0.1

main() {
    set -ex -o pipefail
    echo "Value of bam: '${bam}'"
    echo "Value of bam_index: '${bam_index}'"
    echo "Value of known_sites: '${known_sites[@]}'"
    echo "Value of known_sites: '${known_sites_index[@]}'"
    echo "Value of interval_lists_array: '${interval_lists_array[@]}'"

   
    interval_files=()
    for file in "${interval_lists_array[@]}"; do
        interval_files+=("-iarray_of_scattered_input=${file}")
    done

    map_job=$(dx-jobutil-new-job map \
        "${interval_files[@]}" \
        -i "reference=${reference}" \
        -i "reference_index=${reference_index}" \
        -i "reference_dict=${reference_dict}" \
        -i "bam=${bam}" \
        -i "bam_index=${bam_index}" \
        -i "known_sites:array:file=${known_sites[@]}" \
        -i "known_sites_index:array:file=${known_sites_index[@]}" \
        -i "dockerimage_gatk=${dockerimage_gatk}")

    
    output_name="$bam_prefix.bqsr.bam"
    postprocess_job=$(dx-jobutil-new-job postprocess \
        -i process_outputs:array:jobref=${map_job}:process_outputs \
        -i output_name:string=${output_name} \
        -i bam:file=${bam} \
        -i "bam_index=${bam_index}" \
        -i "reference=${reference}" \
        -i "reference_index=${reference_index}" \
        -i "reference_dict=${reference_dict}" \
        --depends-on $map_job)

    dx-jobutil-add-output vcf --class=jobref "$postprocess_job":vcf 
    dx-jobutil-add-output vcf_index --class=jobref "$postprocess_job":vcf_index
    dx-jobutil-add-output tumor_sample_name --class=string "$tumor_sample_name"
}

map() {
    set -ex -o pipefail
    
    echo "Value of array_of_scattered_input: '${array_of_scattered_input[@]}'" 
    echo "Value of known_sites: '${known_sites[@]}'"
    echo "Value of known_sites_index: '${known_sites_index[@]}'"

    # known_sites_files=()
    # for file in "${known_sites[@]}"; do
    #     known_sites_files+=("-iarray_of_scattered_input=${file}")
    # done

    process_jobs=()
    i=1
    for interval_file in "${array_of_scattered_input[@]}"
    do
        echo "Value of interval_file: '${interval_file}'"
        
        process_jobs+=($(dx-jobutil-new-job process -i "interval_file=${interval_file}" -i "reference=${reference}" -i "reference_index=${reference_index}" -i "reference_dict=${reference_dict}" -i "bam=${bam}" -i "bam_index=${bam_index}" -i shard:string=shard-$i -i "known_sites:array:file=${known_sites[@]}" -i "known_sites_index:array:file=${known_sites_index[@]}" -i "dockerimage_gatk=${dockerimage_gatk}"))
        echo "${process_jobs[@]}"
        i=$((i+1))
        echo "$i.after"
    done

    for process_job in "${process_jobs[@]}"
    do
        dx-jobutil-add-output process_outputs --class=array:jobref "${process_job}:table" --array
    done

}

process() {
    set -ex -o pipefail
    
    dx-download-all-inputs --parallel

    known_sites_files=""
    for file in "${known_sites_path[@]}"; do
        known_sites_filename=$(basename -- "$file")
        mv $file .
        known_sites_files="$known_sites_files --known-sites $known_sites_filename"   
    done

    ## the index all have different paths
    for file in "${known_sites_index_path[@]}"; do
        mv $file .
    done

    mv $bam_index_path ~/in/bam
    mv $reference_index_path ~/in/reference
    mv $reference_dict_path ~/in/reference
    docker load -i ${dockerimage_gatk_path}
  
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -w /home/dnanexus broadinstitute/gatk:4.2.1.0 \
        /bin/bash /usr/local/bin/run_bqsr.sh ${shard}.recal_data.table ${reference_path} ${bam_path} ${interval_file_path} "$known_sites_files"
 

    table=$(dx upload ${shard}.recal_data.table --brief)
    dx-jobutil-add-output table --class=file "$table" 

}

postprocess() {
    set -ex -o pipefail

    #  postprocess_job=$(dx-jobutil-new-job postprocess \
    #     -i process_outputs:array:jobref=${map_job}:process_outputs \
    #     -i output_name:string=${output_name} \
    #     -i bam:file=${bam} \
    #     -i "bam_index=${bam_index}" \
    #     -i "reference=${reference}" \
    #     -i "reference_index=${reference_index}" \
    #     -i "reference_dict=${reference_dict}" \
    #     --depends-on $map_job)

    echo "Value of process_outputs: '${process_outputs[@]}'"
    echo "Value of output_name: '${output_name}'"

    dx download $bam
    dx download $bam_index
    dx download $reference
    dx download $reference_index
    dx download $reference_dict 

    # download vcfs and vcf indices
    mkdir -p /tmp/bqsr_tables
    cd /tmp/bqsr_tables
    for process_output in "${process_outputs[@]}" 
    do
        dx download "$process_output"
        echo "$process_output_path"
    done
    ls -1sh /tmp/bqsr_tables/*.recal_data.table

    inputs=""
    for file in $(ls -1sh /tmp/bqsr_tables/*.recal_data.table); do
        inputs="$inputs -I $file"   
    done

    cd /home/dnanexus

    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -w /home/dnanexus broadinstitute/gatk:4.2.1.0 /usr/local/bin/gather_bqsr.sh "total.recal_data.table" $inputs $bam_name $output_name $reference_name
 
    samtools index ${output_name}

    bam_out=$(dx upload "${output_name}" --brief)
    dx-jobutil-add-output bam_out --class=file "$bam_out" 

    bam_out_index=$(dx upload "${output_name}.bai" --brief)
    dx-jobutil-add-output bam_out_index --class=file "$bam_out_index"

}
