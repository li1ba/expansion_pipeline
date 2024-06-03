#!/bin/bash
# stranger 0.0.1
main() {

    set -e -x -o pipefail

    echo "Value of input_vcf: '$input_vcf'"
    echo "Value of repeats_file: '$repeats_file'"
    echo "Value of table_fields: '$table_fields'"

    #directory for docker
    mkdir stranger

    #load docker image
    docker load -i /docker/stranger_0.8.0_bkus.tar.gz
    image_id=$(docker images --format '{{.ID}}' | awk 'NR==1')

    #download inputs
    dx download "$input_vcf" -o /home/dnanexus/stranger/input.vcf
    input_name=`dx describe "$input_vcf" --name`
    
    dx download "$repeats_file" -o /home/dnanexus/stranger/repeats.json

    docker run -v /home/dnanexus/stranger:/stranger $image_id stranger \
    /stranger/input.vcf -f /stranger/repeats.json > /home/dnanexus/stranger/${input_name%.vcf*}_stranger.vcf

    output_vcf=$(dx upload /home/dnanexus/stranger/${input_name%.vcf*}_stranger.vcf --brief)

    dx-jobutil-add-output output_vcf "$output_vcf" --class=file

    ###
    #Create a table based on chosen VCF fields
    ###
    #read fields into array
    IFS=' ' read -r -a chosen_fields <<< "$table_fields"

    normal_fields=("CHROM" "POS" "ID" "REF" "ALT" "QUAL" "FILTER" "INFO" "FORMAT")
    info_fields=("END" "REPID" "RL" "RU" "SVTYPE" "VARID" "STR_STATUS" "STR_NORMAL_MAX" "STR_PATHOLOGIC_MIN" "SourceDisplay" "Source" "SourceId" "SweGenMean" "SweGenStd" "DisplayRU" "InheritanceMode" "HGNCId" "RankScore" "Disease")
    format_fields=("ADFL" "ADIR" "ADSP" "GT" "LC" "REPCI" "REPCN" "SO")


    sep="\t"
    bcftools_command=""
    fields1=""
    fields2=""
    header1=""
    header2=""

    for i in "${chosen_fields[@]}";
    do
      if [[ " ${normal_fields[*]} " =~ " ${i} " ]]; then
      fields1=${fields1}%${i}${sep}
      header1=${header1}${i}${sep}
      elif [[ " ${info_fields[*]} " =~ " ${i} " ]]; then
      fields2=${fields2}%${i}${sep}
      header2=${header2}INFO/${i}${sep}
      elif [[ " ${format_fields[*]} " =~ " ${i} " ]]; then
      fields2=${fields2}%${i}${sep}
      header2=${header2}FORMAT/${i}${sep}
      fi
    done


    bcftools_command=""
    fields2=${fields2%%'\t'}

    if [[ -z $fields2 ]]; then
      bcftools_command="${fields1%%'\t'}\n"
    else
      bcftools_command="$fields1[${fields2%%'\t'}]\n"
    fi

    docker run -v /home/dnanexus/:/dnanexus $image_id bcftools query -f $bcftools_command /dnanexus/stranger/${input_name%.vcf*}_stranger.vcf -o /dnanexus/stranger/${input_name%.vcf*}_stranger_nh.tsv
    #add header
    full_header=$header1${header2%%'\t'}
    sed "1i $full_header" /home/dnanexus/stranger/${input_name%.vcf*}_stranger_nh.tsv > /home/dnanexus/stranger/${input_name%.vcf*}_stranger.tsv
    output_vcf_as_table=$(dx upload /home/dnanexus/stranger/${input_name%.vcf*}_stranger.tsv --brief)
    dx-jobutil-add-output output_vcf_as_table "$output_vcf_as_table" --class=file

}
