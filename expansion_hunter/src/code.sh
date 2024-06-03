#!/bin/bash
# expansion_hunter 1.0.0

main() {
    set -e -x -o pipefail

    echo "Value of variant_catalog: '$variant_catalog'"
    echo "Value of reads: '$reads'"
    echo "Value of BAM_file_index: '$BAM_file_index'"
    echo "Value of reference: '$reference'"
    echo "Value of fasta_index: '$fasta_index'"
    echo "Value of sex: '$sex'"
    echo "Value of region_extension_length: '$region_extension_length'"
    echo "Output JSON: '$json'"

    #docker environment directory
    mkdir exphunt

    #load docker image
    docker load -i /docker/expansionhunter_4.0.2_bkus.tar.gz
    image_id=$(docker images --format '{{.ID}}' | awk 'NR==1')

    #download inputs
    #bam file
    bamname=`dx describe "${reads}" --name`
    dx download "$reads" -o /home/dnanexus/exphunt/$bamname

    #bai index
    bainame=`dx describe "${BAM_file_index}" --name`
    dx download "$BAM_file_index" -o /home/dnanexus/exphunt/$bainame

    #fasta file with index
    fastaname=`dx describe "${reference}" --name`
    dx download "$reference" -o /home/dnanexus/exphunt/$fastaname
    if [[ $fastaname == *".gz" ]]; then
      zcat /home/dnanexus/exphunt/$fastaname > /home/dnanexus/exphunt/${fastaname%.gz}
      fastaname=${fastaname%.gz}
    fi

    fastaindex=`dx describe "${fasta_index}" --name`
    dx download "$fasta_index" -o /home/dnanexus/exphunt/$fastaindex

    varcat=`dx describe "${variant_catalog}" --name`
    dx download "$variant_catalog" -o /home/dnanexus/exphunt/$varcat
    
    #run ExpansionHunter
    exphunt_command="docker run -v /home/dnanexus/exphunt:/exphunt $image_id ExpansionHunter
--reads /exphunt/$bamname
--reference /exphunt/$fastaname
--variant-catalog /exphunt/$varcat
--region-extension-length $region_extension_length
--output-prefix /exphunt/exphunter_${bamname%.*} "

    exphunt_flags=""
    if [[ $sex == "male" ]]; then
      exphunt_flags+=" --sex $sex"
    fi

    final_exphunt_command=$($exphunt_command$exphunt_flags)
    echo $final_exphunt_command

    #Sort and index ExpansionHunter bamlet
    docker run -v /home/dnanexus/exphunt:/exphunt $image_id samtools sort /exphunt/exphunter_${bamname%.*}_realigned.bam -o /exphunt/exphunter_${bamname%.*}_realigned.sorted.bam
    docker run -v /home/dnanexus/exphunt:/exphunt $image_id samtools index /exphunt/exphunter_${bamname%.*}_realigned.sorted.bam

    #upload exphunter output
    #.bam file
    bamlet=$(dx upload /home/dnanexus/exphunt/exphunter_${bamname%.*}_realigned.sorted.bam -o exphunter_${bamname%.*}_realigned.sorted.bam  --brief)
    dx-jobutil-add-output bamlet "$bamlet" --class=file
    #.bam file
    bamlet_idx=$(dx upload /home/dnanexus/exphunt/exphunter_${bamname%.*}_realigned.sorted.bam.bai -o exphunter_${bamname%.*}_realigned.sorted.bam.bai  --brief)
    dx-jobutil-add-output bamlet_idx "$bamlet_idx" --class=file
    #.json file
    if [ "$json" == "true" ];
    then
      json_file=$(dx upload /home/dnanexus/exphunt/exphunter_${bamname%.*}.json -o exphunter_${bamname%.*}.json  --brief)
      dx-jobutil-add-output json_file "$json_file" --class=file
    fi
     #.vcf file
    variant_vcf=$(dx upload /home/dnanexus/exphunt/exphunter_${bamname%.*}.vcf -o exphunter_${bamname%.*}.vcf  --brief)
    dx-jobutil-add-output variants_vcf "$variant_vcf" --class=file

}
