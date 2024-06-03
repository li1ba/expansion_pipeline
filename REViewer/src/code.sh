#!/bin/bash
# REViewer 0.0.1
main() {
    set -e -x -o pipefail

    echo "Value of reads: '${reads[@]}'"
    echo "Value of reads_index: '${reads_index[@]}'"
    echo "Value of vcf: '$vcf'"
    echo "Value of reference: '$reference'"
    echo "Value of reference_index: '$reference_index'"
    echo "Value of catalog: '$catalog'"
    echo "Value of locus: '$locus'"
    echo "Value of region-extension-length: '$region_extension_length'"
    echo "Value of pdf: '$pdf'"
    echo "Value of metrics: '$metrics'"

    #create docker environment directory
    mkdir reviewer

    #load docker image
    docker load -i /docker/reviewer_0.2.7_bkus.tar.gz
    image_id=$(docker images --format '{{.ID}}' | awk 'NR==1')
    
    #download inputs
    #bam file
    bamname=`dx describe "${reads[0]}" --name`
    dx download "${reads[0]}" -o /home/dnanexus/$bamname

    #bai index
    bainame=`dx describe "${reads_index[0]}" --name`
    dx download "${reads_index[0]}" -o /home/dnanexus/$bainame

    #vcf
    vcfname=`dx describe "$vcf" --name`
    dx download "$vcf" -o /home/dnanexus/$vcfname

    #fasta file with index
    fastaname=`dx describe "${reference}" --name`
    dx download "$reference" -o /home/dnanexus/$fastaname

    fastaindex=`dx describe "${reference_index}" --name`
    dx download "$reference_index" -o /home/dnanexus/$fastaindex

    #variant catalog
    varcat=`dx describe "${catalog}" --name`
    dx download "$catalog" -o /home/dnanexus/$varcat

    #output dir
    mkdir /home/dnanexus/output

    #get all loci from variant catalog
    if [ -z "$locus" ]
    then
        docker run -v  /home/dnanexus:/docker $image_id cat /docker/${varcat} | jq ' .[] | .LocusId'| tr -d '"'  > /home/dnanexus/loci.txt
        IFS=', ' readarray -t array < "/home/dnanexus/loci.txt"
    else
        IFS=', ' read -r -a array <<< "$locus"
    fi
    
    jq
    #get only genotyped loci
    docker run -e vcfname="$vcfname" \
    -v /home/dnanexus:/docker $image_id sh \
    -c 'bcftools view -e'\''FMT/GT="./."'\'' /docker/$vcfname | bcftools query -f "%INFO/REPID\n" > /docker/genotyped.txt'

    for i in "${array[@]}"
    do
      #check if loci is genotyped
      if grep -Fxq "$i" /home/dnanexus/genotyped.txt
      then
        #run REViewer for each gene locus
        docker run -v /home/dnanexus:/docker $image_id REViewer \
            --reads /docker/$bamname \
            --vcf /docker/$vcfname \
            --reference /docker/$fastaname \
            --catalog /docker/$varcat \
            --locus $i \
            --region-extension-length $region_extension_length \
            --output-prefix /docker/output/reviewer_${bamname%.*}
        fi
    done

    #export to PDF
    if [ "$pdf" == "true" ]
    then
        cd output
        for i in reviewer_${bamname%.*}*.svg; do
            docker run -v /home/dnanexus:/docker $image_id inkscape \
                --file="/docker/output/$i" --without-gui \
                --export-pdf="/docker/output/${i%.*}.pdf"
        done
        FILES="/home/dnanexus/output/*.pdf"
    else
        FILES="/home/dnanexus/output/*.svg"
    fi
    
    #upload outputs
    for i in $FILES;
    do
        out=$(dx upload $i -o ${i##*/} --brief)
        dx-jobutil-add-output output_files "$out" --array
    done

    if [ "$metrics" == "true" ]
    then
        METRICS="/home/dnanexus/output/*.tsv"
        for i in $METRICS;
        do
            out=$(dx upload $i -o ${i##*/} --brief)
            dx-jobutil-add-output metrics_files "$out" --array
        done
    fi

}
