
#!/usr/bin/bash

# Function to display usage
usage() {
  echo "Usage: $0 -b <bamname> -f <fastaname> -v <varcat> -s <sex>"
  exit 1
}

# Parse input arguments
while getopts ":b:f:v:s:" opt; do
  case ${opt} in
    b )
      bamname=$OPTARG
      ;;
    f )
      fastaname=$OPTARG
      ;;
    v )
      varcat=$OPTARG
      ;;
    s )
      sex=$OPTARG
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      usage
      ;;
    : )
      echo "Invalid option: -$OPTARG requires an argument" 1>&2
      usage
      ;;
  esac
done
shift $((OPTIND -1))

# Check if all arguments are provided
if [ -z "${bamname}" ] || [ -z "${fastaname}" ] || [ -z "${varcat}" ] || [ -z "${sex}" ]; then
  usage
fi

#Run ExpansionHunter
echo "Running ExpansionHunter.."
docker run -v $PWD:/exphunt li1ba/expansionhunter_bkus:4.0.2 ExpansionHunter \
 --reads /exphunt/$bamname \
 --reference /exphunt/$fastaname \
 --region-extension-length 1000 \
 --variant-catalog /exphunt/$varcat \
 --sex $sex \
 --output-prefix /exphunt/exphunter_${bamname%.*}

#Sort and index ExpansionHunter bamlet
docker run -v $PWD:/exphunt li1ba/expansionhunter_bkus:4.0.2 samtools sort /exphunt/exphunter_${bamname%.*}_realigned.bam -o /exphunt/exphunter_${bamname%.*}_realigned.sorted.bam
docker run -v $PWD:/exphunt li1ba/expansionhunter_bkus:4.0.2 samtools index /exphunt/exphunter_${bamname%.*}_realigned.sorted.bam

vcfname=exphunter_${bamname%.*}.vcf 

#Run Stranger
echo "Running Stranger.."
docker run -v $PWD:/exphunt li1ba/stranger_bkus:0.8.0 stranger /exphunt/exphunter_${bamname%.*}.vcf -f /exphunt/$varcat > "${bamname%.*}_stranger.vcf"

#Run REViewer for each locus
echo "Running REViewer.."

docker run --rm -v $PWD:/exphunt -v /home/dnanexus:/output li1ba/reviewer_bkus:0.2.7 \
sh -c "cat /exphunt/${varcat} | jq '.[] | .LocusId' | tr -d '\"' > /output/loci.txt"


mkdir -p reviewer_output
chmod 777 reviewer_output

while read p;
do
  docker run -v $PWD:/exphunt li1ba/reviewer_bkus:0.2.7 REViewer \
  --reads /exphunt/exphunter_${bamname%.*}_realigned.sorted.bam \
  --output-prefix /exphunt/reviewer_output/reviewer_${bamname%.*} \
  --vcf /exphunt/$vcfname \
  --reference /exphunt/$fastaname \
  --catalog /exphunt/$varcat \
  --locus $p
done < loci.txt

rm loci.txt

echo "Done!"
