##Use this script to pull docker images, add them to app directories and build the workflow###

#Pull & save docker images
docker pull li1ba/expansionhunter_bkus:4.0.2
docker save li1ba/expansionhunter_bkus:4.0.2 | gzip >  expansion_hunter/resources/docker/expansionhunter_4.0.2_bkus.tar.gz

docker pull li1ba/stranger_bkus:0.8.0
docker save li1ba/stranger_bkus:0.8.0 | gzip > stranger/resources/docker/stranger_0.8.0_bkus.tar.gz

docker pull li1ba/reviewer_bkus:0.2.7 
docker save li1ba/reviewer_bkus:0.2.7 | gzip > REViewer/resources/docker/reviewer_0.2.7_bkus.tar.gz

#Build dnanexus applets and save applet IDs
export EH=$(dx build -a expansion_hunter|awk -F'"' '{print $4}')
export stranger=$(dx build -a stranger|awk -F'"' '{print $4}')
export reviewer=$(dx build -a REViewer|awk -F'"' '{print $4}')
export dnanexus_link='$dnanexus_link'

#Create workflow .json file using applet IDs 
(echo "cat <<EOF > dxworkflow.json";
    cat Expansion_detection_pipeline/dxworkflow.json;
    echo "EOF";
    ) > temp.json
. temp.json
cat dxworkflow.json
rm temp.json Expansion_detection_pipeline/dxworkflow.json
mv dxworkflow.json Expansion_detection_pipeline

#Build workflow	
dx build Expansion_detection_pipeline --workflow --keep-open
