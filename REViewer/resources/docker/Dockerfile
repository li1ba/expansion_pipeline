FROM debian:bullseye-slim as jq-installer

# Install jq 1.6
RUN apt update -y && \
    apt install -y --no-install-recommends \
        build-essential \
        cmake \
        libz-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt update -y && \
    apt install -y --no-install-recommends wget && \
    wget https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 -O /usr/local/bin/jq --no-check-certificate && \
    chmod a+x /usr/local/bin/jq && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


FROM jess/inkscape as inkscape-installer

RUN mkdir -p /inkscape-install/usr/bin && \
    mkdir -p /inkscape-install/usr/lib && \
    mkdir -p /inkscape-install/usr/share && \
    mkdir -p /inkscape-install/usr/lib/x86_64-linux-gnu/gdk-pixbuf-2.0/2.10.0

RUN cp /usr/bin/inkscape /inkscape-install/usr/bin/ && \
    ldd /usr/bin/inkscape | grep "=>" | awk '{print $3}' | xargs -I '{}' cp '{}' /inkscape-install/usr/lib/ && \
    cp -r /usr/share/inkscape /inkscape-install/usr/share/ && \
    cp /usr/lib/x86_64-linux-gnu/gdk-pixbuf-2.0/2.10.0/loaders.cache /inkscape-install/usr/lib/x86_64-linux-gnu/gdk-pixbuf-2.0/2.10.0/loaders.cache


FROM staphb/bcftools:1.12 as bcftools-installer
# Create directories to hold bcftools and its dependencies
RUN mkdir -p /bcftools/usr/local/bin && \
    mkdir -p /bcftools/usr/lib

# Copy bcftools binary & dependencies
RUN cp /usr/local/bin/bcftools /bcftools/usr/local/bin/
RUN ldd /usr/local/bin/bcftools | grep "=>" | awk '{print $3}' | xargs -I '{}' cp '{}' /bcftools/usr/lib/

# Build the final image
FROM debian:bullseye-slim
# Copy jq
COPY --from=jq-installer /usr/local/bin/jq /usr/local/bin/jq
# Copy Inkscape & dependencies 
COPY --from=inkscape-installer /inkscape-install/ /
RUN apt-get update && apt-get install -y \
    fontconfig \
    fonts-dejavu \
    fonts-liberation \
    fonts-freefont-ttf
# Copy bcftools       
COPY --from=bcftools-installer /bcftools /

# Install wget & get binary
RUN apt-get update && apt-get install -y wget gzip
RUN wget https://github.com/Illumina/REViewer/releases/download/v0.2.7/REViewer-v0.2.7-linux_x86_64.gz -O /usr/local/bin/REViewer.gz \
    && gzip -d /usr/local/bin/REViewer.gz \
    && chmod +x /usr/local/bin/REViewer

RUN echo $PATH
#working directory
WORKDIR /

