FROM staphb/bcftools:1.12 as bcftools-installer
# Create directories to hold bcftools and its dependencies
RUN mkdir -p /bcftools/usr/local/bin && \
    mkdir -p /bcftools/usr/lib

# Copy bcftools binary & dependencies
RUN cp /usr/local/bin/bcftools /bcftools/usr/local/bin/
RUN ldd /usr/local/bin/bcftools | grep "=>" | awk '{print $3}' | xargs -I '{}' cp '{}' /bcftools/usr/lib/

# Build the final image
FROM clinicalgenomics/stranger:0.8.0

# Copy bcftools       
COPY --from=bcftools-installer /bcftools /

WORKDIR /data
