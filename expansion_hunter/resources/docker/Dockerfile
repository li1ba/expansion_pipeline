FROM staphb/samtools:1.12 as samtools-installer
# Create directories to hold samtools and its dependencies
RUN mkdir -p /samtools/usr/local/bin && \
    mkdir -p /samtools/usr/lib

# Copy samtools binary & dependencies
RUN cp /usr/local/bin/samtools /samtools/usr/local/bin/
RUN ldd /usr/local/bin/samtools | grep "=>" | awk '{print $3}' | xargs -I '{}' cp '{}' /samtools/usr/lib/

# Build the final image
FROM clinicalgenomics/expansionhunter:4.0.2

# Copy samtools       
COPY --from=samtools-installer /samtools /

WORKDIR /data
