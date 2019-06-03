Here is the command to generate the new quality graphs:

```{.bash}
# Generate post-trimmed QC
mkdir postTrimQC/
java -Xmx1G -jar ${BVATOOLS_JAR} readsqc --quality 33 \
  --read1 reads/normal/runD0YR4ACXX_1/normal.t30l50.pair1.fastq.gz \
  --read2 reads/normal/runD0YR4ACXX_1/normal.t30l50.pair2.fastq.gz \
  --threads 2 --regionName normalD0YR4ACXX_1 --output postTrimQC/
```
