



export APP_ROOT=/home/training/Applications/
export PATH=$PATH:$APP_ROOT/IGVTools
export PICARD_JAR=$APP_ROOT/picard-tools/picard.jar
export SNPEFF_HOME=$APP_ROOT/snpEff/
export GATK_JAR=$APP_ROOT/gatk/GenomeAnalysisTK.jar
export BVATOOLS_JAR=$APP_ROOT/bvatools-1.6/bvatools-1.6-full.jar
export TRIMMOMATIC_JAR=$APP_ROOT/Trimmomatic-0.33/trimmomatic-0.33.jar
export STRELKA_HOME=$APP_ROOT/strelka-1.0.14/
export MUTECT_JAR=$APP_ROOT/mutect-src/mutect-1.1.7.jar
export VARSCAN_JAR=$APP_ROOT/varscan2/VarScan.v2.3.9.jar
export REF=/home/training/ebicancerworkshop201507/reference

cd $HOME/ebicancerworkshop201507


zless -S raw_reads/normal/run62DVGAAXX_1/normal.64.pair1.fastq.gz



zcat raw_reads/normal/run62DVGAAXX_1/normal.64.pair1.fastq.gz | head -n4
zcat raw_reads/normal/run62DVGAAXX_1/normal.64.pair2.fastq.gz | head -n4


zgrep -c "^@HWUSI" raw_reads/normal/run62DVGAAXX_1/normal.64.pair1.fastq.gz


zgrep -c "^@" raw_reads/normal/run62DVGAAXX_1/normal.64.pair1.fastq.gz


# Generate original QC
mkdir originalQC/
java -Xmx1G -jar ${BVATOOLS_JAR} readsqc --quality 64 \
  --read1 raw_reads/normal/run62DVGAAXX_1/normal.64.pair1.fastq.gz \
  --read2 raw_reads/normal/run62DVGAAXX_1/normal.64.pair2.fastq.gz \
  --threads 2 --regionName normalrun62DVGAAXX_1 --output originalQC/


cat adapters.fa


# Trim and convert data
for file in raw_reads/*/run*_?/*.pair1.fastq.gz;
do
  FNAME=`basename $file`;
  DIR=`dirname $file`;
  OUTPUT_DIR=`echo $DIR | sed 's/raw_reads/reads/g'`;

  mkdir -p $OUTPUT_DIR;
  java -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred64 \
    $file \
    ${file%.pair1.fastq.gz}.pair2.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.pair1.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.single1.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.pair2.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.single2.fastq.gz \
    TOPHRED33 ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:30 MINLEN:50 \
    2> ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.trim.out ; 
done

cat reads/normal/run62DVGAAXX_1/normal.trim.out


# Align data
for file in reads/*/run62*_4/*.pair1.fastq.gz;
do
  FNAME=`basename $file`;
  DIR=`dirname $file`;
  OUTPUT_DIR=`echo $DIR | sed 's/reads/alignment/g'`;
  SNAME=`echo $file | sed 's/reads\/\([^/]\+\)\/.*/\1/g'`;
  RUNID=`echo $file | sed 's/.*\/run\([^_]\+\)_.*/\1/g'`;
  LANE=`echo $file | sed 's/.*\/run[^_]\+_\(.\).*/\1/g'`;

  mkdir -p $OUTPUT_DIR;

  bwa mem -M -t 3 \
    -R "@RG\\tID:${SNAME}_${RUNID}_${LANE}\\tSM:${SNAME}\\t\
LB:${SNAME}\\tPU:${RUNID}_${LANE}\\tCN:Centre National de Genotypage\\tPL:ILLUMINA" \
    ${REF}/Homo_sapiens.GRCh37.fa \
    $file \
    ${file%.pair1.fastq.gz}.pair2.fastq.gz \
  | java -Xmx2G -jar ${PICARD_JAR}  SortSam \
    INPUT=/dev/stdin \
    OUTPUT=${OUTPUT_DIR}/${SNAME}.sorted.bam \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000
done


# Merge Data
java -Xmx2G -jar ${PICARD_JAR}  MergeSamFiles \
  INPUT=alignment/normal/run62DPDAAXX_8/normal.sorted.bam \
  INPUT=alignment/normal/run62DVGAAXX_1/normal.sorted.bam \
  INPUT=alignment/normal/run62MK3AAXX_5/normal.sorted.bam \
  INPUT=alignment/normal/runA81DF6ABXX_1/normal.sorted.bam \
  INPUT=alignment/normal/runA81DF6ABXX_2/normal.sorted.bam \
  INPUT=alignment/normal/runBC04D4ACXX_2/normal.sorted.bam \
  INPUT=alignment/normal/runBC04D4ACXX_3/normal.sorted.bam \
  INPUT=alignment/normal/runBD06UFACXX_4/normal.sorted.bam \
  INPUT=alignment/normal/runBD06UFACXX_5/normal.sorted.bam \
  OUTPUT=alignment/normal/normal.sorted.bam \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

java -Xmx2G -jar ${PICARD_JAR}  MergeSamFiles \
  INPUT=alignment/tumor/run62DU0AAXX_8/tumor.sorted.bam \
  INPUT=alignment/tumor/run62DUUAAXX_8/tumor.sorted.bam \
  INPUT=alignment/tumor/run62DVMAAXX_4/tumor.sorted.bam \
  INPUT=alignment/tumor/run62DVMAAXX_6/tumor.sorted.bam \
  INPUT=alignment/tumor/run62DVMAAXX_8/tumor.sorted.bam \
  INPUT=alignment/tumor/run62JREAAXX_4/tumor.sorted.bam \
  INPUT=alignment/tumor/run62JREAAXX_6/tumor.sorted.bam \
  INPUT=alignment/tumor/run62JREAAXX_8/tumor.sorted.bam \
  INPUT=alignment/tumor/runAC0756ACXX_5/tumor.sorted.bam \
  INPUT=alignment/tumor/runBD08K8ACXX_1/tumor.sorted.bam \
  INPUT=alignment/tumor/run62DU6AAXX_8/tumor.sorted.bam \
  INPUT=alignment/tumor/run62DUYAAXX_7/tumor.sorted.bam \
  INPUT=alignment/tumor/run62DVMAAXX_5/tumor.sorted.bam \
  INPUT=alignment/tumor/run62DVMAAXX_7/tumor.sorted.bam \
  INPUT=alignment/tumor/run62JREAAXX_3/tumor.sorted.bam \
  INPUT=alignment/tumor/run62JREAAXX_5/tumor.sorted.bam \
  INPUT=alignment/tumor/run62JREAAXX_7/tumor.sorted.bam \
  INPUT=alignment/tumor/runAC0756ACXX_4/tumor.sorted.bam \
  INPUT=alignment/tumor/runAD08C1ACXX_1/tumor.sorted.bam \
  OUTPUT=alignment/tumor/tumor.sorted.bam \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true



ls -l alignment/normal/
samtools view -H alignment/normal/normal.sorted.bam | grep "^@RG"



samtools view alignment/normal/normal.sorted.bam | head -n4


# Realign
java -Xmx2G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/Homo_sapiens.GRCh37.fa \
  -o alignment/normal/realign.intervals \
  -I alignment/normal/normal.sorted.bam \
  -I alignment/tumor/tumor.sorted.bam \
  -L 9

java -Xmx2G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R ${REF}/Homo_sapiens.GRCh37.fa \
  -targetIntervals alignment/normal/realign.intervals \
  --nWayOut .realigned.bam \
  -I alignment/normal/normal.sorted.bam \
  -I alignment/tumor/tumor.sorted.bam

  mv normal.sorted.realigned.ba* alignment/normal/
  mv tumor.sorted.realigned.ba* alignment/tumor/







# Fix Mate
java -Xmx2G -jar ${PICARD_JAR}  FixMateInformation \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000 \
  INPUT=alignment/normal/normal.sorted.realigned.bam \
  OUTPUT=alignment/normal/normal.matefixed.bam
java -Xmx2G -jar ${PICARD_JAR}  FixMateInformation \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000 \
  INPUT=alignment/tumor/tumor.sorted.realigned.bam \
  OUTPUT=alignment/tumor/tumor.matefixed.bam


# Mark Duplicates
java -Xmx2G -jar ${PICARD_JAR}  MarkDuplicates \
  REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/normal/normal.matefixed.bam \
  OUTPUT=alignment/normal/normal.sorted.dup.bam \
  METRICS_FILE=alignment/normal/normal.sorted.dup.metrics

java -Xmx2G -jar ${PICARD_JAR}  MarkDuplicates \
  REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/tumor/tumor.matefixed.bam \
  OUTPUT=alignment/tumor/tumor.sorted.dup.bam \
  METRICS_FILE=alignment/tumor/tumor.sorted.dup.metrics


less alignment/normal/normal.sorted.dup.metrics


# Recalibrate
for i in normal tumor
do
  java -Xmx2G -jar ${GATK_JAR} \
    -T BaseRecalibrator \
    -nct 2 \
    -R ${REF}/Homo_sapiens.GRCh37.fa \
    -knownSites ${REF}/dbSnp-137_chr9.vcf.gz \
    -L 9:130215000-130636000 \
    -o alignment/${i}/${i}.sorted.dup.recalibration_report.grp \
    -I alignment/${i}/${i}.sorted.dup.bam

    java -Xmx2G -jar ${GATK_JAR} \
      -T PrintReads \
      -nct 2 \
      -R ${REF}/Homo_sapiens.GRCh37.fa \
      -BQSR alignment/${i}/${i}.sorted.dup.recalibration_report.grp \
      -o alignment/${i}/${i}.sorted.dup.recal.bam \
      -I alignment/${i}/${i}.sorted.dup.bam
done






# Extract BAM metrics
Once your whole bam is generated, it's always a good thing to check the data again to see if everything makes sens.

**Compute coverage**
If you have data from a capture kit, you should see how well your targets worked

**Insert Size**
It tells you if your library worked

**Alignment metrics**
It tells you if your sample and you reference fit together

## Compute coverage
Both GATK and BVATools have depth of coverage tools. 

Here we'll use the GATK one

```{.bash}
# Get Depth
for i in normal tumor
do
  java  -Xmx2G -jar ${GATK_JAR} \
    -T DepthOfCoverage \
    --omitDepthOutputAtEachBase \
    --summaryCoverageThreshold 10 \
    --summaryCoverageThreshold 25 \
    --summaryCoverageThreshold 50 \
    --summaryCoverageThreshold 100 \
    --start 1 --stop 500 --nBins 499 -dt NONE \
    -R ${REF}/Homo_sapiens.GRCh37.fa \
    -o alignment/${i}/${i}.sorted.dup.recal.coverage \
    -I alignment/${i}/${i}.sorted.dup.recal.bam \
    -L 9:130215000-130636000 
done


less -S alignment/normal/normal.sorted.dup.recal.coverage.sample_interval_summary
less -S alignment/tumor/tumor.sorted.dup.recal.coverage.sample_interval_summary


# Get insert size
for i in normal tumor
do
  java -Xmx2G -jar ${PICARD_JAR}  CollectInsertSizeMetrics \
    VALIDATION_STRINGENCY=SILENT \
    REFERENCE_SEQUENCE=${REF}/Homo_sapiens.GRCh37.fa \
    INPUT=alignment/${i}/${i}.sorted.dup.recal.bam \
    OUTPUT=alignment/${i}/${i}.sorted.dup.recal.metric.insertSize.tsv \
    HISTOGRAM_FILE=alignment/${i}/${i}.sorted.dup.recal.metric.insertSize.histo.pdf \
    METRIC_ACCUMULATION_LEVEL=LIBRARY
done


less -S alignment/normal/normal.sorted.dup.recal.metric.insertSize.tsv
less -S alignment/tumor/tumor.sorted.dup.recal.metric.insertSize.tsv


# Get alignment metrics
for i in normal tumor
do
  java -Xmx2G -jar ${PICARD_JAR}  CollectAlignmentSummaryMetrics \
    VALIDATION_STRINGENCY=SILENT \
    REFERENCE_SEQUENCE=${REF}/Homo_sapiens.GRCh37.fa \
    INPUT=alignment/${i}/${i}.sorted.dup.recal.bam \
    OUTPUT=alignment/${i}/${i}.sorted.dup.recal.metric.alignment.tsv \
    METRIC_ACCUMULATION_LEVEL=LIBRARY
done


less -S alignment/normal/normal.sorted.dup.recal.metric.alignment.tsv
less -S alignment/tumor/tumor.sorted.dup.recal.metric.alignment.tsv



mkdir pairedVariants


# SAMTools mpileup
for i in normal tumor
do
samtools mpileup -L 1000 -B -q 1 \
  -f ${REF}/Homo_sapiens.GRCh37.fa \
  -r 9:130215000-130636000 \
  alignment/${i}/${i}.sorted.dup.recal.bam \
  > pairedVariants/${i}.mpileup
done

# varscan
java -Xmx2G -jar ${VARSCAN_JAR} somatic pairedVariants/normal.mpileup pairedVariants/tumor.mpileup pairedVariants/varscan --output-vcf 1 --strand-filter 1 --somatic-p-value 0.001 


# Variants MuTecT
# Note MuTecT only works with Java 6, 7 will give you an error
# if you get "Comparison method violates its general contract!
# you used java 7"
java -Xmx2G -jar ${MUTECT_JAR} \
  -T MuTect \
  -R ${REF}/Homo_sapiens.GRCh37.fa \
  -dt NONE -baq OFF --validation_strictness LENIENT -nt 2 \
  --dbsnp ${REF}/dbSnp-137_chr9.vcf.gz \
  --cosmic ${REF}/b37_cosmic_v70_140903.vcf.gz \
  --input_file:normal alignment/normal/normal.sorted.dup.recal.bam \
  --input_file:tumor alignment/tumor/tumor.sorted.dup.recal.bam \
  --out pairedVariants/mutect.call_stats.txt \
  --coverage_file pairedVariants/mutect.wig.txt \
  -pow pairedVariants/mutect.power \
  -vcf pairedVariants/mutect.vcf \
  -L 9:130215000-130636000


# Variants Strelka
cp ${STRELKA_HOME}/etc/strelka_config_bwa_default.ini ./
# Fix ini since we subsampled
sed 's/isSkipDepthFilters =.*/isSkipDepthFilters = 1/g' -i strelka_config_bwa_default.ini

${STRELKA_HOME}/bin/configureStrelkaWorkflow.pl \
  --normal=alignment/normal/normal.sorted.dup.recal.bam \
  --tumor=alignment/tumor/tumor.sorted.dup.recal.bam \
  --ref=${REF}/Homo_sapiens.GRCh37.fa \
  --config=$(pwd)/strelka_config_bwa_default.ini \
  --output-dir=pairedVariants/strelka/

  cd pairedVariants/strelka/
  make -j3
  cd ../..

  cp pairedVariants/strelka/results/passed.somatic.snvs.vcf pairedVariants/strelka.vcf


for i in pairedVariants/*.vcf;do bgzip -c $i > $i.gz ; tabix -p vcf $i.gz;done


zless -S pairedVariants/varscan.snp.vcf.gz


# SnpEff
java  -Xmx6G -jar ${SNPEFF_HOME}/snpEff.jar \
  eff -v -c ${SNPEFF_HOME}/snpEff.config \
  -o vcf \
  -i vcf \
  -stats pairedVariants/varscan.snpeff.vcf.stats.html \
  hg19 \
  pairedVariants/paired_varscan.snp.vcf \
  > pairedVariants/varscan.snpeff.vcf


less -S pairedVariants/mpileup.snpeff.vcf


# Coverage Track
for i in normal tumor
do
  igvtools count \
    -f min,max,mean \
    alignment/${i}/${i}.sorted.dup.recal.bam \
    alignment/${i}/${i}.sorted.dup.recal.bam.tdf \
    b37
done


