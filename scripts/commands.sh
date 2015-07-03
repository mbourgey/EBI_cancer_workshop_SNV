



export APP_ROOT=/home/training/Applications/
export PATH=$PATH:$APP_ROOT/IGVTools
export PICARD_HOME=$APP_ROOT/picard-tools-1.115/
export SNPEFF_HOME=$APP_ROOT/snpEff/
export GATK_JAR=$APP_ROOT/gatk/GenomeAnalysisTK.jar
export BVATOOLS_JAR=$APP_ROOT/bvatools-1.3/bvatools-1.3-full.jar
export TRIMMOMATIC_JAR=$APP_ROOT/Trimmomatic-0.32/trimmomatic-0.32.jar
export STRELKA_HOME=$APP_ROOT/strelka-1.0.13/
export MUTECT_JAR=$APP_ROOT/muTect-1.1.4-bin/muTect-1.1.4.jar
export REF=/home/training/Data/DNA_SNV_CNV_bourgey/data/reference

cd $HOME/ebiCancerWorkshop201407


zless -S raw_reads/normal/runD0YR4ACXX_1/normal.64.pair1.fastq.gz



zcat raw_reads/normal/runD0YR4ACXX_1/normal.64.pair1.fastq.gz | head -n4
zcat raw_reads/normal/runD0YR4ACXX_1/normal.64.pair2.fastq.gz | head -n4


zgrep -c "^@HISEQ2" raw_reads/normal/runD0YR4ACXX_1/normal.64.pair1.fastq.gz


zgrep -c "^@" raw_reads/normal/runD0YR4ACXX_1/normal.64.pair1.fastq.gz


# Generate original QC
mkdir originalQC/
java7 -Xmx1G -jar ${BVATOOLS_JAR} readsqc --quality 64 \
  --read1 raw_reads/normal/runD0YR4ACXX_1/normal.64.pair1.fastq.gz \
  --read2 raw_reads/normal/runD0YR4ACXX_1/normal.64.pair2.fastq.gz \
  --threads 2 --regionName normalD0YR4ACXX_1 --output originalQC/


cat adapters.fa


# Trim and convert data
for file in raw_reads/*/run*_?/*.pair1.fastq.gz;
do
  FNAME=`basename $file`;
  DIR=`dirname $file`;
  OUTPUT_DIR=`echo $DIR | sed 's/raw_reads/reads/g'`;

  mkdir -p $OUTPUT_DIR;
  java7 -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred64 \
    $file \
    ${file%.pair1.fastq.gz}.pair2.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.pair1.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.single1.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.pair2.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.single2.fastq.gz \
    TOPHRED33 ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:30 MINLEN:50 \
    2> ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.trim.out ; 
done

cat reads/normal/runD0YR4ACXX_1/normal.trim.out


# Align data
for file in reads/*/run*_?/*.pair1.fastq.gz;
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
    ${REF}/bwa/b37.fasta \
    $file \
    ${file%.pair1.fastq.gz}.pair2.fastq.gz \
  | java7 -Xmx2G -jar ${PICARD_HOME}/SortSam.jar \
    INPUT=/dev/stdin \
    OUTPUT=${OUTPUT_DIR}/${SNAME}.sorted.bam \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000
done


# Merge Data
java7 -Xmx2G -jar ${PICARD_HOME}/MergeSamFiles.jar \
  INPUT=alignment/normal/runC0LWRACXX_1/normal.sorted.bam \
  INPUT=alignment/normal/runC0LWRACXX_6/normal.sorted.bam \
  INPUT=alignment/normal/runC0PTAACXX_6/normal.sorted.bam \
  INPUT=alignment/normal/runC0PTAACXX_7/normal.sorted.bam \
  INPUT=alignment/normal/runC0PTAACXX_8/normal.sorted.bam \
  INPUT=alignment/normal/runC0R2BACXX_6/normal.sorted.bam \
  INPUT=alignment/normal/runC0R2BACXX_7/normal.sorted.bam \
  INPUT=alignment/normal/runC0R2BACXX_8/normal.sorted.bam \
  INPUT=alignment/normal/runD0YR4ACXX_1/normal.sorted.bam \
  INPUT=alignment/normal/runD0YR4ACXX_2/normal.sorted.bam \
  OUTPUT=alignment/normal/normal.sorted.bam \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

java7 -Xmx2G -jar ${PICARD_HOME}/MergeSamFiles.jar \
  INPUT=alignment/tumor/runBC0TV0ACXX_8/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0LVJACXX_6/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0PK4ACXX_7/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0PK4ACXX_8/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0R29ACXX_7/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0R29ACXX_8/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0TTBACXX_3/tumor.sorted.bam \
  INPUT=alignment/tumor/runD114WACXX_8/tumor.sorted.bam \
  OUTPUT=alignment/tumor/tumor.sorted.bam \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true



ls -l alignment/normal/
samtools view -H alignment/normal/normal.sorted.bam | grep "^@RG"



samtools view alignment/normal/normal.sorted.bam | head -n4


# Realign
java7 -Xmx2G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/b37.fasta \
  -o alignment/normal/realign.intervals \
  -I alignment/normal/normal.sorted.bam \
  -I alignment/tumor/tumor.sorted.bam \
  -L 19

java7 -Xmx2G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R ${REF}/b37.fasta \
  -targetIntervals alignment/normal/realign.intervals \
  --nWayOut .realigned.bam \
  -I alignment/normal/normal.sorted.bam \
  -I alignment/tumor/tumor.sorted.bam

  mv normal.sorted.realigned.bam alignment/normal/
  mv tumor.sorted.realigned.bam alignment/tumor/







# Fix Mate
java7 -Xmx2G -jar ${PICARD_HOME}/FixMateInformation.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000 \
  INPUT=alignment/normal/normal.sorted.realigned.bam \
  OUTPUT=alignment/normal/normal.matefixed.bam
java7 -Xmx2G -jar ${PICARD_HOME}/FixMateInformation.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000 \
  INPUT=alignment/tumor/tumor.sorted.realigned.bam \
  OUTPUT=alignment/tumor/tumor.matefixed.bam


# Mark Duplicates
java7 -Xmx2G -jar ${PICARD_HOME}/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/normal/normal.matefixed.bam \
  OUTPUT=alignment/normal/normal.sorted.dup.bam \
  METRICS_FILE=alignment/normal/normal.sorted.dup.metrics

java7 -Xmx2G -jar ${PICARD_HOME}/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/tumor/tumor.matefixed.bam \
  OUTPUT=alignment/tumor/tumor.sorted.dup.bam \
  METRICS_FILE=alignment/tumor/tumor.sorted.dup.metrics


less alignment/normal/normal.sorted.dup.metrics


# Recalibrate
for i in normal tumor
do
  java7 -Xmx2G -jar ${GATK_JAR} \
    -T BaseRecalibrator \
    -nct 2 \
    -R ${REF}/b37.fasta \
    -knownSites ${REF}/dbSnp-137.vcf.gz \
    -L 19:50500000-52502000 \
    -o alignment/${i}/${i}.sorted.dup.recalibration_report.grp \
    -I alignment/${i}/${i}.sorted.dup.bam

    java7 -Xmx2G -jar ${GATK_JAR} \
      -T PrintReads \
      -nct 2 \
      -R ${REF}/b37.fasta \
      -BQSR alignment/${i}/${i}.sorted.dup.recalibration_report.grp \
      -o alignment/${i}/${i}.sorted.dup.recal.bam \
      -I alignment/${i}/${i}.sorted.dup.bam
done




# Check Recalibration
for i in normal tumor
do
  java7 -Xmx2G -jar ${GATK_JAR} \
    -T BaseRecalibrator \
    -nct 2 \
    -R ${REF}/b37.fasta \
    -knownSites ${REF}/dbSnp-137.vcf.gz \
    -L 19:50500000-52502000 \
    -o alignment/${i}/${i}.sorted.dup.recalibration_report.seconnd.grp \
    -I alignment/${i}/${i}.sorted.dup.bam \
    -BQSR alignment/${i}/${i}.sorted.dup.recalibration_report.grp

  java7 -Xmx2G -jar ${GATK_JAR} \
    -T AnalyzeCovariates \
    -R ${REF}/b37.fasta \
    -before alignment/${i}/${i}.sorted.dup.recalibration_report.grp \
    -after alignment/${i}/${i}.sorted.dup.recalibration_report.seconnd.grp \
    -csv alignment/${i}/BQSR.${i}.csv \
    -plots alignment/${i}/BQSR.${i}.pdf
done


# Get Depth
for i in normal tumor
do
  java7  -Xmx2G -jar ${GATK_JAR} \
    -T DepthOfCoverage \
    --omitDepthOutputAtEachBase \
    --summaryCoverageThreshold 10 \
    --summaryCoverageThreshold 25 \
    --summaryCoverageThreshold 50 \
    --summaryCoverageThreshold 100 \
    --start 1 --stop 500 --nBins 499 -dt NONE \
    -R ${REF}/b37.fasta \
    -o alignment/${i}/${i}.sorted.dup.recal.coverage \
    -I alignment/${i}/${i}.sorted.dup.recal.bam \
    -L 19:50500000-52502000 
done


less -S alignment/normal/normal.sorted.dup.recal.coverage.sample_interval_summary
less -S alignment/tumor/tumor.sorted.dup.recal.coverage.sample_interval_summary


# Get insert size
for i in normal tumor
do
  java -Xmx2G -jar ${PICARD_HOME}/CollectInsertSizeMetrics.jar \
    VALIDATION_STRINGENCY=SILENT \
    REFERENCE_SEQUENCE=${REF}/b37.fasta \
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
  java -Xmx2G -jar ${PICARD_HOME}/CollectAlignmentSummaryMetrics.jar \
    VALIDATION_STRINGENCY=SILENT \
    REFERENCE_SEQUENCE=${REF}/b37.fasta \
    INPUT=alignment/${i}/${i}.sorted.dup.recal.bam \
    OUTPUT=alignment/${i}/${i}.sorted.dup.recal.metric.alignment.tsv \
    METRIC_ACCUMULATION_LEVEL=LIBRARY
done


less -S alignment/normal/normal.sorted.dup.recal.metric.alignment.tsv
less -S alignment/tumor/tumor.sorted.dup.recal.metric.alignment.tsv



mkdir pairedVariants


# Variants SAMTools
samtools mpileup -L 1000 -B -q 1 -D -S -g \
  -f ${REF}/b37.fasta \
  -r 19:50500000-52502000 \
  alignment/normal/normal.sorted.dup.recal.bam \
  alignment/tumor/tumor.sorted.dup.recal.bam \
  | bcftools view -vcg -T pair - \
  > pairedVariants/mpileup.vcf


# Variants MuTecT
# Note MuTecT only works with Java 6, 7 will give you an error
# if you get "Comparison method violates its general contract!
# you used java 7"
java -Xmx2G -jar ${MUTECT_JAR} \
  -T MuTect \
  -R ${REF}/b37.fasta \
  -dt NONE -baq OFF --validation_strictness LENIENT -nt 2 \
  --dbsnp ${REF}/dbSnp-137.vcf.gz \
  --cosmic ${REF}/b37_cosmic_v54_120711.vcf \
  --input_file:normal alignment/normal/normal.sorted.dup.recal.bam \
  --input_file:tumor alignment/tumor/tumor.sorted.dup.recal.bam \
  --out pairedVariants/mutect.call_stats.txt \
  --coverage_file pairedVariants/mutect.wig.txt \
  -pow pairedVariants/mutect.power \
  -vcf pairedVariants/mutect.vcf \
  -L 19:50500000-52502000


# Variants Strelka
cp ${STRELKA_HOME}/etc/strelka_config_bwa_default.ini ./
# Fix ini since we subsampled
sed 's/isSkipDepthFilters =.*/isSkipDepthFilters = 1/g' -i strelka_config_bwa_default.ini

${STRELKA_HOME}/bin/configureStrelkaWorkflow.pl \
  --normal=alignment/normal/normal.sorted.dup.recal.bam \
  --tumor=alignment/tumor/tumor.sorted.dup.recal.bam \
  --ref=${REF}/b37.fasta \
  --config=strelka_config_bwa_default.ini \
  --output-dir=pairedVariants/strelka/

  cd pairedVariants/strelka/
  make -j3
  cd $HOME/ebiCancerWorkshop201407

  cp pairedVariants/strelka/results/passed.somatic.snvs.vcf pairedVariants/strelka.vcf


for i in pairedVariants/*.vcf;do bgzip -c $i > $i.gz ; tabix -p vcf $i.gz;done


zless -S variants/mpileup.vcf.gz


# SnpEff
java7  -Xmx6G -jar ${SNPEFF_HOME}/snpEff.jar \
  eff -v -c ${SNPEFF_HOME}/snpEff.config \
  -o vcf \
  -i vcf \
  -stats pairedVariants/mpileup.snpeff.vcf.stats.html \
  hg19 \
  pairedVariants/mpileup.vcf \
  > pairedVariants/mpileup.snpeff.vcf


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


