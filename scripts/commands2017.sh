export APP_ROOT=/usr/local/bin
export APP_EXT=/home/training/tools
export PICARD_JAR=$APP_ROOT/picard.jar
export SNPEFF_HOME=$APP_EXT/snpEff
export GATK_JAR=$APP_ROOT/GenomeAnalysisTK.jar
export BVATOOLS_JAR=$APP_ROOT/bvatools-1.6-full.jar
export TRIMMOMATIC_JAR=$APP_ROOT/trimmomatic-0.36.jar
export VARSCAN_JAR=$APP_EXT/VarScan.v2.3.9.jar

#conpair setup
export CONPAIR_DIR=$APP_EXT/Conpair/
export CONPAIR_DATA=$CONPAIR_DIR/data
export CONPAIR_SCRIPTS=$CONPAIR_DIR/scripts
export PYTHONPATH=$APP_EXT/Conpair/modules:$PYTHONPATH

#VarDict setup
export VARDICT_HOME=$APP_EXT/VarDictJava-1.5.1
export VARDICT_BIN=$VARDICT_HOME/VarDict


#set-up PATH
export PATH=$CONPAIR_DIR/scripts:$APP_EXT/bwa-0.7.15:$VARDICT_HOME/bin:$PATH


export REF=/home/training/ebicancerworkshop2017/reference


cd $HOME/ebicancerworkshop2017/SNV

#zless -S raw_reads/normal/run62DVGAAXX_1/normal.64.pair1.fastq.gz

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
for file in reads/*/run*/*.pair1.fastq.gz;
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
#java -Xmx2G -jar ${PICARD_JAR}  FixMateInformation \
#  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000 \
#  INPUT=alignment/normal/normal.sorted.realigned.bam \
#  OUTPUT=alignment/normal/normal.matefixed.bam
#java -Xmx2G -jar ${PICARD_JAR}  FixMateInformation \
#  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000 \
#  INPUT=alignment/tumor/tumor.sorted.realigned.bam \
#  OUTPUT=alignment/tumor/tumor.matefixed.bam
  
# Mark Duplicates
java -Xmx2G -jar ${PICARD_JAR}  MarkDuplicates \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/normal/normal.sorted.realigned.bam \
  OUTPUT=alignment/normal/normal.sorted.dup.bam \
  METRICS_FILE=alignment/normal/normal.sorted.dup.metrics

java -Xmx2G -jar ${PICARD_JAR}  MarkDuplicates \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/tumor/tumor.sorted.realigned.bam \
  OUTPUT=alignment/tumor/tumor.sorted.dup.bam \
  METRICS_FILE=alignment/tumor/tumor.sorted.dup.metrics
  
#less alignment/normal/normal.sorted.dup.metrics

# Recalibrate
for i in normal tumor
do
  java -Xmx2G -jar ${GATK_JAR} \
    -T BaseRecalibrator \
    -nct 2 \
    -R ${REF}/Homo_sapiens.GRCh37.fa \
    -knownSites ${REF}/dbSnp-137_chr9.vcf \
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

#pileup for the tumor sample
run_gatk_pileup_for_sample.py \
  -m 6G \
  -G $GATK_JAR \
  -D $CONPAIR_DIR \
  -R ${REF}/Homo_sapiens.GRCh37.fa \
  -B alignment/tumor/tumor.sorted.dup.recal.bam \
  -O alignment/tumor/tumor.sorted.dup.recal.gatkPileup

#pileup for the normal sample
run_gatk_pileup_for_sample.py \
  -m 2G \
  -G $GATK_JAR \
  -D $CONPAIR_DIR \
  -R ${REF}/Homo_sapiens.GRCh37.fa \
  -B alignment/normal/normal.sorted.dup.recal.bam \
  -O alignment/normal/normal.sorted.dup.recal.gatkPileup

#Check concordance
verify_concordance.py -H \
  -M  ${CONPAIR_DATA}/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt \
  -N alignment/normal/normal.sorted.dup.recal.gatkPileup \
  -T alignment/tumor/tumor.sorted.dup.recal.gatkPileup \
  > TumorPair.concordance.tsv 

#Esitmate contamination
estimate_tumor_normal_contamination.py  \
  -M ${CONPAIR_DATA}/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt \
  -N alignment/normal/normal.sorted.dup.recal.gatkPileup \
  -T alignment/tumor/tumor.sorted.dup.recal.gatkPileup \
   > TumorPair.contamination.tsv

#less TumorPair.concordance.tsv
#less TumorPair.contamination.tsv

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

#less -S alignment/normal/normal.sorted.dup.recal.coverage.sample_interval_summary
#less -S alignment/tumor/tumor.sorted.dup.recal.coverage.sample_interval_summary

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

#less -S alignment/normal/normal.sorted.dup.recal.metric.insertSize.tsv
#less -S alignment/tumor/tumor.sorted.dup.recal.metric.insertSize.tsv

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

#less -S alignment/normal/normal.sorted.dup.recal.metric.alignment.tsv
#less -S alignment/tumor/tumor.sorted.dup.recal.metric.alignment.tsv

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

# Filtering
grep "^#\|SS=2" pairedVariants/varscan2.snp.vcf > pairedVariants/varscan2.snp.somatic.vcf

# Variants MuTecT2
java -Xmx2G -jar ${GATK_JAR} \
  -T MuTect2 \
  -R ${REF}/Homo_sapiens.GRCh37.fa \
  -dt NONE -baq OFF --validation_strictness LENIENT \
  --dbsnp ${REF}/dbSnp-137_chr9.vcf \
  --cosmic ${REF}/b37_cosmic_v70_140903.vcf.gz \
  --input_file:normal alignment/normal/normal.sorted.dup.recal.bam \
  --input_file:tumor alignment/tumor/tumor.sorted.dup.recal.bam \
  --out pairedVariants/mutect2.vcf \
  -L 9:130215000-130636000
  
# Filtering
vcftools --vcf pairedVariants/mutect2.vcf --stdout --remove-indels --remove-filtered-all --recode --indv NORMAL --indv TUMOR | awk ' BEGIN {OFS="\t"} {if(substr($0,0,2) != "##") {t=$10; $10=$11; $11=t } ;print } ' >  pairedVariants/mutect2.snp.somatic.vcf
  
# Variants Vardict
java -XX:ParallelGCThreads=1 -Xmx4G -classpath $VARDICT_HOME/lib/VarDict-1.5.1.jar:$VARDICT_HOME/lib/commons-cli-1.2.jar:$VARDICT_HOME/lib/jregex-1.2_01.jar:$VARDICT_HOME/lib/htsjdk-2.8.0.jar com.astrazeneca.vardict.Main   -G ${REF}/Homo_sapiens.GRCh37.fa   -N tumor_pair   -b "alignment/tumor/tumor.sorted.dup.recal.bam|alignment/normal/normal.sorted.dup.recal.bam"  -C -f 0.02 -Q 10 -c 1 -S 2 -E 3 -g 4 -th 3 vardict.bed | $VARDICT_BIN/testsomatic.R   | perl $VARDICT_BIN/var2vcf_paired.pl     -N "TUMOR|NORMAL"     -f 0.02 -P 0.9 -m 4.25 -M  > pairedVariants/vardict.vcf

# Filtering
bcftools view -f PASS  -i 'INFO/STATUS ~ ".*Somatic"' pairedVariants/vardict.vcf | awk ' BEGIN {OFS="\t"} { if(substr($0,0,1) == "#" || length($4) == length($5)) {if(substr($0,0,2) != "##") {t=$10; $10=$11; $11=t} ; print}} ' > pairedVariants/vardict.snp.somatic.vcf

#less pairedVariants/varscan2.snp.somatic.vcf
#less pairedVariants/mutect2.snp.somatic.vcf
#less pairedVariants/vardict.snp.somatic.vcf

# Unified callset
bcbio-variation-recall ensemble \
  --cores 2 --numpass 2 --names mutect2,varscan2,vardict \
  pairedVariants/ensemble.snp.somatic.vcf.gz \
  ${REF}/Homo_sapiens.GRCh37.fa \
  pairedVariants/mutect2.snp.somatic.vcf    \
  pairedVariants/varscan2.snp.somatic.vcf    \
  pairedVariants/vardict.snp.somatic.vcf

#zless pairedVariants/ensemble.snp.somatic.vcf.gz

# SnpEff
java  -Xmx6G -jar ${SNPEFF_HOME}/snpEff.jar \
  eff -v -c ${SNPEFF_HOME}/snpEff.config \
  -o vcf \
  -i vcf \
  -stats pairedVariants/ensemble.snp.somatic.snpeff.stats.html \
  GRCh37.75 \
  pairedVariants/ensemble.snp.somatic.vcf.gz \
  > pairedVariants/ensemble.snp.somatic.snpeff.vcf
#less -S pairedVariants/ensemble.snp.somatic.snpeff.vcf
# Coverage Track
for i in normal tumor
do
  igvtools count \
    -f min,max,mean \
    alignment/${i}/${i}.sorted.dup.recal.bam \
    alignment/${i}/${i}.sorted.dup.recal.bam.tdf \
    b37
done
