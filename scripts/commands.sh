
docker run --privileged -v /tmp:/tmp --network host -it -w $PWD -v $HOME:$HOME --user $UID:$GROUPS -v /etc/group:/etc/group  -v /etc/passwd:/etc/passwd  c3genomics/genpipes:0.8


export REF=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/


cd $HOME/ebicancerworkshop2021/SNV



module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/trimmomatic/0.36 mugqic/samtools/1.9 mugqic/bwa/0.7.17 mugqic/GenomeAnalysisTK/4.1.0.0 mugqic/R_Bioconductor/3.5.0_3.7 mugqic/VarScan/2.4.3 mugqic/vcftools/0.1.14 mugqic/bcftools/1.9 mugqic/VarDictJava/1.4.9 mugqic/bcbio.variation.recall/0.1.7 mugqic/snpEff/4.3 mugqic/igvtools/2.3.67 mugqic/perl/5.22.1




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
    ${REF}/genome/bwa_index/Homo_sapiens.GRCh37.fa \
    $file \
    ${file%.pair1.fastq.gz}.pair2.fastq.gz \
  | java -Xmx2G -jar ${GATK_JAR}  SortSam \
    -I /dev/stdin \
    -O ${OUTPUT_DIR}/${SNAME}.sorted.bam \
    --CREATE_INDEX true --SORT_ORDER coordinate --MAX_RECORDS_IN_RAM 500000
done


# Merge Data
java -Xmx2G -jar ${GATK_JAR}  MergeSamFiles \
  -I alignment/normal/run62DPDAAXX_8/normal.sorted.bam \
  -I alignment/normal/run62DVGAAXX_1/normal.sorted.bam \
  -I alignment/normal/run62MK3AAXX_5/normal.sorted.bam \
  -I alignment/normal/runA81DF6ABXX_1/normal.sorted.bam \
  -I alignment/normal/runA81DF6ABXX_2/normal.sorted.bam \
  -I alignment/normal/runBC04D4ACXX_2/normal.sorted.bam \
  -I alignment/normal/runBC04D4ACXX_3/normal.sorted.bam \
  -I alignment/normal/runBD06UFACXX_4/normal.sorted.bam \
  -I alignment/normal/runBD06UFACXX_5/normal.sorted.bam \
  -O alignment/normal/normal.sorted.bam \
  --CREATE_INDEX true

java -Xmx2G -jar ${GATK_JAR}  MergeSamFiles \
  -I alignment/tumor/run62DU0AAXX_8/tumor.sorted.bam \
  -I alignment/tumor/run62DUUAAXX_8/tumor.sorted.bam \
  -I alignment/tumor/run62DVMAAXX_4/tumor.sorted.bam \
  -I alignment/tumor/run62DVMAAXX_6/tumor.sorted.bam \
  -I alignment/tumor/run62DVMAAXX_8/tumor.sorted.bam \
  -I alignment/tumor/run62JREAAXX_4/tumor.sorted.bam \
  -I alignment/tumor/run62JREAAXX_6/tumor.sorted.bam \
  -I alignment/tumor/run62JREAAXX_8/tumor.sorted.bam \
  -I alignment/tumor/runAC0756ACXX_5/tumor.sorted.bam \
  -I alignment/tumor/runBD08K8ACXX_1/tumor.sorted.bam \
  -I alignment/tumor/run62DU6AAXX_8/tumor.sorted.bam \
  -I alignment/tumor/run62DUYAAXX_7/tumor.sorted.bam \
  -I alignment/tumor/run62DVMAAXX_5/tumor.sorted.bam \
  -I alignment/tumor/run62DVMAAXX_7/tumor.sorted.bam \
  -I alignment/tumor/run62JREAAXX_3/tumor.sorted.bam \
  -I alignment/tumor/run62JREAAXX_5/tumor.sorted.bam \
  -I alignment/tumor/run62JREAAXX_7/tumor.sorted.bam \
  -I alignment/tumor/runAC0756ACXX_4/tumor.sorted.bam \
  -I alignment/tumor/runAD08C1ACXX_1/tumor.sorted.bam \
  -O alignment/tumor/tumor.sorted.bam \
  --CREATE_INDEX true


ls -l alignment/normal/
samtools view -H alignment/normal/normal.sorted.bam | grep "^@RG"


samtools view alignment/normal/normal.sorted.bam | head -n4


# Realign
#switch to old GATK 3.8
module unload  mugqic/GenomeAnalysisTK/4.1.0.0
module load mugqic/GenomeAnalysisTK/3.8

java -Xmx2G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -o alignment/normal/realign.intervals \
  -I alignment/normal/normal.sorted.bam \
  -I alignment/tumor/tumor.sorted.bam \
  -L 9

java -Xmx2G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -targetIntervals alignment/normal/realign.intervals \
  --nWayOut .realigned.bam \
  -I alignment/normal/normal.sorted.bam \
  -I alignment/tumor/tumor.sorted.bam

  mv normal.sorted.realigned.ba* alignment/normal/
  mv tumor.sorted.realigned.ba* alignment/tumor/


#return to GATK 4
module unload mugqic/GenomeAnalysisTK/3.8
module load  mugqic/GenomeAnalysisTK/4.1.0.0
  



# Mark Duplicates
java -Xmx2G -jar ${GATK_JAR}  MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  -I alignment/normal/normal.sorted.realigned.bam \
  -O alignment/normal/normal.sorted.dup.bam \
  --METRICS_FILE alignment/normal/normal.sorted.dup.metrics

java -Xmx2G -jar ${GATK_JAR}  MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX=true \
  -I alignment/tumor/tumor.sorted.realigned.bam \
  -O alignment/tumor/tumor.sorted.dup.bam \
  --METRICS_FILE alignment/tumor/tumor.sorted.dup.metrics


^#less alignment/normal/normal.sorted.dup.metrics


# Recalibrate
for i in normal tumor
do
  java -Xmx2G -jar ${GATK_JAR} BaseRecalibrator \
    -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
    --known-sites ${REF}/annotations/Homo_sapiens.GRCh37.dbSNP150.vcf.gz \
    -L 9:130215000-130636000 \
    -O alignment/${i}/${i}.sorted.dup.recalibration_report.grp \
    -I alignment/${i}/${i}.sorted.dup.bam

    java -Xmx2G -jar ${GATK_JAR} ApplyBQSR \
      -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
      -bqsr alignment/${i}/${i}.sorted.dup.recalibration_report.grp \
      -O alignment/${i}/${i}.sorted.dup.recal.bam \
      -I alignment/${i}/${i}.sorted.dup.bam
done


#Pileup table for the tumor sample
java  -Xmx2G -jar ${GATK_JAR} GetPileupSummaries \
   -I alignment/tumor/tumor.sorted.dup.recal.bam \
   -V $REF/annotations/Homo_sapiens.GRCh37.gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz \
   -L 9:130215000-130636000 \
   -O alignment/tumor/tumor.pileups.table

#Pileup table for the normal sample 
java  -Xmx2G -jar ${GATK_JAR} GetPileupSummaries \
   -I alignment/normal/normal.sorted.dup.recal.bam \
   -V $REF/annotations/Homo_sapiens.GRCh37.gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz \
   -L 9:130215000-130636000 \
   -O alignment/normal/normal.pileups.table

#Esitmate contamination
java  -Xmx2G -jar ${GATK_JAR} CalculateContamination \
   -I alignment/tumor/tumor.pileups.table \
   -matched alignment/normal/normal.pileups.table \
   -O contamination.table


^#less TumorPair.concordance.tsv
^#less TumorPair.contamination.tsv


# Get Depth
#switch to old GATK 3.8
module unload  mugqic/GenomeAnalysisTK/4.1.0.0
module load mugqic/GenomeAnalysisTK/3.8

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
    -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
    -o alignment/${i}/${i}.sorted.dup.recal.coverage \
    -I alignment/${i}/${i}.sorted.dup.recal.bam \
    -L 9:130215000-130636000 
done

#return to GATK 4
module unload mugqic/GenomeAnalysisTK/3.8
module load  mugqic/GenomeAnalysisTK/4.1.0.0

^#less -S alignment/normal/normal.sorted.dup.recal.coverage.sample_interval_summary
^#less -S alignment/tumor/tumor.sorted.dup.recal.coverage.sample_interval_summary


# Get insert size
for i in normal tumor
do
  java -Xmx2G -jar ${GATK_JAR}  CollectInsertSizeMetrics \
    -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
    -I alignment/${i}/${i}.sorted.dup.recal.bam \
    -O alignment/${i}/${i}.sorted.dup.recal.metric.insertSize.tsv \
    -H alignment/${i}/${i}.sorted.dup.recal.metric.insertSize.histo.pdf \
    --METRIC_ACCUMULATION_LEVEL LIBRARY
done


^#less -S alignment/normal/normal.sorted.dup.recal.metric.insertSize.tsv
^#less -S alignment/tumor/tumor.sorted.dup.recal.metric.insertSize.tsv


# Get alignment metrics
for i in normal tumor
do
  java -Xmx2G -jar ${GATK_JAR}  CollectAlignmentSummaryMetrics \
    -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
    -I alignment/${i}/${i}.sorted.dup.recal.bam \
    -O alignment/${i}/${i}.sorted.dup.recal.metric.alignment.tsv \
    --METRIC_ACCUMULATION_LEVEL LIBRARY
done


^#less -S alignment/normal/normal.sorted.dup.recal.metric.alignment.tsv
^#less -S alignment/tumor/tumor.sorted.dup.recal.metric.alignment.tsv


mkdir pairedVariants


# SAMTools mpileup
for i in normal tumor
do
samtools mpileup -B -q 1 \
  -f ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -r 9:130215000-130636000 \
  alignment/${i}/${i}.sorted.dup.recal.bam \
  > pairedVariants/${i}.mpileup
done


# varscan
java -Xmx2G -jar ${VARSCAN2_JAR} somatic \
   pairedVariants/normal.mpileup \
   pairedVariants/tumor.mpileup \
   pairedVariants/varscan2 \
   --output-vcf 1 \
   --strand-filter 1 \
   --somatic-p-value 0.001 


# Filtering
grep "^#\|SS=2" pairedVariants/varscan2.snp.vcf > pairedVariants/varscan2.snp.somatic.vcf


# Variants MuTecT2
java -Xmx2G -jar ${GATK_JAR} Mutect2 \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -I alignment/normal/normal.sorted.dup.recal.bam \
  -I alignment/tumor/tumor.sorted.dup.recal.bam \
  -normal normal \
  -tumor tumor \
  --germline-resource $REF/annotations/Homo_sapiens.GRCh37.gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz \
  -O pairedVariants/mutect2.vcf \
  -L 9:130215000-130636000


# Filtering
java -Xmx2G -jar ${GATK_JAR} FilterMutectCalls \
   -V pairedVariants/mutect2.vcf \
   --contamination-table contamination.table \
   -O pairedVariants/mutect2.filtered.vcf

vcftools --vcf pairedVariants/mutect2.vcf --stdout --remove-indels --recode | sed -e "s|normal|NORMAL|g" -e "s|tumor|TUMOR|g"  >  pairedVariants/mutect2.snp.somatic.vcf
  

java -classpath $VARDICT_HOME/lib/VarDict-1.4.9.jar:$VARDICT_HOME/lib/commons-cli-1.2.jar:$VARDICT_HOME/lib/jregex-1.2_01.jar:$VARDICT_HOME/lib/htsjdk-2.8.0.jar com.astrazeneca.vardict.Main \
  -G ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -N tumor_pair \
  -b "alignment/tumor/tumor.sorted.dup.recal.bam|alignment/normal/normal.sorted.dup.recal.bam"  \
  -Q 10 -f 0.05 -c 1 -S 2 -E 3 -g 4 -th 3 \
  -R 9:130215000-130636000 \
  | $VARDICT_BIN/testsomatic.R   \
  | $VARDICT_BIN/var2vcf_paired.pl -N "TUMOR|NORMAL" -f 0.05 > pairedVariants/vardict.vcf
  

bcftools filter \
   -i 'FILTER="PASS"&&TYPE="snp"&&INFO/STATUS="StrongSomatic"' \
   pairedVariants/vardict.vcf \
   | awk ' BEGIN {OFS="\t"} \
   { if(substr($0,0,1) == "#" || length($4) == length($5)) {if(substr($0,0,2) != "##") \
   {t=$10; $10=$11; $11=t} ; print}} ' > pairedVariants/vardict.snp.somatic.vcf

^#less pairedVariants/varscan2.snp.somatic.vcf
^#less pairedVariants/mutect2.snp.somatic.vcf
^#less pairedVariants/vardict.snp.somatic.vcf

# Unified callset
bcbio-variation-recall ensemble \
  --cores 2 --numpass 2 --names mutect2,varscan2,vardict \
  pairedVariants/ensemble.snp.somatic.vcf.gz \
  ${REF}/genome/Homo_sapiens.GRCh37.fa \
  pairedVariants/mutect2.snp.somatic.vcf \
  pairedVariants/varscan2.snp.somatic.vcf \
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

^#less -S pairedVariants/ensemble.snp.somatic.snpeff.vcf

# Coverage Track
for i in normal tumor
do
  igvtools count \
    -f min,max,mean \
    alignment/${i}/${i}.sorted.dup.recal.bam \
    alignment/${i}/${i}.sorted.dup.recal.bam.tdf \
    b37
done

exit


