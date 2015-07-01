cd /lb/project/mugqic/projects/kyotoWorkshop_2015

m mugqic/tabix/0.2.6 mugqic/igvtools/2.3.14 mugqic/snpEff/4.0 mugqic/GenomeAnalysisTK/3.2-2 mugqic/bvatools/1.4 mugqic/trimmomatic/0.32 mugqic/java/openjdk-jdk1.7.0_60 mugqic/bwa/0.7.10 mugqic/picard/1.123

export REF=$(pwd)/references/

#zless -S raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz

zcat raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz | head -n4
zcat raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair2.fastq.gz | head -n4

zgrep -c "^@SRR" raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz

zgrep -c "^@" raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz

mkdir originalQC/

java -Xmx1G -jar ${BVATOOLS_JAR} readsqc   --read1 raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz   --read2 raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair2.fastq.gz   --threads 2 --regionName SRR --output originalQC/

java -Xmx1G -jar ${BVATOOLS_JAR} readsqc   --read1 raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair1.fastq.gz   --read2 raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair2.fastq.gz   --threads 2 --regionName ERR --output originalQC/

cat adapters.fa

mkdir -p reads/NA12878/runSRR_1/

mkdir -p reads/NA12878/runERR_1/

java -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33   raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair1.fastq.gz   raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair2.fastq.gz   reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair1.fastq.gz   reads/NA12878/runERR_1/NA12878.ERR.t20l32.single1.fastq.gz   reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair2.fastq.gz   reads/NA12878/runERR_1/NA12878.ERR.t20l32.single2.fastq.gz   ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:20 MINLEN:32   2> reads/NA12878/runERR_1/NA12878.ERR.trim.out
java -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33   raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz   raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair2.fastq.gz   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair1.fastq.gz   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.single1.fastq.gz   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair2.fastq.gz   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.single2.fastq.gz   ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:20 MINLEN:32   2> reads/NA12878/runSRR_1/NA12878.SRR.trim.out

cat reads/NA12878/runERR_1/NA12878.ERR.trim.out reads/NA12878/runSRR_1/NA12878.SRR.trim.out

mkdir postTrimQC/

java -Xmx1G -jar ${BVATOOLS_JAR} readsqc   --read1 reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair1.fastq.gz   --read2 reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair2.fastq.gz   --threads 2 --regionName ERR --output postTrimQC/

java -Xmx1G -jar ${BVATOOLS_JAR} readsqc   --read1 reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair1.fastq.gz   --read2 reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair2.fastq.gz   --threads 2 --regionName SRR --output postTrimQC/

mkdir -p alignment/NA12878/runERR_1
mkdir -p alignment/NA12878/runSRR_1

bwa mem -M -t 2   -R '@RG\tID:ERR_ERR_1\tSM:NA12878\tLB:ERR\tPU:runERR_1\tCN:Broad Institute\tPL:ILLUMINA'   ${REF}/b37.fasta   reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair1.fastq.gz   reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair2.fastq.gz   | java -Xmx2G -jar ${PICARD_HOME}/SortSam.jar   INPUT=/dev/stdin   OUTPUT=alignment/NA12878/runERR_1/NA12878.ERR.sorted.bam   CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000

bwa mem -M -t 2   -R '@RG\tID:SRR_SRR_1\tSM:NA12878\tLB:SRR\tPU:runSRR_1\tCN:Broad Institute\tPL:ILLUMINA'   ${REF}/b37.fasta   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair1.fastq.gz   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair2.fastq.gz   | java -Xmx2G -jar ${PICARD_HOME}/SortSam.jar   INPUT=/dev/stdin   OUTPUT=alignment/NA12878/runSRR_1/NA12878.SRR.sorted.bam   CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000

java -Xmx2G -jar ${PICARD_HOME}/MergeSamFiles.jar   INPUT=alignment/NA12878/runERR_1/NA12878.ERR.sorted.bam   INPUT=alignment/NA12878/runSRR_1/NA12878.SRR.sorted.bam   OUTPUT=alignment/NA12878/NA12878.sorted.bam   VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

ls -l alignment/NA12878/

samtools view -H alignment/NA12878/NA12878.sorted.bam | grep "^@RG"

samtools view alignment/NA12878/NA12878.sorted.bam | head -n2

samtools view alignment/NA12878/NA12878.sorted.bam | grep ERR001742.6173685

samtools view -c -f4 alignment/NA12878/NA12878.sorted.bam

samtools view -c -F4 alignment/NA12878/NA12878.sorted.bam

java -Xmx2G  -jar ${GATK_JAR}   -T RealignerTargetCreator   -R ${REF}/b37.fasta   -o alnment/NA12878/realign.intervals   -I alignment/NA12878/NA12878.sorted.bam   -L 1

java -Xmx2G -jar ${GATK_JAR}   -T IndelRealigner   -R ${REF}/b37.fasta   -targetIntervals alignment/NA12878/realign.intervals   -o alignment/NA12878/NA12878.realigned.sorted.bam   -I alignment/NA12878/NA12878.sorted.bam

java -Xmx2G -jar ${PICARD_HOME}/FixMateInformation.jar   VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000   INPUT=alignment/NA12878/NA12878.realigned.sorted.bam   OUTPUT=alignment/NA12878/NA12878.matefixed.sorted.bam


java -Xmx2G -jar ${PICARD_HOME}/MarkDuplicates.jar   REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true   INPUT=alignment/NA12878/NA12878.matefixed.sorted.bam   OUTPUT=alignment/NA12878/NA12878.sorted.dup.bam   METRICS_FILE=alignment/NA12878/NA12878.sorted.dup.metrics

java -Xmx2G -jar ${GATK_JAR}   -T BaseRecalibrator   -nct 2   -R ${REF}/b37.fasta   -knownSites ${REF}/dbSnp-137.vcf.gz   -L 1:47000000-47171000   -o alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp   -I alignment/NA12878/NA12878.sorted.dup.bam

java -Xmx2G -jar ${GATK_JAR}   -T PrintReads   -nct 2   -R ${REF}/b37.fasta   -BQSR alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp   -o alignment/NA12878/NA12878.sorted.dup.recal.bam   -I alignment/NA12878/NA12878.sorted.dup.bam

java -Xmx2G -jar ${GATK_JAR}   -T BaseRecalibrator   -nct 2   -R ${REF}/b37.fasta   -knownSites ${REF}/dbSnp-137.vcf.gz   -L 1:47000000-47171000   -o alignment/NA12878/NA12878.sorted.dup.recalibration_report.seconnd.grp   -I alignment/NA12878/NA12878.sorted.dup.bam   -BQSR alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp

java -Xmx2G -jar ${GATK_JAR}   -T AnalyzeCovariates   -R ${REF}/b37.fasta   -before alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp   -after alignment/NA12878/NA12878.sorted.dup.recalibration_report.seconnd.grp   -csv BQSR.csv   -plots BQSR.pdf

java  -Xmx2G -jar ${GATK_JAR}   -T DepthOfCoverage   --omitDepthOutputAtEachBase   --summaryCoverageThreshold 10   --summaryCoverageThreshold 25   --summaryCoverageThreshold 50   --summaryCoverageThreshold 100   --start 1 --stop 500 --nBins 499 -dt NONE   -R ${REF}/b37.fasta   -o alignment/NA12878/NA12878.sorted.dup.recal.coverage   -I alignment/NA12878/NA12878.sorted.dup.recal.bam   -L 1:47000000-47171000

#less -S alignment/NA12878/NA12878.sorted.dup.recal.coverage.sample_interval_summary

java -Xmx2G -jar ${PICARD_HOME}/CollectInsertSizeMetrics.jar   VALIDATION_STRINGENCY=SILENT   REFERENCE_SEQUENCE=${REF}/b37.fasta   INPUT=alignment/NA12878/NA12878.sorted.dup.recal.bam   OUTPUT=alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.tsv   HISTOGRAM_FILE=alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.histo.pdf   METRIC_ACCUMULATION_LEVEL=LIBRARY

#less -S alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.tsv

java -Xmx2G -jar ${PICARD_HOME}/CollectAlignmentSummaryMetrics.jar   VALIDATION_STRINGENCY=SILENT   REFERENCE_SEQUENCE=${REF}/b37.fasta   INPUT=alignment/NA12878/NA12878.sorted.dup.recal.bam   OUTPUT=alignment/NA12878/NA12878.sorted.dup.recal.metric.alignment.tsv   METRIC_ACCUMULATION_LEVEL=LIBRARY

#less -S alignment/NA12878/NA12878.sorted.dup.recal.metric.alignment.tsv

mkdir variants

samtools mpileup -L 1000 -B -q 1 -D -S -g   -f ${REF}/b37.fasta   -r 1:47000000-47171000   alignment/NA12878/NA12878.sorted.dup.recal.bam   | bcftools view -vcg -   > variants/mpileup.vcf

#java -Xmx2G -jar ${GATK_JAR}   -T UnifiedGenotyper   -R ${REF}/b37.fasta   -I alignment/NA12878/NA12878.sorted.dup.recal.bam   -o variants/ug.vcf   --genotype_likelihoods_model BOTH   -dt none   -L 1:46000000-47600000

for i in variants/*.vcf;do bgzip -c $i > $i.gz ; tabix -p vcf $i.gz;done

#zless -S variants/mpileup.vcf.gz

java  -Xmx6G -jar ${SNPEFF_HOME}/snpEff.jar   eff -v -c ${SNPEFF_HOME}/snpEff.config   -o vcf   -i vcf   -stats variants/mpileup.snpeff.vcf.stats.html   GRCh37.74   variants/mpileup.vcf   > variants/mpileup.snpeff.vcf

igvtools count   -f min,max,mean   alignment/NA12878/NA12878.sorted.dup.recal.bam   alignment/NA12878/NA12878.sorted.dup.recal.bam.tdf   b37


