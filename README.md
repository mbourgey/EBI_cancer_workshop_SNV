
# Introduction to DNA-Seq processing for cancer data - SNVs
***By Mathieu Bourgey, Ph.D***  
*https://bitbucket.org/mugqic/mugqic_pipelines*

================================

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

================================

In this workshop, we will present the main steps that are commonly used to process and to analyze cancer sequencing data. We will focus only on whole genome data and provide command lines that allow detecting Single Nucleotide Variants (SNV). This workshop will show you how to launch individual steps of a complete DNA-Seq SNV pipeline using cancer data


## Data Source
We will be working on a CageKid sample pair, patient C0098.
The CageKid project is part of ICGC and is focused on renal cancer in many of it's forms.
The raw data can be found on EGA and calls, RNA and DNA, can be found on the ICGC portal. 
For more details about [CageKid](http://www.cng.fr/cagekid/)

For practical reasons we subsampled the reads from the sample because running the whole dataset would take way too much time and resources.



### Environment setup
```{.bash}

docker run --privileged -v /tmp:/tmp --network host -it -w $PWD -v $HOME:$HOME \
--user $UID:$GROUPS -v /etc/group:/etc/group  -v /etc/passwd:/etc/passwd \
-v /etc/fonts/:/etc/fonts/  c3genomics/genpipes:0.8


export REF=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/


cd $HOME/ebicancerworkshop2019/SNV


```


### Software requirements
These are all already installed, but here are the original links.

  * [BVATools](https://bitbucket.org/mugqic/bvatools/downloads)
  * [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  * [SAMTools](http://sourceforge.net/projects/samtools/)
  * [IGV](http://www.broadinstitute.org/software/igv/download)
  * [BWA](http://bio-bwa.sourceforge.net/)
  * [Genome Analysis Toolkit](http://www.broadinstitute.org/gatk/)
  * [SnpEff](http://snpeff.sourceforge.net/)
  * [Varscan2](http://varscan.sourceforge.net/)
  * [vardict](https://github.com/AstraZeneca-NGS/VarDictJava)
  * [bcbio variation](https://github.com/chapmanb/bcbio.variation)


We should load the corresponding modules 

```{.bash}
module load mugqic/java/openjdk-jdk1.8.0_72 \
   mugqic/bvatools/1.6 \
   mugqic/trimmomatic/0.36 \
   mugqic/samtools/1.9 \
   mugqic/bwa/0.7.17 \
   mugqic/GenomeAnalysisTK/4.1.0.0 \
   mugqic/R_Bioconductor/3.5.0_3.7 \
   mugqic/VarScan/2.4.3 \
   mugqic/vcftools/0.1.14 \
   mugqic/bcftools/1.9 \
   mugqic/VarDictJava/1.4.9 \
   mugqic/bcbio.variation.recall/0.1.7 \
   mugqic/snpEff/4.3 \
   mugqic/igvtools/2.3.67

```


## Original Setup

The initial structure of your folders should look like this:
```
<ROOT>
|-- raw_reads/               # fastqs from the center (down sampled)
    `-- normal               # The blood sample directory
        `-- run*_?           # Lane directory by run number. Contains the fastqs
    `-- tumor                # The tumor sample directory
        `-- run*_?           # Lane directory by run number. Contains the fastqs
|-- savedResults             # Folder containing precomputed results
|-- scripts                  # cheat sheet folder
|-- adapters.fa              # fasta file containing the adapter used for sequencing
```


### Cheat file
* You can find all the unix command lines of this practical in the file: [commands.sh](scripts/commands.sh)



# First data glance
So you've just received an email saying that your data is ready for download from the sequencing center of your choice.

**What should you do ?** [solution](solutions/_data.md)


### Fastq files
Let's first explore the fastq file.

Try these commands

```{.bash}
zless -S raw_reads/normal/run62DVGAAXX_1/normal.64.pair1.fastq.gz

```

**Why was it like that ?** [solution](solutions/_fastq1.md)


Now try these commands:

```{.bash}
zcat raw_reads/normal/run62DVGAAXX_1/normal.64.pair1.fastq.gz | head -n4
zcat raw_reads/normal/run62DVGAAXX_1/normal.64.pair2.fastq.gz | head -n4

```

**What was special about the output ? Why was it like that?** [Solution](solutions/_fastq2.md)

You could also just count the reads

```{.bash}
zgrep -c "^@HWUSI" raw_reads/normal/run62DVGAAXX_1/normal.64.pair1.fastq.gz

```

We should obtain 4003 reads

**Why shouldn't you just do ?** 

```{.bash}
zgrep -c "^@" raw_reads/normal/run62DVGAAXX_1/normal.64.pair1.fastq.gz

```

[Solution](solutions/_fastq3.md)


### Quality
We can't look at all the reads. Especially when working with whole genome 50x data. You could easily have Billions of reads.

Tools like FastQC and BVATools readsqc can be used to plot many metrics from these data sets.

Let's look at the data:

```{.bash}
# Generate original QC
mkdir -p originalQC/
java -Xmx1G -jar ${BVATOOLS_JAR} readsqc --quality 64 \
  --read1 raw_reads/normal/run62DVGAAXX_1/normal.64.pair1.fastq.gz \
  --read2 raw_reads/normal/run62DVGAAXX_1/normal.64.pair2.fastq.gz \
  --threads 2 --regionName normalrun62DVGAAXX_1 --output originalQC/

```

Open the images

All the generated graphics have their uses. But 3 of them are particularly useful to get an overal picture of how good or bad a run went.
        - The Quality box plots 
        - The nucleotide content graphs.
        - The Box plot shows the quality distribution of your data.
 
The quality of a base is computated using the Phread quality score.
[notes](notes/_fastQC1.md) 


The quality of a base is computated using the Phread quality score.
![Phred quality score formula](img/phred_formula.png)

In the case of base quality the probability use represents the probability of base to have been wrongly called
![Base Quality values](img/base_qual_value.png)

The formula outputs an integer that is encoded using an ASCII table. 

The way the lookup is done is by taking the the phred score adding 33 and using this number as a lookup in the table.

Older illumina runs, and the data here, were using phred+64 instead of phred+33 to encode their fastq files.

![ACII table](img/ascii_table.png)


**What stands out in the graphs ?**
[Solution](solutions/_fastqQC1.md)



**Why do we see adapters ?** 
[solution](solutions/_adapter1.md)

Although nowadays this doesn't happen often, it does still happen. In some cases, miRNA, it is expected to have adapters.


### Trimming
Since adapter are not part of the genome they should be removed

To do that we will use Trimmomatic.
 
The adapter file is in your work folder. 

```{.bash}
cat adapters.fa

```

**Why are there 2 different ones ?** [Solution](solutions/_trim1.md)


trimming with trimmomatic:


```{.bash}
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

```

[note on trimmomatic command](notes/_trimmomatic.md)

**What does Trimmomatic says it did ?** [Solution](solutions/_trim2.md)

Exercice: 
**Let's generate the new graphs** [Solution](solutions/_fastqQC2.md)

**How does it look now ?** [Solution](solutions/_trim3.md)



# Alignment
The raw reads are now cleaned up of artefacts we can align each lane separatly.

**Why should this be done separatly?** [Solution](solutions/_aln1.md)

**Why is it important to set Read Group information ?** [Solution](solutions/_aln2.md)

##Alignment with bwa-mem

```{.bash}
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

  bwa mem -M -t 5 \
    -R "@RG\\tID:${SNAME}_${RUNID}_${LANE}\\tSM:${SNAME}\\t\
LB:${SNAME}\\tPU:${RUNID}_${LANE}\\tCN:Centre National de Genotypage\\tPL:ILLUMINA" \
    ${REF}/genome/bwa_index/Homo_sapiens.GRCh37.fa \
    $file \
    ${file%.pair1.fastq.gz}.pair2.fastq.gz \
  | java -Xmx6G -jar ${GATK_JAR}  SortSam \
    -I /dev/stdin \
    -O ${OUTPUT_DIR}/${SNAME}.sorted.bam \
    --CREATE_INDEX true --SORT_ORDER coordinate --MAX_RECORDS_IN_RAM 500000
done

```
 
**Why did we pipe the output of one to the other ?** [Solution](solutions/_aln3.md)

**Could we have done it differently ?** [Solution](solutions/_aln4.md)


## Lane merging
We now have alignments for each of the sequences lanes:
 
   - This is not practical in it's current form. 
   - What we wan't to do now is merge the results into one BAM.

Since we identified the reads in the BAM with read groups, even after the merging, we can still identify the origin of each read.


```{.bash}
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

``` 

You should now have one BAM containing all your data.

Let's double check

```{.bash}
ls -l alignment/normal/
samtools view -H alignment/normal/normal.sorted.bam | grep "^@RG"

```

You should have your 9 read group entries.

**Why did we use the -H switch?** [Solution](solutions/_merge1.md)

**Try without. What happens?** [Solution](solutions/_merge2.md)


## SAM/BAM exploration
Let's spend some time to explore bam files.

```{.bash}
samtools view alignment/normal/normal.sorted.bam | head -n4

```

Here you have examples of alignment results.
A full description of the flags can be found in the SAM specification
http://samtools.sourceforge.net/SAM1.pdf

You can try using picards explain flag site to understand what is going on with your reads
http://broadinstitute.github.io/picard/explain-flags.html

The flag is the 2nd column.

**What do the flags of the first 4th reads mean?** [solutions](solutions/_sambam1.md)

Exercice:
**Let's take the 3rd one, the one that is in proper pair, and find it's mate.** [solutions](solutions/_sambam3.md)

**Why the pairing information is important ?**  [solutions](solutions/_sambam4.md)

## SAM/BAM filtering

You can use samtools to filter reads as well.

Exercice:
**How many reads mapped and unmapped were there?** [solution](solutions/_sambam2.md)


## SAM/BAM CIGAR string
Another useful bit of information in the SAM is the CIGAR string.
It's the 6th column in the file. 

This column explains how the alignment was achieved.
 
        M == base aligns *but doesn't have to be a match*. A SNP will have an M even if it disagrees with the reference.
        I == Insertion
        D == Deletion
        S == soft-clips. These are handy to find un removed adapters, viral insertions, etc.

An in depth explanation of the CIGAR can be found [here](http://genome.sph.umich.edu/wiki/SAM)

The exact details of the cigar string can be found in the SAM spec as well.


We won't go into too much detail at this point since we want to concentrate on cancer specific issues now.


# Cleaning up alignments
We started by cleaning up the raw reads. Now we need to fix some alignments.

The first step for this is to realign around indels and snp dense regions.

The Genome Analysis toolkit has a tool for this called IndelRealigner.

It basically runs in 2 steps:
 
   1. Find the targets
   2. Realign them
	

##GATK IndelRealigner

```{.bash}
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
  
```
**Why did we use both normal and tumor together?** [Solution](solutions/_realign3.md)

**How could we make this go faster ?** [Solution](solutions/_realign1.md)

**How many regions did it think needed cleaning ?** [Solution](solutions/_realign2.md)

Indel Realigner also makes sure the called deletions are left aligned when there is a microsatellite or homopolymer.

```
This
ATCGAAAA-TCG
into
ATCG-AAAATCG

or
ATCGATATATATA--TCG
into
ATCG--ATATATATATCG
```

**Why it is important ?**[Solution](solutions/_realign4.md)


## Mark duplicates
**What are duplicate reads ?** [Solution](solutions/_markdup1.md)

**What are they caused by ?** [Solution](solutions/_markdup2.md)

**What are the ways to detect them ?** [Solution](solutions/_markdup3.md)

Here we will use picards approach:

```{.bash}
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

```

We can look in the metrics output to see what happened.

```{.bash}
less alignment/normal/normal.sorted.dup.metrics

```

**How many duplicates were there ?** [Solution](solutions/_markdup4.md)

Dupliate reads number are estimated separately for each library.
 
**Why is this important to do not combine everything ?** [Solution](solutions/_markdup5.md)

[Note on Duplicate rate](notes/_marduop1.md)

## Base Quality recalibration
**Why do we need to recalibrate base quality scores ?** [Solution](solutions/_recal1.md)


It runs in 2 steps:

    1 - Build covariates based on context and known snp sites
    2 - Correct the reads based on these metrics


GATK BaseRecalibrator:

```{.bash}
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

```



# Extract BAM metrics
Once your whole bam is generated, it's always a good thing to check the data again to see if everything makes sens.

**Contamination**
It tells you if your date are contaminated or if a mix-up had occured

**Compute coverage**
If you have data from a capture kit, you should see how well your targets worked

**Insert Size**
It tells you if your library worked

**Alignment metrics**
It tells you if your sample and you reference fit together

## Estimate Normal/tumor contamination
To estimate these metrics we will use the GATK tool. This run in 2 steps:  

    1 - Generate GATK pileup tables
    2 - Estimate contamination


```{.bash}
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

```

Look at the contamination result file

```{.bash}
less contamination.table

```

**What do you think about this estimation ?** [solution](solutions/_conta.md)

## Compute coverage
Both GATK and BVATools have depth of coverage tools. 

Here we'll use the GATK one

```{.bash}
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
```
[note on DepthOfCoverage command](notes/_DOC.md)

In this project the expected coverages are : ~40x for the normal and ~50x for the tumor

Look at the coverage:

```{.bash}
less -S alignment/normal/normal.sorted.dup.recal.coverage.sample_interval_summary
less -S alignment/tumor/tumor.sorted.dup.recal.coverage.sample_interval_summary

```

**Is the coverage fit with the expectation ?** [solution](solutions/_DOC1.md)

## Insert Size
It corresponds to the size of DNA fragments sequenced.

Different from the gap size (= distance between reads) !

These metrics are computed using Picard:

```{.bash}
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

```

look at the output

```{.bash}
less -S alignment/normal/normal.sorted.dup.recal.metric.insertSize.tsv
less -S alignment/tumor/tumor.sorted.dup.recal.metric.insertSize.tsv

```

There is something interesting going on with our libraries.

**Can you tell what it is?** [Solution](solutions/_insert1.md)

**Which library is the most suitable for cancer analysis ?** [Solution](solutions/_insert2.md)

## Alignment metrics
For the alignment metrics, samtools flagstat is very fast but with bwa-mem since some reads get broken into pieces, the numbers are a bit confusing. 

We prefer the Picard way of computing metrics:

```{.bash}
# Get alignment metrics
for i in normal tumor
do
  java -Xmx2G -jar ${GATK_JAR}  CollectAlignmentSummaryMetrics \
    -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
    -I alignment/${i}/${i}.sorted.dup.recal.bam \
    -O alignment/${i}/${i}.sorted.dup.recal.metric.alignment.tsv \
    --METRIC_ACCUMULATION_LEVEL LIBRARY
done

```

explore the results

```{.bash}
less -S alignment/normal/normal.sorted.dup.recal.metric.alignment.tsv
less -S alignment/tumor/tumor.sorted.dup.recal.metric.alignment.tsv

```

**Do you think the sample and the reference genome fit together ?** [Solution](solutions/_alnMetrics1.md)

# Variant calling

![SNV call summary workflow](img/snv_call.png)

Most of SNV caller use either a Baysian, a threshold or a t-test approach to do the calling

 Here we will try 3 variant callers.
- Varscan 2
- MuTecT2
- Vardict



many, MANY others can be found here:
https://www.biostars.org/p/19104/

In our case, let's start with:

```{.bash}
mkdir pairedVariants

```

## varscan 2

VarScan calls somatic variants (SNPs and indels) using a heuristic method and a statistical test based on the number of aligned reads supporting each allele.


Varscan somatic caller expects both a normal and a tumor file in SAMtools pileup format. From sequence alignments in binary alignment/map (BAM) format. To build a pileup file, you will need:

- A SAM/BAM file ("myData.bam") that has been sorted using the sort command of SAMtools.
- The reference sequence ("reference.fasta") to which reads were aligned, in FASTA format.
- The SAMtools software package.


```{.bash}
# SAMTools mpileup
for i in normal tumor
do
samtools mpileup -B -q 1 \
  -f ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -r 9:130215000-130636000 \
  alignment/${i}/${i}.sorted.dup.recal.bam \
  > pairedVariants/${i}.mpileup
done

```
[note on samtools mpileup command](notes/_mpileup1.md)


Now we can run varscan:

```{.bash}
# varscan
java -Xmx2G -jar ${VARSCAN2_JAR} somatic \
   pairedVariants/normal.mpileup \
   pairedVariants/tumor.mpileup \
   pairedVariants/varscan2 \
   --output-vcf 1 \
   --strand-filter 1 \
   --somatic-p-value 0.001 

```

Then we can extract somatic SNPs:

```{.bash}
# Filtering
grep "^#\|SS=2" pairedVariants/varscan2.snp.vcf > pairedVariants/varscan2.snp.somatic.vcf

```


## GATK MuTecT2

```{.bash}
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

```

Then we can extract somatic SNPs:

```{.bash}
# Filtering
java -Xmx2G -jar ${GATK_JAR} FilterMutectCalls \
   -V pairedVariants/mutect2.vcf \
   --contamination-table contamination.table \
   -O pairedVariants/mutect2.filtered.vcf

vcftools --vcf pairedVariants/mutect2.vcf --stdout --remove-indels --recode | sed -e "s|normal|NORMAL|g" -e "s|tumor|TUMOR|g"  >  pairedVariants/mutect2.snp.somatic.vcf
  
```


## Vardict

```{.bash}
$VARDICT_HOME/bin/VarDict \
  -G ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -N tumor_pair \
  -b "alignment/tumor/tumor.sorted.dup.recal.bam|alignment/normal/normal.sorted.dup.recal.bam"  \
  -C -Q 10 -c 1 -S 2 -E 3 -g 4 -th 3 \
  | $VARDICT_BIN/testsomatic.R   \
  | $VARDICT_BIN/var2vcf_paired.pl -N "TUMOR|NORMAL" -f 0.05 > pairedVariants/vardict.vcf
  
```

Then we can extract somatic SNPs:

```{.bash}
bcftools filter \
   -i 'FILTER="PASS"&&TYPE="snp"&&INFO/STATUS="StrongSomatic"' \
   pairedVariants/vardict.vcf \
   | awk ' BEGIN {OFS="\t"} \
   { if(substr($0,0,1) == "#" || length($4) == length($5)) {if(substr($0,0,2) != "##") \
   {t=$10; $10=$11; $11=t} ; print}} ' > pairedVariants/vardict.snp.somatic.vcf
```

Now we have somatic variants from all two methods. Let's look at the results. 

```{.bash}
less pairedVariants/varscan2.snp.somatic.vcf
less pairedVariants/mutect2.snp.somatic.vcf
less pairedVariants/vardict.snp.somatic.vcf
```

**could you notice something from these vcf files ?** [Solution](solutions/_vcf1.md)


Details on the spec can be found here:
http://vcftools.sourceforge.net/specs.html

Fields vary from caller to caller.
 
Some values are are almost always there: 
 
   - The ref vs alt alleles, 
   - variant quality (QUAL column)
   - The per-sample genotype (GT) values.

[note on the vcf format fields](notes/_vcf1.md)

Choosing the best caller is not an easy task each of them have their pros and cons. Now new methods have been developped to extract the best information from a multiple set of variant caller. These methods refer to the ensemble approach (as developped in bcbio.variation or somaticSeq) and rely on pre-selecting a subset of variants from the interesect of multiple caller and then apply Machine Learning approach to filter the high quality variants.

As we don't have enough variant for the full ensemble approach we will just launch the initial step in order to generate a unifed callset form all the call found in at least 2 different caller:

```{.bash}
# Unified callset
bcbio-variation-recall ensemble \
  --cores 2 --numpass 2 --names mutect2,varscan2,vardict \
  pairedVariants/ensemble.snp.somatic.vcf.gz \
  ${REF}/genome/Homo_sapiens.GRCh37.fa \
  pairedVariants/mutect2.snp.somatic.vcf \
  pairedVariants/varscan2.snp.somatic.vcf \
  pairedVariants/vardict.snp.somatic.vcf

```

look at the unified callset

```{.bash}
zless pairedVariants/ensemble.snp.somatic.vcf.gz

```


# Annotations
The next step in trying to make sense of the variant calls is to assign functional consequence to each variant.

At the most basic level, this involves using gene annotations to determine if variants are sense, missense, or nonsense. 

We typically use SnpEff but many use Annovar and VEP as well.
Let's run snpEff:

```{.bash}
# SnpEff
java  -Xmx6G -jar ${SNPEFF_HOME}/snpEff.jar \
  eff -v -c ${SNPEFF_HOME}/snpEff.config \
  -o vcf \
  -i vcf \
  -stats pairedVariants/ensemble.snp.somatic.snpeff.stats.html \
  GRCh37.75 \
  pairedVariants/ensemble.snp.somatic.vcf.gz \
  > pairedVariants/ensemble.snp.somatic.snpeff.vcf
```

You can learn more about the meaning of snpEff annotations [here](http://snpeff.sourceforge.net/SnpEff_manual.html#output).

Use less to look at the new vcf file: 

```{.bash}
less -S pairedVariants/ensemble.snp.somatic.snpeff.vcf
```

**Can you see the difference with the previous vcf ?** [solution](solutions/_snpeff1.md)


The annotation is presented in the INFO field using the new ANN format. For more information on this field see [here](http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf). Typically, we have: 


`ANN=Allele|Annotation|Putative impact|Gene name|Gene ID|Feature type|Feature ID|Transcript biotype|Rank Total|HGVS.c|...`

Here's an example of a typical annotation: 

`ANN=T|intron_variant|MODIFIER|FAM129B|ENSG00000136830|transcript|ENST00000373312|protein_coding|1/13|c.56-2842T>A|||`

**What does the example annotation actually mean?** [solution](solutions/_snpEff3.md)


Exercice: 
**Find a somatic mutation with a predicted High or Moderate impact** [solution](solutions/_snpeff2.md)


**What effect categories were represented in these variants?** [solution](solutions/_snpEff4.md)




Next, you should view the report generated by snpEff.

Use the procedure described previously to retrieve:

`snpEff_summary.html`



## Data visualisation
The Integrative Genomics Viewer (IGV) is an efficient visualization tool for interactive exploration of large genome datasets. 

![IGV browser presentation](img/igv.png)

Before jumping into IGV, we'll generate a track IGV can use to plot coverage:

```{.bash}
# Coverage Track
for i in normal tumor
do
  igvtools count \
    -f min,max,mean \
    alignment/${i}/${i}.sorted.dup.recal.bam \
    alignment/${i}/${i}.sorted.dup.recal.bam.tdf \
    b37
done
```

Then:
 
   1. Open IGV
   2. Chose the reference genome corresponding to those use for alignment (b37)
   3. Load bam file
   4. Load vcf files

Explore/play with the data: 
 
   - find somatic variants
   - Look around...

[solution](solutions/_igv1.md)

**Open the high impact position in IGV, what do you see?** [solution](solutions/_snpEff5.md)

# Exit the container environment

```{.bash}
exit

```


## Aknowledgments
I would like to thank and acknowledge Louis Letourneau for this help and for sharing his material. The format of the tutorial has been inspired from Mar Gonzalez Porta. I also want to acknowledge Joel Fillon, Louis Letrouneau (again), Robert Eveleigh, Edouard Henrion, Francois Lefebvre, Maxime Caron and Guillaume Bourque for the help in building these pipelines and working with all the various datasets.
