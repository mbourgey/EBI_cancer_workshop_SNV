## SAM/BAM exploration
Let's spend some time to explore bam files.

To have examples of alignment results try:

```
samtools view alignment/NA12878/NA12878.sorted.bam | head -n4
```

A full description of the flags can be found in the SAM specification
http://samtools.sourceforge.net/SAM1.pdf

Try using picards explain flag site to understand what is going on with your reads
http://broadinstitute.github.io/picard/explain-flags.html

The flag is the 2nd column.

**What do the flags of the first 1st and 3rd reads mean?** [solutions](../solutions/_sambam1.md)

Let's take the 2nd one, the one that is in proper pair, and find it's pair.

## SAM/BAM filtering

You can use samtools to filter reads as well.

If you want to count the *un-aligned* reads you can use:

```
samtools view -c -f4 alignment/NA12878/NA12878.sorted.bam
```

Or if you want to count the *aligned* reads you can use:

```
samtools view -c -F4 alignment/NA12878/NA12878.sorted.bam
```

**How many reads mapped and unmapped were there?** [solution](../solutions/_sambam2.md)

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
