If you want to count the *un-aligned* reads you can use:

```{.bash}
samtools view -c -f4 alignment/normal/normal.sorted.bam
```

Number of unmapped reads :

	XXXXX

Or if you want to count the *aligned* reads you can use:

```
samtools view -c -F4 alignment/normal/normal.sorted.bam
```

Number of mapped reads :

	XXXXX
