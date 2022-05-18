There is many tools you can use to  explore vcf (vcftools, GEMIN ...) but for such a simple task the easiest way is tu use the old good grep command:

```
grep "HIGH\|MODERATE"  pairedVariants/ensemble.snp.somatic.snpeff.vcf | less
```

There is actually only 1 somatic variants with high impact and 2 variant with moderate impact.
