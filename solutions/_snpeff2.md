There is many tools you can use to  explore vcf (vcftools, GEMIN ...) but for such a simple task the easiest way is tu use the old good grep command:

```
grep "HIGH\|MODERATE"  pairedVariants/varscan.snpeff.vcf | grep SOMATIC | less
```
