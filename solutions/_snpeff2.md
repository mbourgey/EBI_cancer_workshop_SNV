There is many tools you can use to  explore vcf (vcftools, GEMIN ...) but for such a simple task the easiest way is tu use the old good grep command:

```
grep "HIGH\|MODERATE"  pairedVariants/vmutect2.snpeff.vcf | grep PASS | less
```

there is actually only 3 somatic variants with high or moderate  impact