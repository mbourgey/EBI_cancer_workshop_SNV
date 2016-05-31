this can be done using this command:

```{.bash}
grep "HIGH\|MODERATE" pairedVariants/varscan.snpeff.vcf | grep "SOMATIC" | less
```

A missens mutation in AK1 is predicted to have a moderate impact on the protein