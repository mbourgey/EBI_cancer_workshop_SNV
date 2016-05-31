some commands to find somatic variant in the vcf file

## varscan

```{.bash}
grep SOMATIC pairedVariants/varscan.vcf 
```

## MuTecT
```{.bash}
grep -v REJECT pairedVariants/mutect.vcf | grep -v "^#"
```

## Strelka
```{.bash}
grep -v "^#" pairedVariants/strelka.vcf
```


Then look at some of these position in IGV
