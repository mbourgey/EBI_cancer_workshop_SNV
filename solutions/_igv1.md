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


Here are some variants:   

Looks good but located in intron, so possibly no impact:   

  * chr9 127534620

Called as somatic snp with high/moderate impact and looks real:   

  * chr9 127871831


We can see here, it is usually a good sign when several callers detect the mutation
