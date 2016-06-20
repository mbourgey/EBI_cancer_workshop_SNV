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
 - 9 130296899 

Called as somatic indels with high/moderate impact by mutect2 but looks False Positive calls 

 - 9 130635680 - 
 - 9 130635681

Called as somatic snp with high/moderate impact and looks rea

 - 9 130634110

 