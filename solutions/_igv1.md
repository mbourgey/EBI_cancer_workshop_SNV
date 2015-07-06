some commands to find somatic variant in the vcf file

## SAMtools

```{.bash}
cat pairedVariants/mpileup.vcf | perl -ne 'my @values=split("\t"); my ($clr) = $values[7] =~ /CLR=(\d+)/; if(defined($clr) && $clr >= 90 && $values[5] >= 100) {print "$_"}'
```

## MuTecT
```{.bash}
zgrep -v REJECT pairedVariants/mutect.vcf | grep -v "^#"
```

## Strelka
```{.bash}
 grep -v "^#" pairedVariants/strelka.vcf
```


Then look at some of these position in IGV