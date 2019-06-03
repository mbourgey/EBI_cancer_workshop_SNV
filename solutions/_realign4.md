This makes it easier for down stream analysis tools

For NGS analysis, the convention is to left align indels. 


This is only really needed when calling variants with legacy locus-based tools such as samtools or GATK UnifiedGenotyper. Otherwise you will have worse performance and accuracy.

With more sophisticated tools (like GATK HaplotypeCaller) that involve reconstructing haplotypes (eg through reassembly), the problem of multiple valid representations is handled internally and does not need to be corrected explicitly.