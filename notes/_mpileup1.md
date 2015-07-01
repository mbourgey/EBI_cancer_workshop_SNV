Samtools:

	-L 1000 : max per-sample depth for INDEL calling [1000] ; 
	-B : disable BAQ (per-Base Alignment Quality) ; 
	-q 1 : skip alignments with mapQ smaller than 1 ; 
	-g generate genotype likelihoods in BCF format
Bcftools :

	-v output potential variant sites only
	-c SNP calling (force –e : likelihood based analyses)
	-g call genotypes at variant sites

