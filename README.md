Requirements:
	
	Python 3.9.0
	
	numpy
	
	statsmodels
	
	matplotlib

Description :

	Pyxis2 is a command line utility use for calculating statistical significance of detect clusters

Usages :
	
	Pyxis2.py -v [GFF File] -s -w -i -k -p -q -l
	
	Required arguments:
	
		-v		GFF3 file
	
		-s		Text file containing regulated genes
	
	
	Optional arguments:
	
		-l		Step size
	
		-k		Open reading frames type [0 = Dubious, psuedo and protein coding ORFs ; 1 = Protein coding ORFs]
	
		-i		String to search in GFF3 file
	
		-w		False discovery rate
	
		-p		Lower p-value, most significant value
	
		-q		Maximum p-value, cut off value
	
Example :

	Pyxis2.py -v Example.gff3 -s RegulatedGenes.txt -l 2 -k 0 -i SearchString -w fdr_bh -p 0.0001 -q 0.10
	Pyxis2.py -v Example.gff3 -s RegulatedGenes.txt -l 2 -k 1 -i SearchString -p 0.0001 -q 0.10

False discovery rate string.
	
	bonferroni : one-step correction

	sidak : one-step correction

	holm-sidak : step down method using Sidak adjustments

	holm : step-down method using Bonferroni adjustments

	simes-hochberg : step-up method (independent)

	hommel : closed method based on Simes tests (non-negative)

	fdr_bh : Benjamini/Hochberg (non-negative)

	fdr_by : Benjamini/Yekutieli (negative)

	fdr_tsbh : two stage fdr correction (non-negative)

	fdr_tsbky : two stage fdr correction (non-negative)


