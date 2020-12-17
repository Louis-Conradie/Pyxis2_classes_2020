Requirements:
	Python 3.9.0
	numpy
	statsmodels
	matplotlib
Description :
	Pyxis2 is a command line utility use for calculating something:
Usages :
	Pyxis2.py -v [GFF File] -s -w -i -k -p -q -l
	-v		GFF3 File
	-s		Text File Containing Regulated Genes
	-l		Step Size, Amount of Genes Apart Equaling a Cluster
	-k		Open Reading Frames Type
	-i		Sting to Search For
	-w		False Discovery Rate
	-p		Lower p, most significant value
	-q		Top p, Cut off Value	
Example :
	Pyxis2.py -v Example.gff3 -s RegulatedGenes.txt -l 2 -k 1 -i SearchString -w 2 -p 0 -q 10
	Pyxis2.py -v Example.gff3 -s RegulatedGenes.txt -l 2 -k "" -i SearchString -w 2 -p 0 -q 10

method str
	Method used for testing and adjustment of pvalues. Can be either the full name or initial 	letters. Available methods are:

•	bonferroni : one-step correction
•	sidak : one-step correction
•	holm-sidak : step down method using Sidak adjustments
•	holm : step-down method using Bonferroni adjustments
•	simes-hochberg : step-up method (independent)
•	hommel : closed method based on Simes tests (non-negative)
•	fdr_bh : Benjamini/Hochberg (non-negative)
•	fdr_by : Benjamini/Yekutieli (negative)
•	fdr_tsbh : two stage fdr correction (non-negative)
•	fdr_tsbky : two stage fdr correction (non-negative)
