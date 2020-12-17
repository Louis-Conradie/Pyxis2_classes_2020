**Requirements:**

		Python 3.9.0

**Required Modules:**

		numpy

		statsmodels

		matplotlib

**Description :**

		Pyxis2 is a command line utility use for calculating statistical significance of detect clusters

**Usages :**

		Pyxis2.py -v [GFF3 File] -s [Text File] -w -i -k -p -q -l

*Required arguments:*

		-v		GFF3 file

		-s		Text file containing regulated genes

*Optional arguments:*

		-l		Step size to detect gene clusters (default: 5)

		-k		Open reading frames type [0 = Dubious, psuedo and protein coding ORFs ; 1 = Protein coding ORFs] (default: 1)

		-i		Search string used to identify Genes in GFF file (default: ID=gene:)

		-w		False discovery rate

		-p		Lower P-value used to make statistical test more rigorous (default: 0.01)

		-q		Maximum p-value to make statistical test more rigorous (default:0.05)

**Example :**

		Pyxis2.py -v Example.gff3 -s example_text_file.txt -l 2 -k 0 -i ID=gene: -w fdr_bh -p 0.01 -q 0.05
		Pyxis2.py -v Example.gff3 -s example_text_file.txt -l 2 -k 1 -i ID=gene: -p 0.0001 -q 0.10

**False discovery rate acceptable strings.**

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
