Summary Outline:

1. Data Pre-Processing
	A.) Python 
		i.) Create sets of subtype ids after file cleanup
			a.) file cleanup = removing characters after sample
ID-key
		ii.) ID Cleanup
			a.) Id-keys in the TCGA expression data
(TCGA_489_UE.txt), the clinical data (ov_tcga_pub_clinical_data.tsv), and the
cancer subtype idkey files (fallopian.idkey.txt, etc), do not follow a uniform
naming convention. 
			b.) Id-keys must be truncated to 20 characters from 28,
			c.) the clinical data must me parsed to swap every '-'
in sample id-keys with '.' to match the sample naming format of the expression data file.
		iii.) Phenotype Matrix Generation
			a.) Clinical Data parsed in matrix and partitioned
into resistant/sensitive platinum status based on disease response to
chemotherapy in the absent of explicitly defined platinum status
			b.) clinical data list is partitioned into
sub-matrices for each cancer subtype to be partitioned further into smaller
matrices
			c.) submatrices are combined into phenotype matrix
with 3 columns: "Sample_ID", "Subtype", "Platinum.Status"
	B.) R
		a.) Organizing/Subsetting
		b.) Limma
		c.) output/representation


Input/Output Notes:
1. idkey files are passed via the command line, and must be passed into the
python program execution in the order they appear in subtypeIdClasses.py
script. 
2. data files are hardcoded into the clin_dataParserv1.py script
3. Model 1 -- Uses no contrast matrix: design matrix column 1 is expression
levels of sensitive platinum status and column 2 represents difference in
expression level between expression levels of sensitive and resistant groups. 
   Model 2 -- Uses contrast matrix for different parametrization of data such
that columns in model 2's design matrix represent expression levels in the
platinum status groups. 
4. Model 1 and Model 2 yielded the same results. 
