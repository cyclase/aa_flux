# CAIS

All R scripts used to calculate CAIS slopes across 117 vertebrates, for the proteome, Pfams, and other subsets of the data (Young vs. Old domains, and transmembrane vs. non-transmembrane). 

Steps:
0. Calculate Pfam amino acid frequencies
	Input	AA_for_linearmodeling_species_differences.csv
	Script	AllPfam_AAC.R 
	Output	AllPfam_AAC.csv (amino acid frequencies within Pfam domains)

1. Calculate CAIS slopes for Pfams	
	Inputs	AllPfam_AAC.csv, new_CAIS.csv (CAIS values)
	Script	Vertebrate_PfamAACvsCAIS.R
	Output	PFAM_CAISeffectonAAfreq_SlopesandStderrors_noAntelope_new_CAIS.csv (CAIS slopes for Pfam domains)

2. Calculate CAIS slopes for proteome	
	Inputs	VertebrateAAC_CAIS.csv (amino acid frequencies for whole proteome), new_CAIS.csv
	Script	Vertebrate_PfamAACvsCAIS.R
	Output	CAISeffectonAAfreq_SlopesandStderrors_noAntelope_new_CAIS.csv (CAIS slopes for whole proteome)

3. Calculate amino acid frequencies for domain subsets (Age, Transmembrane)	
	Inputs	Pfamages.csv (Ages), AA_for_linearmodeling_species_differences.csv (Tmhmm results), new_CAIS.csv
	Script	Old_young_trans_nontrans_pfam_AAC.R
	Outputs	Old_Pfam_AAC.csv (amino acid frequencies for old domains), Recent_Pfam_AAC.csv (amino acid frequencies for young domains), NonTransMembrane_Pfam_AAC.csv (amino acid frequencies for non-transmembrane domains), TransMembrane_Pfam_AAC.csv (amino acid frequencies for transmembrane domains)

4. Calculate CAIS slopes for domain subsets (Age, Transmembrane)
	Inputs	Old_Pfam_AAC.csv (amino acid frequencies for old domains), Recent_Pfam_AAC.csv (amino acid frequencies for young domains), NonTransMembrane_Pfam_AAC.csv (amino acid frequencies for non-transmembrane domains), TransMembrane_Pfam_AAC.csv (amino acid frequencies for transmembrane domains), new_CAIS.csv
	Script	Age_Membrane_PfamAAC_VertebrateCAIS.R
	Outputs	CAISeffectonPfamAAC_Trans_NonTrans_Old_Young_newCAIS.csv (CAIS slopes for Old, Young, Transmembrane, and Non-transmembrane domains)

5. Calculate CAIS slopes for center-log ratio transformed amino acid frequencies (Pfam)
	Inputs	AllPfam_AAC.csv, new_CAIS.csv
	Script	clr_pfamAACvsCAIS.R
	Output	clr_Pfam_CAISeffectonAAfreq_SlopesandStderrors_noAntelope_new_CAIS.csv (CAIS slopes for whole proteome, where amino acid frequencies are center-log ratio transformed)

6. Calculate CAIS slopes for center-log ratio transformed amino acid frequencies (whole proteome)
	Inputs	VertebrateAAC_CAIS.csv, new_CAIS.csv
	Script	clr_proteomeAACvsCAIS.R
	Output	clr_proteome_CAISeffectonAAfreq_SlopesandStderrors_noAntelope_new_CAIS.csv (CAIS slopes for whole proteome, where amino acid frequencies are center-log ratio transformed)