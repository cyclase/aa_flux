# Flux

The R script used to calculate amino acid fluxes for 4 mammals: mouse, rat, human, and chimp. 

Calculate fluxes for each amino acid, for each species as well as for rodents and primates. Also calculate shared and species flux values.
	Input	Mouse_substitutions.csv, 
			Rat_substitutions.csv, 
			Human_substitutions.csv, 
			Chimp_substitutions.csv
			
	Script	flux.rmd 
	
	Output	Dall.csv (fluxes for each species), 
			Dprimate.csv (primate fluxes), 
			Drodent.csv (rodent fluxes), 
			Shared.csv (shared fluxes), 
			Species.csv (species fluxes), 
			all_fluxes.csv (all fluxes)
