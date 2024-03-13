library(stringr)
library(protr)
library(ggplot2)
library(seqinr)
library(ggrepel)

setwd('/Users/hanon/Documents/GitHub/aa_flux/CAIS')

# Read input files
# remove antelope ID 442
CAIS_values <- read.csv('CAIS_KLD.txt', header = T)
CAIS_values <-  CAIS_values[-which(CAIS_values$SpeciesUID == '442'),] # remove antelope

#VertebrateAAC_CAIS.csv comes from the protein sequences in vertebrates in sql
Species_ProteinAAC <- read.csv('VertebrateAAC_CAIS.csv', header = T)
Species_ProteinAAC <- Species_ProteinAAC[-which(duplicated(Species_ProteinAAC$SpeciesUID)),]
Species_PFAM_AAC <- read.csv('AllPfam_AAC.csv', header = T)
Species_PFAM_AAC <- Species_PFAM_AAC[-which(Species_PFAM_AAC$SpeciesUID == '442'),]
AA_properties <- read.csv('../Figures/Short_amino_acid_properties_Oct23.csv', header = T)


# choose if the pfam aac df or the protein aac df
SpeciesAAC_wCAIS <- Species_PFAM_AAC
CAIS_values <- CAIS_values[which(CAIS_values$SpeciesUID %in% SpeciesAAC_wCAIS$SpeciesUID),]
match(CAIS_values$SpeciesUID, SpeciesAAC_wCAIS$SpeciesUID) # match

#SpeciesAAC_Acol <- which(colnames(SpeciesAAC_wCAIS) =='A')
#rm(VertebrateCAIS_AAC_df)
VertebrateCAIS_AAC_df <- cbind(CAIS_values, SpeciesAAC_wCAIS[,SpeciesAAC_Acol:dim(SpeciesAAC_wCAIS)[2]])
#write.csv(VertebrateCAIS_AAC_df, 'PfamAAC_CAISvalues_117Vertebrates.csv')

rowSums(VertebrateCAIS_AAC_df[,3:24]) # all species AAC add up to 1

# Intialize dataframe
AAstart <- which(colnames(VertebrateCAIS_AAC_df) == 'A')
AAend <- which(colnames(VertebrateCAIS_AAC_df) == 'V')
model_estimates <- data.frame(matrix(nrow = 2))

# loop model over all amino acids
for (aa in AAstart:AAend) {

model <- lm(VertebrateCAIS_AAC_df[,aa]~VertebrateCAIS_AAC_df$CAIS )
model_summary <- summary(model)
model_summary <- as.data.frame(model_summary$coefficients[2,1:2])
colnames(model_summary) <- colnames(VertebrateCAIS_AAC_df)[aa]
model_estimates <- cbind(model_estimates, model_summary) }


model_estimates <- model_estimates[,-1]
rownames(model_estimates) <- c('Estimates','Std.Error')
model_estimates <- as.data.frame(t(model_estimates))


# Remove O and U
model_estimates <- model_estimates[-which(model_estimates[,1]== 0),]
#model_estimates <- data.frame(model_estimates)
model_estimates <- model_estimates[match(AA_properties$aa,rownames(model_estimates)),]

write.csv(model_estimates, 'PFAM_CAISeffectonAAfreq_SlopesandStderrors_noAntelope_KLD.csv')



