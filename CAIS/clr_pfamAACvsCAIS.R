library(stringr)
library(protr)
library(ggplot2)
library(seqinr)
library(ggrepel)
library(dplyr)
library(compositions)
library(data.table)

setwd("C:/Users/hanon/Documents/Github/aa_flux/CAIS")

#Center log ratio transformation of Pfam amino acid frequencies
pf_freqs <- read.delim(file = "AllPfam_AAC.csv", header=TRUE, stringsAsFactors=FALSE, sep=",")

t_pf_freqs <- transpose(pf_freqs)
rownames(t_pf_freqs) <- colnames(pf_freqs)
colnames(t_pf_freqs) <- rownames(pf_freqs)

t_pf_clr <- apply(t_pf_freqs, 2, clr) %>% as.data.frame()
t_pf_clr2 <- rbind(SpeciesUID = t_pf_freqs[1,], t_pf_clr[2:23,])

pf_clr <- transpose(t_pf_clr2)
rownames(pf_clr) <- colnames(t_pf_clr2)
colnames(pf_clr) <- rownames(t_pf_clr2)

write.csv(pf_clr, file="clr_pfamAAC.csv", quote=F)

# Read input files
# remove antelope ID 442
CAIS_value <- read.csv('CAIS_KLD.csv',header = T)
CAIS_values <-  CAIS_value[-which(CAIS_value$SpeciesUID == '442'),] # remove antelope

clr_pfamAAC <- read.csv('clr_pfamAAC.csv', header = T)
clr_pfamAAC2 <- clr_pfamAAC[-which(clr_pfamAAC$SpeciesUID == '442'),]

# choose if the pfam aac df or the protein aac df
SpeciesAAC_wCAIS <- clr_pfamAAC2
CAIS_values <- CAIS_values[which(CAIS_values$SpeciesUID %in% SpeciesAAC_wCAIS$SpeciesUID),]
match(CAIS_values$SpeciesUID, SpeciesAAC_wCAIS$SpeciesUID) # match

SpeciesAAC_Acol <- which(colnames(SpeciesAAC_wCAIS) =='A')
rm(VertebrateCAIS_AAC_df)
VertebrateCAIS_AAC_df <- cbind(CAIS_values, SpeciesAAC_wCAIS[,SpeciesAAC_Acol:dim(SpeciesAAC_wCAIS)[2]])
#write.csv(VertebrateCAIS_AAC_df, 'PfamAAC_CAISvalues_117Vertebrates.csv')

rowSums(VertebrateCAIS_AAC_df[,3:24]) # all species AAC do not add up to 1 because they are clr transformed

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
model_estimates <- data.frame(model_estimates)
#model_estimates <- model_estimates[match(AA_properties$Letter,rownames(model_estimates)),]

write.csv(model_estimates, 'clr_Pfam_CAISeffectonAAfreq_SlopesandStderrors_noAntelope_KLD.csv')



