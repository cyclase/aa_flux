library(stringr)
library(protr)
library(ggplot2)
library(seqinr)
library(ggrepel)
library(dplyr)
library(compositions)
library(data.table)

setwd("C:/Users/hanon/Documents/Github/aa_flux/CAIS")

#Center log ratio transformation of proteome-wide amino acid frequencies
pr_freqs <- read.delim(file = "VertebrateAAC_CAIS.csv", header=TRUE, stringsAsFactors=FALSE, sep=",")

t_pr_freqs <- transpose(pr_freqs)
rownames(t_pr_freqs) <- colnames(pr_freqs)
colnames(t_pr_freqs) <- rownames(pr_freqs)

t_pr_clr <- apply(t_pr_freqs, 2, clr) %>% as.data.frame()
t_pr_clr2 <- rbind(SpeciesUID = t_pr_freqs[2,], t_pr_clr[3:24,])

pr_clr <- transpose(t_pr_clr2)
rownames(pr_clr) <- colnames(t_pr_clr2)
colnames(pr_clr) <- rownames(t_pr_clr2)

write.csv(pr_clr, file="clr_proteomeAAC.csv", quote=F)

# Read input files
# remove antelope ID 442
CAIS_value <- read.csv('CAIS_KLD.csv',header = T)
CAIS_values <-  CAIS_value[-which(CAIS_value$SpeciesUID == '442'),] # remove antelope

clr_proteomeAAC <- read.csv('clr_proteomeAAC.csv', header = T)
clr_proteomeAAC2 <- clr_proteomeAAC[-which(duplicated(clr_proteomeAAC$SpeciesUID)),]

#clr_proteomeAAC2 <- clr_pfamAAC[-which(clr_proteomeAAC$SpeciesUID == '442'),]

# choose if the pfam aac df or the protein aac df
SpeciesAAC_wCAIS <- clr_proteomeAAC2
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

write.csv(model_estimates, 'clr_proteome_CAISeffectonAAfreq_SlopesandStderrors_noAntelope_KLD.csv')



