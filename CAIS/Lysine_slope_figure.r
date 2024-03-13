setwd('/Users/hanon/Documents/Github/aa_flux/CAIS')


CAIS_values <- read.csv('new_CAIS.csv',header = T) 
CAIS_values <-  CAIS_values[-which(CAIS_values$SpeciesUID == '442'),] # remove antelope

#VertebrateAAC_CAIS.csv comes from the protein sequences in vertebrates in sql
Species_ProteinAAC <- read.csv('VertebrateAAC_CAIS.csv', header = T)
Species_ProteinAAC <- Species_ProteinAAC[-which(duplicated(Species_ProteinAAC$SpeciesUID)),]
Species_PFAM_AAC <- read.csv('AllPfam_AAC.csv', header = T)
Species_PFAM_AAC <- Species_PFAM_AAC[-which(Species_PFAM_AAC$SpeciesUID == '442'),]

SpeciesAAC_wCAIS <- Species_PFAM_AAC
CAIS_values <- CAIS_values[which(CAIS_values$SpeciesUID %in% SpeciesAAC_wCAIS$SpeciesUID),]
match(CAIS_values$SpeciesUID, SpeciesAAC_wCAIS$SpeciesUID) # match

SpeciesAAC_Acol <- which(colnames(SpeciesAAC_wCAIS) =='A')
rm(VertebrateCAIS_AAC_df)
VertebrateCAIS_AAC_df <- cbind(CAIS_values, SpeciesAAC_wCAIS[,SpeciesAAC_Acol:dim(SpeciesAAC_wCAIS)[2]])
#write.csv(VertebrateCAIS_AAC_df, 'PfamAAC_CAISvalues_117Vertebrates.csv')

rowSums(VertebrateCAIS_AAC_df[,3:24]) # all species AAC add up to 1

# Intialize dataframe
AAstart <- which(colnames(VertebrateCAIS_AAC_df) == 'A')
AAend <- which(colnames(VertebrateCAIS_AAC_df) == 'V')
model_estimates <- data.frame(matrix(nrow = 2))

# loop model over all amino acids
  model <- lm(VertebrateCAIS_AAC_df[,"K"]~VertebrateCAIS_AAC_df$CAIS )
  model_summary <- summary(model)

library("ggplot2")
mytheme <- theme(axis.title=element_text(size=16, face="bold"), axis.text=element_text(size=14), panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_line(colour = "grey90"), panel.grid.minor = element_blank(), panel.border = element_rect(color = "grey50", fill=NA), plot.margin = unit(c(0.15, 0.3, 0.1, 0.1), "cm"))

Fig <- ggplot(VertebrateCAIS_AAC_df, aes(x=CAIS, y=K))+
  labs( x = "CAIS", y = "Lysine frequency") + 
  geom_point() +
  geom_smooth(method = "lm", alpha = .15) +
  mytheme
  
  
  geom_abline(data=subset(reg, Method=="MA"), aes(intercept=Intercept, slope=Slope), color="red", size=1.5) + 
  annotate(geom="text", x=0.2, y=ypos_lo(dm$Rnoratio), label=paste("R = ", signif(prp$estimate, digits=3), "\n p = ", signif(prp$p.value, digits=1)), color="black", size=6) + 
  geom_pointrange(aes(xmin=Pnoratio-Perror, xmax=Pnoratio+Perror)) + 
  geom_pointrange(aes(ymin=Rnoratio-Rerror, ymax=Rnoratio+Rerror)) + 
  geom_point(size=3) + 
  geom_text(aes(label=aa),hjust=-0.1, vjust=0, color="blue", size = 6) +
  mytheme
  
myggsave_sq <- function(x){
    ggsave(plot=x, filename=paste0(deparse(substitute(x)), ".svg"), path="C:/Users/hanon/Documents/Github/aa_flux/Figures", device="svg", height=5, width=5, units="in")
}

myggsave_sq(Fig)
