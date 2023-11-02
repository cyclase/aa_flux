library(stringr)
library(protr)
library(ggplot2)
library(seqinr)
library(ggrepel)

setwd('/Users/sawsanwehbi/Desktop/Masel Lab')

# Read input files 
#Species UID 442 is the tibetan antelope; removed bcz its contaminated
CAIS_values <- read.csv('new_CAIS.csv',header = T) 
#CAIS_values <- read.csv('Rotation/CAIS_values.csv', header = T)
CAIS_values <-  CAIS_values[-which(CAIS_values$SpeciesUID == '442'),]
CAIS_values <- CAIS_values[which(CAIS_values$SpeciesUID %in% SpeciesAAC_wCAIS$SpeciesUID),]
NonTrans_file <- read.csv('NonTransMembranePfam_AAC.csv', header = T)
NonTrans_file <- NonTrans_file[-which(NonTrans_file$X == '442'),]
Trans_file <- read.csv('TransMembranePfam_AAC.csv', header = T)
Trans_file <- Trans_file[-which(Trans_file$X == '442'),]
Old_file <- read.csv('Old_Pfam_AAC.csv', header = T)
Old_file <- Old_file[-which(Old_file$X == '442'),]
Young_file <- read.csv('Recent_Pfam_AAC.csv', header = T)
Young_file <- Young_file[-which(Young_file$X == '442'),]

                  ###### CAIS effect on NonTransMembrane #####
# Intialize dataframe
AAstart <- which(colnames(NonTrans_file) == 'A')
AAend <- which(colnames(NonTrans_file) == 'V')
NonTrans_model_estimates <- data.frame(matrix(nrow = 2))


# loop model over all amino acids
for (aa in AAstart:AAend) {
  
  model <- lm(NonTrans_file[,aa]~CAIS_values$CAIS)
  model_summary <- summary(model)
  model_summary <- as.data.frame(model_summary$coefficients[2,1:2])
  colnames(model_summary) <- colnames(NonTrans_file)[aa]
  NonTrans_model_estimates <- cbind(NonTrans_model_estimates, model_summary) }


NonTrans_model_estimates <- NonTrans_model_estimates[,-1]
rownames(NonTrans_model_estimates) <- c('NonTrans_Estimates','NonTrans_Std.Error')
NonTrans_model_estimates <- t(NonTrans_model_estimates)

# Remove O and U
NonTrans_model_estimates <- NonTrans_model_estimates[-which(NonTrans_model_estimates[,1]== 0),]
NonTrans_model_estimates <- data.frame(NonTrans_model_estimates)



                ###### CAIS effect on TransMembrane #####
# Intialize dataframe
AAstart <- which(colnames(Trans_file) == 'A')
AAend <- which(colnames(Trans_file) == 'V')
Trans_model_estimates <- data.frame(matrix(nrow = 2))


# loop model over all amino acids
for (aa in AAstart:AAend) {
  
  model <- lm(Trans_file[,aa]~CAIS_values$CAIS)
  model_summary <- summary(model)
  model_summary <- as.data.frame(model_summary$coefficients[2,1:2])
  colnames(model_summary) <- colnames(Trans_file)[aa]
  Trans_model_estimates <- cbind(Trans_model_estimates, model_summary) }


Trans_model_estimates <- Trans_model_estimates[,-1]
rownames(Trans_model_estimates) <- c('Trans_Estimates','Trans_Std.Error')
Trans_model_estimates <- t(Trans_model_estimates)

# Remove O and U
Trans_model_estimates <- Trans_model_estimates[-which(Trans_model_estimates[,1]== 0),]
Trans_model_estimates <- data.frame(Trans_model_estimates)



                  ###### CAIS effect on Old#####
# Intialize dataframe
AAstart <- which(colnames(Old_file) == 'A')
AAend <- which(colnames(Old_file) == 'V')
Old_model_estimates <- data.frame(matrix(nrow = 2))


# loop model over all amino acids
for (aa in AAstart:AAend) {
  
  model <- lm(Old_file[,aa]~CAIS_values$CAIS)
  model_summary <- summary(model)
  model_summary <- as.data.frame(model_summary$coefficients[2,1:2])
  colnames(model_summary) <- colnames(Old_file )[aa]
  Old_model_estimates <- cbind(Old_model_estimates, model_summary) }


Old_model_estimates <- Old_model_estimates[,-1]
rownames(Old_model_estimates) <- c('Old_Estimates','Old_Std.Error')
Old_model_estimates <- t(Old_model_estimates)

# Remove O and U
Old_model_estimates <- Old_model_estimates[-which(Old_model_estimates[,1]== 0),]
Old_model_estimates <- data.frame(Old_model_estimates)



                     ###### CAIS effect on Young#####
# Intialize dataframe
AAstart <- which(colnames(Young_file) == 'A')
AAend <- which(colnames(Young_file) == 'V')
Young_model_estimates <- data.frame(matrix(nrow = 2))


# loop model over all amino acids
for (aa in AAstart:AAend) {
  
  model <- lm(Young_file[,aa]~CAIS_values$CAIS)
  model_summary <- summary(model)
  model_summary <- as.data.frame(model_summary$coefficients[2,1:2])
  colnames(model_summary) <- colnames(Young_file)[aa]
  Young_model_estimates <- cbind(Young_model_estimates, model_summary) }


Young_model_estimates <- Young_model_estimates[,-1]
rownames(Young_model_estimates) <- c('Young_Estimates','Young_Std.Error')
Young_model_estimates <- t(Young_model_estimates)

# Remove O and U
Young_model_estimates <- Young_model_estimates[-which(Young_model_estimates[,1]== 0),]
Young_model_estimates <- data.frame(Young_model_estimates)


# write all 4 model estimates into csv file
PfamEstimates <- cbind(NonTrans_model_estimates,Trans_model_estimates,
                       Old_model_estimates,Young_model_estimates)
write.csv(PfamEstimates,'CAISeffectonPfamAAC_Trans_NonTrans_Old_Young_newCAIS.csv')


                            ##### FIGURE 6 ######
            ######## Plot TransMembrane vs NonTransmembrane #######
plot_df <- cbind(NonTrans_model_estimates,Trans_model_estimates)

quartz()
plot <- ggplot(data= plot_df, aes(NonTrans_Estimates,Trans_Estimates)) + 
  geom_point(size=2) + labs(x="Effect of CAIS on amino acid frequency\n in Nontransmembrane Pfams", 
                            y="Effect of CAIS on amino acid frequency\n in Transmembrane Pfams")+
  geom_smooth(method = "lm", color="red")
plot <- plot + theme(axis.title=element_text(size=16, face ="bold"),title=element_text(size=8, face ="bold"),
                     panel.background = element_rect(fill = "white", colour = "grey50"),
                     panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))


cor.test(plot_df$Trans_Estimates,plot_df$NonTrans_Estimates)
plot <- plot +
  geom_text(data = plot_df, aes(x = 0.02, y = -0.07, 
                                label = "Pearson's R:0.88\n p-value:1.6e-07"), 
            colour = 'black', size = 6) + 
  geom_errorbar(aes(ymin = Trans_Estimates - Trans_Std.Error, ymax = Trans_Estimates + Trans_Std.Error)) +
  geom_errorbarh(aes(xmin = NonTrans_Estimates - NonTrans_Std.Error, xmax = NonTrans_Estimates + NonTrans_Std.Error)) +
  geom_text_repel(aes(x=NonTrans_Estimates, y=Trans_Estimates, label=rownames(plot_df)),
                  colour = 'blue', max.overlaps = 20)
plot


            ######## Plot Old vs Young #######
plot_df <- cbind(Old_model_estimates,Young_model_estimates)

quartz()
plot <- ggplot(data= plot_df, aes(Old_Estimates,Young_Estimates)) + 
  geom_point(size=2) + labs(x="Effect of CAIS on amino acid frequency\n in Old Pfams", 
                            y="Effect of CAIS on amino acid frequency\n in Young Pfams")+
  geom_smooth(method = "lm", color="red")
plot <- plot + theme(axis.title=element_text(size=16, face ="bold"),title=element_text(size=8, face ="bold"),
                     panel.background = element_rect(fill = "white", colour = "grey50"),
                     panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))


cor.test(plot_df$Old_Estimates,plot_df$Young_Estimates)
plot <- plot +
  geom_text(data = plot_df, aes(x = 0.025, y = -0.15, 
                                label = "Pearson's R:0.84\n p-value:3.2e-06"), 
            colour = 'black', size = 6) + 
  geom_errorbar(aes(ymin = Young_Estimates - Young_Std.Error, ymax = Young_Estimates + Young_Std.Error)) +
  geom_errorbarh(aes(xmin = Old_Estimates - Old_Std.Error, xmax = Old_Estimates + Old_Std.Error)) +
  geom_text_repel(aes(x=Old_Estimates, y=Young_Estimates, label=rownames(plot_df)),
                  colour = 'blue', max.overlaps = 20)

