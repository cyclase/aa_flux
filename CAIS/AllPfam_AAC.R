# set library path to my home directory where packages are installed
myPaths <- .libPaths()
myPaths <- c(myPaths, "/home/u27/sawsanwehbi/R/x86_64-pc-linux-gnu-library/4.0")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)  # add new path

## Read input files
pfam_file <- 	read.csv('AA_for_linearmodeling_species_differences.csv', header = T)
pfam_file <- pfam_file[,-1]
colnames(pfam_file) <- c('SpeciesUID', 'PfamUID','DomainLength','A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
                         'F','P','U','S','T','W','Y','V', 'TmhmmTopology')
head(pfam_file)
names(table(pfam_file$SpeciesUID))

# Vertebrate UID with CAIS values
CAIS_df <- read.csv('new_CAIS.csv', header = T)

##### For All Pfams #####
SpeciesUID_wCAIS <- names(table(pfam_file$SpeciesUID))

# Intialize species AAC % dataframe
Species_pfam_AAC_df <- data.frame(matrix(ncol = 22))
colnames(Species_pfam_AAC_df) <- c('A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
          'F','P','U','S','T','W','Y','V')

for (species in 1:length(SpeciesUID_wCAIS)) {
 SpeciesProteinRows <- which(pfam_file$SpeciesUID == SpeciesUID_wCAIS[species])
 print(SpeciesUID_wCAIS[species])

 #Initialize protein aa count dataframe for each species
 ProteinAAcount_df <- data.frame(matrix(0, ncol=22))
 colnames(ProteinAAcount_df) <- c('A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
                                 'F','P','U','S','T','W','Y','V')

for (proteinrow in 1:length(SpeciesProteinRows)) {
   ProteinAAcount <- as.numeric(pfam_file[SpeciesProteinRows[proteinrow],4:25])*
    as.numeric(pfam_file$DomainLength[SpeciesProteinRows[proteinrow]])
   if (any(is.na(ProteinAAcount))) {next}
  ProteinAAcount_df <- ProteinAAcount_df + as.numeric(ProteinAAcount) }

 ProteomeLength <- sum(pfam_file$DomainLength[SpeciesProteinRows])
 SpeciesAAC <- ProteinAAcount_df/ProteomeLength
 rownames(SpeciesAAC) <- SpeciesUID_wCAIS[species]
 Species_pfam_AAC_df <- rbind(Species_pfam_AAC_df, SpeciesAAC) }

Species_pfam_AAC_df <- Species_pfam_AAC_df[-1,]
print(rowSums(Species_pfam_AAC_df))
Species_pfam_AAC_df <- as.data.frame(Species_pfam_AAC_df)
write.csv(Species_pfam_AAC_df, 'AllPfam_AAC.csv')

