## Read input files
pfam_file <- 	read.csv('AA_for_linearmodeling_species_differences.csv', header = T)
pfam_file <- pfam_file[,-1]
colnames(pfam_file) <- c('SpeciesUID', 'PfamUID','DomainLength','A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
                         'F','P','U','S','T','W','Y','V', 'TmhmmTopology')
head(pfam_file)
names(table(pfam_file$SpeciesUID))

# Vertebrate UID with CAIS values
CAIS_df <- read.csv('CAIS_Values.csv', header = T)

# Pfam ages
pfamages <- read.csv('Pfamages.csv')
# These are the thresholds catherine used
#some pfams in the datafile are considered neither and are thus disregarded 
oldpfams <- pfamages$PfamUID[which(pfamages$Age_oldest >= 2101)]
recentpfams <- pfamages$PfamUID[which(pfamages$Age_oldest <= 1496)]


                       ##### For Non Membrane Pfams #####
#nonTrans_pfam_file <- pfam_file[which(pfam_file$TmhmmTopology<0.5),]
#SpeciesUID_wCAIS <- names(table(nonTrans_pfam_file$SpeciesUID))

# Intialize species AAC % dataframe
#Species_pfam_AAC_df <- data.frame(matrix(ncol = 22))
#colnames(Species_pfam_AAC_df) <- c('A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
                   #          'F','P','U','S','T','W','Y','V')

#for (species in 1:length(SpeciesUID_wCAIS)) {
 # SpeciesProteinRows <- which(nonTrans_pfam_file$SpeciesUID == SpeciesUID_wCAIS[species])
 # print(SpeciesUID_wCAIS[species])
  
  # Initialize protein aa count dataframe for each species
 # ProteinAAcount_df <- data.frame(matrix(0, ncol=22))
 # colnames(ProteinAAcount_df) <- c('A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
  #                                 'F','P','U','S','T','W','Y','V')
  
 # for (proteinrow in 1:length(SpeciesProteinRows)) {
 #   ProteinAAcount <- as.numeric(nonTrans_pfam_file[SpeciesProteinRows[proteinrow],4:25])*
  #    as.numeric(nonTrans_pfam_file$DomainLength[SpeciesProteinRows[proteinrow]])
  #  ProteinAAcount_df <- ProteinAAcount_df + as.numeric(ProteinAAcount)}
  
 # ProteomeLength <- sum(nonTrans_pfam_file$DomainLength[SpeciesProteinRows])
 # SpeciesAAC <- ProteinAAcount_df/ProteomeLength
 # rownames(SpeciesAAC) <- SpeciesUID_wCAIS[species]
 # Species_pfam_AAC_df <- rbind(Species_pfam_AAC_df, SpeciesAAC) }

#Species_pfam_AAC_df <- Species_pfam_AAC_df[-1,]
#print(rowSums(Species_pfam_AAC_df))
#Species_pfam_AAC_df <- as.data.frame(Species_pfam_AAC_df)
#write.csv(Species_pfam_AAC_df, 'NonTransMembranePfam_AAC.csv')



                             ##### For TransMembrane Pfams #####
#Trans_pfam_file <- pfam_file[which(pfam_file$TmhmmTopology>=0.5),]
#SpeciesUID_wCAIS <- names(table(Trans_pfam_file$SpeciesUID))

# Intialize species AAC % dataframe
#Species_pfam_AAC_df <- data.frame(matrix(ncol = 22))
#colnames(Species_pfam_AAC_df) <- c('A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
         #                          'F','P','U','S','T','W','Y','V')

#for (species in 1:length(SpeciesUID_wCAIS)) {
 # SpeciesProteinRows <- which(Trans_pfam_file$SpeciesUID == SpeciesUID_wCAIS[species])
 # print(SpeciesUID_wCAIS[species])
  
  # Initialize protein aa count dataframe for each species
 # ProteinAAcount_df <- data.frame(matrix(0, ncol=22))
 # colnames(ProteinAAcount_df) <- c('A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
                              #     'F','P','U','S','T','W','Y','V')
  
 # for (proteinrow in 1:length(SpeciesProteinRows)) {
  #  ProteinAAcount <- as.numeric(Trans_pfam_file[SpeciesProteinRows[proteinrow],4:25])*
  #    as.numeric(Trans_pfam_file$DomainLength[SpeciesProteinRows[proteinrow]])
  #  ProteinAAcount_df <- ProteinAAcount_df + as.numeric(ProteinAAcount) }
  
 # ProteomeLength <- sum(Trans_pfam_file$DomainLength[SpeciesProteinRows])
  #SpeciesAAC <- ProteinAAcount_df/ProteomeLength
 # rownames(SpeciesAAC) <- SpeciesUID_wCAIS[species]
 # Species_pfam_AAC_df <- rbind(Species_pfam_AAC_df, SpeciesAAC) }

#Species_pfam_AAC_df <- Species_pfam_AAC_df[-1,]
#print(rowSums(Species_pfam_AAC_df))
#Species_pfam_AAC_df <- as.data.frame(Species_pfam_AAC_df)
#write.csv(Species_pfam_AAC_df, 'TransMembranePfam_AAC.csv')


                          ######## Old pfam AAC ########
old_pfam_file <- pfam_file[which(pfam_file$PfamUID %in% oldpfams),]      
SpeciesUID_wCAIS <- names(table(old_pfam_file$SpeciesUID))

# Intialize species AAC % dataframe
Species_pfam_AAC_df <- data.frame(matrix(ncol = 22))
colnames(Species_pfam_AAC_df) <- c('A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
                                   'F','P','U','S','T','W','Y','V')

for (species in 1:length(SpeciesUID_wCAIS)) {
  SpeciesProteinRows <- which(old_pfam_file$SpeciesUID == SpeciesUID_wCAIS[species])
  print(SpeciesUID_wCAIS[species])
  
  # Initialize protein aa count dataframe for each species
  ProteinAAcount_df <- data.frame(matrix(0, ncol=22))
  colnames(ProteinAAcount_df) <- c('A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
                                   'F','P','U','S','T','W','Y','V')
  
  for (proteinrow in 1:length(SpeciesProteinRows)) {
    if(any(is.na(old_pfam_file[SpeciesProteinRows[proteinrow],4:25]))) {next()}
    ProteinAAcount <- as.numeric(old_pfam_file[SpeciesProteinRows[proteinrow],4:25])*
       as.numeric(old_pfam_file$DomainLength[SpeciesProteinRows[proteinrow]])
    ProteinAAcount_df <- ProteinAAcount_df + as.numeric(ProteinAAcount) }
   
  ProteomeLength <- sum(old_pfam_file$DomainLength[SpeciesProteinRows])
  SpeciesAAC <- ProteinAAcount_df/ProteomeLength
  rownames(SpeciesAAC) <- SpeciesUID_wCAIS[species]
  Species_pfam_AAC_df <- rbind(Species_pfam_AAC_df, SpeciesAAC) }

Species_pfam_AAC_df <- Species_pfam_AAC_df[-1,]
print(rowSums(Species_pfam_AAC_df))
Species_pfam_AAC_df <- as.data.frame(Species_pfam_AAC_df)
write.csv(Species_pfam_AAC_df, 'Old_Pfam_AAC.csv')

                    ######## recent pfam AAC ########
recent_pfam_file <- pfam_file[which(pfam_file$PfamUID %in% recentpfams),]      
SpeciesUID_wCAIS <- names(table(recent_pfam_file$SpeciesUID))

# Intialize species AAC % dataframe
Species_pfam_AAC_df <- data.frame(matrix(ncol = 22))
colnames(Species_pfam_AAC_df) <- c('A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
                                   'F','P','U','S','T','W','Y','V')

for (species in 1:length(SpeciesUID_wCAIS)) {
  SpeciesProteinRows <- which(recent_pfam_file$SpeciesUID == SpeciesUID_wCAIS[species])
  print(SpeciesUID_wCAIS[species])
  
  # Initialize protein aa count dataframe for each species
  ProteinAAcount_df <- data.frame(matrix(0, ncol=22))
  colnames(ProteinAAcount_df) <- c('A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
                                   'F','P','U','S','T','W','Y','V')
  
  for (proteinrow in 1:length(SpeciesProteinRows)) {
    if(any(is.na(recent_pfam_file[SpeciesProteinRows[proteinrow],4:25]))){next()}
    ProteinAAcount <- as.numeric(recent_pfam_file[SpeciesProteinRows[proteinrow],4:25])*
      as.numeric(recent_pfam_file$DomainLength[SpeciesProteinRows[proteinrow]])
    ProteinAAcount_df <- ProteinAAcount_df + as.numeric(ProteinAAcount) }
  
  ProteomeLength <- sum(recent_pfam_file$DomainLength[SpeciesProteinRows])
  SpeciesAAC <- ProteinAAcount_df/ProteomeLength
  rownames(SpeciesAAC) <- SpeciesUID_wCAIS[species]
  Species_pfam_AAC_df <- rbind(Species_pfam_AAC_df, SpeciesAAC) }

Species_pfam_AAC_df <- Species_pfam_AAC_df[-1,]
print(rowSums(Species_pfam_AAC_df))
Species_pfam_AAC_df <- as.data.frame(Species_pfam_AAC_df)
write.csv(Species_pfam_AAC_df, 'Recent_Pfam_AAC.csv')


# check why some species have NA AAC
#ProteinAAcount_407 <- data.frame(matrix(0, ncol=22))
#colnames(ProteinAAcount_407) <- c('A','R','N','D','C','E','Q','G','H','O','I','L','K','M',
                     #             'F','P','U','S','T','W','Y','V')
#species407 <- recent_pfam_file[which(recent_pfam_file$SpeciesUID == 407),]

#for (proteinrow in 1:nrow(species407)) {
 # if(any(is.na(species407[proteinrow,4:25]))){next()}
 # ProteinAAcount <- as.numeric(species407[proteinrow,4:25])*
 #   as.numeric(species407$DomainLength[proteinrow])
 # ProteinAAcount_407 <- ProteinAAcount_407 + as.numeric(ProteinAAcount) }
