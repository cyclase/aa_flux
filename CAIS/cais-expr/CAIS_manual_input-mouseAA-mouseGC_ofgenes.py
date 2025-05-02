"""

@author: Andrew
"""

import os, sys, json, csv, mysql.connector,datetime,math
from Bio import SeqIO

"""
The purpose of this file is to determine the codon adaptation index of species (CAIS) for individual species.
This is a metric that describes codon bias patterns. 
CAIS is calculated in the following way:
    
    I) Calculate the probability of a codon given the genomic GC content and amino acid frequencies (null expectation)
    
    II) Count the number of each codon in a sequence
        
    II) Calculate CAIS metric for the species
         
"""

maxInt = sys.maxsize

while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

####################################################################################
#                             Program Executes Below                               #
####################################################################################


#########For Total Dataset metrics
    
#The amino acid frequency values in this table are calculated in another script (WHICH TAKES & HOURS TO RUN)
#USED FOR THE RE_WEIGHTING OF CAIS

Total_AA_freqTable = {'F':0.037743701,
                'L':0.10050291,
                'S':0.084741159,
                'Y':0.027227893,
                '*':0.001834864,
                'C':0.023776481,
                'W':0.012125743,
                'P':0.061358446,
                'H':0.026675899,
                'Q':0.047321992,
                'R':0.055683056,
                'I':0.044144741,
                'M':0.022081637,
                'T':0.054338872,
                'N':0.035681786,
                'K':0.056543419,
                'V':0.060763834,
                'A':0.068168195,
                'D':0.047138947,
                'E':0.068143293,
                'G':0.064003134}
    
#genome wide GC content for each species
#Species_Total_GC_content = {'bin1':0.5169,
#                            'bin2':0.5169,
#                            'bin3':0.5169,
#                            'bin4':0.5169,
#                            'bin5':0.5169,
#                            'bin6':0.5169,
#                            'bin7':0.5169,
#                            'bin8':0.5169,
#                            'bin9':0.5169,
#                            'bin10':0.5169}


############################ACTUAL CAIS CALCULATING CODE#####################################

#input fasta file name
handle = open("mouse_genes_100codons.fasta")
for seq_record in SeqIO.parse(handle, "fasta") :
    GC_total_prob = 0.5169
    notGC_total_prob = 1-GC_total_prob


# We'll keep dictionaries of all the codons that we count, and their expected probabilities,
# sorted by which amino acid they correspond to. This will make calculating the CAIS relatively easy and painless- Sara W

    RawCount = {'F':{'TTT':0,'TTC':0},
                'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
                'Y':{'TAT':0,'TAC':0},
                '*':{'TAA':0,'TAG':0,'TGA':0},
                'C':{'TGT':0,'TGC':0},
                'W':{'TGG':0},
                'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                'H':{'CAT':0,'CAC':0},
                'Q':{'CAA':0,'CAG':0},
                'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
                'I':{'ATT':0,'ATC':0,'ATA':0},
                'M':{'ATG':0},
                'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                'N':{'AAT':0,'AAC':0},
                'K':{'AAA':0,'AAG':0},
                'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                'D':{'GAT':0,'GAC':0},
                'E':{'GAA':0,'GAG':0},
                'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}

    Sum = {'F':0,
           'L':0,
           'S':0,
           'Y':0,
           '*':0,
           'C':0,
           'W':0,
           'P':0,
           'H':0,
           'Q':0,
           'R':0,
           'I':0,
           'M':0,
           'T':0,
           'N':0,
           'K':0,
           'V':0,
           'A':0,
           'D':0,
           'E':0,
           'G':0}

    ProbTable = {'F':{'TTT':0,'TTC':0},
                'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
                'Y':{'TAT':0,'TAC':0},
                '*':{'TAA':0,'TAG':0,'TGA':0},
                'C':{'TGT':0,'TGC':0},
                'W':{'TGG':0},
                'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                'H':{'CAT':0,'CAC':0},
                'Q':{'CAA':0,'CAG':0},
                'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
                'I':{'ATT':0,'ATC':0,'ATA':0},
                'M':{'ATG':0},
                'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                'N':{'AAT':0,'AAC':0},
                'K':{'AAA':0,'AAG':0},
                'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                'D':{'GAT':0,'GAC':0},
                'E':{'GAA':0,'GAG':0},
                'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}
    Sum_ProbTable = {'F':0,
           'L':0,
           'S':0,
           'Y':0,
           '*':0,
           'C':0,
           'W':0,
           'P':0,
           'H':0,
           'Q':0,
           'R':0,
           'I':0,
           'M':0,
           'T':0,
           'N':0,
           'K':0,
           'V':0,
           'A':0,
           'D':0,
           'E':0,
           'G':0}

    Expected = {'F':{'TTT':0,'TTC':0},
                'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
                'Y':{'TAT':0,'TAC':0},
                '*':{'TAA':0,'TAG':0,'TGA':0},
                'C':{'TGT':0,'TGC':0},
                'W':{'TGG':0},
                'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                'H':{'CAT':0,'CAC':0},
                'Q':{'CAA':0,'CAG':0},
                'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
                'I':{'ATT':0,'ATC':0,'ATA':0},
                'M':{'ATG':0},
                'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                'N':{'AAT':0,'AAC':0},
                'K':{'AAA':0,'AAG':0},
                'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                'D':{'GAT':0,'GAC':0},
                'E':{'GAA':0,'GAG':0},
                'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}

    Observed = {'F':{'TTT':0,'TTC':0},
                'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
                'Y':{'TAT':0,'TAC':0},
                '*':{'TAA':0,'TAG':0,'TGA':0},
                'C':{'TGT':0,'TGC':0},
                'W':{'TGG':0},
                'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                'H':{'CAT':0,'CAC':0},
                'Q':{'CAA':0,'CAG':0},
                'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
                'I':{'ATT':0,'ATC':0,'ATA':0},
                'M':{'ATG':0},
                'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                'N':{'AAT':0,'AAC':0},
                'K':{'AAA':0,'AAG':0},
                'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                'D':{'GAT':0,'GAC':0},
                'E':{'GAA':0,'GAG':0},
                'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}


######################################

    CodingSequence = seq_record.seq
    sequenceLength = len(seq_record)
    codonList = [CodingSequence[n:n+3] for n in range(0,sequenceLength,3)]

    for AA in RawCount:
        for Codon in RawCount[AA]:
            CodonCount = codonList.count(Codon)
            RawCount[AA][Codon] += CodonCount

    # Calculate expected distribution given GC content and amino acid frequency
    Summed_codon_instances = 0 

    for AA in RawCount:
        Sum[AA] = sum(RawCount[AA].values())
        for Codon in RawCount[AA]:
            if Codon in ['TTA','TAT','ATT','AAT','ATA','TAA','AAA','TTT']:
                Prob = (notGC_total_prob*notGC_total_prob*notGC_total_prob)*0.125
            elif Codon in ['GGG','CCC','GCG','CGC','GGC','CCG','CGG','GCC']:
                Prob = (GC_total_prob*GC_total_prob*GC_total_prob)*0.125
            elif Codon in ['GAT','CAT','AGT','ACT','ATG','ATC','GTA','CTA','TGA','TCA','TAG','TAC','TTC','AAC','TTG','AAG','ACA','AGA','TCT','TGT','GTT','CTT','GAA','CAA']:
                Prob = (GC_total_prob*notGC_total_prob*notGC_total_prob)*0.125
            elif Codon in ['GGA','GGT','CCA','CCT','GAG','GTG','CAC','CTC','AGC','TGC','ACG','GCT','AGG','TGG','ACC','GAC','CAG','GTC','CTG','TCC','TCG','CGT','CGA','GCA']:
                Prob = (GC_total_prob*GC_total_prob*notGC_total_prob)*0.125
            else:
                print("storm the castle")
            ProbTable[AA][Codon] = Prob

    for AA in ProbTable:
        Sum_ProbTable[AA]= sum(ProbTable[AA].values())
        for Codon in ProbTable[AA]:
            Expected[AA][Codon] = ProbTable[AA][Codon]/Sum_ProbTable[AA]
                    
    #calculate observed codon frequencies
    for AA in RawCount:
        if Sum[AA] == 0:  # Skip amino acids that are no present in this gene
            continue
        for Codon in RawCount[AA]:
            Observed[AA][Codon] = RawCount[AA][Codon] / Sum[AA] if Sum[AA] > 0 else 0


    #CAIS
    CAIS = 0
    for AA in RawCount:
        if Sum[AA] == 0: #skip amino acids that are not present in this gene
            continue
        aa_weight = Total_AA_freqTable[AA]/sum(Observed[AA].values())
        for Codon in RawCount[AA]:
            if Observed[AA][Codon] == 0 or Expected[AA][Codon] == 0: #skip codons that aren't present in this gene
                continue
            CAIS += Observed[AA][Codon]*math.log(Observed[AA][Codon]/Expected[AA][Codon])*aa_weight

    #CAIS_unweighted
    CAIS_unweighted = 0
    for AA in RawCount:
        if Sum[AA] == 0: #skip amino acids that are not present in this gene
            continue
        for Codon in RawCount[AA]:
            if Observed[AA][Codon] == 0 or Expected[AA][Codon] == 0: #skip codons that aren't present in this gene
                continue
            CAIS_unweighted += Observed[AA][Codon]*math.log(Observed[AA][Codon]/Expected[AA][Codon])
       
    #Print results for each species
    print("%s,%s,%s"%(seq_record.id,CAIS,CAIS_unweighted))
    
