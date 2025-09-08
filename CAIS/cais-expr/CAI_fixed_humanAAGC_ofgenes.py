"""

@author: Andrew
"""

import os, sys, json, csv,datetime,math
from Bio import SeqIO
from math import prod

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

Total_AA_freqTable = {'F':0.03619928,
                'L':0.099365156,
                'S':0.083817295,
                'Y':0.026231866,
                '*':0.001714171,
                'C':0.023088572,
                'W':0.012269904,
                'P':0.06389793,
                'H':0.026295268,
                'Q':0.047505033,
                'R':0.056846206,
                'I':0.042770224,
                'M':0.021084171,
                'T':0.05332138,
                'N':0.035485625,
                'K':0.056784013,
                'V':0.059254337,
                'A':0.070226412,
                'D':0.046981582,
                'E':0.070632475,
                'G':0.066229099}
    
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
handle = open("human_genes_100codons.fasta")
for seq_record in SeqIO.parse(handle, "fasta") :
    GC_total_prob = 0.5108
    notGC_total_prob = 1-GC_total_prob


# We'll keep dictionaries of all the codons that we count, and their expected probabilities,
# sorted by which amino acid they correspond to. This will make calculating the CAIS relatively easy and painless- Sara W
    totalCodonCount = 0

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

    RSCUTable = {'F':{'TTT':0,'TTC':0},
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

    RelativeAdaptednessTable = {'F':{'TTT':0,'TTC':0},
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
    codonList = [CodingSequence[n:n+3] for n in range(0, sequenceLength, 3)]

    for AA in RawCount:
        for Codon in RawCount[AA]:
            CodonCount = codonList.count(Codon)
            RawCount[AA][Codon] += CodonCount
            totalCodonCount += CodonCount

    # Once the raw counts are found for each codon, we then start calculating RSCU_i values.
    for AA in RawCount:
        Sum = sum(RawCount[AA].values())
        for Codon in RawCount[AA]:
            # Handling the case where an amino acid never appears
            if Sum != 0:
                RSCU = RawCount[AA][Codon] / Sum
            else:
                RSCU = 1
            RSCUTable[AA][Codon] = RSCU

    # Calculate relative adaptedness values
    for AA in RSCUTable:
        MaxRSCU = 0.3382
        for Codon in RSCUTable[AA]:
            RelativeAdaptedness = RSCUTable[AA][Codon] / MaxRSCU
            RelativeAdaptednessTable[AA][Codon] = RelativeAdaptedness

    # Calculate CAI using a logarithmic transformation
    LogOfCAI = 0
    for AA in RawCount:
        for Codon in RawCount[AA]:
            for k in range(RawCount[AA][Codon]):  # More efficient loop
                LogOfCAI += math.log(RelativeAdaptednessTable[AA][Codon])

    # Compute final CAI value
    LogOfCAI = (1 / totalCodonCount) * LogOfCAI
    CAI = math.exp(LogOfCAI)

    # Print results for each species
    print("%s,%s" % (seq_record.id, CAI))

