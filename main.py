from collections import Counter

from IO import *
from codons import *
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
from scipy.spatial import distance_matrix

mamalian1_list = []
mamalian2_list = []
mamalian3_list = []
mamalian4_list = []

bacterial1_list = []
bacterial2_list = []
bacterial3_list = []
bacterial4_list = []


# function which find all sequences
def findAllSeq(seq, list_seq):
    state = 0
    temp_str = ''

    for i in range(0, len(seq), 3):

        if state == 1:
            temp_str = temp_str + str(seq[i:i + 3])

        if seq[i:i + 3] == 'ATG' and state == 0:
            state = 1
            temp_str = temp_str + str(seq[i:i + 3])

        if ((seq[i:i + 3] == 'TAG' or seq[i:i + 3] == 'TGA' or seq[i:i + 3] == 'TAA') and state == 1):
            state = 0
            list_seq.append(Seq(temp_str))
            temp_str = ''


def findAllSeq1(seq, list_seq):
    state = 0
    temp_str = ''

    for i in range(0, len(seq), 3):

        if state == 1:
            temp_str = temp_str + str(seq[i:i + 3])

        if seq[i:i + 3] == 'ATG' and state == 0:
            state = 1
            temp_str = temp_str + str(seq[i:i + 3])

        if ((seq[i:i + 3] == 'TAG' or seq[i:i + 3] == 'TGA' or seq[i:i + 3] == 'TAA') and state == 1):
            state = 0
            if len(temp_str) > 100:
                list_seq.append(Seq(temp_str))
                temp_str = ''
            else:
                temp_str = ''


mamalian1R = mamalian1.reverse_complement()
mamalian2R = mamalian2.reverse_complement()
mamalian3R = mamalian3.reverse_complement()
mamalian4R = mamalian4.reverse_complement()

bacterial1R = bacterial1.reverse_complement()
bacterial2R = bacterial2.reverse_complement()
bacterial3R = bacterial3.reverse_complement()
bacterial4R = bacterial4.reverse_complement()
# Bacterial sequences
findAllSeq1(bacterial4, bacterial4_list)
findAllSeq1(bacterial4[1:], bacterial4_list)
findAllSeq1(bacterial4[2:], bacterial4_list)

findAllSeq1(bacterial4R, bacterial4_list)
findAllSeq1(bacterial4R[1:], bacterial4_list)
findAllSeq1(bacterial4R[2:], bacterial4_list)

findAllSeq1(bacterial3, bacterial3_list)
findAllSeq1(bacterial3[1:], bacterial3_list)
findAllSeq1(bacterial3[2:0], bacterial3_list)

findAllSeq1(bacterial3R, bacterial3_list)
findAllSeq1(bacterial3R[1:], bacterial3_list)
findAllSeq1(bacterial3R[2:], bacterial3_list)

findAllSeq1(bacterial2, bacterial2_list)
findAllSeq1(bacterial2[1:], bacterial2_list)
findAllSeq1(bacterial2[2:], bacterial2_list)

findAllSeq1(bacterial2R, bacterial2_list)
findAllSeq1(bacterial2R[1:], bacterial2_list)
findAllSeq1(bacterial2R[2:], bacterial2_list)

findAllSeq1(bacterial1, bacterial1_list)
findAllSeq1(bacterial1[1:], bacterial1_list)
findAllSeq1(bacterial1[2:0], bacterial1_list)

findAllSeq1(bacterial1R, bacterial1_list)
findAllSeq1(bacterial1R[1:], bacterial1_list)
findAllSeq1(bacterial1R[2:], bacterial1_list)

findAllSeq1(bacterial4, bacterial4_list)
findAllSeq1(bacterial4[1:], bacterial4_list)
findAllSeq1(bacterial4[2:], bacterial4_list)

findAllSeq1(bacterial4R, bacterial4_list)
findAllSeq1(bacterial4R[1:], bacterial4_list)
findAllSeq1(bacterial4R[2:], bacterial4_list)

findAllSeq1(bacterial3, bacterial3_list)
findAllSeq1(bacterial3[1:], bacterial3_list)
findAllSeq1(bacterial3[2:0], bacterial3_list)

findAllSeq1(bacterial3R, bacterial3_list)
findAllSeq1(bacterial3R[1:], bacterial3_list)
findAllSeq1(bacterial3R[2:], bacterial3_list)

findAllSeq1(bacterial2, bacterial2_list)
findAllSeq1(bacterial2[1:], bacterial2_list)
findAllSeq1(bacterial2[2:], bacterial2_list)

findAllSeq1(bacterial2R, bacterial2_list)
findAllSeq1(bacterial2R[1:], bacterial2_list)
findAllSeq1(bacterial2R[2:], bacterial2_list)

findAllSeq1(bacterial1, bacterial1_list)
findAllSeq1(bacterial1[1:], bacterial1_list)
findAllSeq1(bacterial1[2:0], bacterial1_list)

findAllSeq1(bacterial1R, bacterial1_list)
findAllSeq1(bacterial1R[1:], bacterial1_list)
findAllSeq1(bacterial1R[2:], bacterial1_list)
#Mamalian sequences
#---------------------------------------------------------
findAllSeq1(mamalian1, mamalian1_list)
findAllSeq1(mamalian1[1:], mamalian1_list)
findAllSeq1(mamalian1[2:], mamalian1_list)

findAllSeq1(mamalian1R, mamalian1_list)
findAllSeq1(mamalian1R[1:], mamalian1_list)
findAllSeq1(mamalian1R[2:], mamalian1_list)

findAllSeq1(mamalian3, mamalian3_list)
findAllSeq1(mamalian3[1:], mamalian3_list)
findAllSeq1(mamalian3[2:0], mamalian3_list)

findAllSeq1(mamalian3R, mamalian3_list)
findAllSeq1(mamalian3R[1:], mamalian3_list)
findAllSeq1(mamalian3R[2:], mamalian3_list)

findAllSeq1(mamalian2, mamalian2_list)
findAllSeq1(mamalian2[1:], mamalian2_list)
findAllSeq1(mamalian2[2:], mamalian2_list)

findAllSeq1(mamalian2R, mamalian2_list)
findAllSeq1(mamalian2R[1:], mamalian2_list)
findAllSeq1(mamalian2R[2:], mamalian2_list)

findAllSeq1(mamalian4, mamalian4_list)
findAllSeq1(mamalian4[1:], mamalian4_list)
findAllSeq1(mamalian4[2:0], mamalian4_list)

findAllSeq1(mamalian4R, mamalian4_list)
findAllSeq1(mamalian4R[1:], mamalian4_list)
findAllSeq1(mamalian4R[2:], mamalian4_list)


empty_seq = Seq("")

def concatenate_seq(l, s):
    for item in l:
        s += item
    return s


concatenatedBacterial1 = concatenate_seq(bacterial1_list, empty_seq)
concatenatedBacterial2 = concatenate_seq(bacterial2_list, empty_seq)
concatenatedBacterial3 = concatenate_seq(bacterial3_list, empty_seq)
concatenatedBacterial4 = concatenate_seq(bacterial4_list, empty_seq)
concatenatedMamalian1 = concatenate_seq(mamalian1_list, empty_seq)
concatenatedMamalian2 = concatenate_seq(mamalian2_list, empty_seq)
concatenatedMamalian3 = concatenate_seq(mamalian3_list, empty_seq)
concatenatedMamalian4 = concatenate_seq(mamalian4_list, empty_seq)

def findCodonNumber(list, codon):
    allCodons = 0
    for i in range(0, len(list), 3):
        if (list[i: i + 3] == codon):
            allCodons += 1
    return allCodons


#def find_codons_frequency(seq, dna):
 #   ite = 0
  #  for i in dna:
   #     print(DNA_CODONS[ite] + ' '+ str(
    #        findCodonNumber(seq, DNA_CODONS[ite]) / len(seq) / 3 * 1000))
     #   ite += 1

#find_codons_frequency(concatenatedBacterial1, DNA_CODONS)

aminoBacterial1 = concatenatedBacterial1.translate()
aminoBacterial2 = concatenatedBacterial2.translate()
aminoBacterial3 = concatenatedBacterial3.translate()
aminoBacterial4 = concatenatedBacterial4.translate()
aminoMamalian1 = concatenatedMamalian1.translate()
aminoMamalian2 = concatenatedMamalian2.translate()
aminoMamalian3 = concatenatedMamalian3.translate()
aminoMamalian4 = concatenatedMamalian4.translate()

all_possible_dicodones = []

for x, y in [x + y for x in amino_list for y in amino_list]:
    all_possible_dicodones.append(x+y)

#print(all_possible_dicodones)


split_aminoBacterial1 = [aminoBacterial1[i:i+2] for i in range(0, len(aminoBacterial1), 2)]
split_aminoBacterial2 = [aminoBacterial2[i:i+2] for i in range(0, len(aminoBacterial2), 2)]
split_aminoBacterial3 = [aminoBacterial2[i:i+2] for i in range(0, len(aminoBacterial3), 2)]
split_aminoBacterial4 = [aminoBacterial1[i:i+2] for i in range(0, len(aminoBacterial4), 2)]

split_aminoMamalian1 = [aminoMamalian1[i:i+2] for i in range(0, len(aminoMamalian1), 2)]
split_aminoMamalian2 = [aminoMamalian2[i:i+2] for i in range(0, len(aminoMamalian2), 2)]
split_aminoMamalian3 = [aminoMamalian3[i:i+2] for i in range(0, len(aminoMamalian3), 2)]
split_aminoMamalian4 = [aminoMamalian4[i:i+2] for i in range(0, len(aminoMamalian4), 2)]

count_aminoBacterial1 = Counter(split_aminoBacterial1)
count_aminoBacterial2 = Counter(split_aminoBacterial2)
count_aminoBacterial4 = Counter(split_aminoBacterial4)
count_aminoBacterial3 = Counter(split_aminoBacterial3)
count_aminoMamalian1 = Counter(split_aminoMamalian1)
count_aminoMamalian2 = Counter(split_aminoMamalian2)
count_aminoMamalian3 = Counter(split_aminoMamalian3)
count_aminoMamalian4 = Counter(split_aminoMamalian4)

freq_list_aminoBacterial1 = []
freq_list_aminoBacterial2 = []
freq_list_aminoBacterial3 = []
freq_list_aminoBacterial4 = []

freq_list_aminoMamalian1 = []
freq_list_aminoMamalian2 = []
freq_list_aminoMamalian3 = []
freq_list_aminoMamalian4 = []


def create_list_with_fre(newlist, count, amino):
    for i in range(len(all_possible_dicodones)):
        newlist.append(round((count[str(all_possible_dicodones[i])] / len(amino) * 1000), 4))

create_list_with_fre(freq_list_aminoBacterial1, count_aminoBacterial1, aminoBacterial1)
create_list_with_fre(freq_list_aminoBacterial2, count_aminoBacterial2, aminoBacterial2)
create_list_with_fre(freq_list_aminoBacterial3, count_aminoBacterial3, aminoBacterial3)
create_list_with_fre(freq_list_aminoBacterial4, count_aminoBacterial4, aminoBacterial4)

create_list_with_fre(freq_list_aminoMamalian1, count_aminoMamalian1, aminoMamalian1)
create_list_with_fre(freq_list_aminoMamalian2, count_aminoMamalian2, aminoMamalian2)
create_list_with_fre(freq_list_aminoMamalian3, count_aminoMamalian3, aminoMamalian3)
create_list_with_fre(freq_list_aminoMamalian4, count_aminoMamalian4, aminoMamalian4)

print(freq_list_aminoBacterial1)


freq_of_aminoB12 = []
freq_of_aminoB13 = []
freq_of_aminoB14 = []
freq_of_aminoB23 = []
freq_of_aminoB24 = []
freq_of_aminoB34 = []


freq_of_aminoM1M2 = []
freq_of_aminoM1M3 = []
freq_of_aminoM1M4 = []
freq_of_aminoM2M3 = []
freq_of_aminoM2M4 = []
freq_of_aminoM3M4 = []

freq_of_aminoB1M1 = []
freq_of_aminoB1M2 = []
freq_of_aminoB1M3 = []
freq_of_aminoB1M4 = []

freq_of_aminoB2M1 = []
freq_of_aminoB2M2 = []
freq_of_aminoB2M3 = []
freq_of_aminoB2M4 = []

freq_of_aminoB3M1 = []
freq_of_aminoB3M2 = []
freq_of_aminoB3M3 = []
freq_of_aminoB3M4 = []

freq_of_aminoB4M1 = []
freq_of_aminoB4M2 = []
freq_of_aminoB4M3 = []
freq_of_aminoB4M4 = []

# Create seq list for 2 amino acid frequency
def create_two_seq_list(empty_list, first_seq, second_seq):
    for i in range(len(first_seq)):
        empty_list.append([first_seq[i], second_seq[i]])



create_two_seq_list(freq_of_aminoB3M1, freq_list_aminoBacterial3, freq_list_aminoMamalian1)
create_two_seq_list(freq_of_aminoB3M2, freq_list_aminoBacterial3, freq_list_aminoMamalian2)
create_two_seq_list(freq_of_aminoB3M3, freq_list_aminoBacterial3, freq_list_aminoMamalian3)
create_two_seq_list(freq_of_aminoB3M4, freq_list_aminoBacterial3, freq_list_aminoMamalian4)

met_listy = ''

# Function which create line for each amino acid
def create_line_for_fre(dicodons, lisyt, which):
    for i in range(20):
        tik = abs(dicodons[i+which][0] - dicodons[i+which][1])
        lisyt = lisyt + str(round(tik, 4)) + ' '
    return lisyt

# Function which is for creating Phylip matrix
def function_all(amino, a, w):
    listy = ''
    inter = 0
    inter1 = 0
    while inter != 400:
        at = create_line_for_fre(a, listy, inter)
        print(amino[inter1] + ' ' + at)
        listy = ''
        inter += 20
        inter1 += 1

# Print Phylip format
#B3M1 mean Bacterial3.fasta and Mamalian1.fasta data
print(20)
print(function_all(amino_list_3letters, freq_of_aminoB3M1, 0))
print(20)
print(function_all(amino_list_3letters, freq_of_aminoB3M2, 0))
print(20)
print(function_all(amino_list_3letters, freq_of_aminoB3M3, 0))
print(20)
print(function_all(amino_list_3letters, freq_of_aminoB3M4, 0))


