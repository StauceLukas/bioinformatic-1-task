from Bio import SeqIO
from Bio.Seq import Seq

for record in SeqIO.parse("data/bacterial1.fasta", "fasta"):
    bacterial1 = record.seq

for record in SeqIO.parse("data/bacterial2.fasta", "fasta"):
    bacterial2 = record.seq

for record in SeqIO.parse("data/bacterial3.fasta", "fasta"):
    bacterial3 = record.seq

for record in SeqIO.parse("data/bacterial4.fasta", "fasta"):
    bacterial4 = record.seq

for record in SeqIO.parse("data/mamalian1.fasta", "fasta"):
    mamalian1 = record.seq

for record in SeqIO.parse("data/mamalian2.fasta", "fasta"):
    mamalian2 = record.seq

for record in SeqIO.parse("data/mamalian3.fasta", "fasta"):
    mamalian3 = record.seq

for record in SeqIO.parse("data/mamalian4.fasta", "fasta"):
    mamalian4 = record.seq
