from Bio import SeqIO

record, = SeqIO.parse("staphylococcus_aureus.fasta", "fasta")

seq = record.seq
print(record)
print(len(seq))