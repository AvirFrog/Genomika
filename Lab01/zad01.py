import random
nucleotides = ['A', 'T', 'C', 'G']
length = 1000


def random_seq(nucleotide, n):
	random_set = [random.choice(nucleotide) for i in range(n)]
	return ''.join(random_set)


seq = random_seq(nucleotides, length)
print(len(seq))

seq_file = open("zad01.fasta", "w")
seq_file.write("> Seq zad01 \n")
seq_file.write(seq)

seq_file.close()
