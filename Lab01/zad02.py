import random
nucleotides = {'A': 0.1, 'C': 0.3, 'G': 0.4, 'T': 0.2}
length = 1000


def random_seq(nucleotides, n):
    random_seq = random.choices(list(nucleotides.keys()), list(nucleotides.values()), k=n)
    return "".join(random_seq)


seq = random_seq(nucleotides, length)
print(len(seq))

seq_file = open("zad02.fasta", "w")
seq_file.write(f"> Seq zad02 \n{seq}")
seq_file.close()


