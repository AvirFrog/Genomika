import random

le = 1000
tm= {'A': {'A':0.6, 'C':0.1, 'T':0.2, 'G':0.1},
     'C': {'A':0.1, 'C':0.5, 'T':0.1, 'G':0.3},
     'T': {'A':0.4, 'C':0.05, 'T':0.5, 'G':0.05},
     'G': {'A':0.05, 'C':0.2, 'T':0.05, 'G':0.7}}


def random_seq(tm, n):
    seq = [random.choice(list(tm.keys()))]
    for i in range(n-1):
        nuc_out = random.choices(list(tm[seq[i]].keys()), 
                                 list(tm[seq[i]].values()))
        seq.extend(nuc_out)
    return "".join(seq)


seq = random_seq(tm, le)
print(len(seq))

seq_file = open("zad03.fasta", "w")
seq_file.write(f"> Seq zad03 \n{seq}")
seq_file.close()
