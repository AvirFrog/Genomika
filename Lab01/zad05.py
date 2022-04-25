import random

le = 1000
tm = {'A': {'A': 0.95, 'B': 0.05},
      'B': {'A': 0.1, 'B': 0.9}}

em = {'A': {'A': 0.25, 'C': 0.25, 'T': 0.25, 'G': 0.25},
      'B': {'A': 0.1, 'C': 0.4, 'T': 0.1, 'G': 0.4}}


def random_seq(tm, em, n):
    state = random.choice(list(tm.keys()))
    seq = []
    hidden = [state]
    for i in range(n):
        nuc_out = random.choices(list(em[state].keys()), list(em[state].values()))
        seq.extend(nuc_out)
        state = random.choices(list(tm[state].keys()), list(tm[state].values()))[0]
        hidden.extend(state)
    return [''.join(seq), ''.join(hidden)]


seq = random_seq(tm, em, le)
print(len(seq[0]))

seq_file = open("zad05.fasta", "w")
seq_file.write(f"> Seq zad05 \n{seq[0]}")
seq_file.close()
