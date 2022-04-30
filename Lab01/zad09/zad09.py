import random
import numpy as np
import matplotlib.pyplot as plt


def random_seq(tm, n):
      seq = [random.choice(list(tm.keys()))]

      for i in range(n - 1):
            seq.extend(random.choices(list(tm[seq[i]].keys()), list(tm[seq[i]].values())))
      return ''.join(seq)


def kmer_freq(seq, k=2):
      length = len(seq)
      start = 0
      stop = k
      kmer_dic = {}
      while(stop <= length):        ### right border exception
            subseq = seq[start:stop]
            if subseq in kmer_dic:
                  kmer_dic[subseq] += 1
            else:
                  kmer_dic[subseq] = 1
            start += 1
            stop += 1
      return kmer_dic


def bootstrap(seq, reps, kmer):
      kmer_out = {kmer:[]}
      for rep in range(1,reps):
            out = ""
            for i in(range(1,len(seq))):
                  idx = random.randint(0,len(seq))
                  out += seq[idx:idx+1]
            kmer_out[kmer].append(kmer_freq(out)[kmer])
      return kmer_out


def perm(seq, reps, kmer):
      kmer_out = {kmer:[]}
      print(seq)
      for rep in range(1,reps):
            out = "".join(random.sample(seq, len(seq)))
            kmer_out[kmer].append(kmer_freq(out)[kmer])
      return kmer_out


random.seed(0)

length = 1000

tm = {'A': {'A': 0.6, 'C': 0.1, 'T': 0.2, 'G': 0.1},
      'C': {'A': 0.1, 'C': 0.5, 'T': 0.1, 'G': 0.3},
      'T': {'A': 0.4, 'C': 0.05, 'T': 0.5, 'G': 0.05},
      'G': {'A': 0.05, 'C': 0.2, 'T': 0.05, 'G': 0.7}}

seq = random_seq(tm, length)
seq = seq.lower()

kmers = kmer_freq(seq, k=2)
max_kmer = sorted(kmers.items(), key=lambda n: n[1], reverse=True)[0]

bs = perm(seq, 1000, max_kmer[0])
print(bs)

l_ci = np.quantile(bs[max_kmer[0]], 0.025)
l_hi = np.quantile(bs[max_kmer[0]], 0.975)

print((l_ci, l_hi))
print(kmer_freq(seq, k=2)[max_kmer[0]])

plt.hist(bs[max_kmer[0]])
plt.axvline(x=l_ci, color = "black", linestyle = ':')
plt.axvline(x=l_hi, color = "black", linestyle = ':')
plt.axvline(x=kmer_freq(seq, k=2)[max_kmer[0]], color = "orange", linestyle = '--')
plt.show()
