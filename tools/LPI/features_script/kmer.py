from functools import reduce
from utils import convert_list_to_str
import os


def get_rna_kmer_fea(ids, seqs, k, save_dir, freq=True):
    feas = []
    elements = reduce(lambda x, y: [i + j for i in x for j in y], [['A', 'T', 'C', 'G']] * k)
    elements_dict = {}

    for seq in seqs:
        for e in elements:
            elements_dict[e] = 0
        for start in range(k):
            for i in range(len(seq)//k):
                tmp = seq[start+k*i: start+k*(i+1)]
                if(len(tmp) == k):
                    elements_dict[seq[start+k*i: start+k*(i+1)]] += 1
        fea = []
        if freq == True:
            counts = sum(elements_dict.values())
            for e in elements:
                fea.append(elements_dict[e]/counts)
        else:
            for e in elements:
                fea.append(elements_dict[e])

        feas.append(fea)

    fea_str = convert_list_to_str(feas)
    fea_file = os.path.join(save_dir, "lncRNA_kmer_features")
    f = open(fea_file, "w")
    for k, v in zip(ids, fea_str):
        f.write(k + " " + v + "\n")
    f.close()


def get_pro_kmer_fea(ids, seqs, k, save_dir, freq=True):
    feas = []
    elements = reduce(lambda x, y: [i + j for i in x for j in y], [['D', 'H', 'C', 'A']] * k)
    elements_dict = {}

    for seq in seqs:
        for e in elements:
            elements_dict[e] = 0
        seq = seq.replace('R', 'H')
        seq = seq.replace('N', 'C')
        seq = seq.replace('Q', 'C')
        seq = seq.replace('E', 'D')
        seq = seq.replace('G', 'C')
        seq = seq.replace('I', 'A')
        seq = seq.replace('L', 'A')
        seq = seq.replace('K', 'H')
        seq = seq.replace('M', 'A')
        seq = seq.replace('F', 'A')
        seq = seq.replace('P', 'A')
        seq = seq.replace('S', 'C')
        seq = seq.replace('T', 'C')
        seq = seq.replace('W', 'A')
        seq = seq.replace('Y', 'C')
        seq = seq.replace('V', 'A')
        for start in range(k):
            for i in range(len(seq)//k):
                tmp = seq[start+k*i: start+k*(i+1)]
                if(len(tmp) == k):
                    elements_dict[seq[start+k*i: start+k*(i+1)]] += 1
        fea = []
        if freq == True:
            counts = sum(elements_dict.values())
            for e in elements:
                fea.append(elements_dict[e]/counts)
        else:
            for e in elements:
                fea.append(elements_dict[e])

        feas.append(fea)

    fea_str = convert_list_to_str(feas)
    fea_file = os.path.join(save_dir, "protein_kmer_features")
    f = open(fea_file, "w")
    for k, v in zip(ids, fea_str):
        f.write(k + " " + v + "\n")
    f.close()