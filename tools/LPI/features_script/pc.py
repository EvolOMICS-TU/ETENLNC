import math
import os
from utils import convert_list_to_str


def fourier(numseq_list, fourier_len=10):
    numseq_len = len(numseq_list)
    out_put = []
    for k in range(fourier_len):
        tmp = 0
        for n in range(numseq_len):
            tmp += numseq_list[n]*math.cos(math.pi*(n+0.5)*(k+0.5)/numseq_len)
        tmp *= math.sqrt(2/numseq_len)
        out_put.append(tmp)
    return out_put


def get_hydrogenbonding_fea(seq, fourier_len=10):
    seq = seq.replace('A', ' 24 17 40 10 ')
    seq = seq.replace('C', ' 49 21 26 ')
    seq = seq.replace('G', ' 21 86 17 41 29 ')
    seq = seq.replace('T', ' 24 17 22 ')
    numseq_list = [float(x) for x in seq.split()]
    out_put = fourier(numseq_list, fourier_len)
    return out_put


def get_vanderwall_fea(seq, fourier_len=10):
    seq = seq.replace('A', ' 79 98 69 40 53 37 84 62 49 28 ')
    seq = seq.replace('C', ' 14 44 98 42 30 50 39 19 ')
    seq = seq.replace('G', ' 26 74 24 37 22 21 19 67 48 44 21 ')
    seq = seq.replace('T', ' 25 42 74 53 43 67 44 24 ')
    numseq_list = [float(x) for x in seq.split()]
    out_put = fourier(numseq_list, fourier_len)
    return out_put


def get_HphobBullbreese_fea(seq, fourier_len=10):
    seq = seq.replace('A', ' 0.61 ')
    seq = seq.replace('R', ' 0.69 ')
    seq = seq.replace('N', ' 0.89 ')
    seq = seq.replace('D', ' 0.61 ')
    seq = seq.replace('C', ' 0.36 ')
    seq = seq.replace('Q', ' 0.97 ')
    seq = seq.replace('E', ' 0.51 ')
    seq = seq.replace('G', ' 0.81 ')
    seq = seq.replace('H', ' 0.69 ')
    seq = seq.replace('I', ' -1.45 ')
    seq = seq.replace('L', ' -1.65 ')
    seq = seq.replace('K', ' 0.46 ')
    seq = seq.replace('M', ' -0.66 ')
    seq = seq.replace('F', ' -1.52 ')
    seq = seq.replace('P', ' -0.17 ')
    seq = seq.replace('S', ' 0.42 ')
    seq = seq.replace('T', ' 0.29 ')
    seq = seq.replace('W', ' -1.20 ')
    seq = seq.replace('Y', ' -1.43 ')
    seq = seq.replace('V', ' -0.75 ')
    numseq_list = [float(x) for x in seq.split()]
    out_put = fourier(numseq_list, fourier_len)
    return out_put


def get_PolarityGrantham_fea(seq, fourier_len=10):
    seq = seq.replace('A', ' 8.1 ')
    seq = seq.replace('R', ' 10.5 ')
    seq = seq.replace('N', ' 11.6 ')
    seq = seq.replace('D', ' 13 ')
    seq = seq.replace('C', ' 5.5 ')
    seq = seq.replace('Q', ' 10.5 ')
    seq = seq.replace('E', ' 12.3 ')
    seq = seq.replace('G', ' 9 ')
    seq = seq.replace('H', ' 10.4 ')
    seq = seq.replace('I', ' 5.2 ')
    seq = seq.replace('L', ' 4.9 ')
    seq = seq.replace('K', ' 11.3 ')
    seq = seq.replace('M', ' 5.7 ')
    seq = seq.replace('F', ' 5.2 ')
    seq = seq.replace('P', ' 8 ')
    seq = seq.replace('S', ' 9.2 ')
    seq = seq.replace('T', ' 8.6 ')
    seq = seq.replace('W', ' 5.4 ')
    seq = seq.replace('Y', ' 6.2 ')
    seq = seq.replace('V', ' 5.9 ')
    numseq_list = [float(x) for x in seq.split()]
    out_put = fourier(numseq_list, fourier_len)
    return out_put


def get_PolarityZimmerman_fea(seq, fourier_len=10):
    seq = seq.replace('A', ' 0 ')
    seq = seq.replace('R', ' 52 ')
    seq = seq.replace('N', ' 3.38 ')
    seq = seq.replace('D', ' 49.7 ')
    seq = seq.replace('C', ' 1.48 ')
    seq = seq.replace('Q', ' 3.53 ')
    seq = seq.replace('E', ' 49.9 ')
    seq = seq.replace('G', ' 0 ')
    seq = seq.replace('H', ' 51.6 ')
    seq = seq.replace('I', ' 0.13 ')
    seq = seq.replace('L', ' 0.13 ')
    seq = seq.replace('K', ' 49.5 ')
    seq = seq.replace('M', ' 1.43 ')
    seq = seq.replace('F', ' 0.35 ')
    seq = seq.replace('P', ' 1.58 ')
    seq = seq.replace('S', ' 1.67 ')
    seq = seq.replace('T', ' 1.66 ')
    seq = seq.replace('W', ' 2.1 ')
    seq = seq.replace('Y', ' 1.61 ')
    seq = seq.replace('V', ' 0.13 ')
    numseq_list = [float(x) for x in seq.split()]
    out_put = fourier(numseq_list, fourier_len)
    return out_put


def get_BulkinessZimmerman_fea(seq, fourier_len=10):
    seq = seq.replace('A', ' 11.5 ')
    seq = seq.replace('R', ' 14.28 ')
    seq = seq.replace('N', ' 12.82 ')
    seq = seq.replace('D', ' 11.68 ')
    seq = seq.replace('C', ' 13.46 ')
    seq = seq.replace('Q', ' 14.45 ')
    seq = seq.replace('E', ' 13.57 ')
    seq = seq.replace('G', ' 3.4 ')
    seq = seq.replace('H', ' 13.69 ')
    seq = seq.replace('I', ' 21.4 ')
    seq = seq.replace('L', ' 21.4 ')
    seq = seq.replace('K', ' 15.71 ')
    seq = seq.replace('M', ' 16.25 ')
    seq = seq.replace('F', ' 19.8 ')
    seq = seq.replace('P', ' 17.43 ')
    seq = seq.replace('S', ' 9.47 ')
    seq = seq.replace('T', ' 15.77 ')
    seq = seq.replace('W', ' 21.67 ')
    seq = seq.replace('Y', ' 18.03 ')
    seq = seq.replace('V', ' 21.57 ')
    numseq_list = [float(x) for x in seq.split()]
    out_put = fourier(numseq_list, fourier_len)
    return out_put


def get_IsoelectricPointZimmerman_fea(seq, fourier_len=10):
    seq = seq.replace('A', ' 6 ')
    seq = seq.replace('R', ' 10.76 ')
    seq = seq.replace('N', ' 5.41 ')
    seq = seq.replace('D', ' 2.77 ')
    seq = seq.replace('C', ' 5.05 ')
    seq = seq.replace('Q', ' 5.65 ')
    seq = seq.replace('E', ' 3.22 ')
    seq = seq.replace('G', ' 5.97 ')
    seq = seq.replace('H', ' 7.59 ')
    seq = seq.replace('I', ' 6.02 ')
    seq = seq.replace('L', ' 5.98 ')
    seq = seq.replace('K', ' 9.74 ')
    seq = seq.replace('M', ' 5.74 ')
    seq = seq.replace('F', ' 5.48 ')
    seq = seq.replace('P', ' 6.3 ')
    seq = seq.replace('S', ' 5.68 ')
    seq = seq.replace('T', ' 5.66 ')
    seq = seq.replace('W', ' 5.89 ')
    seq = seq.replace('Y', ' 5.66 ')
    seq = seq.replace('V', ' 5.96 ')
    numseq_list = [float(x) for x in seq.split()]
    out_put = fourier(numseq_list, fourier_len)
    return out_put


def get_HphobKyteDoolottle_fea(seq, fourier_len=10):
    seq = seq.replace('A', ' 1.8 ')
    seq = seq.replace('R', ' -4.5 ')
    seq = seq.replace('N', ' -3.5 ')
    seq = seq.replace('D', ' -3.5 ')
    seq = seq.replace('C', ' 2.5 ')
    seq = seq.replace('Q', ' -3.5 ')
    seq = seq.replace('E', ' -3.5 ')
    seq = seq.replace('G', ' -0.4 ')
    seq = seq.replace('H', ' -3.2 ')
    seq = seq.replace('I', ' 4.5 ')
    seq = seq.replace('L', ' 3.8 ')
    seq = seq.replace('K', ' -3.9 ')
    seq = seq.replace('M', ' 1.9 ')
    seq = seq.replace('F', ' 2.8 ')
    seq = seq.replace('P', ' -1.6 ')
    seq = seq.replace('S', ' -0.8 ')
    seq = seq.replace('T', ' -0.7 ')
    seq = seq.replace('W', ' -0.9 ')
    seq = seq.replace('Y', ' -1.3 ')
    seq = seq.replace('V', ' 4.2 ')
    numseq_list = [float(x) for x in seq.split()]
    out_put = fourier(numseq_list, fourier_len)
    return out_put


def get_HphobEisenberg_fea(seq, fourier_len=10):
    seq = seq.replace('A', ' 0.25 ')
    seq = seq.replace('R', ' -1.76 ')
    seq = seq.replace('N', ' -0.64 ')
    seq = seq.replace('D', ' -0.72 ')
    seq = seq.replace('C', ' 0.04 ')
    seq = seq.replace('Q', ' -0.69 ')
    seq = seq.replace('E', ' -0.62 ')
    seq = seq.replace('G', ' 0.16 ')
    seq = seq.replace('H', ' -0.4 ')
    seq = seq.replace('I', ' 0.73 ')
    seq = seq.replace('L', ' 0.53 ')
    seq = seq.replace('K', ' -1.1 ')
    seq = seq.replace('M', ' 0.26 ')
    seq = seq.replace('F', ' 0.61 ')
    seq = seq.replace('P', ' -0.07 ')
    seq = seq.replace('S', ' -0.26 ')
    seq = seq.replace('T', ' -0.18 ')
    seq = seq.replace('W', ' 0.37 ')
    seq = seq.replace('Y', ' 0.02 ')
    seq = seq.replace('V', ' 0.54 ')
    numseq_list = [float(x) for x in seq.split()]
    out_put = fourier(numseq_list, fourier_len)
    return out_put


def get_HphobHoppWoods_fea(seq, fourier_len=10):
    seq = seq.replace('A', ' -0.5 ')
    seq = seq.replace('R', ' 3 ')
    seq = seq.replace('N', ' 0.2 ')
    seq = seq.replace('D', ' 3 ')
    seq = seq.replace('C', ' -1 ')
    seq = seq.replace('Q', ' 0.2 ')
    seq = seq.replace('E', ' 3 ')
    seq = seq.replace('G', ' 0 ')
    seq = seq.replace('H', ' -0.5 ')
    seq = seq.replace('I', ' -1.8 ')
    seq = seq.replace('L', ' -1.8 ')
    seq = seq.replace('K', ' 3 ')
    seq = seq.replace('M', ' -1.3 ')
    seq = seq.replace('F', ' -2.5 ')
    seq = seq.replace('P', ' 0 ')
    seq = seq.replace('S', ' 0.3 ')
    seq = seq.replace('T', ' -0.4 ')
    seq = seq.replace('W', ' -3.4 ')
    seq = seq.replace('Y', ' -2.3 ')
    seq = seq.replace('V', ' -1.5 ')
    numseq_list = [float(x) for x in seq.split()]
    out_put = fourier(numseq_list, fourier_len)
    return out_put


def get_rna_pc_fea(ids, seqs, save_dir, fourier_len=10):
    fea = []
    for seq in seqs:
        hydrogenbonding_fea = get_hydrogenbonding_fea(seq, fourier_len)
        vanderwall_fea = get_vanderwall_fea(seq, fourier_len)
        fea.append(hydrogenbonding_fea + vanderwall_fea)

    fea_str = convert_list_to_str(fea)
    rna_pc_fea_file = os.path.join(save_dir, "lncRNA_pc_features")
    f = open(rna_pc_fea_file, "w")
    for k, v in zip(ids, fea_str):
        f.write(k + " " + v + "\n")
    f.close()

    # return fea


def get_pro_pc_fea(ids, seqs, save_dir, fourier_len=10):
    fea = []
    for seq in seqs:
        HphobBullbreese_fea = get_HphobBullbreese_fea(seq, fourier_len)
        PolarityGrantham_fea = get_PolarityGrantham_fea(seq, fourier_len)
        PolarityZimmerman_fea = get_PolarityZimmerman_fea(seq, fourier_len)
        BulkinessZimmerman_fea = get_BulkinessZimmerman_fea(seq, fourier_len)
        IsoelectricPointZimmerman_fea = get_IsoelectricPointZimmerman_fea(seq, fourier_len)
        HphobKyteDoolottle_fea = get_HphobKyteDoolottle_fea(seq, fourier_len)
        HphobEisenberg_fea = get_HphobEisenberg_fea(seq, fourier_len)
        HphobHoppWoods_fea = get_HphobHoppWoods_fea(seq, fourier_len)
        fea.append(HphobBullbreese_fea + PolarityGrantham_fea
                   + PolarityZimmerman_fea + BulkinessZimmerman_fea + IsoelectricPointZimmerman_fea
                   + HphobKyteDoolottle_fea + HphobEisenberg_fea + HphobHoppWoods_fea)

    fea_str = convert_list_to_str(fea)
    pro_pc_fea_file = os.path.join(save_dir, "protein_pc_features")
    f = open(pro_pc_fea_file, "w")
    for k, v in zip(ids, fea_str):
        f.write(k + " " + v + "\n")
    f.close()

    # return fea
