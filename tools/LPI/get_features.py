import os

from utils import GetFasta
from features_script.ss import get_rna_ss_fea
from features_script.ss import get_pro_ss_fea
from features_script.kmer import get_rna_kmer_fea
from features_script.kmer import get_pro_kmer_fea
from features_script.motif import get_rna_motif_fea
from features_script.motif import get_pro_motif_fea
from features_script.pc import get_rna_pc_fea
from features_script.pc import get_pro_pc_fea


def get_pairs_file(rna_seq_file, pro_seq_file):
    output_pairs_file = "/LPI/data/user_data/predict_pairs_data/created_predict_pairs"
    f = open(output_pairs_file, 'w')
    rna_ids, rna_seqs = GetFasta(rna_seq_file)
    pro_ids, pro_seqs = GetFasta(pro_seq_file)
    for r in rna_ids:
        for p in pro_ids:
            f.write(r + ' ' + p + '\n')
    f.close()
    return output_pairs_file


def get_features_files(rna_seq_file, pro_seq_file):
    features_dir = "/LPI/data/user_data/features_data"
    rna_ids, rna_seqs = GetFasta(rna_seq_file)
    pro_ids, pro_seqs = GetFasta(pro_seq_file)
    get_rna_ss_fea(rna_ids, rna_seqs, features_dir)
    get_pro_ss_fea(pro_ids, pro_seqs, features_dir)
    get_rna_kmer_fea(rna_ids, rna_seqs, 4, features_dir, freq=True)
    get_pro_kmer_fea(pro_ids, pro_seqs, 3, features_dir, freq=True)
    get_rna_motif_fea(rna_ids, rna_seqs, features_dir)
    get_pro_motif_fea(pro_ids, pro_seqs, features_dir)
    get_rna_pc_fea(rna_ids, rna_seqs, features_dir, fourier_len=10)
    get_pro_pc_fea(pro_ids, pro_seqs, features_dir, fourier_len=10)
    print("Extract features has all finished.")


def read_fea_file(fea_file):
    """read features"""
    fea_dict = {}

    f = open(fea_file, 'r')
    for line in f.readlines():
        line = line.strip()
        if len(line.split()) < 5:
            continue
        fea_dict[line.split()[0]] = [float(x) for x in line.split()[1:]]
    f.close()

    return fea_dict


def read_features_files(pairs_file):
    features_dir = "/LPI/data/user_data/features_data"
    rna_stu_fea_file = os.path.join(features_dir, "lncRNA_Structure_features")
    pro_stu_fea_file = os.path.join(features_dir, "protein_Structure_features")
    lncRNA_kmer_fea_file = os.path.join(features_dir, "lncRNA_kmer_features")
    protein_kmer_fea_file = os.path.join(features_dir, "protein_kmer_features")
    lncRNA_motif_fea_file = os.path.join(features_dir, "lncRNA_motif_features")
    protein_motif_fea_file = os.path.join(features_dir, "protein_motif_features")
    lncRNA_pc_fea_file = os.path.join(features_dir, "lncRNA_pc_features")
    protein_pc_fea_file = os.path.join(features_dir, "protein_pc_features")

    rna_stru_fea = read_fea_file(rna_stu_fea_file)
    pro_stru_fea = read_fea_file(pro_stu_fea_file)
    lncRNA_kmer_fea = read_fea_file(lncRNA_kmer_fea_file)
    protein_kmer_fea = read_fea_file(protein_kmer_fea_file)
    lncRNA_motif_fea = read_fea_file(lncRNA_motif_fea_file)
    protein_motif_fea = read_fea_file(protein_motif_fea_file)
    lncRNA_pc_fea = read_fea_file(lncRNA_pc_fea_file)
    protein_pc_fea = read_fea_file(protein_pc_fea_file)

    # print(len(rna_stru_fea['n1114']))  # 30
    # print(len(pro_stru_fea['Q15717']))  # 50
    # print(len(lncRNA_kmer_fea['n1114']))  # 256
    # print(len(protein_kmer_fea['Q15717']))  # 64
    # print(len(lncRNA_motif_fea['n1114']))  # 18
    # print(len(protein_motif_fea['Q15717']))  # 11
    # print(len(lncRNA_pc_fea['n1114']))  # 20
    # print(len(protein_pc_fea['Q15717']))  # 80

    rna_id = []
    pro_id = []

    f = open(pairs_file, 'r')
    for line in f.readlines():
        line = line.strip()
        rna_id.append(line.split()[0])
        pro_id.append(line.split()[1])
    f.close()

    features = []

    for r, p in zip(rna_id, pro_id):
        features.append(rna_stru_fea[r] + lncRNA_kmer_fea[r]
                        + lncRNA_motif_fea[r] + lncRNA_pc_fea[r]
                        + pro_stru_fea[p] + protein_kmer_fea[p]
                        + protein_motif_fea[p] + protein_pc_fea[p])
    return features
