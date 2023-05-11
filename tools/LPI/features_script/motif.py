import os
from utils import convert_list_to_str


def get_rna_motif_fea(ids, seqs, save_dir):
    """NOTE:Since lncRNA mostly uses T instead of U, this function uses T.
    If used this function alone, the function or sequence may need to be modified."""
    motif_fea = []
    Fox1 = ("TGCATGT", )
    Nova = ("TCATTTCAC", "TCATTTCAT", "CCATTTCAC", "CCATTTCAT")
    Slm2 = ("TAAAC", "TAAAA", "TAATC", "TAATA")
    Fusip1 = ("AAAGA", "AAAGG", "AGAGA", "AGAGG", "CAAGA", "CAAGG", "CGAGA", "CGAGG")
    PTB = ("TTTTT", "TTTCT", "TCTTT", "TCTCT")
    ARE = ("TATTTATT", )
    hnRNPA1 = ("TAGGGT", "TAGGGA")
    PUM = ("TGTAAATA", "TGTAGATA", "TGTATATA", "TGTACATA")
    U1A = ("ATTGCAC", )
    HuD = ("TTATTT", )
    QKI = ("ATTAAT", "ATTAAC", "ACTAAT", "ACTAAC")
    U2B = ("ATTGCAG", )
    SF1 = ("TACTAAC", )
    HuR = ("TTTATTT", "TTTGTTT", "TTTCTTT", "TTTTTTT")
    YB1 = ("CCTGCG", "TCTGCG")
    AU = ("AT", )
    UG = ("TG", )
    FIVE_OP = ""
    motifs_tuple = (Fox1, Nova, Slm2, Fusip1, PTB, ARE, hnRNPA1, PUM, U1A, HuD, QKI, U2B, SF1, HuR, YB1, AU, UG, FIVE_OP)
    for seq in seqs:
        fea_tmp = []
        for motif in motifs_tuple:
            cont = 0
            if motif == "":
                cont = fea_tmp[0] + fea_tmp[1] + fea_tmp[5] + fea_tmp[7] + fea_tmp[8]
                fea_tmp.append(cont)
                break
            for m in motif:
                cont += seq.count(m)
            fea_tmp.append(cont)

        motif_fea.append(fea_tmp)

    motif_fea_str = convert_list_to_str(motif_fea)
    rna_motif_fea_file = os.path.join(save_dir, "lncRNA_motif_features")
    f = open(rna_motif_fea_file, "w")
    for k, v in zip(ids, motif_fea_str):
        f.write(k + " " + v + "\n")
    f.close()
    # return motif_fea


def get_pro_motif_fea(ids, seqs, save_dir):
    motif_fea = []
    E = ("E",)
    K = ("K",)
    EE = ("EE",)
    KK = ("KK",)
    H_R = ("H", "R")
    RS_SR = ("RS", "SR")
    RGG = ("RGG",)
    YGG = ("YGG",)
    R = ("R",)
    H = ("H",)
    HR_RH = ("HR", "RH")
    motifs_tuple = (E, K, EE, KK, H_R, RS_SR, RGG, YGG, R, H, HR_RH)
    for seq in seqs:
        fea_tmp = []
        for motif in motifs_tuple:
            cont = 0
            for m in motif:
                cont += seq.count(m)
            fea_tmp.append(cont)

        motif_fea.append(fea_tmp)

    motif_fea_str = convert_list_to_str(motif_fea)
    pro_motif_fea_file = os.path.join(save_dir, "protein_motif_features")
    f = open(pro_motif_fea_file, "w")
    for k, v in zip(ids, motif_fea_str):
        f.write(k + " " + v + "\n")
    f.close()

    # return motif_fea


# def get_features(rna_fasta_file, pro_fasta_file, save_dir):
#     rna_ids, rna_Seqs = GetFasta(rna_fasta_file)
#     pro_ids, pro_seqs = GetFasta(pro_fasta_file)
#     assert len(rna_ids) == len(pro_ids)
#
#     rna_seqs = []
#     for seq in rna_Seqs:
#         seq = seq.replace('U', 'T')
#         rna_seqs.append(seq)
#
#
#     # 还有个检验rna,pro规不规范的
#     # rna_kmer_fea = get_rna_kmer_fea(rna_ids, rna_seqs, 4, save_dir, freq=True)
#     pro_kmer_fea = get_pro_kmer_fea(pro_ids, pro_seqs, 3, save_dir, freq=True)
#
#     # rna_seq_fea = RNA_SequenceFeatures(rna_ids, rna_seqs, save_dir)
#     # pro_seq_fea = Protein_SequenceFeatures(pro_ids, pro_seqs, save_dir)
#     # # rna_seq_fea_np = np.array(rna_seq_fea)
#     # # pro_seq_fea_np = np.array(pro_seq_fea)
#
#     # rna_motif_fea = get_rna_motif_fea(rna_ids, rna_seqs, save_dir)
#     # pro_motif_fea = get_pro_motif_fea(pro_ids, pro_seqs, save_dir)
#     # # rna_motif_fea_np = np.array(rna_motif_fea)
#     # # pro_motif_fea_np = np.array(pro_motif_fea)
#
#     # rna_pc_fea = get_rna_pc_fea(rna_ids, rna_seqs, save_dir)
#     # pro_pc_fea = get_pro_pc_fea(pro_ids, pro_seqs, save_dir)
#     # rna_pc_fea_np = np.array(rna_pc_fea)
#     # pro_pc_fea_np = np.array(pro_pc_fea)
#
#     # rna_ss_fea = RNA_StructureFeatures(rna_ids, rna_seqs, save_dir)
#     # pro_ss_fea = Protein_StructureFeatures(pro_ids, pro_seqs, save_dir)
#     # rna_ss_fea_np = np.array(rna_ss_fea)
#     # pro_ss_fea_np = np.array(pro_ss_fea)
#
#     # fea = np.hstack(rna_seq_fea_np, pro_seq_fea_np)
#     # fea = np.hstack(rna_motif_fea_np, pro_motif_fea_np)
#     # fea = np.hstack(rna_pc_fea_np, pro_pc_fea_np)
#
#     # return motif_fea


# def CountOccurrences(string, substring):
#     # 初始化count和start为0
#     count = 0
#     start = 0
#
#     # 搜索字符串直到终点
#     while start < len(string):
#
#         # 检查从"start"位置到结束是否存在子字符串
#         flag = string.find(substring, start)
#
#         if flag != -1:
#             # 如果存在子字符串，将"start"从子字符串的开始移动到下一个位置
#             start = flag + 1
#
#             # count自增
#             count += 1
#         else:
#             # 如果没有其他子字符串，返回count的值
#             return count


if __name__ == "__main__":
    pos_rna_fasta_file = "D:/program/get_features2/data/wash_dir/NPI_6204_pos_seq_pairs_rna_w"
    pos_pro_fasta_file = "D:/program/get_features2/data/wash_dir/NPI_6204_pos_seq_pairs_pro_w"
    fea_save_dir = "D:/program/get_features2/data/features_data"

    get_features(pos_rna_fasta_file, pos_pro_fasta_file, fea_save_dir)

