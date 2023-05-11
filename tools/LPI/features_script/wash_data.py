from utils import GetFasta


def wash_rna_data(rna_fasta_file, out_file):
    ids, Seqs = GetFasta(rna_fasta_file)
    seqs = []
    for seq in Seqs:
        seq = seq.replace('U', 'T')
        seqs.append(seq)
    element = ['A', 'C', 'G', 'T']
    for i in range(len(seqs)):
        for j in range(len(seqs[i])):
            if seqs[i][j] not in element:
                s = seqs[i][j]
                print(s)
                print(seqs[i])
                seqs[i] = seqs[i].replace(s, 'A')

    f = open(out_file, 'w')
    for i, s in zip(ids, seqs):
        f.write('>'+i+'\n'+s+'\n')
    f.close()


def wash_pro_data(pro_fasta_file, out_file):
    ids, seqs = GetFasta(pro_fasta_file)

    element = ['G', 'A', 'V', 'L', 'I', 'F', 'W', 'Y', 'D', 'N', 'E', 'K', 'Q', 'M', 'S', 'T', 'C', 'P', 'H', 'R']
    for i in range(len(seqs)):
        for j in range(len(seqs[i])):
            if seqs[i][j] not in element:
                s = seqs[i][j]
                print(s)
                print(seqs[i])
                seqs[i] = seqs[i].replace(s, 'A')
        while len(seqs[i]) < 15:
            seqs[i] += 'A'
        if len(seqs[i]) >= 10000:
            seqs[i] = seqs[i][:10000]
    f = open(out_file, 'w')
    for i, s in zip(ids, seqs):
        f.write('>'+i+'\n'+s+'\n')
    f.close()
