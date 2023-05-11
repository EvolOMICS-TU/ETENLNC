import os
import subprocess
from tqdm import tqdm
import shutil


def get_rna_ss_fea(rnaID, rna_seq, out_prefix):
    '''to compute structure features of rnas and proteins'''

    rna_out = os.path.join(out_prefix, "lncRNA_Structure_features")
    if os.path.exists(rna_out):
        os.remove(rna_out)

    features_workdir = os.path.join(out_prefix, "features_workdir")
    if os.path.exists(features_workdir):
        shutil.rmtree(features_workdir)
    os.mkdir(features_workdir)
    #####################
    # lncRNA structure
    rna_file_part = os.path.join(features_workdir, "tmp.rna.file.")
    rna_file_list = os.path.join(features_workdir, "tmp.filelist")

    rna_Seq = []
    for seq in rna_seq:
        seq = seq.replace('U', 'T')
        rna_Seq.append(seq)
    i = 0
    for rnaid, rnaseq in zip(rnaID, rna_Seq):
        f_tmp = open(rna_file_part + str(i), "w")
        i += 1
        f_tmp.write(">" + rnaid + "\n")
        f_tmp.write(rnaseq + "\n")
        f_tmp.close()

    file_list_cmd = "ls " + features_workdir + " |grep tmp.rna.file > " + rna_file_list
    #print file_list_cmd
    subprocess.call(file_list_cmd, shell=True)

    RNAScore2 = "./data/tools/RNAScore2"

    with open(rna_file_list, "r") as fr:
        bar = tqdm(fr.readlines())
        for tmp in bar:
            tmp = tmp.strip()
            tmpfile = os.path.join(features_workdir, tmp)
            tmpout = os.path.join(features_workdir, tmp + ".r_score")
            rna_cmd = RNAScore2 + " -i " + tmpfile + " -o " + tmpout + " -l 250 -r -noPS"
            #print rna_cmd
            subprocess.call(rna_cmd, shell=True)

            combine_cmd = "cat " + tmpout + " >> " + rna_out
            #print combine_cmd
            subprocess.call(combine_cmd, shell=True)
            os.remove(tmpfile)
            os.remove(tmpout)
            bar.set_description("Extract lncRNA Struct Features:")
        fr.close()
    #####################

    os.remove(rna_file_list)
    os.rmdir(features_workdir)
    print("Extract lncRNA struct features has finished.")
    return rna_out


def get_pro_ss_fea(proID, proSeq, out_prefix):
    '''to compute structure features of rnas and proteins'''

    protein_out = os.path.join(out_prefix, "protein_Structure_features")
    if os.path.exists(protein_out):
        os.remove(protein_out)

    #####################
    # protein structure
    features_workdir = os.path.join(out_prefix, "features_workdir")
    if os.path.exists(features_workdir):
        shutil.rmtree(features_workdir)
    os.mkdir(features_workdir)
    protein_file_part = os.path.join(features_workdir, "tmp.protein.file.")
    protein_file_list = os.path.join(features_workdir, "tmp.filelist")

    stride_dat = "./data/tools/stride.dat"
    stride_cmd = "cp " + stride_dat + " " + os.path.abspath('.')
    tmp_stride_dat = os.path.join(os.path.abspath('.'), "stride.dat")
    subprocess.call(stride_cmd, shell=True)

    i = 0
    for proid, proseq in zip(proID, proSeq):
        f_tmp = open(protein_file_part + str(i), "w")
        i += 1
        f_tmp.write(">" + proid + "\n")
        f_tmp.write(proseq + "\n")
        f_tmp.close()

    file_list_cmd = "ls " + features_workdir + \
        " |grep tmp.protein.file > " + protein_file_list
    #print file_list_cmd
    subprocess.call(file_list_cmd, shell=True)

    RNAScore2 = "./data/tools/RNAScore2"

    with open(protein_file_list, "r") as fp:
        bar = tqdm(fp.readlines())
        for tmp in bar:
            tmp = tmp.strip()
            tmpfile = os.path.join(features_workdir, tmp)
            tmpout = os.path.join(features_workdir, tmp + ".pro_score")
            protein_cmd = RNAScore2 + " -i " + tmpfile + " -o " + tmpout + " -p"
            #print protein_cmd
            subprocess.call(protein_cmd, shell=True)

            combine_cmd = "cat " + tmpout + " >> " + protein_out
            #print combine_cmd
            subprocess.call(combine_cmd, shell=True)
            os.remove(tmpfile)
            os.remove(tmpout)
            bar.set_description("Extract protein Struct Features:")

        fp.close()
    #####################

    os.remove(protein_file_list)
    os.remove(tmp_stride_dat)
    os.rmdir(features_workdir)
    print("Extract protein struct features has finished.")

    return protein_out