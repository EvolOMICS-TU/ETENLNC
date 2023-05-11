import os

import torch
import torch.nn as nn
import torch.utils.data as dataset

import pandas as pd
import numpy as np
import random
import copy

from features_script.wash_data import wash_rna_data
from features_script.wash_data import wash_pro_data
from get_features import get_features_files
from get_features import read_features_files
from get_features import get_pairs_file
from capsnet import NET
from utils import normalize_save

torch.manual_seed(1)
torch.cuda.manual_seed(1)

# Hyper parameters
EPOCH = 30
BATCH_SIZE = 100


def retrain(p_pairs_file, n_pairs_file, lr):

    p_data = read_features_files(p_pairs_file)
    n_data = read_features_files(n_pairs_file)

    random.seed(1)
    random.shuffle(p_data)  # 防止相似数据挨在一起
    random.shuffle(n_data)

    train_data = p_data + n_data

    train_label = torch.cat((torch.ones(len(p_data)), torch.zeros(len(n_data))))

    mean_save_file = "/LPI/data/user_data/model/normalize_mean"
    std_save_file = "/LPI/data/user_data/model/normalize_std"
    train_data = normalize_save(train_data, mean_save_file, std_save_file)  # normalize data, and save mean and std in a file

    train_data = torch.tensor(train_data)

    model_save_file = "/LPI/data/user_data/model/parameter.pkl"

    torch_dataset = dataset.TensorDataset(train_data, train_label)
    loader = dataset.DataLoader(
        dataset=torch_dataset,
        batch_size=BATCH_SIZE,
        shuffle=True,
    )

    # get net,optimizer and loss_fun
    net = NET()

    net = net
    optimizer = torch.optim.Adam(net.parameters(), lr=lr)
    loss_fn = nn.BCELoss()

    # train
    for epoch in range(EPOCH):
        for batch_idx, (train_data, train_label) in enumerate(loader):
            train_probs = net(train_data)
            train_loss = loss_fn(train_probs, train_label)
            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()

    torch.save(net.state_dict(), model_save_file)

    print("Retrain Finshed.")


def predict(predict_pairs_file, output_file, user_model=False):
    if user_model:
        parameter_path = '/LPI/data/user_data/model/parameter.pkl'
        normalize_mean = np.load('/LPI/data/user_data/model/normalize_mean.npy')
        normalize_std = np.load('/LPI/data/user_data/model/normalize_std.npy')
        net = NET()
        net.load_state_dict(torch.load(parameter_path))

        predict_data = read_features_files(predict_pairs_file)

        predict_data_np = np.array(predict_data, dtype=float)
        predict_data_np = (predict_data_np - normalize_mean) / normalize_std  # normalize
        predict_data = predict_data_np.tolist()
        predict_data = torch.tensor(predict_data)

        net.eval()
        with torch.no_grad():
            predict_probs = net(predict_data)
            one = torch.ones_like(predict_probs)
            zero = torch.zeros_like(predict_probs)
            predict_pred = torch.where(predict_probs > 0.5, one, zero)

    else:
        predict_data_save = read_features_files(predict_pairs_file)
        flag = 0
        for i in range(33):
            parameter_path = '/LPI/data/default_data/default_model/parameter'+str(i)+'.pkl'
            normalize_mean = np.load('/LPI/data/default_data/default_model/normalize_mean'+str(i)+'.npy')
            normalize_std = np.load('/LPI/data/default_data/default_model/normalize_std'+str(i)+'.npy')
            net = NET()
            net.load_state_dict(torch.load(parameter_path))
            predict_data_np = np.array(predict_data_save, dtype=float)
            predict_data_np = (predict_data_np - normalize_mean) / normalize_std  # normalize
            predict_data = predict_data_np.tolist()
            predict_data = torch.tensor(predict_data)

            net.eval()
            with torch.no_grad():
                predict_probs = net(predict_data)
                one = torch.ones_like(predict_probs)
                zero = torch.zeros_like(predict_probs)
                predict_pred = torch.where(predict_probs > 0.5, one, zero)

            predict_pred = predict_pred.cpu().numpy().tolist()

            if flag == 0:
                all_pred = copy.copy(predict_pred)
                flag = 1
            else:
                for j in range(len(predict_pred)):
                    if predict_pred[j] == 1:
                        all_pred[j] += 1

        all_pred = torch.tensor(all_pred)
        one = torch.ones_like(all_pred)
        zero = torch.zeros_like(all_pred)
        predict_pred = torch.where(all_pred > 16, one, zero)

    predict_pred = predict_pred.cpu().numpy()
    rna_id = []
    pro_id = []
    f = open(predict_pairs_file, 'r')
    for line in f.readlines():
        line = line.strip()
        if len(line) == 0:
            continue
        rna_id.append(line.split()[0])
        pro_id.append(line.split()[1])
    f.close()

    result = pd.DataFrame({'lncRNA-ID': rna_id, 'protein-ID': pro_id,
                           "predict('1' is interact,'0' is not interact)": predict_pred})
    result.to_csv(output_file + '.csv')
    print("\nPrediction is finished, check the result at:" + output_file + '.csv')


if __name__ == '__main__':
    # Prepare fasta files of lncRNAs and proteins that contain all the lncRNAs and proteins used in prediction and
    # retraining.(See README)
    rna_fasta_file_name = "Novel_LncRNAs.fa"  # Enter the name of your file
    pro_fasta_file_name = "proteins.fa"  # Enter the name of your file


    # If you want to retrain the model, you need to prepare retrain pairs files. (See README)
    RETRAIN = False  # choose: True, False
    retrain_interacting_pairs_file_name = "example_interacting_pairs"  # Enter the name of your file
    retrain_uninteracting_pairs_file_name = "example_negative_pairs"  # Enter the name of your file


    # If you want to make prediction, set "PREDICT = True". And the lncRNAs and proteins in "rna_fasta_file_name" and
    # "pro_fasta_file_name", respectively, will form lncRNA-protein pairs one by one.
    PREDICT = True  # choose: True, False
    predict_output_file_name = "lnc_protein"  # Enter the name of your file
    USE_RETRAIN_MODEL = False   # choose: True, False  (If you have retrained the model, you can choose True.)
    # You can also assign specific pairs to predict. In this way, you need to prepare predict pairs file. (See README)
    predict_designated_pairs = False
    predict_pairs = "example_predict_pairs"  # Enter the name of your file


    # When a complete feature extraction is performed once, the feature file is saved. This option can be set to False
    # to accelerate prediction or retrain.
    EXTRACT_FEATURES = True  # choose: True, False

    if EXTRACT_FEATURES:
        fasta_dir = "/home/cluster/nath/ETENLNC_cluster/tools/fasta_data/"
        wash_dir = "/LPI/data/user_data/wash_dir"
        rna_fasta = os.path.join(fasta_dir, rna_fasta_file_name)
        pro_fasta = os.path.join(fasta_dir, pro_fasta_file_name)

        rna_fasta_washed = os.path.join(wash_dir, "lncRNA_washed.fa")  # Correcting illegal sequences
        pro_fasta_washed = os.path.join(wash_dir, "protein_washed.fa")
        wash_rna_data(rna_fasta, rna_fasta_washed)
        wash_pro_data(pro_fasta, pro_fasta_washed)

        get_features_files(rna_fasta_washed, pro_fasta_washed)

    if RETRAIN:
        retrain_pairs_dir = '/LPI/data/user_data/retrain_pairs_data'
        p_pairs_file = os.path.join(retrain_pairs_dir, retrain_interacting_pairs_file_name)
        n_pairs_file = os.path.join(retrain_pairs_dir, retrain_uninteracting_pairs_file_name)

        retrain(p_pairs_file, n_pairs_file, lr=0.001)

    if PREDICT:
        predict_pairs_dir = '/LPI/data/user_data/predict_pairs_data'
        if predict_designated_pairs:
            predict_pairs_file = os.path.join(predict_pairs_dir, predict_pairs)
        else:
            fasta_dir = "/home/cluster/nath/ETENLNC_cluster/tools/fasta_data/"
            rna_fasta = os.path.join(fasta_dir, rna_fasta_file_name)
            pro_fasta = os.path.join(fasta_dir, pro_fasta_file_name)
            output_pairs_file = get_pairs_file(rna_fasta, pro_fasta)
            predict_pairs_file = os.path.join(output_pairs_file)

        output_dir = '/home/cluster/nath/ETENLNC_cluster/tools/fasta_data/'
        output_file = os.path.join(output_dir, predict_output_file_name)
        predict(predict_pairs_file, output_file, user_model=USE_RETRAIN_MODEL)


