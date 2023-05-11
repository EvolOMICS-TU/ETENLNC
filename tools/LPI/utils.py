import numpy as np
import torch
import sys


def GetFasta(inputfile):
    '''Get sequence from input
    '''
    try:
        f = open(inputfile, 'r')    # input file
    except (IOError,ValueError) as e:
        print(sys.stderr, str(e))
        sys.exit(1)

    tmpseq = ''
    seqlist = []
    seqID = []
    #seq_dict = {}

    for line in f.readlines():
        line = line.strip()
        if not len(line):
            continue
        elif line[0] == '>':
            seqID.append(line.split()[0][1:])
            #seq_dict[line.split()[1][1:]] = ''
            if tmpseq != '':
                seqlist.append(tmpseq)
            tmpseq = ''
        else:
            tmpseq += line.upper()
    seqlist.append(tmpseq)      ## append the last sequence
    f.close()

    return [seqID, seqlist]


def convert_list_to_str(lists):
    list_str = []
    for l in lists:
        s = ""
        for n in l:
            s += (str(n)+' ')
        list_str.append(s)
    return list_str


def get_k_fold_data(k, test_idx, x):
    """用于产生交叉验证数据,返回第test_idx组数据做test时的训练和验证数据（输入为list,输出也为list）"""
    assert k > 1
    each_fold_size = len(x)//k
    x_train, x_test = None, None

    for j in range(k):
        idx = slice(j * each_fold_size, (j + 1) * each_fold_size)
        x_part = x[idx]
        if j == test_idx:
            x_test = x_part
        elif x_train is None:
            x_train = x_part
        else:
            x_train += x_part  # 要是x_train为一个样本 这就有bug

    return x_train, x_test


def normalize(train_data, test_data):
    """对数据进行归一化操作（输入为list，输出也为list）"""
    train_data_np = np.array(train_data, dtype=float)
    test_data_np = np.array(test_data, dtype=float)

    mean = np.mean(train_data_np, axis=0, keepdims=True)
    std = np.std(train_data_np, axis=0, ddof=1, keepdims=True)
    index = np.where(std == 0)  # 防止除数为零
    std[index] = 1e-7
    train_data_np = (train_data_np - mean) / std
    test_data_np = (test_data_np - mean) / std

    train_data = train_data_np.tolist()
    test_data = test_data_np.tolist()

    return train_data, test_data


def normalize_save(train_data, mean_save_file, std_save_file):
    """对数据进行标准化化操作（输入为list，输出也为list）"""
    train_data_np = np.array(train_data, dtype=float)

    mean = np.mean(train_data_np, axis=0, keepdims=True)
    std = np.std(train_data_np, axis=0, ddof=1, keepdims=True)
    index = np.where(std == 0)  # 防止除数为零
    std[index] = 1e-7
    train_data_np = (train_data_np - mean) / std
    train_data = train_data_np.tolist()

    print(mean.size)
    print(std.size)
    np.save(mean_save_file, mean)
    np.save(std_save_file, std)
    return train_data


def find_miss_data(net, train_id, train_data, train_label, test_id, test_data, test_label, p_pairs_file, n_pairs_file, output_file):
    """首先先整合出所有存在的rna或pro，整合为两个dict，并统计在正负样本中出现的次数。
    dict中含list
    dict()[0]为正样本中次数，
    dict()[1]为负样本中次数
    由于正样本包含所有用到的rna与pro,因此用正样本为模板创造dict"""
    rid = {}
    pid = {}
    f = open(p_pairs_file, 'r')
    for line in f.readlines():
        line = line.strip()
        if line.split()[0] in rid.keys():
            rid[line.split()[0]][0] += 1
        else:
            rid[line.split()[0]] = [1, 0, 0, [], 0, [], 0, [], 0, []]
            # rna[正样本中出现次数，负样本中次数，训练集miss的p次数，训练集miss的p目标...(后面等同)训n,测p,测n]
        if line.split()[1] in pid.keys():
            pid[line.split()[1]][0] += 1
        else:
            pid[line.split()[1]] = [1, 0, 0, [], 0, [], 0, [], 0, []]
    f.close()

    """再将负样本个数信息添加到dict()[1]"""
    f = open(n_pairs_file, 'r')
    for line in f.readlines():
        line = line.strip()
        rid[line.split()[0]][1] += 1
        pid[line.split()[1]][1] += 1
    f.close()

    """再将train_p train_n test_p test_n添加到dict"""

    net.eval()
    with torch.no_grad():
        train_probs = net(train_data)
        one = torch.ones_like(train_probs)
        zero = torch.zeros_like(train_probs)
        train_pred = torch.where(train_probs > 0.5, one, zero)
        train_compare = train_pred == train_label
        miss_idx = np.argwhere(train_compare.cpu().data.numpy() == 0).reshape(-1)
        train_miss_num = miss_idx.shape

        all_train_miss_target = train_id[miss_idx]
        all_train_miss_target_label = train_label.cpu().numpy()[miss_idx]

        p_train_miss_target_idx = np.argwhere(all_train_miss_target_label == 1).reshape(-1)
        p_train_miss_target = all_train_miss_target[p_train_miss_target_idx].tolist()
        for tp in p_train_miss_target:
            tp = tp.strip()
            rid[tp.split()[0]][2] += 1
            rid[tp.split()[0]][3].append(tp.split()[1])

            pid[tp.split()[1]][2] += 1
            pid[tp.split()[1]][3].append(tp.split()[0])

        n_train_miss_target_idx = np.argwhere(all_train_miss_target_label == 0).reshape(-1)
        n_train_miss_target = all_train_miss_target[n_train_miss_target_idx].tolist()
        for tn in n_train_miss_target:
            tn = tn.strip()
            rid[tn.split()[0]][4] += 1
            rid[tn.split()[0]][5].append(tn.split()[1])

            pid[tn.split()[1]][4] += 1
            pid[tn.split()[1]][5].append(tn.split()[0])


        test_probs = net(test_data)
        one = torch.ones_like(test_probs)
        zero = torch.zeros_like(test_probs)
        test_pred = torch.where(test_probs > 0.5, one, zero)
        test_compare = test_pred == test_label
        miss_idx = np.argwhere(test_compare.cpu().data.numpy() == 0).reshape(-1)
        test_miss_num = miss_idx.shape

        all_test_miss_target = test_id[miss_idx]
        all_test_miss_target_label = test_label.cpu().numpy()[miss_idx]

        p_test_miss_target_idx = np.argwhere(all_test_miss_target_label == 1).reshape(-1)
        p_test_miss_target = all_test_miss_target[p_test_miss_target_idx].tolist()
        for tp in p_test_miss_target:
            tp = tp.strip()
            rid[tp.split()[0]][6] += 1
            rid[tp.split()[0]][7].append(tp.split()[1])

            pid[tp.split()[1]][6] += 1
            pid[tp.split()[1]][7].append(tp.split()[0])

        n_test_miss_target_idx = np.argwhere(all_test_miss_target_label == 0).reshape(-1)
        n_test_miss_target = all_test_miss_target[n_test_miss_target_idx].tolist()
        for tn in n_test_miss_target:
            tn = tn.strip()
            rid[tn.split()[0]][8] += 1
            rid[tn.split()[0]][9].append(tn.split()[1])

            pid[tn.split()[1]][8] += 1
            pid[tn.split()[1]][9].append(tn.split()[0])

    net.train()

    rid = sorted(rid.items(), key=lambda x: x[1][6]+x[1][8], reverse=True)
    pid = sorted(pid.items(), key=lambda x: x[1][6]+x[1][8], reverse=True)

    f = open(output_file, 'w')
    for p in pid:
        f.write(p[0]+'\tp_num:'+str(p[1][0])+'\tn_num:'+str(p[1][1])+'\ttrain_miss:'+str(p[1][2]+p[1][4])+\
                '\ttest_miss:'+str(p[1][6]+p[1][8])+'\ntrain_p_miss:'+str(p[1][2])+'\t'+str(p[1][3]) +\
                '\ntrain_n_miss:'+str(p[1][4])+'\t'+str(p[1][5])+'\ntest_p_miss:'+str(p[1][6])+'\t'+str(p[1][7]) +\
                '\ntest_n_miss:'+str(p[1][8])+'\t'+str(p[1][9])+'\n')
    f.write('\n\n\n\n\n')
    for r in rid:
        f.write(r[0]+'\tp_num:'+str(r[1][0])+'\tn_num:'+str(r[1][1])+'\ttrain_miss:'+str(r[1][2]+r[1][4])+\
                '\ttest_miss:'+str(r[1][6]+r[1][8])+'\ntrain_p_miss:'+str(r[1][2])+'\t'+str(r[1][3]) +\
                '\ntrain_n_miss:'+str(r[1][4])+'\t'+str(r[1][5])+'\ntest_p_miss:'+str(r[1][6])+'\t'+str(r[1][7]) +\
                '\ntest_n_miss:'+str(r[1][8])+'\t'+str(r[1][9])+'\n')
    f.close()
    # print(train_miss_num)
    # print(test_miss_num)
    # print(rid)
    # print(pid)
    return train_miss_num, test_miss_num, rid, pid
