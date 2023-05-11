import torch
import torch.nn as nn
import torch.nn.functional as func
import math


def squash(x):
    length2 = x.pow(2).sum(dim=2)+1e-7
    length = length2.sqrt()
    x = x*(length2/(length2+1)/length).view(x.size(0), x.size(1), -1)
    return x


class CapsLayer(nn.Module):
    def __init__(self, input_caps, input_dim, output_caps, output_dim):
        super().__init__()
        self.input_caps = input_caps
        self.input_dim = input_dim
        self.output_caps = output_caps
        self.output_dim = output_dim
        self.weights = nn.Parameter(torch.Tensor(self.input_caps, self.input_dim, self.output_caps * self.output_dim))
        # self.routing_module = AgreementRouting(self.input_caps, self.output_caps)
        self.reset_parameters()

    def reset_parameters(self):
        stdv = 1. / math.sqrt(self.input_caps)
        self.weights.data.uniform_(-stdv, stdv)

    # 输入的x格式为：(batch_m, input_caps, input_dim)
    def forward(self, u):
        u = u.unsqueeze(2)
        u_predict = u.matmul(self.weights)
        u_predict = u_predict.view(u_predict.size(0), self.input_caps, self.output_caps, self.output_dim)
        s = u_predict.sum(dim=1)
        v = squash(s)
        probs = v.pow(2).sum(dim=2).sqrt()
        return v, probs


class NET(nn.Module):
    def __init__(self):
        super().__init__()
        self.seq_togather1 = nn.Linear(320, 200)
        self.seq_togather2 = nn.Linear(200, 100)
        self.seq_togather3 = nn.Linear(100, 50)
        self.seq_togather4 = nn.Linear(50, 3)

        self.mot_togather1 = nn.Linear(29, 30)
        self.mot_togather2 = nn.Linear(30, 30)
        self.mot_togather3 = nn.Linear(30, 30)
        self.mot_togather4 = nn.Linear(30, 3)

        self.pc_togather1 = nn.Linear(100, 100)
        self.pc_togather2 = nn.Linear(100, 50)
        self.pc_togather3 = nn.Linear(50, 20)
        self.pc_togather4 = nn.Linear(20, 3)

        self.ss_togather1 = nn.Linear(20, 30)
        self.ss_togather2 = nn.Linear(30, 30)
        self.ss_togather3 = nn.Linear(30, 30)
        self.ss_togather4 = nn.Linear(30, 3)

        self.caps = CapsLayer(input_caps=4, input_dim=3, output_caps=1, output_dim=3)

        self.prelu = nn.PReLU()
        self.dropout = nn.Dropout(p=0.5)

    def forward(self, x):
        rna_ss = x[:, :10]
        pro_ss = x[:, 324:334]
        ss = torch.cat((rna_ss, pro_ss), dim=1)
        rna_seq = x[:, 30:286]
        pro_seq = x[:, 374:438]
        seq = torch.cat((rna_seq, pro_seq), dim=1)
        rna_mot = x[:, 286:304]
        pro_mot = x[:, 438:449]
        mot = torch.cat((rna_mot, pro_mot), dim=1)
        rna_pc = x[:, 304:324]
        pro_pc = x[:, 449:529]
        pc = torch.cat((rna_pc, pro_pc), dim=1)

        seq_togater = self.prelu(self.dropout(self.seq_togather1(seq)))
        seq_togater = self.prelu(self.dropout(self.seq_togather2(seq_togater)))
        seq_togater = self.prelu(self.dropout(self.seq_togather3(seq_togater)))
        seq_togater = torch.tanh(self.seq_togather4(seq_togater)).view(seq_togater.shape[0], 1, -1)
        # probs = torch.sigmoid(self.seq_togather4(seq_togater))

        mot_togater = self.prelu(self.dropout(self.mot_togather1(mot)))
        mot_togater = self.prelu(self.dropout(self.mot_togather2(mot_togater)))
        mot_togater = self.prelu(self.dropout(self.mot_togather3(mot_togater)))
        mot_togater = torch.tanh(self.mot_togather4(mot_togater)).view(mot_togater.shape[0], 1, -1)
        # probs = torch.sigmoid(self.mot_togather4(mot_togater))

        pc_togater = self.prelu(self.dropout(self.pc_togather1(pc)))
        pc_togater = self.prelu(self.dropout(self.pc_togather2(pc_togater)))
        pc_togater = self.prelu(self.dropout(self.pc_togather3(pc_togater)))
        pc_togater = torch.tanh(self.pc_togather4(pc_togater)).view(pc_togater.shape[0], 1, -1)
        # probs = torch.sigmoid(self.pc_togather4(pc_togater))

        ss_togater = self.prelu(self.dropout(self.ss_togather1(ss)))
        ss_togater = self.prelu(self.dropout(self.ss_togather2(ss_togater)))
        ss_togater = self.prelu(self.dropout(self.ss_togather3(ss_togater)))
        ss_togater = torch.tanh(self.ss_togather4(ss_togater)).view(ss_togater.shape[0], 1, -1)
        # probs = torch.sigmoid(self.ss_togather4(ss_togater))

        togater = torch.cat((seq_togater, pc_togater, mot_togater, ss_togater), dim=1)

        v, probs = self.caps(togater)

        return probs.squeeze(1)


# class AgreementRouting(nn.Module):
#     def __init__(self, input_caps, output_caps, n_iterations=3):
#         super().__init__()
#         self.n_iterations = n_iterations
#         self.b = torch.zeros((input_caps, output_caps)).cuda()
#
#     def forward(self, u_predict):
#         batch_size, input_caps, output_caps, output_dim = u_predict.size()
#         self.b.zero_()
#         c = func.softmax(self.b, dim=1)
#         s = (c.unsqueeze(2) * u_predict).sum(dim=1)
#         v = squash(s)
#
#         if self.n_iterations > 0:
#             b_batch = self.b.expand((batch_size, input_caps, output_caps))  # 这块带上batch是因为每个样本的c都不一样
#             for r in range(self.n_iterations):
#                 v = v.unsqueeze(1)
#                 b_batch = b_batch + (u_predict * v).sum(-1)
#
#                 c = func.softmax(b_batch.view(-1, output_caps), dim=1).view(-1, input_caps, output_caps, 1)
#                 s = (c * u_predict).sum(dim=1)
#                 v = squash(s)
#
#         return v


# # loss_fun
# class MarginLoss(nn.Module):
#     def __init__(self, m_pos, m_neg, lambda_):
#         super(MarginLoss, self).__init__()
#         self.m_pos = m_pos
#         self.m_neg = m_neg
#         self.lambda_ = lambda_
#
#     def forward(self, lengths, targets, size_average=True):
#         t = torch.zeros(lengths.size()).long()
#         if targets.is_cuda:
#             t = t.cuda()
#         t = t.scatter_(1, targets.data.view(-1, 1), 1)  # 按行传播 相当于one-hot
#         losses = t.float() * func.relu(self.m_pos - lengths).pow(2) + \
#             self.lambda_ * (1. - t.float()) * func.relu(lengths - self.m_neg).pow(2)
#         return losses.mean() if size_average else losses.sum()

