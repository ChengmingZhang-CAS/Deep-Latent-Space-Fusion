"""
By Chengming Zhang (zhangchengming2017@sibcb.ac.cn), July 7, 2020.
All rights reserved.

architecture:
    X_1 -> Z_1 -> X_1-recon -> Z_1-recon
    X_2 -> Z_2 -> X_2-recon -> Z_2-recon
    X_3 -> Z_3 -> X_3-recon -> Z_3-recon

"""
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn import preprocessing
from sklearn.metrics import silhouette_score
from collections import Counter
from sklearn.decomposition import PCA
import os


class AutoEncoder(nn.Module):
    def __init__(self, input_dims):
        super(AutoEncoder, self).__init__()

        assert isinstance(input_dims, list)
        self.encoder = nn.Sequential(
            nn.Linear(input_dims[0], input_dims[1]),
            nn.LeakyReLU(),
            nn.Linear(input_dims[1], input_dims[2]),
            nn.Sigmoid()
        )
        self.decoder = nn.Sequential(
            nn.Linear(input_dims[2], input_dims[1]),
            nn.LeakyReLU(),
            nn.Linear(input_dims[1], input_dims[0]),
            nn.Sigmoid()
        )

    def forward(self, x):
        z = self.encoder(x)
        x_recon = self.decoder(z)
        return x_recon


class build_AE(nn.Module):
    def __init__(self, input_dims):
        """
        :param input_dims: network architecture neurons
        """
        super(build_AE, self).__init__()
        assert isinstance(input_dims, list)
        self.encoder = nn.Sequential()
        for i in range(0, len(input_dims) - 1):
            #  Each layer will divide the size of feature map by 2
            self.encoder.add_module('Encoder%d' % i, nn.Linear(input_dims[i], input_dims[i + 1]))
            self.encoder.add_module('LeakyRelu%d' % i, nn.LeakyReLU())
        input_dims = list(reversed(input_dims))
        self.decoder = nn.Sequential()
        for i in range(0, len(input_dims) - 1):
            # W'Z + LeakyRelu
            self.decoder.add_module('Decoder%d' % i, nn.Linear(input_dims[i], input_dims[i + 1]))
            self.decoder.add_module('LeakyRelu%d' % i, nn.LeakyReLU())

    def forward(self, x):
        z = self.encoder(x)
        y = self.decoder(z)
        return y


def train_AE(model, X, epochs, lr=1e-3, device='cpu', show_freq=10):
    # ensure input in the form of Tensor
    if not isinstance(X, torch.Tensor):
        X = torch.tensor(np.array(X), dtype=torch.float32, device=device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    loss_record = []
    for epoch in range(epochs):
        X_recon = model.forward(X)
        loss = F.mse_loss(X_recon, X, reduction='mean')
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        loss_record.append(loss)
        if epoch % show_freq == 0 or epoch == epochs - 1:
            # loss X-to-X reconstruction
            print('Epoch {}: Loss_XX: {}'.format(epoch, loss))
            print(X.numpy()[100, 10:14], X_recon.detach().numpy()[100, 10:14])

    return loss_record


class CycleAE(nn.Module):
    def __init__(self, input_dims):
        super(CycleAE, self).__init__()

        assert isinstance(input_dims, list)
        self.encoder = nn.Sequential(
            nn.Linear(input_dims[0], input_dims[1]),
            nn.LeakyReLU(),
            nn.Linear(input_dims[1], input_dims[2]),
            nn.Sigmoid()
        )
        self.decoder = nn.Sequential(
            nn.Linear(input_dims[2], input_dims[1]),
            nn.LeakyReLU(),
            nn.Linear(input_dims[1], input_dims[0]),
            nn.Sigmoid()
        )

    def forward(self, x):
        z = self.encoder(x)
        x_recon = self.decoder(z)
        z_recon = self.encoder(x_recon)
        return z, x_recon, z_recon


def train_CycleAE(model, X, epochs, lr=1e-3, weight_xx=1.0, weight_zz=1.0, device='cpu', show_freq=10):
    # ensure input in the form of Tensor
    if not isinstance(X, torch.Tensor):
        X = torch.tensor(np.array(X), dtype=torch.float32, device=device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    loss_record = []
    for epoch in range(epochs):
        Z, X_recon, Z_recon = model.forward(X)
        loss1 = F.mse_loss(X_recon, X, reduction='mean')
        loss2 = F.mse_loss(Z_recon, Z, reduction='mean')
        loss = weight_xx * loss1 + weight_zz * loss2
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        loss_record.append([loss1, loss2])
        if epoch % show_freq == 0 or epoch == epochs - 1:
            # loss: X-X, Z-Z
            print('Epoch {}: Loss_XX: {}'.format(epoch, loss1))
            print('Epoch {}: Loss_ZZ: {}'.format(epoch, loss2))
            print('Epoch {}: total loss: {}'.format(epoch, loss))
            print(X.numpy()[100, 10:14], X_recon.detach().numpy()[100, 10:14])
            print(Z.detach().numpy()[100, 10:14], Z_recon.detach().numpy()[100, 10:14])

    return loss_record


class MultiAE(nn.Module):
    def __init__(self, num_omics, input_dims, num_sample):
        super(MultiAE, self).__init__()
        self.K = num_omics
        self.n = num_sample
        self.AE = [build_AE(input_dims[k]) for k in range(num_omics)]

    def forward(self, X):  # shape=[n, d]
        Z = []
        X_recon = []
        Z_recon = []
        for k in range(self.K):
            x_input = X[k]
            z = self.AE[k].encoder(x_input)
            x_recon = self.AE[k].decoder(z)
            z_recon = self.AE[k].encoder(x_recon)

            Z.append(z)
            X_recon.append(x_recon)
            Z_recon.append(z_recon)

        return Z, X_recon, Z_recon

    def loss_fn(self, X, Z, X_recon, Z_recon, weight_xx, weight_zz):
        Loss_XX = []
        Loss_ZZ = []
        for k in range(self.K):
            loss_xx = F.mse_loss(X_recon[k], X[k], reduction='sum')
            loss_zz = F.mse_loss(Z_recon[k], Z[k], reduction='sum')
            Loss_XX.append(loss_xx)
            Loss_ZZ.append(loss_zz)
        loss = weight_xx * np.sum(Loss_XX) + weight_zz * np.sum(Loss_ZZ)
        print(loss)
        return loss, [Loss_XX, Loss_ZZ]


def train_MultiAE(model,  # type: MultiAE
                  X, epochs, num_omics, num_sample, lr=1e-3,
                  weight_xx=1.0, weight_zz=1.0, device='cpu', show_freq=10):
    for k in range(num_omics):
        if not isinstance(X[k], torch.Tensor):
            X[k] = torch.tensor(np.array(X[k]), dtype=torch.float32, device=device)

    optimizer = optim.Adam(model.parameters(), lr=lr)
    loss_record = []
    for epoch in range(epochs):
        Z, X_recon, Z_recon = model(X)
        loss, loss_item = model.loss_fn(X, Z, X_recon, Z_recon, weight_xx, weight_zz)
        loss_record.append(loss_item)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        if epoch % show_freq == 0 or epoch == epochs - 1:
            # loss X-to-X reconstruction
            print('Epoch {}: Loss_XX: {}'.format(epoch,
                                                 [np.round(l.detach().numpy() / num_sample, 3) for l in loss_item[0]]))
            # loss Z-to-Z reconstruction
            print('Epoch {}: Loss_ZZ: {}'.format(epoch,
                                                 [np.round(l.detach().numpy() / num_sample, 3) for l in loss_item[1]]))

            k = 0  # show layer index
            # row = np.random.permutation(num_sample)
            row = range(num_sample)
            col = 10
            print(X[k].numpy()[row[0:3], col], X_recon[k].detach().numpy()[row[0:3], col])
            print(Z[k].detach().numpy()[row[0:3], col], Z_recon[k].detach().numpy()[row[0:3], col])

    return loss_record


def load_omics_data(omics_files):
    omics_list = []
    for file in omics_files:
        temp_omics = pd.read_csv(file, header=0, index_col=0).T
        omics_list.append(temp_omics)
    return omics_list


def keep_high_var_features(omics_list, num_features=2000):
    retained_omics_list = []
    for i in range(len(omics_list)):
        temp_omics = omics_list[i]
        if temp_omics.shape[1] > num_features:
            features_vars = temp_omics.var(axis=0)
            threshold = sorted(features_vars, reverse=True)[num_features]
            new_omics = temp_omics.loc[:, features_vars > threshold]
            retained_omics_list.append(new_omics)
        else:
            retained_omics_list.append(temp_omics)

    return retained_omics_list


def normalize_matrix(omics_list, type='z-score'):
    retained_omics_list = []
    for i in range(len(omics_list)):
        temp_omics = omics_list[i]
        if type == 'z-score':
            new_omics = preprocessing.scale(temp_omics, axis=0)
            retained_omics_list.append(new_omics)
        elif type == 'min-max':
            new_omics = preprocessing.minmax_scale(temp_omics, axis=0)
            retained_omics_list.append(new_omics)
        else:
            print("Error! required z-score or min-max")

    return retained_omics_list


def PCA_preprocessing(omics_list, explain_variance=0.6):
    retained_omics_list = []
    for i in range(len(omics_list)):
        temp_omics = omics_list[i]
        pca = PCA(n_components=explain_variance)
        new_omics = pca.fit_transform(temp_omics)
        retained_omics_list.append(new_omics)

    return retained_omics_list


if __name__ == "__main__":

    # ========================== Set file dir and device ==========================
    pretrained_model_dir = 'Pretrained Models'
    if not os.path.exists(pretrained_model_dir):
        os.makedirs(pretrained_model_dir)
    # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # ========================== load data ==========================
    omics_names = ['gene', 'methy', 'mirna']
    omics_files = ['SimData1_gene_expression.csv', 'SimData1_DNA_methylation.csv', 'SimData1_mirna_expression.csv']
    omics_list = load_omics_data(omics_files=omics_files)
    # keep features with high variance
    omics_list = keep_high_var_features(omics_list, num_features=2000)
    # normalization feature(z-score or min-max, default: z-score)
    omics_list = normalize_matrix(omics_list, 'min-max')
    # PCA dimension reduction
    # omics_list = PCA_preprocessing(omics_list)

    omics_var = pd.DataFrame(index=['gene', 'methyl', 'mirna'], columns=['min', 'max', 'mean', 'median'])
    for i in range(len(omics_list)):
        omics_var.iloc[i, 0] = np.min(omics_list[i].var(axis=0))
        omics_var.iloc[i, 1] = np.max(omics_list[i].var(axis=0))
        omics_var.iloc[i, 2] = np.mean(omics_list[i].var(axis=0))
        omics_var.iloc[i, 3] = np.median(omics_list[i].var(axis=0))
        print(omics_var.index[i], omics_list[i].shape)
    print('omics var view: \n', omics_var)

    # ========================== Network and optimization parameters ==========================
    num_sample = omics_list[0].shape[0]
    input_dims = []
    omics_weights = []  # Contribution to the Coefficient of the different omics
    for k in range(len(omics_list)):
        z_dim = min(1024, omics_list[k].shape[1])
        input_dims.append([omics_list[k].shape[1], z_dim, 512])
        omics_weights.append(1 / len(omics_list))
    layer_weights = [1.0, 1.0, 1.0]
    epochs = 1000
    weight_xx = 1.0  # X reconstruction weight
    weight_zz = 0.5  # Z reconstruction weight
    num_omics = len(omics_list)
    show_freq = 10

    for k in range(3):
        omics = omics_names[k]
        print("the {}th omics: {}".format(k, omics_names[k]))  # the kth omics data

        # ========================== Cycle AutoEncoder ==========================
        # k = 0  # the kth omics data
        # omics = omics_names[k]
        print('cycle autoencoder')
        cycle_ae = CycleAE(input_dims=input_dims[k])
        cycle_ae.to(device)

        loss_record = train_CycleAE(cycle_ae, omics_list[k], epochs, lr=1e-3,
                                    weight_xx=weight_xx, weight_zz=weight_zz, device=device, show_freq=100)
        # torch.save(omics_ae.state_dict(), save_dir + '/%s-model.ckp' % db)
        torch.save(cycle_ae.state_dict(), pretrained_model_dir + '/%s-omics cae.pkl' % omics)
        # loss_xx = [loss_record[i][0].detach().numpy() for i in range(epochs)]
        # loss_zz = [loss_record[i][1].detach().numpy() for i in range(epochs)]
        # total_loss = [weight_xx * loss_xx[i] + weight_zz * loss_zz[i] for i in range(epochs)]
        # plt.figure()
        # plt.subplot(221)
        # plt.plot(range(epochs), np.clip(loss_xx, 0, 10), label='x-x')
        # plt.legend()
        # plt.subplot(222)
        # plt.plot(range(epochs), np.clip(loss_zz, 0, 10), label='z-z')
        # plt.legend()
        # plt.subplot(212)
        # plt.plot(range(epochs), np.clip(total_loss, 0, 10), label='total_loss')
        # plt.legend()
        # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        # plt.show()
        # plt.savefig('{} omics cae training loss.png'.format(omics))
        #
        # X = torch.tensor(np.array(omics_list[k]), dtype=torch.float32, device=device)
        # Z, X_recon, Z_recon = cycle_ae.forward(X)
        # loss1 = F.mse_loss(X_recon, X, reduction='mean')
        # loss2 = F.mse_loss(Z_recon, Z, reduction='mean')
        # loss = weight_xx * loss1 + weight_zz * loss2
        # print('the {} omics xx loss: {}'.format(omics, loss1))
        # print('the {} omics zz loss: {}'.format(omics, loss2))
        # print('the {} omics total loss: {}'.format(omics, loss))
        # Z0 = cycle_ae.encoder(X)
        # sns.clustermap(Z0.detach().numpy())
        # plt.savefig('{} omics cae Z.png'.format(omics))

