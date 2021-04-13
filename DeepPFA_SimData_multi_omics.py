"""
Deep Pattern Fusion Analysis for multi-omics data (Deep multi-modal subspace clustering)
By Chengming Zhang (zhangchengming2017@sibcb.ac.cn), July 7, 2020.
All rights reserved.

architecture:
    X_1 -> Z_1 -> C*Z_1 -> X_1-recon -> Z_1
    X_2 -> Z_2 -> C*Z_2 -> X_2-recon -> Z_2
    X_3 -> Z_3 -> C*Z_3 -> X_3-recon -> Z_3

"""
import os
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
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from sklearn.decomposition import PCA
import scipy.io as sio
import math
from snf import compute
from sklearn import cluster
from scipy.spatial import distance
from sklearn.utils.validation import check_array


def silhouette_plot(codes, max_k=15):
    # plt.rcParams['font.sans-serif'] = ['SimHei']
    Scores = []  # 存放轮廓系数
    for k in range(2, max_k):
        estimator = KMeans(n_clusters=k, random_state=555)  # 构造聚类器
        estimator.fit(codes)
        Scores.append(
            silhouette_score(codes, estimator.labels_, metric='euclidean'))
    X = range(2, max_k)
    plt.figure()
    plt.xlabel('Cluster num K', fontsize=15)
    plt.ylabel('Silhouette Coefficient', fontsize=15)
    plt.plot(X, Scores, 'o-')
    plt.show()
    return Scores


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


class SelfExpression(nn.Module):
    def __init__(self, n, knn):
        super(SelfExpression, self).__init__()
        self.knn = knn
        # self.Coefficient = nn.Parameter((1 / n) * torch.ones(n, n, dtype=torch.float32), requires_grad=True)
        self.Coefficient = nn.Parameter(1.0e-8 * torch.ones((n, n), dtype=torch.float32), requires_grad=True)

    def forward(self, z):  # shape=[n, d]
        coef = self.Coefficient
        coef = coef * torch.tensor(self.knn > 0, dtype=torch.float32)
        # coef = coef * torch.tensor(self.knn > 1, dtype=torch.float32)
        coef = coef * (coef > 0)
        self.Coefficient.data = coef
        # coef = coef - torch.diag(torch.diag(coef))
        # coef = torch.abs(coef)
        # coef = 1 / 2 * (coef + coef.T)
        zs = torch.matmul(coef, z)  # z self-expression
        # y = torch.matmul(self.Coefficient, x)
        return zs


class DMSCNet(nn.Module):
    def __init__(self, num_omics, input_dims, num_samples, kernel, w):
        super(DMSCNet, self).__init__()
        self.K = num_omics
        self.n = num_samples
        self.knn = kernel
        self.w = w
        self.CAE = [CycleAE(input_dims[k]) for k in range(num_omics)]
        self.self_expression = SelfExpression(self.n, self.knn)

    def forward(self, X):  # shape=[n, d]
        Z = []
        Z_selfExp = []
        X_recon = []
        Z_recon = []
        for k in range(self.K):
            x_input = X[k]
            z = self.CAE[k].encoder(x_input)
            z_selfExp = self.self_expression(z)
            x_recon = self.CAE[k].decoder(z_selfExp)
            z_recon = self.CAE[k].encoder(x_recon)

            Z.append(z)
            Z_selfExp.append(z_selfExp)
            X_recon.append(x_recon)
            Z_recon.append(z_recon)

        return Z, Z_selfExp, X_recon, Z_recon

    def loss_fn(self, X, X_recon, Z, Z_selfExp, Z_recon, weight_xx, weight_zz, weight_selfExp, weight_coef):
        Loss_XX = []
        Loss_ZZ = []
        Loss_SelfExp = []
        for k in range(self.K):
            x, x_recon, z, z_selfExp, z_recon = X[k], X_recon[k], Z[k], Z_selfExp[k], Z_recon[k]
            loss_xx = F.mse_loss(x_recon, x, reduction='mean')
            loss_zz = F.mse_loss(z_recon, z, reduction='mean')
            loss_selfExp = F.mse_loss(z_selfExp, z, reduction='mean') * self.w[k]
            Loss_XX.append(loss_xx)
            Loss_ZZ.append(loss_zz)
            Loss_SelfExp.append(loss_selfExp)

        # loss_coef = torch.sum(torch.pow(self.self_expression.Coefficient, 2))
        loss_coef = F.mse_loss(torch.sum(self.self_expression.Coefficient, dim=1), torch.ones(self.n), reduction='mean')
        total_loss = weight_xx * np.sum(Loss_XX) + weight_zz * np.sum(Loss_ZZ) + weight_selfExp * np.sum(Loss_SelfExp) \
                     + weight_coef * loss_coef
        # print(loss_ae, loss_coef, loss_selfExp)

        return total_loss, [Loss_XX, Loss_ZZ, Loss_SelfExp, loss_coef]


def train_DMSCNet(model,  # type: DMSCNet
                  X, epochs, num_omics, lr=1e-3, weight_xx=1.0, weight_zz=0.2, weight_selfExp=1.0,
                  weight_coef=1.0, device='cpu', show_freq=10):
    for k in range(num_omics):
        if not isinstance(X[k], torch.Tensor):
            X[k] = torch.tensor(X[k], dtype=torch.float32, device=device)

    optimizer = optim.Adam(model.parameters(), lr=lr)
    losses = []
    losses_item = []
    for epoch in range(epochs):
        Z, Z_selfExp, X_recon, Z_recon = model(X)
        loss, loss_item = model.loss_fn(X, X_recon, Z, Z_selfExp, Z_recon, weight_xx, weight_zz, weight_selfExp,
                                        weight_coef)
        losses.append(loss)
        losses_item.append(loss_item)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        if epoch % show_freq == 0 or epoch == epochs - 1:
            # loss_item: Loss_XX, Loss_ZZ, Loss_SelfExp, loss_coef
            print('Epoch {}: loss_xx: {}'.format(epoch, [loss_xx.detach().numpy() for loss_xx in loss_item[0]]))
            print('Epoch {}: loss_zz: {}'.format(epoch, [loss_zz.detach().numpy() for loss_zz in loss_item[1]]))
            print('Epoch {}: loss_selfExp: {}'.format(epoch, [loss_selfExp.detach().numpy() for loss_selfExp in loss_item[2]]))
            print('Epoch {}: loss_coefficient: {}'.format(epoch, loss_item[3]))
            print('Epoch {}: total loss: {}'.format(epoch, loss))

            sample_id = 10
            for k in range(num_omics):
                print('the {}th omics loss'.format(k+1))
                print(X[k].numpy()[sample_id, 10:14], X_recon[k].detach().numpy()[sample_id, 10:14])
                print(Z[k].detach().numpy()[sample_id, 10:14], Z_recon[k].detach().numpy()[sample_id, 10:14])
                print(Z[k].detach().numpy()[sample_id, 10:14], Z_selfExp[k].detach().numpy()[sample_id, 10:14])
            C = model.self_expression.Coefficient.detach().to('cpu').numpy()
            print('C.min(), C.max(), C.mean()', C.min(), C.max(), C.mean(), )
            print('np.sum(C, axis=1).min(), np.sum(C, axis=1).max()', np.sum(C, axis=1).min(), np.sum(C, axis=1).max())

    return losses, losses_item


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


def knn_kernel(data, k=10):
    # knn sample weight set to 1 for each row
    num_samples = data.shape[0]
    Dis_Mat = distance.cdist(data, data)
    kernel = np.ones_like(Dis_Mat)
    sort_dist = np.sort(Dis_Mat, axis=1)
    threshold = sort_dist[:, k].reshape(-1, 1)
    sig = (Dis_Mat <= np.repeat(threshold, num_samples, axis=1))
    # kernel = sig * kernel
    kernel = sig * kernel - np.identity(num_samples)

    return kernel


def get_n_clusters(arr, n_clusters=range(2, 6)):
    """
    Finds optimal number of clusters in `arr` via eigengap method

    Parameters
    ----------
    arr : (N, N) array_like
        Input array (e.g., the output of :py:func`snf.compute.snf`)
    n_clusters : array_like
        Numbers of clusters to choose between

    Returns
    -------
    opt_cluster : int
        Optimal number of clusters
    second_opt_cluster : int
        Second best number of clusters
    """

    # confirm inputs are appropriate
    n_clusters = check_array(n_clusters, ensure_2d=False)
    n_clusters = n_clusters[n_clusters > 1]

    # don't overwrite provided array!
    graph = arr.copy()

    graph = (graph + graph.T) / 2
    graph[np.diag_indices_from(graph)] = 0
    degree = graph.sum(axis=1)
    degree[np.isclose(degree, 0)] += np.spacing(1)
    di = np.diag(1 / np.sqrt(degree))
    laplacian = di @ (np.diag(degree) - graph) @ di

    # perform eigendecomposition and find eigengap
    eigs = np.sort(np.linalg.eig(laplacian)[0])
    eigengap = np.abs(np.diff(eigs))
    eigengap = eigengap * (1 - eigs[:-1]) / (1 - eigs[1:])
    n = eigengap[n_clusters - 1].argsort()[::-1]

    return n_clusters[n[:2]]


def draw_eigen_plot(arr, max_cluster_num=15):
    """
    Draw the eigenvalues plot to finds optimal number of clusters in `arr` via eigengap method

    """

    # don't overwrite provided array!
    graph = arr.copy()

    graph = (graph + graph.T) / 2
    graph[np.diag_indices_from(graph)] = 0
    degree = graph.sum(axis=1)
    degree[np.isclose(degree, 0)] += np.spacing(1)
    di = np.diag(1 / np.sqrt(degree))
    laplacian = di @ (np.diag(degree) - graph) @ di

    # perform eigendecomposition and find eigengap
    eigs = np.sort(np.linalg.eig(laplacian)[0])
    plt.plot(eigs[0:max_cluster_num])
    plt.scatter(range(max_cluster_num), eigs[0:max_cluster_num])
    plt.xticks(range(max_cluster_num), np.arange(1, max_cluster_num+1))
    plt.title('get the best cluster num by eigengap')
    plt.savefig('get_the_best_cluster_num_by_eigengap.png')


if __name__ == "__main__":

    pretrained_model_dir = 'Pretrained Models'
    if not os.path.exists(pretrained_model_dir):
        os.makedirs(pretrained_model_dir)
    save_model_dir = 'Final Models'
    if not os.path.exists(save_model_dir):
        os.makedirs(save_model_dir)

    # ================================= load data ====================================
    omics_names = ['gene', 'methy', 'mirna']
    omics_files = ['SimData1_gene_expression.csv', 'SimData1_DNA_methylation.csv', 'SimData1_mirna_expression.csv']
    omics_list = load_omics_data(omics_files=omics_files)
    samples_id = omics_list[0].index

    # keep features with high variance
    omics_list = keep_high_var_features(omics_list, num_features=2000)
    # normalization feature(z-score or min-max, default: z-score)
    omics_list = normalize_matrix(omics_list, 'min-max')

    omics_var = pd.DataFrame(index=['gene', 'methyl', 'mirna'], columns=['min', 'max', 'mean', 'median'])
    for i in range(len(omics_list)):
        omics_var.iloc[i, 0] = np.min(omics_list[i].var(axis=0))
        omics_var.iloc[i, 1] = np.max(omics_list[i].var(axis=0))
        omics_var.iloc[i, 2] = np.mean(omics_list[i].var(axis=0))
        omics_var.iloc[i, 3] = np.median(omics_list[i].var(axis=0))
        print(omics_var.index[i], omics_list[i].shape)
    print('omics var view: \n', omics_var)

    kernel_list = []
    num_samples = omics_list[0].shape[0]
    fused_kernel = np.zeros((num_samples, num_samples))
    neighbor_num = round(num_samples/10)
    if neighbor_num < 25:
        neighbor_num = 25
    elif neighbor_num > 50:
        neighbor_num = 50
    for i in range(len(omics_list)):
        omics_kernel = knn_kernel(omics_list[i], k=neighbor_num)
        # fused_kernel += omics_kernel / len(omics_list)
        fused_kernel += omics_kernel
        kernel_list.append(omics_kernel)
    print(Counter(fused_kernel.ravel()))
    # fused_kernel = fused_kernel + fused_kernel
    # load the label data y
    labels_true = np.zeros(400)
    labels_true[0:100] = 1
    labels_true[100:200] = 2
    labels_true[200:300] = 3
    labels_true[300:400] = 4

    # ========================== Network and optimization parameters ==========================
    num_omics = len(omics_files)
    input_dims = []
    omics_weights = []  # Contribution to the Coefficient of the different omics
    for k in range(len(omics_list)):
        z_dim = min(1024, omics_list[k].shape[1])
        input_dims.append([omics_list[k].shape[1], z_dim, 512])
        omics_weights.append(1 / len(omics_list))
    w = [0.5, 0.5, 0.5]
    # w = [0.2, 0.2, 0.2]
    epochs = 100
    lr = 1e-3
    weight_xx = 1.0  # X reconstruction weight
    weight_zz = 0.05  # Z reconstruction weight
    weight_selfExp = 0.2  # self-expression loss weight
    weight_coef = 0.01  # coefficient regularization
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    show_freq = 10

    # ========================== Deep subspace clustering ==========================
    mode = 'fusion'
    print('preprocessing omics: ', mode)
    dmscnet = DMSCNet(num_omics, input_dims, num_samples, fused_kernel, w)
    dmscnet.to(device)

    # load the pretrained weights which are provided by the cycle autoencoder
    for k in range(num_omics):
        cae_state_dict = torch.load('Pretrained Models/%s-omics cae.pkl' % omics_names[k])
        # cae_state_dict = torch.load('pretrained_cae_dir/%s.pkl' % omics)
        dmscnet.CAE[k].load_state_dict(cae_state_dict)
    print("Pretrained cae weights are loaded successfully.")

    losses, losses_item = train_DMSCNet(dmscnet, omics_list, epochs, num_omics, lr=1e-3, weight_xx=weight_xx,
                                        weight_zz=weight_zz, weight_coef=weight_coef, weight_selfExp=weight_selfExp,
                                        device='cpu', show_freq=10)

    torch.save(dmscnet.state_dict(), save_model_dir + '/%s-dmscnet-model.pkl' % mode)

    Loss_XX = [losses_item[i][0] for i in range(epochs)]
    Loss_ZZ = [losses_item[i][1] for i in range(epochs)]
    Loss_SelfExp = [losses_item[i][2] for i in range(epochs)]
    loss_coef = [losses_item[i][3].detach().numpy() for i in range(epochs)]
    total_loss = [losses[i].detach().numpy() for i in range(epochs)]

    # fig, ax = plt.subplots(nrows=5, ncols=num_omics, sharex=True, sharey=False, figsize=(16, 8))

    plt.figure()
    for k in range(num_omics):
        plt.subplot(321)
        plt.plot(range(epochs), [l[k].detach().numpy() for l in Loss_XX], label='x-x-{}'.format(k+1))
        plt.legend()
        plt.subplot(322)
        plt.plot(range(epochs), [l[k].detach().numpy() for l in Loss_ZZ], label='z-z-{}'.format(k+1))
        plt.legend()
        plt.subplot(323)
        plt.plot(range(epochs), [l[k].detach().numpy()/w[k] for l in Loss_SelfExp], label='SelfExp-{}'.format(k+1))
        plt.legend()
    plt.subplot(324)
    plt.plot(range(epochs), np.clip(loss_coef, 0, 10), label='coef-reg')
    plt.legend()
    plt.subplot(313)
    plt.plot(range(epochs), np.clip(total_loss, 0, 10), label='total loss')
    plt.legend()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.show()
    plt.savefig('{} omics deepPFA training loss.png'.format(mode))

    for k in range(num_omics):
        if not isinstance(omics_list[k], torch.Tensor):
            omics_list[k] = torch.tensor(omics_list[k], dtype=torch.float32, device=device)
    Z, Z_selfExp, X_recon, Z_recon = dmscnet.forward(omics_list)

    Z = [z.detach().numpy() for z in Z]

    C = dmscnet.self_expression.Coefficient.detach().to('cpu').numpy()
    C = preprocessing.normalize(C, norm='l1', axis=1)
    plt.figure()
    sns.heatmap(C)
    plt.savefig('{} omics deepPFA heatmap C.png'.format(mode))

    # ========================== Spectral clustering ==========================

    # get the processed C
    alpha = 1
    fused = (C + C.T)/2
    plt.figure()
    sns.heatmap(fused)
    plt.savefig('{} omics deepPFA heatmap sym-C.png'.format(mode))

    # estimate the number of clusters present in the fused matrix, derived via
    # an "eigengap" method (i.e., largest difference in eigenvalues of the
    # laplacian of the graph). note this function returns the top two options;
    # we'll only use the first
    first, second = compute.get_n_clusters(fused)
    print(first, second)
    # first = max(first, second)

    # draw the eigenvalue plot
    plt.figure()
    draw_eigen_plot(fused)

    # apply clustering procedure
    # you can use any clustering method here, but since DSCNet returns an affinity
    # matrix (i.e., all entries are positively-valued and indicate similarity)
    # spectral clustering makes a lot of sense
    labels_pred = cluster.spectral_clustering(fused, n_clusters=first) + 1
    fused_cluster = pd.DataFrame(labels_pred, index=samples_id, columns=['cluster'])
    fused_cluster.to_csv("fused_cluster_first.csv")

    print("spectral-first:", Counter(labels_pred))
    coef_ARI = metrics.adjusted_rand_score(labels_true, labels_pred)
    coef_AMI = metrics.adjusted_mutual_info_score(labels_true, labels_pred)
    coef_NMI = metrics.normalized_mutual_info_score(labels_true, labels_pred)
    print("spectral coef_ARI", coef_ARI)  # 0.5287292479428645 (-1, 1)
    print("spectral coef_AMI", coef_AMI)  # 0.5287292479428645 (-1, 1)
    print("spectral coef_NMI", coef_NMI)  # 0.5287292479428645 (0, 1)


