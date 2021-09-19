# Construct simulation data from TCGA Melanoma data by SVD
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.mstats import winsorize


def load_omics_data(omics_files):
    omics_list = []
    for file in omics_files:
        temp_omics = pd.read_csv(file, header=0, index_col=0).T
        omics_list.append(temp_omics)
    return omics_list


# load the three omics datasets
omics_names = ['gene', 'methy', 'mirna']
omics_files = ['aligned_exp.csv', 'aligned_methy.csv', 'aligned_mirna.csv']
data_gene = pd.read_csv('aligned_exp.csv', header=0, index_col=0)
data_methy = pd.read_csv('aligned_methy.csv', header=0, index_col=0)
data_mirna = pd.read_csv('aligned_mirna.csv', header=0, index_col=0)
print('data_gene.shape: ', data_gene.shape)
print('data_methy.shape: ', data_methy.shape)
print('data_mirna.shape: ', data_mirna.shape)


def construct_simData(i):

    np.random.seed(i)
    num_sample = data_gene.shape[1]
    sample_id = np.random.permutation(num_sample)[0:400]

    idx1 = np.random.permutation(data_gene.shape[0])[0:2000]
    gene_mat = data_gene.iloc[idx1, sample_id]

    idx2 = np.random.permutation(data_methy.shape[0])[0:2000]
    methy_mat = data_methy.iloc[idx2, sample_id]

    idx3 = np.random.permutation(data_mirna.shape[0])
    mirna_mat = data_mirna.iloc[idx3, sample_id]

    # apply SVD to origin omics data
    U1, S1, VT1 = np.linalg.svd(gene_mat, full_matrices=False)
    U2, S2, VT2 = np.linalg.svd(methy_mat, full_matrices=False)
    U3, S3, VT3 = np.linalg.svd(mirna_mat, full_matrices=False)

    # reconstruct the sample matrix pattern
    X1 = np.zeros_like(gene_mat)
    X2 = np.zeros_like(methy_mat)
    X3 = np.zeros_like(mirna_mat)

    # subgroup(1~100; 101~200; 201~300; 301~400)
    mean_set1 = np.random.permutation([0, 2, 4])
    X1[:, 0:100] += mean_set1[0]
    X1[:, 100:300] += mean_set1[1]
    X1[:, 300:400] += mean_set1[2]

    mean_set2 = np.random.permutation([0, 2, 4])
    X2[:, 0:100] += mean_set2[0]
    X2[:, 100:200] += mean_set2[1]
    X2[:, 200:300] += mean_set2[0]
    X2[:, 300:400] += mean_set2[2]

    mean_set3 = np.random.permutation([0, 2, 4])
    X3[:, 0:200] += mean_set3[0]
    X3[:, 200:300] += mean_set3[1]
    X3[:, 300:400] += mean_set3[2]

    # add noise
    X1 += np.random.normal(0, 2, size=X1.shape)
    X2 += np.random.normal(0, 2, size=X2.shape)
    X3 += np.random.normal(0, 2, size=X3.shape)

    # apply SVD to artificial data
    U_sim1, S_sim1, VT_sim1 = np.linalg.svd(X1, full_matrices=False)
    U_sim2, S_sim2, VT_sim2 = np.linalg.svd(X2, full_matrices=False)
    U_sim3, S_sim3, VT_sim3 = np.linalg.svd(X3, full_matrices=False)

    # reconstruct the omics data
    # the number of retained singular value (60% explained info)
    k1 = np.argmax(np.cumsum(S1)/S1.sum() > 0.6) + 1
    k2 = np.argmax(np.cumsum(S2)/S2.sum() > 0.6) + 1
    k3 = np.argmax(np.cumsum(S3)/S3.sum() > 0.6) + 1
    gene_sim = (U1[:, 0:k1].dot(np.diag(S1[0:k1]))).dot(VT_sim1[0:k1, :])
    gene_sim += np.random.normal(0, 1e-3*gene_sim.std(), size=gene_sim.shape)
    gene_sim -= np.min(gene_sim, axis=1).reshape(gene_sim.shape[0], -1)
    values = winsorize(gene_sim, limits=[0.0, 0.2], axis=1)  # 同时做下0%和上20%的处理

    methy_sim = (U2[:, 0:k2].dot(np.diag(S2[0:k2]))).dot(VT_sim2[0:k2, :])
    methy_sim += np.random.normal(0, 1e-3*methy_sim.std(), size=methy_sim.shape)
    methy_sim -= np.min(methy_sim, axis=1).reshape(methy_sim.shape[0], -1)

    mirna_sim = (U3[:, 0:k3].dot(np.diag(S3[0:k3]))).dot(VT_sim3[0:k3, :])
    mirna_sim += np.random.normal(0, 1e-3*mirna_sim.std(), size=mirna_sim.shape)
    mirna_sim -= np.min(mirna_sim, axis=1).reshape(mirna_sim.shape[0], -1)

    # sns.heatmap(np.corrcoef(mirna_sim.T) - np.diag(np.ones(400)), cmap='Reds')
    # sns.heatmap(np.corrcoef(mirna_sim.T))
    save_dir = "SimData" + str(i)
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    sns.clustermap(gene_sim, standard_scale=0, col_cluster=False)
    plt.savefig(os.path.join(save_dir, 'SimData'+str(i)+'_gene_expression.png'))
    plt.close()
    sns.clustermap(methy_sim, standard_scale=0, col_cluster=False)
    plt.savefig(os.path.join(save_dir, 'SimData'+str(i)+'_DNA_methylation.png'))
    plt.close()
    sns.clustermap(mirna_sim, standard_scale=0, col_cluster=False)
    plt.savefig(os.path.join(save_dir, 'SimData'+str(i)+'_mirna_expression.png'))
    plt.close()

    gene_file_name = os.path.join(save_dir, 'SimData'+str(i)+'_gene_expression.csv')
    methy_file_name = os.path.join(save_dir, 'SimData'+str(i)+'_DNA_methylation.csv')
    mirna_file_name = os.path.join(save_dir, 'SimData'+str(i)+'_mirna_expression.csv')
    pd.DataFrame(gene_sim).to_csv(gene_file_name)
    pd.DataFrame(methy_sim).to_csv(methy_file_name)
    pd.DataFrame(mirna_sim).to_csv(mirna_file_name)


for i in range(1, 11):
    print("Simdata" + str(i))
    construct_simData(i)
