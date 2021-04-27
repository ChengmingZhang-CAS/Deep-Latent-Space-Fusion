# Deep-Pattern-Fusion-Analysis
In this work, we propose a novel deep neural network architecture for integrating multi-omics data by learning the consistent manifold in a sample space, termed Deep Pattern Fusion Analysis (DeepPFA). DeepPFA is built upon a cycle deep autoencoder with a shared self-expressive layer, which can naturally and adaptively merge nonlinear features at each omics level into one unified sample manifold and produce adaptive representation at multi-omics level, thus alleviating the ineffective problems.

Firstly, run the Pretrained_CAE_SimData_multi_omics.py to pretrain the cycle-autoeocders. Then run the DeepPFA_SimData_multi_omics.py to obtain global affininity matrix and the clusterings of multi-omics data.

All the processed TCGA raw data are available at http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html.
