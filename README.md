# Deep-Latent-Space-Fusion
In this work, we propose a novel deep neural network architecture for integrating multi-omics data by learning the consistent manifold in a sample space, termed Deep Latent Fusion Analysis (DLSF). DLSF is built upon a cycle deep autoencoder with a shared self-expressive layer, which can naturally and adaptively merge nonlinear features at each omics level into one unified sample manifold and produce adaptive representation at multi-omics level, thus alleviating the ineffective problems.

In DLSF-main.py, we first pretrain the cycle-autoeocders and save the parameters of the cycle-autoencoder network, then we train the cycle-autoencoders with a self-expressive layer to obtain global affininity matrix and the clusterings of multi-omics data.

All the processed raw data are downloaded from http://acgt.cs.tau.ac.il/multi omic benchmark/download.html. The CLL data are available at the European GenomenPhenome Archive under accession EGAS00001001746 and data tables as R objects can be downloaded from http://pace.embl.de/
