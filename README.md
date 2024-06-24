# Entropic optimal transport eigenmaps for nonlinear alignment and joint embedding of high-dimensional datasets

Embedding high-dimensional data into a low-dimensional space is an indispensable component in data analysis, commonly used for tasks such as clustering, visualization, and manifold learning. In numerous applications, it is often necessary to align and jointly embed multiple datasets obtained under different experimental conditions. Such datasets may share underlying structures of interest but exhibit individual distortions, resulting in misaligned embeddings using traditional techniques. In this work, we propose *Entropic Optimal Transport (EOT) eigenmaps*, a principled approach for aligning and jointly embedding a pair of datasets with theoretical guarantees. Our approach leverages the leading singular vectors of the EOT plan matrix between two datasets to extract the underlying structure shared between the datasets and align them accordingly in a common embedding space. We interpret our approach as an inter-data variant of the classical Laplacian eigenmaps and diffusion maps embeddings, showing that our approach enjoys many favorable properties analogous to these methods. Under a high-dimensional model where two datasets contain a shared low-dimensional manifold structure, but each dataset is subject to possibly data-specific translation, scaling, nuisance structures, and noise, we show that the EOT plan recovers the shared manifold structure by approximating a kernel function evaluated at the locations of the latent variables. The proposed EOT eigenmaps can be related to the eigenfunctions of certain population-level operators that encode the density and geometry of the shared manifold, admitting a geometric interpretation for the obtained low-dimensional embedding.

The method is based on the paper:

Landa, B., Kluger, Y., and Ma, R. (2023) Entropic Optimal Transport Eigenmaps for Nonlinear Alignment and Joint Embedding of High-Dimensional Datasets


# Content

The R script `main.R` includes the R function for the proposed EOT eigenmap algorithm.

The R script `simulation_align.R` includes the R codes for the noisy manifold alignment simulations studies in the manuscript.

The R script `simulation_cluster.R` includes the R codes for the joint clustering simulations studies in the manuscript.

# System Requirements

The method package requires only a standard computer with enough RAM to support the operations defined by a user. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

The R implementation of the method is tested under R version 4.2.3. Reproducing our simulation results requires the R packages: `Rfast`,`RSpectra`,`BiocNeighbors`,`pcaPP`,`uwot`,`clusterSim`,`ggplot2`,`scatterplot3d`,`fossil`.

