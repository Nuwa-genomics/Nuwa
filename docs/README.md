# Nuwa

[![Stars](https://img.shields.io/github/stars/ch1ru/nuwa?logo=GitHub&color=yellow)](https://github.com/ch1ru/nuwa/stargazers)
![CI](https://github.com/ch1ru/nuwa/actions/workflows/run_tests.yml/badge.svg?branch=main)

A bioinformatics web tool built with scanpy for genomics data processing and analysis.

## What it does?

Nuwa aims to integrate several deep learning models in a visual, easy to use interface with other filtering and data analysis familiar to most scanpy users.

## Quick start

```bash
#clone repo
git clone https://github.com/ch1ru/Nuwa.git && cd Nuwa
#If you have a Nvidia GPU
docker-compose -f cuda.docker-compose.yml up -d --build
#Or if you have a CPU
docker-compose -f cpu.docker-compose.yml up -d --build
```

## Usage

Navigate to http://localhost to use the app. Files are mounted to the streamlit-volume directory in the installation path.

## Features

- Upload scRNA expression matrices in multiple formats (h5ad, mtx, loom supported) including built-in datasets and online datasets using their EBI accession key
- Preprocess data (includes mito, ribo and haem annotation, cell & gene filtering, highly variable genes, data scaling and normalization, cell-cycle scoring, batch effect correction, doublet detection, downsampling, subsampling and automated preprocessing recipes included in scanpy.)
- Integrate different dataset types using BBKNN, scanorama and Ingest
- Visualise differential gene expression using multiple statistical tests across clusters or individual genes
- A variety of deep learning models such as doublet detection with Solo and cluster analysis using Cite-seq autoencoder. Adjustable hyperparameters and hardware devices
- Trajectory inference on data and easily adding PAGA embeddings to data
- Spatial transcriptomics with metrics such as Ripley and Centrality scoring
- 3D interactive chart to view generated clusters
- Terminal to interact with the docker container running a bioconda environment


## License

Nuwa is under MIT license, a short and simple permissive license with conditions only requiring preservation of copyright and license notices. Licensed works, modifications, and larger works may be distributed under different terms and without source code. [See licence](https://github.com/ch1ru/Nuwa/blob/main/LICENSE)
