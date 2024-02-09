# Nuwa

[![](https://dcbadge.vercel.app/api/server/wBDavdWp8n)](https://discord.gg/wBDavdWp8n)
[![Stars](https://img.shields.io/github/stars/ch1ru/nuwa?logo=GitHub&color=yellow)](https://github.com/ch1ru/nuwa/stargazers)
[![GitHub last commit](https://img.shields.io/github/last-commit/ch1ru/Nuwa)](https://github.com/ch1ru/Nuwa/pulse)
![CI tests](https://github.com/ch1ru/nuwa/actions/workflows/run_tests.yml/badge.svg?branch=main)
![docs](https://github.com/ch1ru/nuwa/actions/workflows/jekyll-gh-pages.yml/badge.svg?branch=main)

A bioinformatics web tool built with scanpy for genomics data processing and analysis.

```warning
Currently in beta version, **not** recommended for use in research or commercial use. 
```

## What it does?

Nuwa is an open-source, graphical web tool for single cell RNA seq analysis. It primarily focuses on functions available within the Scanpy library, but also incorporates other software packages and deep learning models. Nuwa runs on a local web server which users access in their browser. 

## Quick start

**Clone the repo:**

```bash
#clone repo
git clone https://github.com/ch1ru/Nuwa.git && cd Nuwa
```

```note
## If using a GPU
If you are planning to use a GPU for faster training:
- Make sure cuda drivers are installed on the host machine.
- Install and configure [Nvidia container toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) for docker
```

**Then bring up containers:**

for GPU:
```bash
docker-compose -f cuda.docker-compose.yml up -d --build
```

Or for CPU only:
```bash
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
