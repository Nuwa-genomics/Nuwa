# Nuwa 🧬🐍
A bioinformatics web tool built with scanpy for genomics data processing and analysis. 

### \*\*Work in progress!\*\*

Deep neural networks have many potential use cases for genomic analyses including quality control, dimensionality reduction and spatial transcriptomics. Nuwa aims to integrate some of these models in a visual, easy to use interface with other filtering and data analysis familiar to most scanpy users. 

## Getting Started

First, clone the repo:
```bash
git clone https://github.com/ch1ru/Nuwa.git && cd Nuwa
```

The easiest way to get started is using docker compose:
```bash
docker-compose up -d --build
```
You can also use [docker desktop](https://www.docker.com/products/docker-desktop/)

Then visit http://localhost in your browser.

## Preprocess

Filter genes and cell metrics, find mitochrondrial and ribosomal genes, look at variability in gene expression.

![preprocess](screenshots/Preprocess.png "Preprocess data")

## Build model

Build an deep autoencoder based on [Cite-seq model](https://github.com/naity/citeseq_autoencoder) for cluster analysis. Automatically selects a Cuda capable GPU for faster training if one is available.

![build model](screenshots/model.png "Build Model")

## Analysis

Currently analysis consists of:
- Autoencoder cluster plot
- Principal Component Analysis of selected genes
- Variance ratio of principal components
- Neighbourhood graph

![Analysis](screenshots/analysis.png "Analysis")

## Future work

- Integrate postgresql database for persistent record keeping
- Add logging, make visible on ui
- Add other autoencoder models
- support other files types
- Add other analysis scores/graphs