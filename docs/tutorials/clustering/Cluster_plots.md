---
sort: 4
---

# Cluster plots

We can now go ahead and visualise our data in plots. If you trained a model in the previous chapter, your results will show up here. The learned embeddings will be added to your datasets which will allow you to see clusters representing cell types or other individual gene expression within the sample. We can also use other dimensionality reduction algorithms popular in machine leanring.

```note
### What is dimensionality reduction?

Dimensionality reduction is the transformation of data from a high-dimensional space into a low-dimensional space so that the low-dimensional representation retains some meaningful properties of the original data. In our case we have many variables (genes) which cannot be meaningfully plotting into a chart. We can reduce the thousands of dimensions into 2 or 3 to be represented by our graph axis. 
```

## Autoencoder plots

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cluster_plots_autoencoder_1.png'>

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cluster_plots_autoencoder_2.png'>

## PCA plot

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cluster_plot_pca_1.png'>

## tSNE plot

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cluster_plot_tsne_1.png'>

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cluster_plot_tsne_2.png'>

## UMAP plot

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cluster_plot_umap_1.png'>

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cluster_plot_umap_2.png'>

## PCA variance elbo plot

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cluster_plot_elbo_variance_1.png'>

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cluster_plot_elbo_variance_2.png'>