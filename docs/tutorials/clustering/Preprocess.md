---
sort: 2
---

# Preprocessing

We will now preprocess our raw data which will reduce noise from non-biological signal as well as filter out low quality cells. Part of the preprocessing stage is also tranforming the data such as log-normalizing or scaling data to a given distribution. This will help give more meaningful results and help remove outliers.

## Exploratory data analysis

In order to better understand our data, we will perform the EDA part of our analysis. This will involve using statistical methods including visualizations to provide a summary of different aspects or characteristics.

### Highest expressed genes

First let's see which genes are the most expressed across our dataset:

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/highest_expr.png'>

```tip
Hover over the box plots to see their median, q1 & q3, lowest and highest counts.
```

We can see the most expressed gene by far with a median count of around 3.75, followed by some mitochondrial genes (denoted by the 'MT-' prefix).

```note
### Malat1
Malat1 gene can be a result of unwanted technical noise so can often be removed from the dataset.
```

### Sample sex

If the sex of the single cell donors is unknown, mislabeled, or we wish to verify the sex from the supplementary materials, we can do this by measuring genes found contained in the sex chromosomes. In this example, we measure the levels of the Xist gene across samples (select the subsample tab). Xist is a non-coding RNA found exclusively in female somatic cells (although can also be found in male germline cells) so will work well for our blood samples. It should be noted that while the coding sequence for Xist is on the X chromosome (so is present in both males and females), the RNA transcript is not produced in males and therefore won't appear in our dataset. We can see very clearly the sex of each sample:

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/xist_counts.png'>

Sex of the donors (note that we don't have any male control groups in this particular dataset):

| Sample   | Sex     |
| -------- | ------- |
| covid_1  | male    |
| covid_15 | female  |
| covid_17 | male    |
| ctr_13   | female  |
| ctr_14   | female  |
| ctr_5    | female  |


### Cell cycle states 

It would be useful to know which stages of the cell cycle our samples are in, as this leads to potentially unwanted variation between cells. Here's a brief reminder of the cell cycle states:

<img alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cell_cycle.png'>

*Adapted from Wikipedia (Image License is CC BY-SA 3.0)*

- **G0**: Some cells enter G0 once they are fully differentiated. These would include neuron cells and red blood cells and are permanently arrested in this distinct phase
- **G1**: Gap 1 phase is the beginning of interphase. No DNA synthesis happens here and includes growth of non-chromosomal components of the cell
- **S**: Synthesis of DNA through replication of the chromosomes
- **G2**: The cell enlarges for the preparation of mitosis and the spindle fibres start to form
- **M**: Mitosis phase is the nuclear division of the cell (consisting of prophase, metaphase, anaphase and telophase).

We will be looking at cells in the S and G2 or M phase. To do this we will first need a list of marker genes. To score a gene list, the algorithm calculates the difference of mean expression of the given list and the mean expression of reference genes. To build the reference, the function randomly chooses a bunch of genes matching the distribution of the expression of the given list.

First we will need the [marker gene list](https://raw.githubusercontent.com/Nuwa-genomics/Nuwa/main/app/data/cell_cycle_marker_genes.csv)

Once our marker genes file has been loaded in, ensure the gene and phase columns correspond to the correct column name and select 'Sample' as our group by ID.

```tip
Make sure the gene format of your imported marker genes is the same as the one in your dataset. In our case we are using the gene names (not ensembl IDs). If your marker genes are in the wrong format, either try to convert them, for example with the biomart API. Alternatively convert the dataset genes using the gene format toggle switch in the sidebar.
```

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cell_cycle_panel.png'>

When we run our analysis we see the results as violin plots and scatter plots of the predicted scores:

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cell_cycle_violins.png'>

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/cell_cycle_scores_plot.png'>

You can also see a PCA map with the highlighted scores.


## Quality control

An important part of preprocessing is assessing the quality of our dataset. This means removing low quality cells, or cells which have too low or high counts which may interfere with our analysis. 

### Annotating mitochondrial genes

A popular and simple method of predicting low quality cells is by looking at the percentage of mitochondrial genes in each library/cell sample. A high proportion of mitochondrial genes are indicative of low quality cells (Islam et al. 2014; Ilicic et al. 2016). If cytoplasmic RNA is lost due to perforated cells, this leads to an artificial relative increase in mitochondrial trascripts. Therefore, we will first annotate whether the genes are mitochondrial (denoted by the 'MT-' prefix in gene names) which we can display as scatter or violin plots. 

```warning
Feature names must be in the correct format for detecting mitochondrial/ribosomal/haemoglobin genes (with the MT-, RB- and HB- prefixes respectively). If gene names are in the ensembl format (e.g. ENSG00000139618) you can convert to gene symbols using the gene format toggle on the sidebar.
```

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/pct_mt_scatter.png'>

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/pct_mt_violin.png'>

```tip
Select a color key to compare observations across the dataset. In the above example we can compare individual samples. If no color key is selected the data will be in aggregate.
```

Next, we will remove cells that containing more than 10% mitochondrial genes. This is an example of filtering based on a fixed threshold (since the 10% doesn't take into account the distribution of mitochondrial reads). In the future we will likely implement adaptive thresholds also. 

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/pct_mt_filter_fixed.png'>

### Filtering based on number of genes and cells

Another simple way to filter our data is to only keep cells which contain a minimum number of genes, as well as genes which are expressed in a minimum number of cells. In our case we will only keep:

- cells with at least 200 genes
- genes which appear in at least 3 cells

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/filter_cells.png'>

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/filter_genes.png'>

This helps to remove low quality cells which may contain a lower reading of transcripts due to technical error. Cells with a lower number of features are also less useful in analysis, as are genes which only appear only a few times in the dataset. In general we are trying to make our data more meaningful and less prone to containing technical noise leading to false conclusions.