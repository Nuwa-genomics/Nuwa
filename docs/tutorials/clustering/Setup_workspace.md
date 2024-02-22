---
sort: 1
---

# Setting up workspace

To start our analysis we first need to set up a workspace environment which will load in our datasets and create other files necessary for storing analysis results. Creating a workspace also allows us to store information relating to our investigation in a database. This includes file locations where results are stored, the different datasets used in the investigation and scripts which allow us to reproduce any analysis in the future.

Creating a workspace just requires a name and description of the investigation:

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/covid_workspace.png'>

Next, in the upload page select the covid dataset. You may wish to learn more about the dataset at [ncbi](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689).

## Obs and var

In the file info in the sidebar we can see there are 9000 obs and 33538 vars. 

```note
## What are obs?
Obs (observations) include things like cell annotations (e.g. tumour or normal cell) or other information relating to the dataset. 
```

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/covid_upload.png'>

In our case we have 3 obs columns: 
- Type: Whether the patient has covid or not
- Sample: Identifier for sample donor
- Batch: Identifyer for sample batch

```note
## What are vars?
Var (variables) refers to annotation of gene/feature metadata. This is a dataframe indexed by unique gene names (or other gene identifier such as ensembl ID).
```
 In our case we have 3 columns:
- gene_ids: Ensembl ID of gene
- feature_types: what the type of data represents (i.e "gene_expression" to signify values are gene counts)
- genome: Reference genome which these genes were mapped to (in our case human reference genome GRCh38, or the 38th build of Genome Reference Consortium human)

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/upload_sidebar.png'>

To see more about the structure of Scanpy's anndata oject see [here](https://cellgeni.readthedocs.io/en/latest/visualisations.html#anndata)

## Other useful functions

- Gene format: the format of gene symbols (the way the .var dataframe is indexed) whether they are ensembl IDs or gene symbols. It is recommended not to change these unless necessary. 

- Obs/var make names unique: Make obs/var names unique

## Other dataset sources

Generally you will be uploading your own datasets. Currently H5AD, H5, loom and mtx are supported. You can also import a dataset from the EBI expression atlas using an accession key.