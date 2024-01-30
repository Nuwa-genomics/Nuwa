## filter_highest_expr_genes
Fraction of counts assigned to each gene over all cells. Computes, for each gene, the fraction of counts assigned to that gene within a cell. The n_top genes with the highest mean fraction over all cells are plotted as boxplots.
![screenshot](https://raw.githubusercontent.com/ch1ru/Nuwa/main/docs/assets/images/screenshots/highest_expr_genes.png)
## Parameters
```n_top_genes: int``` Number of top gene symbols to show.