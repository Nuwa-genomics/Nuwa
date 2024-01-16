import scanpy as sc
import pandas as pd
import numpy as np


class Gene_info:
    def __init__(self, adata):
        self.adata = adata
        self.adata_copy = adata.copy()

        self.fetch_annots()


    def fetch_annots(self):
        self.annot = sc.queries.biomart_annotations(
            "mmusculus",
            ['ensembl_gene_id', 'external_gene_name', 'gene_biotype'],
        )

        self.adata_copy.var["gene_name"] = np.nan
        self.adata_copy.var = self.adata_copy.var.reset_index()
        self.adata_copy.var = self.adata_copy.var.rename(columns={"index": "ensembl_gene_id"})

        for i, row in self.adata_copy.var.iterrows():
            self.adata_copy.var["gene_name"][i] = str(self.annot[self.annot.ensembl_gene_id == row.ensembl_gene_id].external_gene_name.values.astype(str))[2:-2]


    def convert_enseml_to_symbols(self):
        self.adata.var = self.adata.var.set_index(self.adata_copy.var.gene_name)
        self.adata.var.index.name = None

    def convert_symbols_to_ensembl(self):
        self.adata.var = self.adata.var.set_index(self.adata_copy.var.ensembl_gene_id)
        self.adata.var.index.name = None

