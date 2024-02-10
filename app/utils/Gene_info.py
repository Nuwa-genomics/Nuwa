import scanpy as sc
import pandas as pd
import numpy as np
from utils.species import *

class Gene_info:
    def __init__(self, adata):
        self.adata = adata
        self.adata_copy = adata.copy()
        self.fail = 0 # no fail status

        self.fetch_annots()


    def fetch_annots(self):
        species_selected = st.session_state.sb_sidebar_species
        if species_selected == "None selected":
            st.toast("Please select a species", icon="⚠️")
            self.fail = -1
            return self.fail
        
        short_name = get_short_species_name_from_long(name=species_selected)
        if short_name == -1:
            st.toast("Something went wrong", icon="❌")
            self.fail = -1
            return self.fail

        self.annot = sc.queries.biomart_annotations(
            short_name,
            ['ensembl_gene_id', 'external_gene_name', 'gene_biotype'],
        )

        self.adata_copy.var = self.adata_copy.var.reset_index()
        
        

    def convert_enseml_to_symbols(self):
        self.adata_copy.var["gene_name"] = np.nan
        self.adata_copy.var = self.adata_copy.var.rename(columns={"index": "ensembl_gene_id"})

        for i, row in self.adata_copy.var.iterrows():
            self.adata_copy.var["gene_name"][i] = str(self.annot[self.annot.ensembl_gene_id == row.ensembl_gene_id].external_gene_name.values.astype(str))[2:-2]

        self.adata.var = self.adata.var.set_index(self.adata_copy.var.gene_name)
        self.adata.var.index.name = None

    def convert_symbols_to_ensembl(self):
        self.adata_copy.var["ensembl_gene_id"] = np.nan
        self.adata_copy.var = self.adata_copy.var.rename(columns={"index": "gene_name"})

        for i, row in self.adata_copy.var.iterrows():
            self.adata_copy.var["ensembl_gene_id"][i] = str(self.annot[self.annot.external_gene_name == row.gene_name].ensembl_gene_id.values.astype(str))[2:-2]

        self.adata.var = self.adata.var.set_index(self.adata_copy.var.ensembl_gene_id)
        self.adata.var.index.name = None

