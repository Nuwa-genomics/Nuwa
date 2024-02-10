import streamlit as st

species = dict(
    ENSG={'short': 'hsapiens', 'long': 'H. Sapiens (human)'},
    ENSMUSG={'short': 'mmusculus', 'long': 'M. musculus (mouse)'},
    ENSDARG={'short': 'drerio', 'long': 'D. Rerio (zebrafish)'},
)

def get_species_names_short():
    species_list = []
    for sp in species.items():
        species_list.append(sp[1]["short"])
    return species_list

def get_species_names_long():
    species_list = []
    for sp in species.items():
        species_list.append(sp[1]["long"])
    return species_list

def infer_species():
    sample_gene = st.session_state.adata_state.current.adata.var_names[0]
    for i, k in enumerate(species.keys()):
        if sample_gene.startswith(k):
            return i, species[k]
    return -1, -1 # not found

def get_short_species_name_from_long(name: str):
    for item in species.items():
        if item[1]["long"].__contains__(name):
            return item[1]["short"]
    return -1


