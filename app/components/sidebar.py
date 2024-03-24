import os
import streamlit as st
import scanpy as sc
from database.schemas import schemas
from scipy import io

from models.AdataModel import AdataModel
from database.database import SessionLocal
from sqlalchemy.orm import Session
from state.AdataState import AdataState
from pathlib import Path
from utils.Gene_info import Gene_info
from utils.species import *
import numpy as np
from state.StateManager import StateManager
from enums.ErrorMessage import ErrorMessage, WarningMessage


class Sidebar:
    """
    Contains utility functions for dataset including saving and loading, exporting python scripts of analysis and adding notes.
    """
    def __init__(self):
        self.conn: Session = SessionLocal()
        self.state_manager = StateManager()

    def show_preview(self, integrate=False):
        """
        Show preview of anndata object attributes.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/show_preview.png

        Example
        -------
        # Same as __repr__ attribute of anndata
        print(adata)
        """
        with st.sidebar:
            with st.expander(label="Show Preview"):
                st.subheader("Anndata preview")
                if integrate:
                    st.markdown(f"""<h3 style='font-size: 16px; color: rgba(255, 255, 255, 0.75)'>{st.session_state.adata_ref.adata_name}</h3>""", unsafe_allow_html=True)
                    st.markdown(f"""<p style='font-size: 14px; color: rgba(255, 255, 255, 0.75)'>{st.session_state.adata_ref.adata}</p>""", unsafe_allow_html=True)
                    st.markdown(f"""<h3 style='font-size: 16px; color: rgba(255, 255, 255, 0.75)'>{st.session_state.adata_target.adata_name}</h3>""", unsafe_allow_html=True)
                    st.markdown(f"""<p style='font-size: 14px; color: rgba(255, 255, 255, 0.75)'>{st.session_state.adata_target.adata}</p>""", unsafe_allow_html=True)
                else:
                    st.markdown(f"""<p style='font-size: 14px; color: rgba(255, 255, 255, 0.75)'>{self.state_manager.get_current_adata()}</p>""", unsafe_allow_html=True)

    def delete_experiment_btn(self):
        with st.sidebar:
            delete_btn = st.button(label="🗑️ Delete Experiment", use_container_width=True, key="btn_delete_adata", type='primary')
            if delete_btn:
                self.state_manager.adata_state() \
                .delete_record(adata_name=st.session_state.sb_adata_selection)
                
                
    def export_script(self):
        """
        Export a python or R script containing any processing and analysis done in the web app. This is useful for reproducibility and documenting analysis done on data.
        
        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/export_script.png
        """
        with st.sidebar:
            with st.expander(label="Export script"):
                python, r = st.tabs(['Python', 'R'])
                with python:
                    python_scripts = st.session_state.script_state.load_script(language = "python")
                    if python_scripts:
                        full_python_script = ""
                        for script in python_scripts:
                            full_python_script += script + '\n'
                        st.code(full_python_script, language="python", line_numbers=True)
                    else:
                        st.warning("Python script not available")
                with r:
                    r_scripts = st.session_state.script_state.load_script(language = "R")
                    if r_scripts:
                        full_r_script = ""
                        for script in r_scripts:
                            full_r_script += script + '\n'
                        st.code(full_r_script, language="r", line_numbers=True)
                    else:
                        st.warning("R script not available")


    def notes(self):
        """
        Add notes to current experiment.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/notes.png
        """
        with st.sidebar:
            try:
                load_adata = self.state_manager.adata_state() \
                .load_adata(workspace_id=st.session_state.current_workspace.id, adata_name=st.session_state.sb_adata_selection)
            
                notes = load_adata.notes
                display_notes = notes if notes != None else ""
                notes_ta = st.text_area(label="Notes", placeholder="Notes", value=display_notes, key="sidebar_notes")

                self.state_manager.adata_state().current.notes = notes_ta
                self.state_manager.adata_state().update_record()
            except Exception as e:
                st.toast(ErrorMessage.NOTES_FAILED_TO_SAVE.value, icon="❌")
                print("Error: ", e)

    def get_adata(self, adataList, name) -> AdataModel:
        for adata in adataList:
            if adata.adata_name == name:
                return adata
        return 0
    
    
    def add_experiment(self):
        """
        Add a new experiment using a pre-existing dataset.

        Parameters
        ----------
        name: str
            Name of experiment.

        dataset: Anndata
            Dataset to use in new experiment.

        notes: str
            Add notes to new experiment if desired.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/add_experiment.png
        """
        with st.sidebar:
            with st.expander(label="Add experiment", expanded=False):
                try:
                    st.subheader("Create New Adata")
                    st.text_input(label="Name", key="ti_new_adata_name")
                    adata_options = self.state_manager.adata_state().get_adata_options()
                    if isinstance(adata_options, ErrorMessage):
                        st.toast(adata_options.value, icon="❌")
                        return
                    st.selectbox(label="Dataset", key="sb_new_experiment_dataset", options=adata_options)
                    st.text_input(label="Notes", key="ti_new_adata_notes")
                    st.button(label="Save", on_click=self.write_adata, key="btn_add_adata")
                except Exception as e:
                    print("Error: ", e)
                    st.error(e)


    def steps(self):
        with st.sidebar:
            with st.form("form_session_steps"):
                current_adata_id = self.state_manager.adata_state().current.id
                sessions = self.conn.query(schemas.Session).filter(schemas.Session.adata_id == current_adata_id).all()
                if sessions:
                    session_names = [session.description for session in sessions]
                
                    steps = st.select_slider(label="Pipeline steps", options=session_names)
                    submit_btn = st.form_submit_button(label="🔧 Undo", use_container_width=True)

                    if submit_btn:
                        print("hi")
                    
                
    def download_adata(self):
        """
        Download adata in h5ad or seurat format.

        Parameters
        ----------
        seurat_format: bool
            Save in seurat format.

        download_directory: str
            File path where the data will be saved. This must be in the 'streamlit-volume' directory to be accessible on the host's filesystem.

        Notes
        -----
        .. image:: https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/download_adata.png
        """
        with st.sidebar:
            with st.expander(label="Download adata file", expanded=False):
                try:
                    st.checkbox(label="Use Seurat format", value=False, key="cb_seurat_format")
                    st.text_input(label="Download directory", value=os.path.join(os.getenv('WORKDIR'), 'downloads', st.session_state.sb_adata_selection), key="ti_save_adata_dir")
                    save_adata_btn = st.button(label="Download", key="btn_save_adata")
                    if save_adata_btn:
                        selected_adata = self.state_manager.adata_state() \
                        .load_adata(workspace_id=st.session_state.current_workspace.id, adata_name=st.session_state.sb_adata_selection)

                        if isinstance(selected_adata, ErrorMessage):
                            st.toast(selected_adata.value, icon="❌")
                        else:
                            download_path = os.path.join(os.getenv('WORKDIR'), 'downloads', st.session_state.sb_adata_selection)
                            if not os.path.exists(download_path):
                                os.mkdir(download_path)
                            if st.session_state.ti_save_adata_dir.find('streamlit-volume') == -1:
                                raise Exception("Download filename must be within the 'streamlit-volume' directory")
                            
                            if st.session_state.cb_seurat_format:
                                if not os.path.isdir(os.path.join(st.session_state.ti_save_adata_dir, 'seurat')):
                                    os.mkdir(os.path.join(st.session_state.ti_save_adata_dir, 'seurat'))
                                with st.spinner(text="Converting adata into seurat"):
                                    #matrix
                                    io.mmwrite(os.path.join(st.session_state.ti_save_adata_dir, 'seurat', 'matrix'), selected_adata.adata.X.T)
                                    #barcodes
                                    with open(os.path.join(st.session_state.ti_save_adata_dir, 'seurat', 'barcodes.tsv'), 'w') as f:
                                        for item in selected_adata.adata.obs_names:
                                            f.write(item + '\n')
                                    #features
                                    with open(os.path.join(st.session_state.ti_save_adata_dir, 'seurat', 'features.tsv'), 'w') as f:
                                        for item in selected_adata.adata.var_names:
                                            f.write(item + '\n')
                                    #metadata
                                    selected_adata.adata.obs.to_csv(os.path.join(st.session_state.ti_save_adata_dir, 'seurat', 'metadata.csv'))
                                    st.toast("Downloaded file", icon='✅')
                                        
                            else:
                                with st.spinner(text="Saving adata"):
                                    sc.write(filename=os.path.join(st.session_state.ti_save_adata_dir, f"{selected_adata.adata_name}.h5ad"), adata=selected_adata.adata)
                                st.toast("Downloaded file", icon='✅')
                except Exception as e:
                    print("Error ", e)
                    st.toast(e, icon="❌")


    def write_adata(self):
        try:
            name = st.session_state.ti_new_adata_name
            notes = st.session_state.ti_new_adata_notes

            self.state_manager.adata_state().insert_record(AdataModel(
                work_id=st.session_state.current_workspace.id,
                adata_name=name,
                filename=os.path.join(os.getenv('WORKDIR'), "adata", f"{name}.h5ad"),
                notes=notes,
            ))
            
        except Exception as e:
            st.error(e)
            print("Error: ", e)


    @staticmethod
    def species():
        """
        Select species which gene profile belongs to. This is needed when changing gene format if spacies cannot be inferred from gene names. Available species are H. sapiens, M. musculus and D. rerio.

        Parameters
        ----------
        species: str
            Species name in format <genus abbr> <species>. 
        """
        st.selectbox(label="Species", key="sb_sidebar_species", options=np.append('None selected', get_species_names_long()), index=infer_species()[0]+1)
        


    @staticmethod
    def gene_format():
        """
        Change the format from gene symbols to ensembl ID or vice versa. This uses the biomart external api in scanpy to fetch gene info. Available for H. sapiens, M. musculus and D. rerio.

        Parameters
        ----------
        ensembl_id: bool
            Toggle to switch between gene symbols and ensembl IDs.
        """
        try:
            format = ""
            for var in StateManager().get_current_adata().var_names[:5]:
                if not var.startswith('ENS'):
                    format = "gene_symbol"
                else:
                    format = "ensembl"

            def change_gene_format():
                gene_info = Gene_info(StateManager().get_current_adata())
                if gene_info.fail == -1:
                    st.session_state["toggle_gene_format"] = not st.session_state["toggle_gene_format"]
                    return -1
                if st.session_state.toggle_gene_format: #change to gene ensembl format
                    gene_info.convert_symbols_to_ensembl()
                else:
                    gene_info.convert_enseml_to_symbols() #change format to gene symbols
                    
                StateManager().get_current_adata().var = gene_info.adata.var

                if st.session_state.toggle_gene_format:
                    st.toast("Changed format to Ensembl IDs", icon='✅')
                else:
                    st.toast("Changed format to gene symbols", icon='✅')

            with st.sidebar:
                st.subheader("Gene format")
                st.toggle(label="Ensembl ID", value=(format == "ensembl"), key="toggle_gene_format", on_change=change_gene_format)

        except Exception as e:
            st.toast(e, icon="❌")
            
            
    def show_version(self):
        with st.sidebar:
            st.markdown(f"""<div style='margin-top: 10px; margin-left: 5px;'>
                        <div style='font-size: 16px; color: rgba(255, 255, 255, 0.4)'>Nuwa v{os.getenv('NUWA_VERSION')}</div>
                        </div>""", unsafe_allow_html=True)



    def show(self, integrate = False):
        with st.sidebar:
            def set_adata():
                switch_adata_result = self.state_manager.adata_state().switch_adata(st.session_state.sb_adata_selection)
                if isinstance(switch_adata_result, ErrorMessage):
                    st.toast(switch_adata_result.value, icon="❌")
                
            
            if integrate:
                options = self.state_manager.adata_state().get_adata_options()
                index = self.state_manager.adata_state().get_index_of_current()
                if isinstance(options, ErrorMessage):
                    st.toast(options.value, icon="❌")
                    return
                if isinstance(index, ErrorMessage):
                    st.toast(options.value, icon="❌")
                    return
                st.selectbox(
                    label="Current Experiment (reference adata):", 
                    options=options, 
                    key="sb_adata_selection", 
                    on_change=set_adata, 
                    index=index
                )
                st.selectbox(
                    label="Integrate into:", 
                    options=options, 
                    key="sb_adata_selection_target"
                )
                    
                #set integrate adata
                
                st.session_state['adata_ref'] = self.state_manager.adata_state().current
                st.session_state['adata_target'] = self.state_manager.adata_state().load_adata(
                    workspace_id=st.session_state.current_workspace.id, 
                    adata_name=st.session_state.sb_adata_selection_target
                )
            
                #sync genes
                #test var names
                if len(st.session_state.adata_ref.adata.var_names) == len(st.session_state.adata_target.adata.var_names):
                    matching = True
                    for i, gene in enumerate(st.session_state.adata_ref.adata.var_names):
                        if gene != st.session_state.adata_target.adata.var_names[i]:
                            matching = False
                            break
                else:
                    matching = False
                    
                empty_msg = st.empty()
                
                if st.session_state.adata_ref.adata_name == st.session_state.adata_target.adata_name:
                    st.warning("Datasets can't be the same")
                    st.session_state["sync_genes"] = False
                elif matching:
                    st.session_state["sync_genes"] = True
                    empty_msg.success("Gene names match, ready to proceed to integrating datasets.")
                else:
                    st.session_state["sync_genes"] = False
                    empty_msg.warning("Datasets have mismatched var names. Toggle 'Sync genes' to reduce the genes to only those present in both datasets.")
                    
                toggle = st.toggle(label="Sync genes", key="toggle_sync_genes", value=False, disabled=(st.session_state.adata_ref.adata_name == st.session_state.adata_target.adata_name))
                
                empty = st.empty()
                if toggle:
                    #use intersecting var names
                    var_names = st.session_state.adata_ref.adata.var_names.intersection(st.session_state.adata_target.adata.var_names)
                    st.session_state['adata_ref'].adata = st.session_state.adata_ref.adata[:, var_names]
                    st.session_state['adata_target'].adata = st.session_state.adata_target.adata[:, var_names]
                    st.session_state["sync_genes"] = True
                    empty_msg.success("Gene names match, ready to proceed to integrating datasets.")
                    empty.markdown(f"""<p style='font-size: 15px; color: rgba(255, 255, 255, 0.75)'>Current experiment {len(st.session_state.adata_ref.adata.var_names)} genes</p><p style='font-size: 15px; color: rgba(255, 255, 255, 0.75)'>Target experiment {len(st.session_state.adata_target.adata.var_names)} genes</p>""", unsafe_allow_html=True)
                else:
                    #set integrate adata
                    st.session_state['adata_ref'] = self.state_manager.adata_state().current
                    st.session_state['adata_target'] = self.state_manager.adata_state().load_adata(
                        workspace_id=st.session_state.current_workspace.id, 
                        adata_name=st.session_state.sb_adata_selection_target
                    )
                    empty.markdown(f"""<p style='font-size: 15px; color: rgba(255, 255, 255, 0.75)'>Current experiment {len(st.session_state.adata_ref.adata.var_names)} genes</p><p style='font-size: 15px; color: rgba(255, 255, 255, 0.75)'>Target experiment {len(st.session_state.adata_target.adata.var_names)} genes</p>""", unsafe_allow_html=True)
                    
                st.divider()

            else:
                options = self.state_manager.adata_state().get_adata_options()
                if isinstance(options, ErrorMessage):
                    st.toast(options.value, icon="❌")
                    return
                st.selectbox(label="Current Experiment:", options=options, key="sb_adata_selection", on_change=set_adata, index=self.state_manager.adata_state().get_index_of_current())

                Sidebar.species()
                
                Sidebar.gene_format()
                
            
            self.download_adata()
            self.add_experiment()
            self.notes()
            self.show_preview()
            self.export_script()
            

            
            
