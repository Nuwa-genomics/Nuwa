import os 
from ml.DeepST.deepst.DeepST import run
import matplotlib.pyplot as plt
from pathlib import Path
import scanpy as sc
import pandas as pd
import streamlit as st

class DeepSTModel:
	def __init__(self):
		self.data_path = "ml/DeepST/data/DLPFC"
		self.data_name = '151673' #### project name
		self.save_path = "ml/DeepST/Results" #### save path
		self.n_domains = 7 ###### the number of spatial domains.

		self.deepen = run(save_path = self.save_path,
			task = "Identify_Domain", #### DeepST includes two tasks, one is "Identify_Domain" and the other is "Integration"
			pre_epochs = 800, ####  choose the number of training
			epochs = 1000, #### choose the number of training
			use_gpu = True)
		
		###### Read in 10x Visium data, or user can read in themselves.
		self.adata = self.deepen._get_adata(platform="Visium", data_path=self.data_path, data_name=self.data_name)

		###### Segment the Morphological Image
		self.adata = self.deepen._get_image_crop(self.adata, data_name=self.data_name) 


		###### Data augmentation. spatial_type includes three kinds of "KDTree", "BallTree" and "LinearRegress", among which "LinearRegress"
		###### is only applicable to 10x visium and the remaining omics selects the other two.
		###### "use_morphological" defines whether to use morphological images.
		self.adata = self.deepen._get_augment(self.adata, spatial_type="LinearRegress", use_morphological=True)

	def run(self):
		self.build_graph()
		self.data_process()
		self.train()
                
	def build_graph(self):
		###### Build graphs. "distType" includes "KDTree", "BallTree", "kneighbors_graph", "Radius", etc., see adj.py
		self.graph_dict = self.deepen._get_graph(self.adata.obsm["spatial"], distType = "BallTree")

	def data_process(self):
		###### Enhanced data preprocessing
		self.data = self.deepen._data_process(self.adata, pca_n_comps = 200)

	def train(self):
		###### Training models
		deepst_embed = self.deepen._fit(
				data = self.data,
				graph_dict = self.graph_dict,)
		###### DeepST outputs
		self.adata.obsm["DeepST_embed"] = deepst_embed

		self.adata = self.deepen._get_cluster_data(self.adata, n_domains=self.n_domains, priori = True)

	def get_adata_df(self):
		ax_df = pd.DataFrame({'fr1': self.adata.obsm['spatial'][:,0], 'fr2': self.adata.obsm['spatial'][:,1], 'color': self.adata.obs['DeepST_refine_domain']})
		st.scatter_chart(ax_df, x='fr1', y='fr2', color='color', height=600)

