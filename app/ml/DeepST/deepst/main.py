import os 
from ml.DeepST.deepst.DeepST import run
import matplotlib.pyplot as plt
from pathlib import Path
import scanpy as sc
import pandas as pd
import streamlit as st

class DeepSTModel:
	def __init__(self, task = "Identify_Domain", n_domains=7, pre_epochs=800, epochs=1000, use_gpu=True, platform="Visium", dist_type = "BallTree", n_pca_comps = 200):
		self.data_path = "ml/DeepST/data/DLPFC"
		self.data_name = '151673' #### project name
		self.save_path = "ml/DeepST/Results" #### save path
		self.n_domains = n_domains
		self.task=task
		self.pre_epochs = pre_epochs
		self.epochs = epochs
		self.platform = platform
		self.use_gpu = use_gpu
		self.dist_type = dist_type
		self.n_pca_comps = n_pca_comps

	def extract_image_features(self):

		self.deepen = run(save_path = self.save_path,
			task = self.task, #### DeepST includes two tasks, one is "Identify_Domain" and the other is "Integration"
			pre_epochs = self.pre_epochs, ####  choose the number of training
			epochs = self.epochs, #### choose the number of training
			use_gpu = self.use_gpu)

		###### Read in 10x Visium data, or user can read in themselves.
		self.adata = self.deepen._get_adata(platform=self.platform, data_path=self.data_path, data_name=self.data_name)

		###### Segment the Morphological Image
		self.adata = self.deepen._get_image_crop(self.adata, data_name=self.data_name, callback=self.callback) 

	def augment(self):
		###### Data augmentation. spatial_type includes three kinds of "KDTree", "BallTree" and "LinearRegress", among which "LinearRegress"
		###### is only applicable to 10x visium and the remaining omics selects the other two.
		###### "use_morphological" defines whether to use morphological images.
		self.adata = self.deepen._get_augment(self.adata, spatial_type="LinearRegress", use_morphological=True)

	def run(self, callback):
		self.callback = callback
		callback(percent=1, text="Tiling image")
		self.extract_image_features()
		callback(percent=26, text="Augmenting data")
		self.augment()
		callback(percent=35, text="Building graph")
		self.build_graph()
		callback(percent=42, text="Preprocessing data")
		self.data_process()
		callback(percent=50, text="Training model")
		self.train(callback=callback)
		self.callback(text=f"Your task has been completed", percent=100)
		return self
                
	def build_graph(self):
		###### Build graphs. "distType" includes "KDTree", "BallTree", "kneighbors_graph", "Radius", etc., see adj.py
		self.graph_dict = self.deepen._get_graph(self.adata.obsm["spatial"], distType = self.dist_type)

	def data_process(self):
		###### Enhanced data preprocessing
		self.data = self.deepen._data_process(self.adata, pca_n_comps = self.n_pca_comps)

	def train(self, callback):
		###### Training models
		deepst_embed = self.deepen._fit(
				data = self.data,
				graph_dict = self.graph_dict, callback=callback)
		###### DeepST outputs
		self.adata.obsm["DeepST_embed"] = deepst_embed

		self.adata = self.deepen._get_cluster_data(self.adata, n_domains=self.n_domains, priori = True)

	def get_adata_df(self):
		ax_df = pd.DataFrame({'fr1': self.adata.obsm['spatial'][:,0], 'fr2': self.adata.obsm['spatial'][:,1], 'color': self.adata.obs['DeepST_refine_domain']})
		return ax_df

