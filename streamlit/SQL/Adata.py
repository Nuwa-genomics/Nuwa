SAVE_ADATA = lambda work_id, adata_name, filename, notes: \
    f"INSERT INTO adata (work_id, adata_name, filename, notes) VALUES({work_id}, '{adata_name}', '{filename}', '{notes}')"