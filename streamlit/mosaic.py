import streamlit as st
import defaults as DFT
import render
import interface
import compute

from missionbio.h5.constants import DNA_ASSAY, PROTEIN_ASSAY

st.set_page_config(page_title='Mosaic', layout='wide')
interface.init()

# ----------- Loading
file, load_raw, apply_filter, kind, should_save, save_name, info = render.files()
if kind == DFT.S3:
    file = compute.download(file)

if file is None:
    interface.error('Please use the options available in the sidebar to load a sample.')

sample = compute.load(file, load_raw, apply_filter)
interface.status('Done.')
interface.info(f'Currently loaded {sample.name}', info)

available_assays = compute.init_assays(sample)

# # ----------- Pre processing
assay_type, clicked, assay_args = render.preprocess(sample, available_assays)

preprocess_funcs = {DNA_ASSAY: compute.preprocess_dna,
                    PROTEIN_ASSAY: compute.preprocess_protein}

for a in available_assays:
    preprocess_funcs[a.name](sample, clicked, *assay_args[a.name])

if clicked:
    interface.rerun()

for a in available_assays:
    if a.name == assay_type:
        assay = a
        break

interface.subheader(f'Analysing {assay_type} | {assay.shape[0]} cells | {assay.shape[1]} ids | {len(set(assay.get_labels()))} clusters')

# # ----------- Data preparation
clicked, scale_attribute, pca_attribute, umap_attribute, pca_comps, info = render.prepare(assay)
if clicked:
    compute.prepare(assay, scale_attribute, pca_attribute, umap_attribute, pca_comps)
    interface.rerun()
interface.info(f'Current transformations are:<br>'
               f'Scale on {assay.metadata[DFT.SCALE_ATTR]}<br>'
               f'PCA on {assay.metadata[DFT.PCA_ATTR]}<br>'
               f'UMAP on {assay.metadata[DFT.UMAP_ATTR]}', info)

# ----------- Clustering
clicked, method, description, cluster_kwargs, info = render.cluster(assay)
if clicked:
    compute.cluster(assay, method, description, **cluster_kwargs)
    interface.rerun()
interface.info(f'Currently clustered using {assay.metadata[DFT.CLUSTER_DESCRIPTION]}', info)

# # # ----------- Customizations
lab_map, pal = render.customize(assay)
compute.customize(assay, lab_map, pal)

# # # ----------- Saving
if should_save:
    compute.save(sample, save_name)

# # # ----------- Visualizations
plot_columns, kind, visualization_kwargs = render.visual(assay, sample)
compute.visual(sample, assay, kind, plot_columns, visualization_kwargs)

interface.status('Done.')
compute.state.init_dna = False
compute.state.init_protein = False
