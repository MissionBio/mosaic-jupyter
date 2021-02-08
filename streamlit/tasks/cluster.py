import streamlit as st

import interface
import defaults as DFT

from missionbio.h5.constants import DNA_ASSAY
from missionbio.mosaic.constants import (
    PCA_LABEL,
    UMAP_LABEL,
    AF_MISSING
)


def run(assay, available_assays):
    clicked, method, description, cluster_kwargs, info = render(assay)

    first_pass_cluster(available_assays)

    if clicked:
        cluster(assay, method, description, **cluster_kwargs)
        interface.rerun()
    interface.info(f'Currently clustered using {assay.metadata[DFT.CLUSTER_DESCRIPTION]}', info)


def render(assay):
    title = 'Clustering'
    if not assay.metadata[DFT.CLUSTERED]:
        title += ' *'

    with st.sidebar.beta_expander(title):
        info = st.empty()

        all_methods = ['dbscan', 'hdbscan', 'graph-community', 'kmeans']
        if assay.name == DNA_ASSAY:
            all_methods.append('count')
        method = st.selectbox('Method', all_methods)

        cluster_kwargs = {}
        if method == 'count':
            layer = st.selectbox('Layer', ['NGT', 'NGT_FILTERED'])
            min_clone_size = st.slider('Minimum clone size (%)', 0.0, 10.0, 1.0)
            group_missing = st.checkbox('Group missing', True)
            ignore_zygosity = st.checkbox('Ignore Zygosity', False)
            features = st.multiselect('Variants', list(assay.ids()))

            description = f'{layer} counts on {len(features)} variants with {min_clone_size}% minimum clone size'
            if ignore_zygosity:
                description += ', ignoring zygosity'
            if group_missing:
                description += ', and grouped missing clones'

            cluster_kwargs = {
                'features': features,
                'layer': layer,
                'group_missing': group_missing,
                'min_clone_size': min_clone_size,
                'ignore_zygosity': ignore_zygosity
            }

            method = assay.count
        else:
            attribute = st.selectbox('Attribute', [UMAP_LABEL, PCA_LABEL], key='Prepare Attribute')
            cluster_arg = st.slider(*DFT.CLUSTER_OPTIONS[method][:-1])

            description = f'{method} on {attribute} with {DFT.CLUSTER_OPTIONS[method][0].lower()} set to {cluster_arg}'
            cluster_kwargs = {
                'attribute': attribute,
                'method': method,
                DFT.CLUSTER_OPTIONS[method][-1]: cluster_arg
            }

            if assay.name == DNA_ASSAY:
                cluster_kwargs['similarity'] = st.slider('Similarity', 0.0, 1.0, 0.8)
                description = description + f" with {cluster_kwargs['similarity']} similarity"

            method = assay.cluster

        clicked = st.button('Cluster')

    return clicked, method, description, cluster_kwargs, info


def first_pass_cluster(available_assays):
    for prep_assay in available_assays:
        if prep_assay.metadata[DFT.INITIALIZE]:
            attribute, method, cluster_arg = DFT.CLUSTER_METHOD[prep_assay.name]
            description = f'{method} on {attribute} with {DFT.CLUSTER_OPTIONS[method][0]} set to {cluster_arg}'
            cluster_kwargs = {
                'attribute': attribute,
                'method': method,
                DFT.CLUSTER_OPTIONS[method][-1]: cluster_arg
            }
            if prep_assay.name == DNA_ASSAY:
                cluster_kwargs['similarity'] = 0.8
                description = description + f" with {cluster_kwargs['similarity']} similarity"
            cluster(prep_assay, prep_assay.cluster, description, **cluster_kwargs)


@st.cache(max_entries=1, hash_funcs=DFT.MOHASH, show_spinner=False)
def cluster(assay, method_func, description, **kwargs):
    similarity = None
    if 'similarity' in kwargs:
        similarity = kwargs['similarity']
        del kwargs['similarity']

    if DFT.CLUSTER_DESCRIPTION not in assay.metadata or assay.metadata[DFT.CLUSTER_DESCRIPTION] != description or not assay.metadata[DFT.CLUSTERED]:
        interface.status(f'Clustering {assay.name.replace("_", " ")}')
        method_func(**kwargs)
        if similarity is not None:
            assay.cluster_cleanup(AF_MISSING, similarity)

        assay.add_metadata(DFT.CLUSTER_DESCRIPTION, description)

    assay.add_metadata(DFT.CLUSTERED, True)
