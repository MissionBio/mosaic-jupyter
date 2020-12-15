import streamlit as st
import os
import numpy as np

import defaults as DFT
import interface
from compute import state

from missionbio.h5.constants import DNA_ASSAY
from missionbio.mosaic.constants import (
    PCA_LABEL,
    UMAP_LABEL
)


def files():
    with st.sidebar.beta_expander('Files', expanded=True):
        interface.info('Load or download a file from s3')
        info = st.empty()

        col1, col2 = st.beta_columns([0.3, 1])
        with col1:
            st.markdown(f"<sup><p style='margin-bottom:22px'></p></sup>", unsafe_allow_html=True)
            load_raw = st.checkbox('Raw')
        with col2:
            st.markdown(f"<sup><p style='margin-bottom:22px'></p></sup>", unsafe_allow_html=True)
            apply_filter = st.checkbox('Filter', True)

        link = st.text_input('Load from s3', value='')

        if not os.path.exists('./h5'):
            os.mkdir('./h5/')

        if not os.path.exists('./h5/downloads'):
            os.mkdir('./h5/downloads/')

        if not os.path.exists('./h5/analyzed'):
            os.mkdir('./h5/analyzed/')

        downloaded_files = np.array(os.listdir('./h5/downloads/'))
        analyzed_files = np.array(os.listdir('./h5/analyzed/'))
        filenames = list(analyzed_files[analyzed_files.argsort()]) + list(downloaded_files[downloaded_files.argsort()])
        filenames = [None] + [f for f in filenames if f[-3:] == '.h5']

        def shownames(name):
            nonlocal analyzed_files
            if name in analyzed_files:
                return '* ' + name
            else:
                return name

        kind = None
        index = filenames.index(state.prev_file) if link == '' else 0
        selector = st.empty()
        file = selector.selectbox('Load existing file', filenames, format_func=shownames, index=index)

        if state.prev_file != file and link == '':
            index = filenames.index(file)
            file = selector.selectbox('Load existing file', filenames, format_func=shownames, index=index)

        state.prev_file = file

        if link != '':
            kind = DFT.S3
            file = link
        elif file is not None:
            if file in downloaded_files:
                file = f'./h5/downloads/{file}'
            else:
                file = f'./h5/analyzed/{file}'
            kind = DFT.LOCAL

        typed_name = st.text_input('Save, download or delete the given file', value='')

        def _get_file_from_name(typed_name):
            if typed_name[-3:] == '.h5':
                typed_name = typed_name[:-3]

            if typed_name + '.h5' in analyzed_files:
                typed_name = './h5/analyzed/' + typed_name + '.h5'
            elif typed_name + '.h5' in downloaded_files:
                typed_name = './h5/downloads/' + typed_name + '.h5'
            else:
                interface.error(f'Cannot find "{typed_name}" in the available files')

            return typed_name

        col1, col2, col3 = st.beta_columns([0.25, 0.4, 0.4])
        with col1:
            st.markdown('')
            should_save = st.button('Save')
        with col3:
            st.markdown('')
            if st.button('Delete'):
                typed_name = _get_file_from_name(typed_name)
                if file is not None and typed_name == file:
                    interface.error('Cannot delete the file used in the current analysis.')
                os.remove(typed_name)
                interface.rerun()
        with col2:
            st.markdown('')
            if st.button('Download'):
                download_path = _get_file_from_name(typed_name)
                interface.download(download_path)

    return file, load_raw, apply_filter, kind, should_save, typed_name, info


def preprocess(sample):
    drop_vars, keep_vars, dp, gq, af, std = state.dna_preprocess
    drop_abs = state.protein_preprocess[0]

    with st.sidebar.beta_expander('Prepocessing'):
        info = st.empty()
        assay_type = st.selectbox('Assay', ['DNA', 'Protein'])

        if assay_type == 'DNA':
            dp = st.slider('Minimum read depth (DP)', min_value=0, max_value=100, value=int(dp))
            gq = st.slider('Minimum genotype quality (GQ)', min_value=0, max_value=100, value=int(gq))
            af = st.slider('Minimum allele frequency (VAF)', min_value=0, max_value=100, value=int(af))
            std = st.slider('Minimum standard deviation of AF', min_value=0, max_value=100, value=int(std))
            ids = sample.dna.ids()
            ids = list(ids[ids.argsort()])
            drop_vars = st.multiselect('Variants to discard', ids, default=drop_vars)
            keep_vars = st.multiselect('Variants to keep', ids, default=keep_vars)

            if len(keep_vars) != 0 and len(drop_vars) != 0:
                interface.error('Cannot keep and drop variants both. Choose only one of the options')

        elif assay_type == 'Protein':
            ids = sample.protein.ids()
            ids = list(ids[ids.argsort()])
            drop_abs = st.multiselect('Antibodies to discard', ids, default=drop_abs)

        interface.info(f'{assay_type} currently loaded', info)

        clicked = st.button('Process')

    return assay_type, clicked, [drop_vars, keep_vars, dp, gq, af, std], [drop_abs]


def prepare(assay):
    title = 'Data preparation'
    if not state.prepped:
        title += ' *'

    with st.sidebar.beta_expander(title):
        info = st.empty()

        scale_attribute = st.selectbox('Scale - Attribute', DFT.LAYERS[assay.name])
        pca_attribute = st.selectbox('PCA - Attribute', ['scaled_counts'] + DFT.LAYERS[assay.name])

        pca_comps = 6
        pca_comps = st.slider('Number of components', 3, 20, 6)

        umap_attribute = st.selectbox('UMAP - Attribute', ['pca', 'scaled_counts'] + DFT.LAYERS[assay.name])

        clicked = st.button('Prepare')

    return clicked, scale_attribute, pca_attribute, umap_attribute, pca_comps, info


def cluster(assay, dna):
    title = 'Clustering'
    if not state.clustered:
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
            features = st.multiselect('Variants', list(dna.ids()))

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


def customize(assay):
    with st.sidebar.beta_expander('Customizations'):
        interface.info('Rename the labels.<br>Merge by giving the same name.')

        lab_map = {}
        pal = assay.get_palette()
        for lab in np.unique(assay.get_labels()):
            col1, col2 = st.beta_columns([1, 0.15])
            with col1:
                new_name = st.text_input(f'Give a new name to {lab}', lab)
            with col2:
                st.markdown(f"<p style='margin-bottom:34px'></p>", unsafe_allow_html=True)
                pal[lab] = st.color_picker('', pal[lab], key=f'colorpicker-{lab}')

            if new_name != lab:
                lab_map[lab] = new_name
                pal[new_name] = pal[lab]
                del pal[lab]

    return lab_map, pal


def visual(assay, dna, protein, sample):
    interface.status('Creating visuals.')

    col1, col2 = st.beta_columns([0.07, 1])
    with col1:
        if st.button(DFT.PRIMARY):
            state.analysis_type = DFT.PRIMARY
    with col2:
        if st.button(DFT.ADVANCED):
            state.analysis_type = DFT.ADVANCED

    if state.analysis_type == DFT.PRIMARY:
        types_container, args_conatiner, plot_columns = st.beta_columns([0.5, 0.75, 2])
        with types_container:
            kind = st.radio('', [DFT.SIGNATURES,
                                 DFT.HEATMAP,
                                 DFT.SCATTERPLOT,
                                 DFT.FEATURE_SCATTER,
                                 DFT.VIOLINPLOT,
                                 DFT.RIDGEPLOT,
                                 DFT.DNA_PROTEIN_PLOT,
                                 DFT.DNA_PROTEIN_HEATMAP], index=1)
    elif state.analysis_type == DFT.ADVANCED:
        columns = st.beta_columns(DFT.LAYOUT[state.prev_plot])
        types_container = columns[0]
        args_conatiner = columns[1]
        plot_columns = columns[2:]
        with types_container:
            plot_types = [DFT.COLORS,
                          DFT.METRICS,
                          DFT.READ_DEPTH,
                          DFT.ASSAY_SCATTER,
                          DFT.DOWNLOAD]
            kind = st.radio('', plot_types, index=plot_types.index(state.prev_plot))

            if state.prev_plot != kind:
                state.prev_plot = kind
                interface.rerun()

    with args_conatiner:
        kwargs = {}
        analyte_map = {'protein': 'Protein', 'dna': 'DNA'}

        if kind == DFT.SIGNATURES:
            kwargs['layer'] = st.selectbox('Layer', DFT.LAYERS[assay.name])
            kwargs['attribute'] = st.selectbox('Signature', ['Median', 'Standard deviation', 'p-value'])
        elif kind == DFT.HEATMAP:
            kwargs['attribute'] = st.selectbox('Attribute', DFT.LAYERS[assay.name], key='Visualization Attribute')
            kwargs['splitby'] = st.selectbox('Split by', DFT.SPLITBY[assay.name])
            kwargs['cluster'] = st.checkbox('Hierarchical cluster', True)
            kwargs['convolve'] = st.slider('Smoothing', 0, 100)
        elif kind == DFT.SCATTERPLOT:
            kwargs['attribute'] = st.selectbox('Attribute', DFT.ATTRS_2D)
            kwargs['colorby'] = st.selectbox('Color by', DFT.COLORBY[assay.name])
            if kwargs['colorby'] not in DFT.SPLITBY[assay.name] + ['density']:
                features = st.multiselect('Features', list(assay.ids()), list(assay.ids())[:min(len(assay.ids()), 4)])
                if len(features) != 0:
                    kwargs['features'] = features
        elif kind == DFT.FEATURE_SCATTER:
            kwargs['layer'] = st.selectbox('Layer', DFT.LAYERS[assay.name])
            feature1 = st.selectbox('Feature 1', list(assay.ids()), index=0)
            feature2 = st.selectbox('Feature 1', list(assay.ids()), index=2)
            kwargs['ids'] = [feature1, feature2]
            kwargs['colorby'] = st.selectbox('Color by', DFT.COLORBY[assay.name])
        elif kind == DFT.VIOLINPLOT:
            kwargs['attribute'] = st.selectbox('Attribute', DFT.LAYERS[assay.name])
            kwargs['splitby'] = st.selectbox('Split by', DFT.SPLITBY[assay.name])
            features = st.multiselect('Features', list(assay.ids()), list(assay.ids())[:min(len(assay.ids()), 4)])
            if len(features) != 0:
                kwargs['features'] = features
        elif kind == DFT.RIDGEPLOT:
            kwargs['attribute'] = st.selectbox('Attribute', DFT.LAYERS[assay.name])
            kwargs['splitby'] = st.selectbox('Split by', DFT.SPLITBY[assay.name])
            features = st.multiselect('Features', list(assay.ids()), list(assay.ids())[:min(len(assay.ids()), 4)])
            if len(features) != 0:
                kwargs['features'] = features
        elif kind == DFT.DNA_PROTEIN_PLOT:
            kwargs['analyte'] = st.selectbox('Analyte', ['protein'], format_func=lambda a: analyte_map[a])
            kwargs['dna_features'] = st.multiselect('DNA features', list(dna.ids()), dna.ids()[:4])
            kwargs['protein_features'] = st.multiselect('Protein features', list(protein.ids()), protein.ids()[:4])
        elif kind == DFT.DNA_PROTEIN_HEATMAP:
            kwargs['clusterby'] = st.selectbox('Cluster by', ['dna', 'protein'], format_func=lambda a: analyte_map[a])
            kwargs['sortby'] = st.selectbox('Sort by', ['dna', 'protein'], format_func=lambda a: analyte_map[a])
            kwargs['dna_features'] = st.multiselect('DNA features', list(dna.ids()), dna.ids())
            kwargs['protein_features'] = st.multiselect('Protein features', list(protein.ids()), protein.ids())

        elif kind == DFT.METRICS:
            st.header('')
            interface.info('<b>Some values might be missing in case the raw<br> files are not loaded.</b> These metrics can be<br> pasted into the metrics sheet as is.')
        elif kind == DFT.READ_DEPTH:
            if assay.name == protein.name:
                kwargs['layer'] = st.selectbox('Layer', DFT.LAYERS[assay.name])
                kwargs['colorby'] = st.selectbox('Color by', DFT.COLORBY[assay.name])
                kwargs['features'] = st.multiselect('Features', list(assay.ids()), list(assay.ids())[:min(len(assay.ids()), 4)])
            else:
                st.header('')
                interface.info('<b>Only applicable for the protein assay</b>')
        elif kind == DFT.ASSAY_SCATTER:
            kwargs['draw'] = sample.protein_raw is not None
            if not kwargs['draw']:
                interface.info('<b>Raw files needed for this plot.</b>')
        elif kind == DFT.DOWNLOAD:
            kwargs['item'] = st.selectbox('Object to Download', DFT.DOWNLOAD_ITEMS)
            kwargs['download'] = st.button('Download', key='download_button')

    return plot_columns, kind, kwargs
