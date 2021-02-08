import boto3
import numpy as np
import os

import streamlit as st
import defaults as DFT
import interface

import missionbio.mosaic.io as mio


def run():
    file, load_raw, apply_filter, kind, should_save, save_name, info = render()
    if kind == DFT.S3:
        file = download(file)

    if file is None:
        interface.error('Please use the options available in the sidebar to load a sample.<br>'
                        'New h5 files should be copied to the <i>/h5/downloads/</i> folder where the app is stored.')

    sample = load(file, load_raw, apply_filter)
    interface.info(f'Currently loaded {sample.name}', info)

    return sample, should_save, save_name


def render():
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

        if not os.path.exists(DFT.ROOT / 'h5'):
            os.mkdir(DFT.ROOT / 'h5')

        if not os.path.exists(DFT.ROOT / 'h5/downloads'):
            os.mkdir(DFT.ROOT / 'h5/downloads/')

        if not os.path.exists(DFT.ROOT / 'h5/analyzed'):
            os.mkdir(DFT.ROOT / 'h5/analyzed/')

        downloaded_files = np.array(os.listdir(DFT.ROOT / 'h5/downloads/'))
        analyzed_files = np.array(os.listdir(DFT.ROOT / 'h5/analyzed/'))
        filenames = list(analyzed_files[analyzed_files.argsort()]) + list(downloaded_files[downloaded_files.argsort()])
        filenames = [None] + [f for f in filenames if f[-3:] == '.h5']

        def shownames(name):
            nonlocal analyzed_files
            if name in analyzed_files:
                return '* ' + name
            else:
                return name

        kind = None
        selector = st.empty()
        file = selector.selectbox('Load existing file', filenames, format_func=shownames)

        if link != '':
            kind = DFT.S3
            file = link
        elif file is not None:
            if file in downloaded_files:
                file = DFT.ROOT / f'h5/downloads/{file}'
            else:
                file = DFT.ROOT / f'h5/analyzed/{file}'
            kind = DFT.LOCAL

        typed_name = st.text_input('Save, download or delete the given file', value='')

        def _get_file_from_name(typed_name):
            if typed_name[-3:] == '.h5':
                typed_name = typed_name[:-3]

            if typed_name + '.h5' in analyzed_files:
                typed_name = DFT.ROOT / f'h5/analyzed/{typed_name}.h5'
            elif typed_name + '.h5' in downloaded_files:
                typed_name = DFT.ROOT / f'h5/downloads/{typed_name}.h5'
            else:
                interface.error(f'Cannot find "{typed_name}" in the available files')

            return typed_name

        col1, col2, col3 = st.beta_columns([0.25, 0.4, 0.4])
        with col1:
            st.markdown('')
            should_save = st.button('Save')
        with col2:
            st.markdown('')
            if st.button('Download'):
                download_path = _get_file_from_name(typed_name)
                interface.download(download_path)
        with col3:
            st.markdown('')
            if st.button('Delete'):
                typed_name = _get_file_from_name(typed_name)
                if file is not None and typed_name == file:
                    interface.error('Cannot delete the file used in the current analysis.')
                os.remove(typed_name)
                interface.rerun()

    return file, load_raw, apply_filter, kind, should_save, typed_name, info


@st.cache(max_entries=1, show_spinner=False)
def download(link):
    interface.status('Downloading from s3.')

    s3 = boto3.client('s3')
    link = link.replace('s3://', '')
    link = link.split('/')
    bucket, file = link[0], '/'.join(link[1:])
    filename = file.split('/')[-1]
    filename = DFT.ROOT / f'h5/downloads/{filename}'
    filename = str(filename)
    try:
        s3.download_file(bucket, file, filename)
    except Exception as e:
        interface.status('Done.')
        interface.error(f'Could not find the given h5 file. {e}')

    return filename


@st.cache(max_entries=1, show_spinner=False, allow_output_mutation=True)
def load(path, load_raw, apply_filter):
    interface.status('Reading h5 file.')

    sample = mio.load(path, apply_filter=apply_filter, raw=load_raw)

    if sample.protein is not None:
        try:
            new_ids = np.array([ab.split(' ')[2] for ab in sample.protein.col_attrs['id']])
        except IndexError:
            new_ids = sample.protein.ids()

        sample.protein.add_col_attr('id', new_ids)
        if sample.protein_raw is not None:
            sample.protein_raw.add_col_attr('id', new_ids)

    init_defaults(sample)

    return sample


def init_defaults(sample):
    def add_arg(assay, key, val):
        if key not in assay.metadata:
            assay.add_metadata(key, val)

    add_arg(sample.dna, DFT.PREPROCESS_ARGS, [DFT.MIN_DP, DFT.MIN_GQ, DFT.MIN_VAF, DFT.MIN_STD])
    add_arg(sample.dna, DFT.DROP_IDS, [])
    add_arg(sample.dna, DFT.KEEP_IDS, [])
    add_arg(sample.dna, DFT.ALL_IDS, sample.dna.ids())
    add_arg(sample.dna, DFT.INITIALIZE, True)
    add_arg(sample.dna, DFT.PREPPED, True)
    add_arg(sample.dna, DFT.CLUSTERED, True)

    add_arg(sample.protein, DFT.DROP_IDS, [])
    add_arg(sample.protein, DFT.ALL_IDS, sample.protein.ids())
    add_arg(sample.protein, DFT.INITIALIZE, True)
    add_arg(sample.protein, DFT.PREPPED, True)
    add_arg(sample.protein, DFT.CLUSTERED, True)
