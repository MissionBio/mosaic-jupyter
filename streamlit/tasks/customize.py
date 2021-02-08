import numpy as np
import streamlit as st

import interface
import defaults as DFT


def run(assay):
    lab_map, pal = render(assay)
    customize(assay, lab_map, pal)


def render(assay):
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


@st.cache(max_entries=1, hash_funcs=DFT.MOHASH, show_spinner=False)
def customize(assay, lab_map, pal):
    old_labs = set(assay.get_labels())
    old_pal = assay.get_palette()
    assay.rename_labels(lab_map)
    assay.set_palette(pal)
    new_labs = set(assay.get_labels())

    if new_labs != old_labs or old_pal != pal:
        interface.rerun()
