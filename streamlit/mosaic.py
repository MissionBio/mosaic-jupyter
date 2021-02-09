import streamlit as st

import interface
import defaults as DFT
from tasks import (
    load,
    preprocess,
    prepare,
    cluster,
    customize,
    save,
    visual
)

st.set_page_config(page_title='Mosaic', layout='wide')
interface.init()
interface.subheader('---')
interface.status('Processing')

sample, should_save, save_name = load.run()

current_assay, available_assays = preprocess.run(sample)

prepare.run(current_assay, available_assays)
cluster.run(current_assay, available_assays)
customize.run(current_assay)

if should_save:
    save.run(sample, save_name)

visual.run(sample, current_assay)

for a in available_assays:
    a.add_metadata(DFT.INITIALIZE, False)
