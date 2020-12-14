import os
import boto3
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import missionbio.mosaic.io as mio
from missionbio.mosaic.dna import Dna as modna
from missionbio.mosaic.protein import Protein as moprotein
from missionbio.mosaic.sample import Sample as mosample
from missionbio.mosaic.constants import (
    AF_MISSING,
    NORMALIZED_READS,
    PCA_LABEL,
    SCALED_LABEL,
    UMAP_LABEL,
    READS,
    COLORS
)

import defaults as DFT
import interface
import session


MOHASH = {moprotein: lambda a: a.name + a.title + str(a.shape),
          modna: lambda a: a.name + a.title + str(a.shape),
          mosample: lambda a: a.name + str(a.dna.shape),
          type(lambda _: None): lambda _: None}

MOHASH_BOOL = {**MOHASH, int: lambda _: None, list: lambda _: None, type(''): lambda _: None}

state = session.get(prev_file=None,
                    dna=None,
                    protein=None,
                    assays=None,
                    init_dna=True,
                    init_protein=True,
                    dna_preprocess=None,
                    protein_preprocess=None,
                    prepped=True,
                    clustered=True,
                    analysis_type=DFT.PRIMARY,
                    prev_plot=DFT.COLORS)


@st.cache(max_entries=1, show_spinner=False)
def download(link):
    interface.status('Downloading from s3.')

    s3 = boto3.client('s3')
    link = link.replace('s3://', '')
    link = link.split('/')
    bucket, file = link[0], '/'.join(link[1:])
    filename = file.split('/')[-1]
    filename = f'./h5/downloads/{filename}'
    try:
        s3.download_file(bucket, file, filename)
    except Exception:
        interface.status('Done.')
        interface.error('Could not find the given h5 file.')

    return filename


@st.cache(max_entries=1, show_spinner=False)
def load(path, load_raw, apply_filter):
    interface.status('Reading h5 file.')

    sample = mio.load(path, apply_filter=apply_filter, raw=load_raw)
    state.init_dna = True

    if sample.protein is not None:
        try:
            new_ids = np.array([ab.split(' ')[2] for ab in sample.protein.col_attrs['id']])
        except IndexError:
            new_ids = sample.protein.ids()

        sample.protein.add_col_attr('id', new_ids)
        if sample.protein_raw is not None:
            sample.protein_raw.add_col_attr('id', new_ids)

        state.init_protein = True

    return sample


def save(sample, dna, protein, should_save, name):
    interface.status('Saving h5 file.')
    if name == '':
        interface.error('Please provide a name to save by.')
    elif name[-3:] == '.h5':
        name = name[:-3]

    prot = sample.protein
    if prot is not None:
        prot = protein

    samp = mosample(protein=prot, dna=dna, cnv=sample.cnv)
    try:
        os.remove(f'./h5/analyzed/{name}.h5')
    except FileNotFoundError:
        pass
    mio.save(samp, f'./h5/analyzed/{name}.h5')

    interface.status('Done.')
    interface.rerun()


def init_states(sample):
    if state.init_dna:
        state.dna = sample.dna[:, :]

        if DFT.FILTERS in sample.dna.metadata:
            state.dna_preprocess = [[], *sample.dna.metadata[DFT.FILTERS]]
        else:
            state.dna_preprocess = [[], DFT.MIN_DP, DFT.MIN_GQ, DFT.MIN_VAF, DFT.MIN_STD]

        if DFT.UMAP_LABEL in sample.dna.row_attrs:
            state.init_dna = False

    if state.init_protein:
        state.protein = sample.protein[:, :]
        state.protein_preprocess = [[], DFT.NSP in sample.protein.layers]
        if DFT.UMAP_LABEL in sample.protein.row_attrs:
            state.init_protein = False


@st.cache(max_entries=1, hash_funcs=MOHASH_BOOL, show_spinner=False, allow_output_mutation=True)
def preprocess_dna(sample, clicked, drop_vars, dp, gq, af, std):
    if state.init_dna or (state.dna_preprocess != [drop_vars, dp, gq, af, std] and clicked):
        interface.status('Processing DNA assay.')

        dna = sample.dna.drop(drop_vars) if len(drop_vars) > 0 else sample.dna[:, :]
        dna_vars = dna.filter_variants(min_dp=dp, min_gq=gq, min_vaf=af, min_std=std)

        if len(dna_vars) == 0:
            interface.status('Done.')
            interface.error('No variants found. Adjust the filters and process again.')

        dna = dna[:, dna_vars]

        if state.dna is not None and state.dna.shape[0] == dna.shape[0]:
            for key, val in state.dna.row_attrs.items():
                dna.add_row_attr(key, val)

            for key, val in state.dna.metadata.items():
                dna.add_metadata(key, val)

        state.dna_preprocess = [drop_vars, dp, gq, af, std]
        state.dna = dna

        if clicked:
            state.prepped = False
            state.clustered = False

    return state.dna


@st.cache(max_entries=1, hash_funcs=MOHASH_BOOL, show_spinner=False, allow_output_mutation=True)
def preprocess_protein(sample, clicked, drop_abs):
    if state.init_protein or (state.protein_preprocess != [drop_abs, True] and clicked):
        interface.status('Processing protein assay.')

        protein = sample.protein.drop(drop_abs) if len(drop_abs) > 0 else sample.protein[:, :]

        for norm in [DFT.CLR, DFT.ASINH, DFT.NSP]:
            protein.normalize_reads(norm)
            protein.add_layer(norm, protein.layers[NORMALIZED_READS])

        if state.protein is not None and state.protein.shape[0] == protein.shape[0]:
            for key, val in state.protein.row_attrs.items():
                protein.add_row_attr(key, val)

            for key, val in state.protein.metadata.items():
                protein.add_metadata(key, val)

        state.protein_preprocess = [drop_abs, True]
        state.protein = protein

        if clicked:
            state.prepped = False
            state.clustered = False

    return state.protein


@st.cache(max_entries=1, hash_funcs=MOHASH, show_spinner=False)
def prepare(assay, scale_attribute, pca_attribute, umap_attribute, pca_comps):
    interface.status(f'Preparing {assay.name.replace("_", " ")} data.')

    attr = scale_attribute
    if SCALED_LABEL not in assay.layers or assay.metadata[DFT.SCALE_ATTR] != attr or not state.prepped:
        assay.scale_data(scale_attribute)
        assay.add_metadata(DFT.SCALE_ATTR, attr)

    attr = f'{pca_attribute}'
    if pca_attribute == SCALED_LABEL:
        attr = f'scaled {scale_attribute}'

    if PCA_LABEL not in assay.row_attrs or assay.metadata[DFT.PCA_ATTR] != attr or not state.prepped:
        assay.run_pca(pca_attribute, components=pca_comps)
        assay.add_metadata(DFT.PCA_ATTR, attr)

    attr = f'{umap_attribute}'
    if umap_attribute == SCALED_LABEL:
        attr = f'scaled {scale_attribute}'
    if umap_attribute == PCA_LABEL:
        if pca_attribute == SCALED_LABEL:
            attr = f'PCA of scaled {scale_attribute}'
        else:
            attr = f'PCA of {pca_attribute}'

    if UMAP_LABEL not in assay.row_attrs or assay.metadata[DFT.UMAP_ATTR] != attr or not state.prepped:
        assay.run_umap(attribute=umap_attribute, random_state=42)
        assay.add_metadata(DFT.UMAP_ATTR, attr)

    state.prepped = True


def preliminary_prepare(dna, protein):
    for prep_assay in [dna, protein]:
        if DFT.UMAP_LABEL not in prep_assay.row_attrs:
            if prep_assay.name == dna.name:
                prepare(prep_assay, AF_MISSING, SCALED_LABEL, AF_MISSING, min(len(prep_assay.ids()), 8))
            elif prep_assay.name == protein.name:
                prepare(prep_assay, DFT.LAYERS[prep_assay.name][0], SCALED_LABEL, PCA_LABEL, min(len(prep_assay.ids()), 8))


@st.cache(max_entries=1, hash_funcs=MOHASH, show_spinner=False)
def cluster(assay, method_func, description, random_state=42, **kwargs):
    similarity = None
    if 'similarity' in kwargs:
        similarity = kwargs['similarity']
        del kwargs['similarity']

    if DFT.CLUSTER_DESCRIPTION not in assay.metadata or assay.metadata[DFT.CLUSTER_DESCRIPTION] != description or not state.clustered:
        interface.status(f'Clustering {assay.name.replace("_", " ")}')
        method_func(**kwargs)
        if similarity is not None:
            assay.cluster_cleanup(AF_MISSING, similarity)

        assay.add_metadata(DFT.CLUSTER_DESCRIPTION, description)

    state.clustered = True


def preliminary_cluster(dna, protein):
    for prep_assay in [dna, protein]:
        if DFT.CLUSTER_DESCRIPTION not in prep_assay.metadata:
            attribute, method, cluster_arg = DFT.CLUSTER_METHOD[prep_assay.name]
            description = f'{method} on {attribute} with {DFT.CLUSTER_OPTIONS[method][0]} set to {cluster_arg}'
            cluster_kwargs = {
                'attribute': attribute,
                'method': method,
                DFT.CLUSTER_OPTIONS[method][-1]: cluster_arg
            }
            if prep_assay.name == dna.name:
                cluster_kwargs['similarity'] = 0.8
                description = description + f" with {cluster_kwargs['similarity']} similarity"
            cluster(prep_assay, prep_assay.cluster, description, random_state=np.random.random(), **cluster_kwargs)


def metrics(sample):
    prot_metrics = sample.protein.metadata
    dna_metrics = sample.dna.metadata

    prot_total_reads = prot_metrics['n_read_pairs']

    st.text(f'name {sample.name}')
    st.text(f'total-reads: {prot_total_reads}')
    st.text(f"reads-after-cutadapt: {prot_metrics['n_read_pairs_trimmed']}")
    st.text(f"%reads-trimmed: {prot_metrics['n_read_pairs_trimmed'] / prot_total_reads}")
    st.text(f"reads-barcode-correction: {prot_metrics['n_read_pairs_valid_cell_barcodes']}")
    st.text(f"%reads-barcode: {prot_metrics[f'n_read_pairs_valid_cell_barcodes'] / prot_total_reads}")
    st.text(f"reads-ab-correction: {prot_metrics['n_read_pairs_valid_ab_barcodes']}")
    st.text(f"%reads-ab: {prot_metrics['n_read_pairs_valid_ab_barcodes'] / prot_total_reads}")

    if sample.protein_raw is None:
        prot_read_count = 0
        read_per_bar = 0
        num_candidate_prot_bars = 0
        candidate_reads = 0
        percentage_aggregates = 0
    else:
        prot_read_count = sample.protein_raw.layers['read_counts']
        read_per_bar = prot_read_count.sum(axis=1)
        num_candidate_prot_bars = (read_per_bar > 10).sum()
        candidate_reads = read_per_bar[read_per_bar > 10].sum()
        percentage_aggregates = get_aggregates(sample.protein_raw)['per_total'].sum() / 100

    cells_after_merge = len(sample.protein.row_attrs['barcode'])

    if sample.cnv_raw is None:
        dna_reads_to_cells = 0
        uniformity = 0
    else:
        dna_reads_to_cells = sample.cnv_raw.layers['read_counts'][sample.cnv_raw.row_attrs['cellfinder_pass'].astype(bool), :].sum() / dna_metrics['n_read_pairs']
        cell_df = pd.DataFrame(sample.cnv_raw.layers['read_counts'][sample.cnv_raw.row_attrs['cellfinder_pass'].astype(bool), :])
        uniformity = cell_df.mean().between(left=(0.2 * cell_df.mean().mean()), right=(5 * cell_df.mean().mean()), inclusive=True).sum() / cell_df.shape[1]

    dna_reads_to_merged_cells = sample.cnv.layers['read_counts'].sum() / dna_metrics['n_read_pairs']
    prot_reads_to_cells = sample.protein.layers['read_counts'].sum() / candidate_reads
    prot_reads_to_cells_total_reads = sample.protein.layers['read_counts'].sum() / prot_total_reads
    prot_bar_per_cell = cells_after_merge / len(sample.dna.row_attrs['barcode'])

    prot_cell_reads = pd.DataFrame(sample.protein.layers['read_counts'], columns=sample.protein.ids())
    avg_reads_per_ab_cell = prot_cell_reads.mean().mean()
    med_reads_per_ab_cell = prot_cell_reads[prot_cell_reads > 1].median().median()
    total_ab = (prot_cell_reads.mean(axis=0) > 0).sum()

    iso_reads = sample.protein.layers['read_counts'][:, np.where(np.isin(sample.protein.ids(), ['IgG1', 'IgG2a', 'IgG2b']))[0]]
    noise_reads = iso_reads.sum(axis=1)
    signal_reads = sample.protein.layers['read_counts'].sum(axis=1)
    nsb = 100 * noise_reads / signal_reads
    nsb_ratio_cell = (nsb > 0.1).sum() / len(nsb)
    nsb_ratio_all = noise_reads.sum() / signal_reads.sum()

    st.text(f"candidate-barcodes: {num_candidate_prot_bars}")
    st.text(f"candidate-reads: {candidate_reads}")
    st.text(f"%candidate-reads {candidate_reads / prot_total_reads:0.3f}")
    st.text(f"cells-after-merge: {cells_after_merge}")

    st.text(f"%reads-to-cells-of-valid-barcodes {prot_reads_to_cells:0.3f}")
    st.text(f"%reads-to-cells-of-total: {prot_reads_to_cells_total_reads:0.3f}")
    st.text(f"%valid-ab-barcodes-that-are-cells: {prot_bar_per_cell:0.3f}")
    st.text(f"%reads-by-aggregates: {percentage_aggregates:.3f}")
    st.text(f"avg-ab-reads: {avg_reads_per_ab_cell}")
    st.text(f"median-ab-reads: {med_reads_per_ab_cell}")
    st.text(f"NSB: {nsb_ratio_all}")
    st.text(f"NSB-proportion: {nsb_ratio_cell}")
    st.text(f"total-ab: {total_ab}")

    st.text(f"%dna-reads-to-cells: {dna_reads_to_cells:0.3f}")
    st.text(f"%dna-reads-to-merged-cells: {dna_reads_to_merged_cells:0.3f}")
    st.text(f"dna-reads-per-cell-per-amp: {sample.cnv.layers['read_counts'].mean().mean():.2f}")
    st.text(f"ado-rate: {sample.dna.metadata['ado_rate']}")
    st.text(f"avg-dna-panel-uniformity {sample.dna.metadata['avg_panel_uniformity']}")
    st.text(f"0.2x-dna-panel-uniformity {uniformity}")


def get_aggregates(assay_object, threshold=10):
    # Takes protein assay with all barcodes as argument
    # Percentage of the total reads to that AB in one cell required to be considered an aggreagte

    # Reads taken by the barcode-antibody as a percentage of total reads to that antibody
    read_counts = assay_object.layers['read_counts']
    percent_reads = 100 * read_counts / read_counts.sum(axis=0)
    df_per_ab = pd.DataFrame(percent_reads,
                             index=assay_object.row_attrs['barcode'],
                             columns=assay_object.col_attrs['id'])

    # Reads taken by the barcode-antibody as a percentage of total overall reads
    read_counts = assay_object.layers['read_counts']
    percent_reads = 100 * read_counts / read_counts.sum().sum()
    df_per_total = pd.DataFrame(percent_reads,
                                index=assay_object.row_attrs['barcode'],
                                columns=assay_object.col_attrs['id'])

    # Apply the threshold on the per antibody percentage
    agg_bars = df_per_ab.loc[(df_per_ab > threshold).any(axis=1), :]
    agg_percent_ab = agg_bars[agg_bars > threshold]
    agg_percent_ab = agg_percent_ab.reset_index().melt(id_vars='index',
                                                       var_name='antibody',
                                                       value_name='per_ab').dropna()

    # Apply the threshold on the per antibody percentage but get values from the DataFrame with the total percentages
    agg_bars_per_total = df_per_total.loc[(df_per_ab > threshold).any(axis=1), :]
    agg_percent_total = agg_bars_per_total[agg_bars > threshold]
    agg_percent_total = agg_percent_total.reset_index().melt(id_vars='index',
                                                             var_name='antibody',
                                                             value_name='per_total').dropna()

    # Merge the two calculated values and return as a DataFrame
    agg_percent = pd.merge(agg_percent_ab, agg_percent_total)

    return agg_percent


def visual(sample, assay, dna, protein, kind, plot_columns, kwargs):
    if kind == DFT.SIGNATURES:
        with plot_columns:
            med, std, pval, _ = assay.feature_signature(layer=kwargs['layer'])
            if kwargs['attribute'] == 'Median':
                df = med
            elif kwargs['attribute'] == 'Standard deviation':
                df = std
            elif kwargs['attribute'] == 'p-value':
                df = pval.applymap("{0:.2E}".format)
            st.write(kwargs['attribute'])
            st.dataframe(df)
    elif kind in [DFT.HEATMAP, DFT.SCATTERPLOT, DFT.FEATURE_SCATTER, DFT.VIOLINPLOT, DFT.RIDGEPLOT]:
        with plot_columns:
            plot_funcs = {
                DFT.HEATMAP: assay.heatmap,
                DFT.SCATTERPLOT: assay.scatterplot,
                DFT.FEATURE_SCATTER: assay.feature_scatter,
                DFT.VIOLINPLOT: assay.violinplot,
                DFT.RIDGEPLOT: assay.ridgeplot
            }

            labelby = None
            if 'splitby' in kwargs:
                labelby = 'splitby'
            elif 'colorby' in kwargs:
                labelby = 'colorby'

            org_lab = assay.get_labels().copy()
            org_pal = assay.get_palette()
            new_lab = org_lab
            new_pal = org_pal

            if kwargs[labelby] == DFT.PROTEIN_LABEL:
                new_pal = protein.get_palette()
                new_lab = protein.get_labels().copy()
                kwargs[labelby] = 'label'

            if kwargs[labelby] == DFT.DNA_LABEL:
                new_pal = dna.get_palette()
                new_lab = dna.get_labels().copy()
                kwargs[labelby] = 'label'

            assay.set_labels(new_lab)
            assay.set_palette(new_pal)

            if 'cluster' in kwargs:
                bars_ordered = assay.clustered_barcodes(orderby=kwargs['attribute'])
                if not kwargs['cluster']:
                    labels = assay.get_labels()[[np.where(assay.barcodes() == b)[0][0] for b in bars_ordered]]
                    bars_ordered = []
                    for lab in pd.unique(labels):
                        bars_ordered.extend(assay.barcodes(lab))
                    bars_ordered = np.array(bars_ordered)

                kwargs['bars_order'] = bars_ordered
                del kwargs['cluster']

            fig = plot_funcs[kind](**kwargs)
            st.plotly_chart(fig)

            assay.set_labels(org_lab)
            assay.set_palette(org_pal)

    elif kind == DFT.DNA_PROTEIN_PLOT:
        with plot_columns:
            samp = mosample(protein=protein[:, kwargs['protein_features']], dna=dna[:, kwargs['dna_features']])
            samp.clone_vs_analyte(kwargs['analyte'])
            st.pyplot(plt.gcf())
    elif kind == DFT.DNA_PROTEIN_HEATMAP:
        with plot_columns:
            samp = mosample(protein=protein[:, kwargs['protein_features']], dna=dna[:, kwargs['dna_features']])
            fig = samp.heatmap(clusterby=kwargs['clusterby'], sortby=kwargs['sortby'], drop='cnv', flatten=False)
            st.plotly_chart(fig)

    elif kind == DFT.METRICS:
        with plot_columns[0]:
            st.header('')
            metrics(sample)
    elif kind == DFT.READ_DEPTH:
        with plot_columns[0]:
            if assay.name == protein.name:
                total_reads = assay.layers[READS].sum(axis=1)
                layer = assay.layers[kwargs['layer']]

                x = np.log10(total_reads)
                st.write(assay.title)

                feats = kwargs['features']
                if len(feats) == 0:
                    feats = assay.ids()

                for feat in feats:
                    y = layer[:, np.where(assay.ids() == feat)[0]].flatten()
                    data = np.array([x, y]).T
                    fig = assay.scatterplot(attribute=data, colorby=kwargs['colorby'])
                    fig.layout.title = ''
                    fig.layout.xaxis.title.text = 'log10(Total reads)'
                    fig.layout.yaxis.title.text = f'Expression of {feat}'
                    fig.layout.width = 400
                    fig.layout.height = 400
                    st.plotly_chart(fig)
    elif kind == DFT.ASSAY_SCATTER:
        if kwargs['draw']:
            with plot_columns[0]:
                # sample.cnv_raw.normalize_barcodes()
                # sample.protein_raw.normalize_barcodes()
                sample.assay_scatter()
                st.pyplot(plt.gcf())
    elif kind == DFT.COLORS:
        colors = COLORS.copy()
        del colors[20]

        for i in range(len(plot_columns)):
            plot_columns[i].header('')

        for i in range(len(colors)):
            plot_columns[i % len(plot_columns)].color_picker(colors[i], colors[i], key=f'constant-colors-{colors[i]}-{i}')
    elif kind == DFT.DOWNLOAD:
        with plot_columns[0]:
            if kwargs['download']:
                if kwargs['item'] == DFT.ANNOTATION:
                    data = np.array([dna.barcodes(), dna.get_labels(), protein.get_labels()]).T
                    df = pd.DataFrame(data, columns=['barcode', 'dna', 'protein'])
                    name = f'./{sample.name}.annotation.csv'
                    df.to_csv(name, index=None)
                    interface.download(name)
                    os.remove(name)


@st.cache(max_entries=1, hash_funcs=MOHASH, show_spinner=False)
def customize(assay, lab_map, pal):
    old_labs = set(assay.get_labels())
    old_pal = assay.get_palette()
    assay.rename_labels(lab_map)
    assay.set_palette(pal)
    new_labs = set(assay.get_labels())

    if new_labs != old_labs or old_pal != pal:
        interface.rerun()
