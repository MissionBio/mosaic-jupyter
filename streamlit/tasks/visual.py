import os
import pandas as pd
import numpy as np
import streamlit as st

import requests

import matplotlib.pyplot as plt
from plotly.subplots import make_subplots

import interface
import defaults as DFT

from missionbio.h5.constants import PROTEIN_ASSAY
from missionbio.mosaic.constants import READS, COLORS
from missionbio.mosaic.sample import Sample as mosample


def run(sample, assay):
    plot_columns, kind, visualization_kwargs = render(sample, assay)
    visual(sample, assay, kind, plot_columns, visualization_kwargs)

    interface.status('Done.')


def render(sample, assay):
    interface.status('Creating visuals.')

    category, kind = assay.metadata[DFT.VISUAL_TYPE]
    options = DFT.VISUALS[category][1]
    column_sizes = DFT.VISUALS[category][0]
    columns = st.beta_columns(column_sizes)
    with columns[0]:
        new_category = st.selectbox("", list(DFT.VISUALS.keys()))
        if new_category != category:
            assay.add_metadata(DFT.VISUAL_TYPE, [new_category, DFT.VISUALS[new_category][1][0]])
            interface.rerun()

    for i in range(len(options)):
        with columns[i + 1]:
            st.markdown(f"<p style='margin-bottom:33px'></p>", unsafe_allow_html=True)
            clicked = st.button(options[i], key=f'visual-{options[i]}')
            if clicked:
                kind = options[i]
                assay.add_metadata(DFT.VISUAL_TYPE, [category, kind])

    if kind in DFT.LAYOUT:
        columns = st.beta_columns(DFT.LAYOUT[kind])
        args_conatiner = columns[0]
        plot_columns = columns[1:]
    else:
        columns = st.beta_columns([0.75, 0.1, 2])
        args_conatiner = columns[0]
        plot_columns = columns[2]

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
            kwargs['points'] = st.checkbox('Box and points', False)
            features = st.multiselect('Features', list(assay.ids()), list(assay.ids())[:min(len(assay.ids()), 4)])
            if len(features) != 0:
                kwargs['features'] = features
        elif kind == DFT.RIDGEPLOT:
            kwargs['attribute'] = st.selectbox('Attribute', DFT.LAYERS[assay.name])
            kwargs['splitby'] = st.selectbox('Split by', DFT.SPLITBY[assay.name])
            features = st.multiselect('Features', list(assay.ids()), list(assay.ids())[:min(len(assay.ids()), 4)])
            if len(features) != 0:
                kwargs['features'] = features
        elif kind == DFT.STRIPPLOT:
            kwargs['attribute'] = st.selectbox('Attribute', DFT.LAYERS[assay.name])
            kwargs['colorby'] = st.selectbox('Colorby', DFT.LAYERS[assay.name])
            features = st.multiselect('Features', list(assay.ids()), list(assay.ids())[:min(len(assay.ids()), 4)])
            if len(features) != 0:
                kwargs['features'] = features
        elif kind == DFT.DNA_PROTEIN_PLOT:
            kwargs['analyte'] = st.selectbox('Analyte', ['protein'], format_func=lambda a: analyte_map[a])
            kwargs['dna_features'] = st.multiselect('DNA features', list(sample.dna.ids()), sample.dna.ids()[:4])
            kwargs['protein_features'] = st.multiselect('Protein features', list(sample.protein.ids()), sample.protein.ids()[:4])
        elif kind == DFT.DNA_PROTEIN_HEATMAP:
            kwargs['clusterby'] = st.selectbox('Cluster by', ['dna', 'protein'], format_func=lambda a: analyte_map[a])
            kwargs['sortby'] = st.selectbox('Sort by', ['dna', 'protein'], format_func=lambda a: analyte_map[a])
            kwargs['dna_features'] = st.multiselect('DNA features', list(sample.dna.ids()), sample.dna.ids())
            kwargs['protein_features'] = st.multiselect('Protein features', list(sample.protein.ids()), sample.protein.ids())

        elif kind == DFT.METRICS:
            st.header('')
            interface.info('<b>Some values might be missing in case the raw<br> files are not loaded.</b> These metrics can be<br> pasted into the metrics sheet as is.')
        elif kind == DFT.READ_DEPTH:
            if assay.name == PROTEIN_ASSAY:
                kwargs['layer'] = st.selectbox('Layer', DFT.LAYERS[assay.name])
                kwargs['colorby'] = st.selectbox('Color by', ['density', None])
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


def visual(sample, assay, kind, plot_columns, kwargs):
    if kind in [DFT.HEATMAP, DFT.SCATTERPLOT, DFT.FEATURE_SCATTER, DFT.VIOLINPLOT, DFT.RIDGEPLOT, DFT.STRIPPLOT, DFT.DNA_PROTEIN_HEATMAP]:
        with plot_columns:
            fig = draw_plots(sample, assay, kind, kwargs)
            st.plotly_chart(fig)

    elif kind == DFT.VAR_ANNOTATIONS:
        df = get_annotations(sample.dna.ids())
        with st.beta_columns([0.2, 10, 1])[1]:
            st.dataframe(df, height=1080 * 4)
    elif kind == DFT.SIGNATURES:
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
    elif kind == DFT.METRICS:
        with plot_columns:
            st.header('')
            metrics(sample)
    elif kind == DFT.COLORS:
        colors = COLORS.copy()
        del colors[20]

        for i in range(len(plot_columns)):
            plot_columns[i].header('')

        for i in range(len(colors)):
            plot_columns[i % len(plot_columns)].color_picker(colors[i], colors[i], key=f'constant-colors-{colors[i]}-{i}')
    elif kind == DFT.DOWNLOAD:
        with plot_columns:
            if kwargs['download']:
                if kwargs['item'] == DFT.ANNOTATION:
                    data = np.array([sample.dna.barcodes(), sample.dna.get_labels(), sample.protein.get_labels()]).T
                    df = pd.DataFrame(data, columns=['barcode', 'dna', 'protein'])
                    name = f'./{sample.name}.annotation.csv'
                    df.to_csv(name, index=None)
                    interface.download(name)
                    os.remove(name)

    elif kind == DFT.DNA_PROTEIN_PLOT:
        with plot_columns:
            samp = mosample(protein=sample.protein[:, kwargs['protein_features']], dna=sample.dna[:, kwargs['dna_features']])
            samp.clone_vs_analyte(kwargs['analyte'])
            st.pyplot(plt.gcf())
    elif kind == DFT.READ_DEPTH:
        with plot_columns:
            if assay.name == PROTEIN_ASSAY:
                fig = draw_plots(sample, assay, kind, kwargs)
                st.plotly_chart(fig)
    elif kind == DFT.ASSAY_SCATTER:
        if kwargs['draw']:
            with plot_columns:
                samp = sample[:]
                samp.cnv_raw.normalize_barcodes()
                samp.protein_raw.normalize_barcodes()
                samp.assay_scatter()
                st.pyplot(plt.gcf())


@st.cache(max_entries=50, hash_funcs=DFT.MOHASH, show_spinner=False, allow_output_mutation=True, ttl=3600)
def draw_plots(sample, assay, kind, kwargs):
    if kind in [DFT.HEATMAP, DFT.SCATTERPLOT, DFT.FEATURE_SCATTER, DFT.VIOLINPLOT, DFT.RIDGEPLOT, DFT.STRIPPLOT]:
        plot_funcs = {
            DFT.HEATMAP: assay.heatmap,
            DFT.SCATTERPLOT: assay.scatterplot,
            DFT.FEATURE_SCATTER: assay.feature_scatter,
            DFT.VIOLINPLOT: assay.violinplot,
            DFT.RIDGEPLOT: assay.ridgeplot,
            DFT.STRIPPLOT: assay.stripplot
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
            new_pal = sample.protein.get_palette()
            new_lab = sample.protein.get_labels().copy()
            kwargs[labelby] = 'label'

        if kwargs[labelby] == DFT.DNA_LABEL:
            new_pal = sample.dna.get_palette()
            new_lab = sample.dna.get_labels().copy()
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

        if kind == DFT.VIOLINPLOT:
            update = kwargs['points']
            del kwargs['points']

        fig = plot_funcs[kind](**kwargs)

        if kind == DFT.VIOLINPLOT and update:
            fig.update_traces(points='all', pointpos=-0.5, box_width=0.6, side='positive', marker_size=3)

        assay.set_labels(org_lab)
        assay.set_palette(org_pal)

    elif kind == DFT.DNA_PROTEIN_HEATMAP:
        samp = mosample(protein=sample.protein[:, kwargs['protein_features']], dna=sample.dna[:, kwargs['dna_features']])
        fig = samp.heatmap(clusterby=kwargs['clusterby'], sortby=kwargs['sortby'], drop='cnv', flatten=False)

    elif kind == DFT.READ_DEPTH:
        total_reads = assay.layers[READS].sum(axis=1)
        layer = assay.layers[kwargs['layer']]

        x = np.log10(total_reads)

        feats = kwargs['features']
        if len(feats) == 0:
            feats = assay.ids()

        nplots = len(feats)
        nrows = round(nplots ** 0.5)
        ncols = nplots // nrows + min(1, nplots % nrows)

        fig = make_subplots(rows=nrows, cols=ncols,
                            x_title='log10(Total reads)',
                            y_title=kwargs['layer'],
                            subplot_titles=feats)

        for i in range(len(feats)):
            row_num = i // ncols + 1
            col_num = i % ncols + 1
            feat = feats[i]

            y = layer[:, np.where(assay.ids() == feat)[0]].flatten()
            data = np.array([x, y]).T
            scatter = assay.scatterplot(attribute=data, colorby=kwargs['colorby'])
            fig.add_trace(scatter.data[0], row=row_num, col=col_num)

        fig.update_yaxes(title='')
        fig.update_xaxes(title='')
        layout = scatter.layout
        layout.update(coloraxis=dict(colorbar=dict(thickness=25,
                                                   len=1 / nrows,
                                                   yanchor="top",
                                                   y=1.035,
                                                   x=1.05,
                                                   ticks="outside")))
        fig.update_layout(layout,
                          width=max(500, 300 * ncols),
                          height=max(500, max(300 * nrows, 30 * len(feats) * nrows)))

    return fig


@st.cache(show_spinner=False)
def get_annotations(variants):
    interface.status('Fetching DNA annotations')
    renamed_variants = np.array([var.replace(':', '-').replace('/', '-') for var in variants], dtype='str')

    url = 'https://api.missionbio.io/annotations/v1/variants?ids=' + ','.join(renamed_variants)
    r = requests.get(url=url)

    data = r.json()
    data = [d['annotations'] for d in data]

    function = [', '.join(d['function']['value']) for d in data]
    gene = [d['gene']['value'] for d in data]
    protein = [d['protein']['value'] for d in data]
    coding_impact = [d['protein_coding_impact']['value'] for d in data]
    clinvar = [', '.join(d['clinvar']['value']) for d in data]
    dann = np.array([d['impact']['value'] for d in data])
    dann[dann == ''] = 0
    dann = np.round(dann.astype(float), 2)

    annot_types = ['Gene', 'Function', 'Protein', 'Coding Impact', 'ClinVar', 'DANN']
    df = pd.DataFrame([gene, function, protein, coding_impact, clinvar, dann], index=annot_types).T
    df['Variant'] = variants

    df = df[['Variant'] + annot_types]

    return df


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
