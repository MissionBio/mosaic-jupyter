import streamlit as st
import pathlib

from missionbio.h5.constants import (
    DNA_ASSAY,
    PROTEIN_ASSAY,
    NGT,
    AF,
    GQ,
    DP
)

from missionbio.mosaic.constants import (
    PCA_LABEL,
    UMAP_LABEL,
    SCALED_LABEL,
    AF_MISSING,
    NGT_FILTERED
)


STREAMLIT_STATIC_PATH = pathlib.Path(st.__path__[0]) / 'static'
# We create a downloads directory within the streamlit static asset directory
# and we write output files to it
DOWNLOADS_PATH = (STREAMLIT_STATIC_PATH / "downloads")
if not DOWNLOADS_PATH.is_dir():
    DOWNLOADS_PATH.mkdir()

# ----------- Interface
RED = '#BF3035'
BLUE = '#1C5C6C'

# ----------- Loading
S3 = 's3'
LOCAL = 'local'

RUN_NUMBER = '__mosaic_run_number'
PREPPED = '__mosaic_prep_arguments'
CLUSTER_DESCRIPTION = '__mosaic_cluster_description'

SCALE_ATTR = '__mosaic_data_prep_scale'
PCA_ATTR = '__mosaic_data_prep_pca'
UMAP_ATTR = '__mosaic_data_prep_umap'

# ----------- Preprocessing
MIN_DP = 10
MIN_GQ = 30
MIN_VAF = 20
MIN_STD = 20

FILTERS = '__mosaic_filters'


CLUSTER_METHOD = {
    DNA_ASSAY: [UMAP_LABEL, 'dbscan', 0.2],
    PROTEIN_ASSAY: [PCA_LABEL, 'graph-community', 30]
}

CLUSTER_OPTIONS = {
    'dbscan': ('Proximity', 0.05, 2.0, 0.2, 'eps'),
    'hdbscan': ('Cluster size', 10, 500, 100, 'min_cluster_size'),
    'kmeans': ('Neighbours', 2, 30, 5, 'n_clusters'),
    'graph-community': ('Neighbours', 10, 500, 100, 'k')
}

DNA_LABEL = 'dna_label'
PROTEIN_LABEL = 'protein_label'

CLR = 'CLR'
ASINH = 'asinh'
NSP = 'NSP'
TOTAL_READS = 'total_reads'

# ----------- Primary visuals
PRIMARY = 'Primary'

SIGNATURES = 'Signatures'
HEATMAP = 'Heatmap'
SCATTERPLOT = 'Scatterplot'
FEATURE_SCATTER = 'Feature Scatter'
VIOLINPLOT = 'Violinplot'
RIDGEPLOT = 'Ridgeplot'
DNA_PROTEIN_PLOT = 'DNA vs Protein'
DNA_PROTEIN_HEATMAP = 'Sample Heatmap'

SPLITBY = {
    DNA_ASSAY: [DNA_LABEL, PROTEIN_LABEL, 'sample_name', None],
    PROTEIN_ASSAY: [PROTEIN_LABEL, DNA_LABEL, 'sample_name', None]
}
ATTRS_2D = [UMAP_LABEL]
LAYERS = {
    DNA_ASSAY: [AF_MISSING, AF, NGT, NGT_FILTERED, GQ, DP],
    PROTEIN_ASSAY: [CLR, ASINH, NSP, SCALED_LABEL]
}
COLORBY = {
    DNA_ASSAY: [DNA_LABEL, PROTEIN_LABEL, 'sample_name', AF, AF_MISSING, NGT, NGT_FILTERED, 'density', None],
    PROTEIN_ASSAY: [PROTEIN_LABEL, DNA_LABEL, 'sample_name', CLR, ASINH, NSP, SCALED_LABEL, 'density', None]
}

# ----------- Advanced visuals
ADVANCED = 'Advanced'
QC = 'QC'
DOWNLOAD = 'Download'
METRICS = 'Metrics'
READ_DEPTH = 'Read Depth'
ASSAY_SCATTER = 'Assay Scatter'
COLORS = 'Colors'

LAYOUT = {
    COLORS: [1.5, 0.5] + [1] * 10,
    METRICS: [0.5, 1, 2],
    READ_DEPTH: [0.5, 0.75, 2],
    ASSAY_SCATTER: [0.5, 0.75, 2],
    QC: [0.5, 0.75, 2],
    DOWNLOAD: [0.5, 0.75, 2]
}

ANNOTATION = 'Annotation'
DOWNLOAD_ITEMS = [ANNOTATION]
