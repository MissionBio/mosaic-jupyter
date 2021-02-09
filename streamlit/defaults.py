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

from missionbio.mosaic.dna import Dna as modna
from missionbio.mosaic.protein import Protein as moprotein
from missionbio.mosaic.sample import Sample as mosample


MOHASH = {moprotein: lambda a: a.name + a.title + str(a.shape),
          modna: lambda a: a.name + a.title + str(a.shape),
          mosample: lambda a: a.name + str(a.dna.shape),
          type(lambda _: None): lambda _: None}

MOHASH_BOOL = {**MOHASH, int: lambda _: None, list: lambda _: None, type(''): lambda _: None}

ROOT = pathlib.Path(__file__).parent

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
INITIALIZE = '__mosaic_initialize'

# ----------- Preprocessing
MIN_DP = 10
MIN_GQ = 30
MIN_VAF = 20
MIN_STD = 20

PREPROCESS_ARGS = '__mosaic_preprocess_args'
DROP_IDS = '__mosaic_drop_ids'
KEEP_IDS = '__mosaic_keep_ids'
ALL_IDS = '__mosaic_all_ids'

# ----------- Preparing
PREPPED = '__mosaic_prepped'

# ----------- Clustering
CLUSTERED = '__mosaic_clustered'
CLUSTER_DESCRIPTION = '__mosaic_cluster_description'

SCALE_ATTR = '__mosaic_data_prep_scale'
PCA_ATTR = '__mosaic_data_prep_pca'
UMAP_ATTR = '__mosaic_data_prep_umap'


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
VISUAL_TYPE = '__mosaic_visual_type'

VISUAL_CATEGORY_1 = 'Plots'
VISUAL_CATEGORY_2 = 'Tables'
VISUAL_CATEGORY_3 = 'Multiomic plots'

HEATMAP = 'Heatmap'
SCATTERPLOT = 'Scatterplot'

FEATURE_SCATTER = 'Feature Scatter'
VIOLINPLOT = 'Violinplot'
RIDGEPLOT = 'Ridgeplot'
STRIPPLOT = 'Stripplot'

VAR_ANNOTATIONS = 'Annotations'
SIGNATURES = 'Signatures'
COLORS = 'Colors'
METRICS = 'Metrics'
DOWNLOAD = 'Download'

DNA_PROTEIN_PLOT = 'DNA vs Protein'
DNA_PROTEIN_HEATMAP = 'Sample Heatmap'
READ_DEPTH = 'Read Depth'
ASSAY_SCATTER = 'Assay Scatter'

LAYOUT = {
    COLORS: [0.5] + [1] * 10,
}

VISUALS = {VISUAL_CATEGORY_1: [[2, 1, 1.1, 1.45, 1, 1, 1, 5],
                               [HEATMAP,
                                SCATTERPLOT,
                                FEATURE_SCATTER,
                                VIOLINPLOT,
                                RIDGEPLOT,
                                STRIPPLOT]],

           VISUAL_CATEGORY_2: [[2, 1, 1.15, 0.75, 0.75, 1, 6.8],
                               [SIGNATURES,
                                VAR_ANNOTATIONS,
                                COLORS,
                                METRICS,
                                DOWNLOAD]],

           VISUAL_CATEGORY_3: [[1.8, 1.25, 1.35, 1, 1.5, 5],
                               [DNA_PROTEIN_PLOT,
                                DNA_PROTEIN_HEATMAP,
                                READ_DEPTH,
                                ASSAY_SCATTER]]}

# ----------- Visuals options
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
    DNA_ASSAY: [DNA_LABEL, PROTEIN_LABEL, 'sample_name', AF, AF_MISSING, NGT, NGT_FILTERED, GQ, DP, 'density', None],
    PROTEIN_ASSAY: [PROTEIN_LABEL, DNA_LABEL, 'sample_name', CLR, ASINH, NSP, SCALED_LABEL, 'density', None]
}

ANNOTATION = 'Annotation'
DOWNLOAD_ITEMS = [ANNOTATION]
