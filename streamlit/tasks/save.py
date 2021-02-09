import os

import interface
import defaults as DFT

import missionbio.mosaic.io as mio


def run(sample, name, should_save):
    for assay, og_assay in zip([sample.dna, sample.protein], [sample._original_dna, sample._original_protein]):
        if assay is not None:
            for key in assay.metadata:
                og_assay.add_metadata(key, assay.metadata[key])

            for key in assay.row_attrs:
                og_assay.add_row_attr(key, assay.row_attrs[key])

    if should_save:
        interface.status('Saving h5 file.')
        if name == '':
            interface.error('Please provide a name to save by.')
        elif name[-3:] == '.h5':
            name = name[:-3]

        try:
            os.remove(DFT.ROOT / f'h5/analyzed/{name}.h5')
        except FileNotFoundError:
            pass

        samp = sample[:]
        set_defaults(samp)
        mio.save(samp, DFT.ROOT / f'h5/analyzed/{name}.h5')

        interface.status('Saved.')
        interface.rerun()


def set_defaults(sample):
    def del_arg(assay, key):
        if assay is not None:
            assay.del_metadata(key)

    del_arg(sample.dna, DFT.ALL_IDS)
    del_arg(sample.dna, DFT.DROP_IDS)
    del_arg(sample.dna, DFT.KEEP_IDS)

    del_arg(sample.protein, DFT.ALL_IDS)
    del_arg(sample.protein, DFT.DROP_IDS)
