import os

import interface
import defaults as DFT

import missionbio.mosaic.io as mio


def run(sample, name):
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
