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
    sample.dna.del_metadata(DFT.ALL_IDS)
    sample.dna.del_metadata(DFT.DROP_IDS)
    sample.dna.del_metadata(DFT.KEEP_IDS)

    sample.protein.del_metadata(DFT.ALL_IDS)
    sample.protein.del_metadata(DFT.DROP_IDS)
