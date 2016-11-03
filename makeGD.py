#!/usr/bin/env python
"""
This script can be run from the command line and it creates geodata h5s
based on which satellites I want and how high I want the IPP
(ionospheric peirce point)
NOTE: this only works on my computer, all the directories are specified
only on my computer. if someone else wants to use this, they must change
the file paths
"""

import GeoData.GeoData as GD
from pathlib import Path
from gps import GDfromRinex

dpath = 'Examples/data'
#%%
dpath = Path(dpath).expanduser()
files = dpath.glob('mah*.h5')
navfile = dpath/'brdc2800.15n'
satfile = dpath/'jplg2800.15i'
C1BiasFile = dpath/'P1C11510.DCB'


if __name__== '__main__':

    from argparse import ArgumentParser
    descr = "this will create GeoData h5s from rinex/h5 at a certain pph"

    p = ArgumentParser(description=descr)
    p.add_argument('-s','--sats',help='satellites to include',default='9 23')
    p.add_argument("-a", "--height",help='pierce point height',default=130)
    p = p.parse_args()

    satlist = [int(a) for a in p.sats.split()]
    height = int(p.height)

    files = sorted(files)
    if not files:
        raise FileNotFoundError('No {}/mah*.h5 files found'.format(dpath))

    for f in files:
        site = f.stem[:4]
        rinexfile = sorted(dpath.glob('{}*.15o'.format(site)))[0]
        a = GD.GeoData(GDfromRinex,(rinexfile,navfile,satfile,C1BiasFile,f,False,height,satlist))

        tail = ''
        for i in satlist:
            tail += 'sat'+str(i)

        fn =  dpath/('h5files/{}km'.format(height))
        fn.mkdir(parents=True,exist_ok=True)

        fn = fn / '{}GD{}.h5'.format(site, tail)

        print('writing {}'.format(fn))

        a.write_h5(str(fn))

