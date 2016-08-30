# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 09:33:02 2016

@author: greg
"""

import GeoData.GeoData as GD
from glob import glob
from gps import GDfromRinex

files = glob("/home/greg/Documents/greg/rinex/mah*.h5")
navfile = '/home/greg/Documents/greg/brdc2800.15n'
satfile = '/home/greg/Documents/greg/jplg2800.15i'
C1BiasFile = '/home/greg/Documents/greg/P1C11510.DCB'


if __name__== '__main__':
    
    from argparse import ArgumentParser
    descr = "this will create GeoData h5s from rinex/h5 at a certain pph"
    
    p = ArgumentParser(description=descr)
    p.add_argument('-s','--sats',help='satellites to include',default=None)
    p.add_argument("-a", "--height",help='pierce point height',default=130)
    p = p.parse_args()
    
    satlist = [int(a) for a in p.sats.split()]
    height = int(p.height)

    for h5file in files:
        site = h5file.split('/')[-1][:4]
        rinexfile = glob('/home/greg/Documents/greg/rinex/{}*.15o'.format(site))[0]
        a = GD.GeoData(GDfromRinex,(rinexfile,navfile,satfile,C1BiasFile,h5file,False,height,satlist))
        fn = '/home/greg/Documents/greg/h5files/'+str(height)+'km/'+site+'GD'
        for i in satlist:
            fn += 'sat'+str(i)+'.h5'
        a.write_h5(fn)
    
    