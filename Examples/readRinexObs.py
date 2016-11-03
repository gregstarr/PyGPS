#!/usr/bin/env python3
"""
RINEX 2 OBS reader
under testing
Michael Hirsch, Greg Starr
MIT License

Program overviw:
1) read the OBS file header:   readHead()
2) parse the OBS file header, obtaining the times, satellites, and data measurement types:   makeSvset()
3) read the OBS files in blocks, where each block is one time interval of data (all sats, all measurements):  makeBlock()

makeBlock dumps the results into a preallocated pandas 3-D Panel with axes:
items / page: time
rows / major_axis: SV
column / minor_axis: data type P1,P2, etc.
"""
from __future__ import division,absolute_import
import numpy as np
from itertools import chain
from datetime import datetime, timedelta
from pandas import DataFrame,Panel,Series,Panel4D
from pandas.io.pytables import read_hdf
from os.path import splitext,expanduser,getsize
from io import BytesIO
import time

def rinexobs(obsfn,writeh5=None,maxtimes=None):
    stem,ext = splitext(expanduser(obsfn))
    if ext[-1].lower() == 'o': #raw text file
        with open(obsfn,'r') as f:
            t=time.time()
            lines = f.read().splitlines(True)
            lines.append('')
            header,version,headlines,obstimes,sats,svset = scan(lines)
            print('{} is a RINEX {} file, {} kB.'.format(obsfn,version,getsize(obsfn)/1000.0))
            data = processBlocks(lines,header,obstimes,svset,headlines,sats)
            print("finished in {0:.2f} seconds".format(time.time()-t))
    #%% save to disk (optional)
        if writeh5:
            h5fn = stem + '.h5'
            print('saving OBS data to {}'.format(h5fn))
            data.to_hdf(h5fn,key='OBS',mode='a',complevel=6,append=False)
    elif ext.lower() == '.h5':
        data = read_hdf(obsfn,key='OBS')
        print('loaded OBS data from {} to {}'.format(blocks.items[0],blocks.items[-1]))
    return data


# this will scan the document for the header info and for the line on
# which each block starts
def scan(lines):
    header={}        
    eoh=0
    for i,line in enumerate(lines):
        if "END OF HEADER" in line:
            eoh=i
            break
        if line[60:].strip() not in header:
            header[line[60:].strip()] = line[:60].strip()
        else:
            header[line[60:].strip()] += " "+line[:60].strip()
    verRinex = float(header['RINEX VERSION / TYPE'].split()[0])
    header['APPROX POSITION XYZ'] = [float(i) for i in header['APPROX POSITION XYZ'].split()]
    header['# / TYPES OF OBSERV'] = header['# / TYPES OF OBSERV'].split()
    header['# / TYPES OF OBSERV'][0] = int(header['# / TYPES OF OBSERV'][0])
    header['INTERVAL'] = float(header['INTERVAL'])
        
    headlines=[]
    obstimes=[]
    sats=[]
    svset=set()
    i=eoh+1
    while True:
        if not lines[i]: break
        if not int(lines[i][28]):
            #no flag or flag=0
            headlines.append(i)
            obstimes.append(_obstime([lines[i][1:3],lines[i][4:6],
                                   lines[i][7:9],lines[i][10:12],
                                   lines[i][13:15],lines[i][16:26]]))
            numsvs = int(lines[i][30:32])
            if(numsvs>12):
                sats.append([int(lines[i][33+s*3:35+s*3]) for s in range(12)])
                i += 1
                sats.append([int(lines[i][33+s*3:35+s*3]) for s in range(numsvs-12)])
            else:
                sats.append([int(lines[i][33+s*3:35+s*3]) for s in range(numsvs)])
        
            i+=numsvs*int(np.ceil(header['# / TYPES OF OBSERV'][0]/5))+1
        else:
            #there was a comment or some header info
            flag=int(lines[i][28])
            if(flag!=4):
                print(flag)
            skip=int(lines[i][30:32])
            i+=skip+1
    for s in sats:
        svset = svset.union(set(s))

    return header,verRinex,headlines,obstimes,sats,svset



def processBlocks(lines,header,obstimes,svset,headlines,sats):
    obstypes = header['# / TYPES OF OBSERV'][1:]
    blocks = Panel4D(labels=obstimes,
                     items=list(svset),
                     major_axis=obstypes,
                     minor_axis=['data','lli','ssi'])
    ttime1 = 0
    ttime2 = 0
    for i in range(len(headlines)):
        linesinblock = len(sats[i])*int(np.ceil(header['# / TYPES OF OBSERV'][0]/5))
        block = ''.join(lines[headlines[i]+1:headlines[i]+linesinblock+1])
        t1 = time.time()
        bdf = _block2df(block,obstypes,sats[i],len(sats[i]))
        ttime1 += (time.time()-t1)
        t2 = time.time()
        blocks.loc[obstimes[i],sats[i]] = bdf
        ttime2 += (time.time()-t2)            
    print("{0:.2f} seconds for _block2df".format(ttime1))
    print("{0:.2f} seconds for panel assignments".format(ttime2))
    return blocks       
        

def _obstime(fol):
    year = int(fol[0])
    if 80<= year <=99:
        year+=1900
    elif year<80: #because we might pass in four-digit year
        year+=2000
    return datetime(year=year, month=int(fol[1]), day= int(fol[2]),
                    hour= int(fol[3]), minute=int(fol[4]),
                    second=int(float(fol[5])),
                    microsecond=int(float(fol[5]) % 1) *100000
                    )

def _block2df(block,obstypes,svnames,svnum):
    """
    input: block of text corresponding to one time increment INTERVAL of RINEX file
    output: 2-D array of float64 data from block. Future: consider whether best to use Numpy, Pandas, or Xray.
    """
    nobs = len(obstypes)
    stride=3

    strio = BytesIO(block.encode())
    barr = np.genfromtxt(strio, delimiter=(14,1,1)*5).reshape((svnum,-1), order='C')

    data = barr[:,0:nobs*stride:stride]
    lli  = barr[:,1:nobs*stride:stride]
    ssi  = barr[:,2:nobs*stride:stride]

    data = np.vstack(([data.T],[lli.T],[ssi.T])).T

    return data

