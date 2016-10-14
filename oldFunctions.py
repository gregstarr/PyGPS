def getRanges(data,svn,maxgap=3,maxjump=1.0): 
    if c2p2(data,svn):
        nans = np.logical_or.reduce((np.isnan(data[['L1','L2','C1','C2'],svn,:,'data']).T))
    else:
        nans = np.logical_or.reduce((np.isnan(data[['L1','L2','C1','P2'],svn,:,'data']).T))
    inarc=False
    start=[]
    end=[]
    phase=2.85E9*(data['L1',svn,:,'data']/f1-data['L2',svn,:,'data']/f2)
    lgi=0
    for i in range(len(nans)):
        if inarc:
            if nans[i]:
                if i-lgi>maxgap:
                    inarc=False
                    end.append(lgi)
            else:
                if abs(phase[i]-phase[lgi])>maxjump:
                    end.append(lgi)
                    start.append(i)
                lgi=i

        else:
            if not nans[i]:
                inarc=True
                start.append(i)
                lgi=i
    if len(start)!=len(end): end.append(i)
    ranges = [(data.major_axis[a],data.major_axis[b]) for a,b in zip(start,end)]

    return ranges
    
    
###############
#Shouldn't work
###############
def c2p2alt(data,site,drange=(None,None)):
    return (np.sum(~np.isnan(
        data[site,drange[0]:drange[1],'C2','data']))>
            np.sum(~np.isnan(
                data[site,drange[0]:drange[1],'P2','data'])))
###############
def getTecalt(data,site,drange,satbias=None):
    if c2p2alt(data,site,drange):
        diffRange = (2.85E9/3.0E8)*(
            data[site,drange[0]:drange[1],'C2','data']
            -data[site,drange[0]:drange[1],'C1','data'])
    else:
        diffRange = (2.85E9/3.0E8)*(
            data[site,drange[0]:drange[1],'P2','data']
            -data[site,drange[0]:drange[1],'C1','data'])
        
    phase=2.85E9*(data[site,drange[0]:drange[1],'L1','data']/f1
                  -data[site,drange[0]:drange[1],'L2','data']/f2)        
        
    diffList = sorted(phase-diffRange)
    medianDiff = diffList[int(len(diffList)/2)]
    distWidth = diffList[int(len(diffList)*.75)]-diffList[int(len(diffList)*.25)]
    medianErr = distWidth/np.sqrt(len(diffList))
    tec = phase - medianDiff
    if satbias!=None:
        tec-=satbias
    return tec,medianErr
####################
def getRangesalt(data,site,maxgap=3,maxjump=2.0): 
    if c2p2alt(data,site):
        nans = np.logical_or.reduce((np.isnan(data[site,:,'L1','data']),
                                     np.isnan(data[site,:,'L2','data']),
                                     np.isnan(data[site,:,'C1','data']),
                                     np.isnan(data[site,:,'C2','data'])))
    else:
        nans = np.logical_or.reduce((np.isnan(data[site,:,'L1','data']),
                                     np.isnan(data[site,:,'L2','data']),
                                     np.isnan(data[site,:,'C1','data']),
                                     np.isnan(data[site,:,'P2','data'])))
    inarc=False
    start=[]
    end=[]
    phase=2.85E9*(data[site,:,'L1','data']/f1-data[site,:,'L2','data']/f2)
    lgi=0

    for i in range(len(nans)):
        if inarc:
            if nans[i]:
                if i-lgi>maxgap:
                    inarc=False
                    end.append(lgi)
            else:
                if abs(phase[i]-phase[lgi])>maxjump:
                    end.append(lgi)
                    start.append(i)
                lgi=i

        else:
            if not nans[i]:
                inarc=True
                start.append(i)
                lgi=i
    
    if len(start)!=len(end): end.append(i)
    ranges = [(data.items[a],data.items[b]) for a,b in zip(start,end)]

    return ranges
####################
