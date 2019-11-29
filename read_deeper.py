import numpy as np
import pandas as pd
import codecs
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

from increase_accuracy import kalmanfilter_2d, judge_stability
from local_map import LocalMap
from hexagonal_tiling import Honeycomb
from transform_coordination import latlon2enu, enu2latlon 

__all__ = ['read_deeper','deeper_csv2localdata','deeper_csv2localdata','delete_landing']

def read_deeper(input_file:str ,*,drop:bool = True):
    d = pd.read_csv(input_file,names = ['lat','lon','depth','t'])
    if drop:
        d = d[d.lat != 0]
        d = d[d.lon != 0]
        d = d[d.t != 0]
        d.index = range(len(d))
        unique_t,unique_index = np.unique(d.t,return_index=True)
        d = d.loc[unique_index]
        d = d.sort_values(by='t')
        d.index = range(len(d))
    return d

def deeper_csv2LocalMap(input_files,*,localmap=LocalMap()):
    lm = localmap
    files = np.array(input_files)
    files = files.reshape(len(files),1)
    for i in range(len(files[:])):
        d = read_deeper(files[i,0])
        e,n,u = latlon2enu(d.loc[:,'lat'],d.loc[:,'lon'],0,lm.origin[0],lm.origin[1],0)
        e,n = kalmanfilter_2d(e,n,d.loc[:,'t'])
        depth0 = d.loc[:,'depth']
        judge = judge_stability(d.loc[:,'t'])
        e = e[judge]
        n = n[judge]
        depth0 = depth0[judge]
        if i == 0:
            x = e
            y = n
            depth =  depth0
        else :
            x = np.hstack((x,e))
            y = np.hstack((y,n))
            depth = np.hstack((depth,depth0))
    hc = Honeycomb(1,x=[lm.xrange[0],lm.xrange[1]],y=[lm.yrange[0],lm.yrange[1]])
    depth = hc.depth_weighted_mean(x,y,depth)
    x = hc.voronoi_point[:,0]
    y = hc.voronoi_point[:,1]
    judge = ~np.isnan(depth)
    x = x[judge]
    y = y[judge]
    depth = depth[judge]
    lm.fill_grid(x,y,depth)
    return lm,x,y,depth

def deeper_csv2localdata(input_files,*,localmap=LocalMap()):
    files = np.array(input_files)
    files = files.reshape(len(files),1)
    for i in range(len(files[:])):
        d = read_deeper(files[i,0])
        e,n,u = latlon2enu(d.loc[:,'lat'],d.loc[:,'lon'],0,localmap.origin[0],localmap.origin[1],0)
        e,n = kalmanfilter_2d(e,n,d.loc[:,'t'])
        depth0 = d.loc[:,'depth']
        judge = judge_stability(d.loc[:,'t'])
        e = e[judge]
        n = n[judge]
        depth0 = depth0[judge]
        if i == 0:
            x = e
            y = n
            depth =  depth0
        else :
            x = np.hstack((x,e))
            y = np.hstack((y,n))
            depth = np.hstack((depth,depth0))
    hc = Honeycomb(1,x=[localmap.xrange[0],localmap.xrange[1]],y=[localmap.yrange[0],localmap.yrange[1]])
    depth = hc.depth_weighted_mean(x,y,depth)
    x = hc.voronoi_point[:,0]
    y = hc.voronoi_point[:,1]
    judge = ~np.isnan(depth)
    x = x[judge]
    y = y[judge]
    lat,lon = enu2latlon(x,y,0,localmap.origin[0],localmap.origin[1],0)
    depth = depth[judge]
    data = np.hstack((lat,lon,x,y,depth))
    data = data.reshape(len(depth),5,order='F')
    localdata = pd.DataFrame(data,columns=['lat','lon','x','y','depth'])
    return localdata

def delete_landing(localmap,coastlines):
    dele = np.array([35.628,139.780])
    for i in range(len(coastlines[0,0,:])):
        for j in range(len(coastlines[:,0,i])):
            lat = coastlines[j,0,i]
            lon = coastlines[j,1,i]
            if lat >= 35.6238 and lat <= 35.6286 and lon >= 139.780 and lon <= 139.78501:
                dele = np.vstack((dele,[lat,lon]))
    tri = Delaunay(dele)
    for i in range(len(localmap.latlabel)):
        for j in range(len(localmap.lonlabel)):
            sim = Delaunay.find_simplex(tri,[localmap.latlabel[i],localmap.lonlabel[j]])
            if sim != -1:
                localmap.values[i,j] = np.nan
    return localmap

