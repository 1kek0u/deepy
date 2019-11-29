import numpy as np
from scipy.spatial import Delaunay
from transform_coordination import Ellipsoid, latlon2enu, enu2latlon 

__all__ = ['LocaMap']

class LocalMap:

    def __init__(self,*,lat=[35.6242,35.626612],lon=[139.7825,139.785514],grid_interval=2,ell=Ellipsoid()):
        Ell = ell
        self.interval = grid_interval
        self.origin = [sum(lat)/2,sum(lon)/2]
        e0,n0,u0 = latlon2enu(lat[0],lon[0],0,self.origin[0],self.origin[1],0,ell=Ell)
        e1,n1,u1 = latlon2enu(lat[1],lon[1],0,self.origin[0],self.origin[1],0,ell=Ell)
        self.xrange = np.array([e0,e1])
        self.yrange = np.array([n0,n1])
        self.xlabel = np.arange(self.xrange[0],(self.xrange[1]+self.interval),self.interval)
        self.ylabel = np.arange(self.yrange[0],(self.yrange[1]+self.interval),self.interval)
        latlabel = enu2latlon(0,self.ylabel,0,self.origin[0],self.origin[1],0)
        lonlabel = enu2latlon(self.xlabel,0,0,self.origin[0],self.origin[1],0)
        self.latlabel = latlabel[0]
        self.lonlabel = lonlabel[1]
        self.latrange = [self.latlabel[0],self.latlabel[-1]]
        self.lonrange = [self.lonlabel[0],self.lonlabel[-1]]
        a = np.zeros([len(self.ylabel),len(self.xlabel)])
        a[:,:] = np.nan
        self.values = a

    def fill_grid(self,x,y,depth):
        x = np.array(x)
        y = np.array(y)
        x = x.reshape(len(x),1)
        y = y.reshape(len(y),1)
        xy = np.hstack((x,y))
        tri = Delaunay(xy)
        print(xy)
        for i in range(len(self.xlabel)):
            for j in range(len(self.ylabel)):
                sim = tri.find_simplex([self.xlabel[i],self.ylabel[j]])
                if sim != -1:
                    p = xy[tri.simplices[sim]] 
                    a = 1/(np.linalg.norm(p[0] - [self.xlabel[i],self.ylabel[j]]) + self.interval/1000)
                    b = 1/(np.linalg.norm(p[1] - [self.xlabel[i],self.ylabel[j]]) + self.interval/1000)
                    c = 1/(np.linalg.norm(p[2] - [self.xlabel[i],self.ylabel[j]]) + self.interval/1000)
                    depth_a = depth[np.arange(len(xy))[(xy[:,0] == p[0][0]) & (xy[:,1] == p[0][1])]]
                    depth_b = depth[np.arange(len(xy))[(xy[:,0] == p[1][0]) & (xy[:,1] == p[1][1])]]
                    depth_c = depth[np.arange(len(xy))[(xy[:,0] == p[2][0]) & (xy[:,1] == p[2][1])]]
                    s = a + b + c
                    self.values[j,i] = depth_a*a/s + depth_b*b/s + depth_c*c/s
        return self.values

