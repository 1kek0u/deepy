import numpy as np
from scipy.spatial import Delaunay, delaunay_plot_2d, Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import pandas as pd

__all__ = ['read_deeper','kalmanfilter_2d','judge_stability','Honeycomb']

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

def kalmanfilter_2d(x,y,t,*,return_MSE:bool = False,
                    R = np.matrix([[0.420905,0,0,0],[0,1.24298,0,0],[0,0,0.420905,0],[0,0,0,1.24298]]),
                    Q = np.matrix([[0.21,0,0,0],[0,0.31,0,0],[0,0,0.21,0],[0,0,0,0.31]]),
                    P0 = np.matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])):
    G = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    I = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    H = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    xyMSE = 0
    P = P0
    f_x = np.zeros(len(x))
    f_y = np.zeros(len(y))
    f_x[0] = x[0]
    f_y[0] = y[0]
    Xhat = np.matrix([[f_x[0]],[f_y[0]],[0],[0]])
    for i in range(1,len(t)):
        dt =(t[i] - t[i-1])/1000
        if dt < 15 :
            F = np.matrix([[1,0,dt,0],[0,1,0,dt],[0,0,1,0],[0,0,0,1]])
            Xhat_k1 = np.dot(F,Xhat)
            P = np.dot(np.dot(F,P),F.T) + np.dot(np.dot(G,Q),G.T)
            Zvx = (x[i]-x[i-1])/dt
            Zvy = (y[i]-y[i-1])/dt
            Z = np.matrix([[x[i]],[y[i]],[Zvx],[Zvy]])
            e = Z - np.dot(H,Xhat_k1)
            S = R + np.dot(np.dot(H,P),H.T)
            K = np.dot(np.dot(P,H),np.linalg.inv(S))
            Xhat = Xhat_k1 + np.dot(K,e)
            P = np.dot((I-np.dot(K,H)),P)
            f_x[i] = Xhat[0,0]
            f_y[i] = Xhat[1,0]
            if return_MSE:
                xyMSE = xyMSE+((Xhat[0,0]-Xhat_k1[0,0])**2)+((Xhat[1,0]-Xhat_k1[1,0])**2)+((Xhat[0,0]-Z[0,0])**2)+((Xhat[1,0]-Z[1,0])**2)
        else:
            P = P0
            f_x[i] = x[i]
            f_y[i] = y[i]
            Xhat = np.matrix([[f_x[i]],[f_y[i]],[0],[0]])
    if return_MSE:
        xyMSE = xyMSE / len(t)
        return f_x,f_y,xyMSE
    else:
        return f_x,f_y

def judge_stability(t,*,casting_interval=10,time2stable=15):
    """
    Parameters
    ----------
        t: numpy.array 
        casting_interval:float . 
        time2stable: float

    Return
    ------
        bool. 
            False : not stable .
            True  : stable.
    """
    judge = np.array([False])
    i = 1
    while True:
        dt = (t[i] - t[0])/1000
        if dt >= time2stable:
            break
        judge = np.hstack((judge,[False]))
        i = i + 1
    while i < len(t):
        dt = (t[i]-t[i-1])/1000
        if dt > casting_interval:
            judge = np.hstack((judge,[False]))
            t_l = t[i]
            while True:
                i = i +1
                if i >= len(t):
                    break
                dt = (t[i] - t_l)/1000
                if dt >= time2stable:
                    break
                judge = np.hstack((judge,[False]))
        else:
            judge = np.hstack((judge,[True]))
            i = i+1
    return judge

class Honeycomb:

    def __init__(self,r,*,x=[0,0],y=[0,0]):
        self.r = r
        x1 = np.arange(x[0],x[1]+3*r*np.sqrt(3)/2,r*np.sqrt(3))
        y1 = np.arange(y[0],y[1]+3*r*3/2,r*3)
        x2 = np.arange(x[0]-r*np.sqrt(3)/2,x[1]+3*r*np.sqrt(3)/2,r*np.sqrt(3))
        y2 = np.arange(y[0]-r*3/2,y[1]+3*r*3/2,r*3)
        x1x1,y1y1 = np.meshgrid(x1,y1)
        x2x2,y2y2 = np.meshgrid(x2,y2)
        p1 = np.c_[x1x1.ravel(),y1y1.ravel()]
        p2 = np.c_[x2x2.ravel(),y2y2.ravel()]
        p = np.vstack((p1,p2))
        self.Voronoi = Voronoi(p)
        self.xrange = [x[0],x[1]+r*np.sqrt(3)/2]
        self.yrange = [y[0],y[1]+r/2]
        p = p[p[:,0] >= (x[0] - r/10000)]
        p = p[p[:,0] <= (x[1] + r*np.sqrt(3)/2 + r/10000)]
        p = p[p[:,1] >= (y[0] - r/10000)]
        p = p[p[:,1] <= (y[1] + r/2 + r/10000)]
        self.voronoi_points = p
        ver = self.Voronoi.vertices
        ver = ver[ver[:,0] >= (x[0] - r/10000)]
        ver = ver[ver[:,0] <= (x[1] + r*np.sqrt(3)/2 + r/10000)]
        ver = ver[ver[:,1] >= (y[0] - r/2 - r/10000)]
        ver = ver[ver[:,1] <= (y[1] + r/2 + r/10000)]
        self.delaunay_points = np.vstack((p,ver))
        self.Delaunay = Delaunay(self.delaunay_points)

    def hexagon_plot(self,*,tri=False):
        fig,ax = plt.subplots(figsize=(4,4))
        if tri :
            delaunay_plot_2d(self.Delaunay,ax)
        voronoi_plot_2d(self.Voronoi,ax)
        plt.xlim(self.xrange)
        plt.ylim(self.yrange)
        return fig     

    def point_plot(self):
        fig,ax = plt.subplots(figsize=(4,4))
        ax.scatter(self.delaunay_points[:,0],self.delaunay_points[:,1])
        ax.scatter(self.voronoi_points[:,0],self.voronoi_points[:,1],marker='x')
        plt.xlim(self.xrange)
        plt.ylim(self.yrange)
        return fig

    def find_hexagon(self,x,y,*,return_coordinates=True):
        x = np.array(x)
        y = np.array(y)       
        x = x.reshape(len(x),1)
        y = y.reshape(len(y),1)
        xy = np.hstack((x,y))
        if return_coordinates:
            p = np.array([[np.nan,np.nan]])
        else:
            p = np.array([[np.nan]])
        for i in range(len(xy)):
            sim = self.Delaunay.find_simplex(xy[i])
            if sim != -1:
                sim = self.delaunay_points[self.Delaunay.simplices[sim]]
                for j in range(4):
                    n = np.arange(len(self.voronoi_points))[(abs(self.voronoi_points[:,0] - sim[j][0]) < self.r/10000) & (abs(self.voronoi_points[:,1] - sim[j][1]) < self.r/10000)]
                    if len(n) == 1:
                        if return_coordinates :
                            p = np.vstack((p,self.voronoi_points[n,:]))
                            break
                        else:
                            p = np.vstack((p,n))
                            break
            else:
                if return_coordinates:
                    p = np.vstack((p,[[np.nan,np.nan]]))
                else:
                    p = np.vstack((p,[[np.nan]]))
        p = p[1:]
        return p

    def depth_sum(self,x,y,depth,*,depth_sum=0,N=0):
        depth = np.array(depth)
        x = np.array(x)
        y = np.array(y)       
        x = x.reshape(len(x),1)
        y = y.reshape(len(y),1)
        xy = np.hstack((x,y))
        if N == 0:
            N = np.zeros([len(self.voronoi_points)])
        if depth_sum == 0:
            depth_sum = np.zeros([len(self.voronoi_points)])
        for i in range(len(xy)):
            sim = self.Delaunay.find_simplex(xy[i])
            if sim != -1:
                sim = self.delaunay_points[self.Delaunay.simplices[sim]]
                for j in range(4):
                    n = np.arange(len(self.voronoi_points))[(abs(self.voronoi_points[:,0] - sim[j][0]) < self.r/10000) & (abs(self.voronoi_points[:,1] - sim[j][1]) < self.r/10000)]
                    if len(n) == 1:
                        N[n] = N[n] + 1
                        depth_sum[n] = depth_sum[n] + depth[i]
                        break
        return depth_sum,N

    def depth_weightedmean(self,x,y,depth,*,return_coordinates=True):
        p = self.find_hexagon(x,y) 
        v_p = self.voronoi_points
        depth = np.array(depth)
        x = np.array(x)
        y = np.array(y)       
        x = x.reshape(len(x),1)
        y = y.reshape(len(y),1)
        xy = np.hstack((x,y))
        depth_v = np.array([np.nan])
        for i in range(len(v_p)):
            n = np.arange(len(p))[(abs(p[:,0] -  v_p[i,0]) <  self.r/10000) & (abs(p[:,1] - v_p[i,1]) < self.r/10000)]
            if len(n) == 1:
                depth_v = np.hstack((depth_v,depth[n]))
            elif len(n) > 1:
                dep = 0
                pxy = xy[n]
                pp = p[n]
                pdepth = depth[n]
                dist = np.sqrt((pxy[:,0]-pp[:,0])**2 + (pxy[:,1]-pp[:,1])**2) + self.r/10
                dist = 1/dist
                for j in range(len(n)):
                    dep = dep + pdepth[j]*dist[j]/dist.sum() 
                depth_v = np.hstack((depth_v,[dep]))
            else:
                depth_v = np.hstack((depth_v,[np.nan]))
        depth_v = depth_v[1:]
        if return_coordinates:
            vp_x = self.voronoi_points[:,0]
            vp_y = self.voronoi_points[:,1]
            judge = ~np.isnan(depth_v)
            vp_x = vp_x[judge]
            vp_y = vp_y[judge]
            depth_v = depth_v[judge]
            return vp_x,vp_y,depth_v
        else:
            return depth_v
        

