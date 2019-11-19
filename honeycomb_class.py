from scipy.spatial import Delaunay, delaunay_plot_2d, Voronoi, voronoi_plot_2d
import numpy as np
import matplotlib.pyplot as plt

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
        self.voronoi_point = p
        ver = self.Voronoi.vertices
        ver = ver[ver[:,0] >= (x[0] - r/10000)]
        ver = ver[ver[:,0] <= (x[1] + r*np.sqrt(3)/2 + r/10000)]
        ver = ver[ver[:,1] >= (y[0] - r/2 - r/10000)]
        ver = ver[ver[:,1] <= (y[1] + r/2 + r/10000)]
        self.delaunay_point = np.vstack((p,ver))
        self.Delaunay = Delaunay(self.delaunay_point)

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
        ax.scatter(self.delaunay_point[:,0],self.delaunay_point[:,1])
        ax.scatter(self.voronoi_point[:,0],self.voronoi_point[:,1],marker='x')
        plt.xlim(self.xrange)
        plt.ylim(self.yrange)
        return fig

    def find_hexagon(self,x,y,*,coordinate=True):
        x = np.array(x)
        y = np.array(y)       
        x = x.reshape(len(x),1)
        y = y.reshape(len(y),1)
        xy = np.hstack((x,y))
        if coordinate:
            v_p = np.array([[np.nan,np.nan]])
        else:
            v_p = np.array([[np.nan]])
        for i in range(len(xy)):
            sim = self.Delaunay.find_simplex(xy[i])
            if sim != -1:
                sim = self.delaunay_point[self.Delaunay.simplices[sim]]
                for j in range(4):
                    n = np.arange(len(self.voronoi_point))[(abs(self.voronoi_point[:,0] - sim[j][0]) < self.r/10000) & (abs(self.voronoi_point[:,1] - sim[j][1]) < self.r/10000)]
                    if len(n) == 1:
                        if coordinate :
                            v_p = np.vstack((v_p,self.voronoi_point[n,:]))
                            break
                        else:
                            v_p = np.vstack((v_p,n))
                            break
            else:
                if coordinate:
                    v_p = np.vstack((v_p,[[np.nan,np.nan]]))
                else:
                    v_p = np.vstack((v_p,[[np.nan]]))
        v_p = v_p[1:,:]
        return v_p

    def fill_in_cell(self,x,y,depth,*,N=0,depth_sum=0):
        x = np.array(x)
        y = np.array(y)       
        x = x.reshape(len(x),1)
        y = y.reshape(len(y),1)
        xy = np.hstack((x,y))
        if N == 0:
            N = np.zeros([len(self.voronoi_point)])
        if depth_sum == 0:
            depth_sum = np.zeros([len(voronoi_point)])
        for i in range(len(xy)):
            sim = self.Delaunay.find_simplex(xy[i])
            if sim != -1:
                sim = self.delaunay_point[self.Delaunay.simplices[sim]]
                for j in range(4):
                    n = np.arange(len(self.voronoi_point))[(abs(self.voronoi_point[:,0] - sim[j][0]) < self.r/10000) & (abs(self.voronoi_point[:,1] - sim[j][1]) < self.r/10000)]
                    if len(n) == 1:
                        N[n] = N[n] + 1
                        depth_sum[n] == depth_sum[n] + depth[i]
        return depth_sum,N


