import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append('../../../')
import deepy as dp2

def make_map_graduation2019(input_files:str,coastline_npy:np.ndarray,*,show=False,chart_title=' ',localmesh=dp2.LocalMesh()):
    lm = dp2.LocalMesh()
    files = np.array(input_files)
    files = files.reshape(len(files),1)
    for i in range(len(files[:])):
        d = dp2.read_deeper(files[i,0])
        e,n,u = dp2.latlonalt2enu(d.loc[:,'lat'],d.loc[:,'lon'],0,lm.origin[0],lm.origin[1],0)
        e,n = dp2.kalmanfilter_2d(e,n,d.loc[:,'t'])
        depth0 = d.loc[:,'depth']
        judge = dp2.judge_stability(d.loc[:,'t'])
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
    hc = dp2.Honeycomb(1,x=[lm.xrange[0],lm.xrange[1]],y=[lm.yrange[0],lm.yrange[1]])
    x,y,depth = hc.depth_weightedmean(x,y,depth)
    lm.values = lm.fill_grid(x,y,depth)
    lm.values = dp2.rm_landingvalues(lm,coastline_npy)
    if show:
        lm.values[lm.values >= 13.5] = 13.5
        plt.contourf(lm.lonlabel,lm.latlabel,lm.values,cmap='viridis_r',levels=np.arange(0,14,0.5))
        dp2.draw_coastline(coastline_npy)
        plt.xlim(lm.lonrange)
        plt.ylim(lm.latrange)
        plt.yticks(ticks=(35.62483,35.626),labels=['35°37.49\'N','35°37.56\'N'],fontsize=8)
        plt.xticks(ticks=(139.78333,139.7852),labels=['139°47.0\'E','139°47.11\'E'],fontsize=8)
        plt.title(chart_title,fontsize=17)
        cbar = plt.colorbar()
        cbar.set_label('depth [m]')
        plt.show()
    return lm 

def draw_ax(lm,coastline_npy,ax,*,chart_title='',latlabel=True,lonlabel=True):
    ax.contourf(lm.lonlabel,lm.latlabel,lm.values,cmap='viridis_r',levels=np.arange(0,14,0.5))
    dp2.draw_coastline(coastline_npy,fig=ax)
    ax.set_xlim(lm.lonrange)
    ax.set_ylim(lm.latrange)
    ax.set_yticks((35.62483,35.626))
    ax.set_xticks((139.78333,139.7852))
    if latlabel:
        ax.set_yticklabels(['35°37.49\'N','35°37.56\'N'],fontsize=6)
    else :
        ax.set_yticklabels([])
    if lonlabel:
        ax.set_xticklabels(['139°47.0\'E','139°47.11\'E'],fontsize=6)
    else :
        ax.set_xticklabels([])
    ax.set_title(chart_title,fontsize=10)

def make_deviationmap(lm,chart,coastline_npy,*,chart_title=' '):
    d=np.array(lm.values)-np.array(chart.mean(axis=2))
    plt.contourf(lm.lonlabel,lm.latlabel,d,cmap='bwr',levels=np.arange(-5,5.1,0.5))
    dp2.draw_coastline(coastline_npy)
    plt.xlim(lm.lonrange)
    plt.ylim(lm.latrange)
    plt.yticks(ticks=(35.62483,35.626),labels=['35°37.49\'N','35°37.56\'N'],fontsize=8)
    plt.xticks(ticks=(139.78333,139.7852),labels=['139°47.0\'E','139°47.11\'E'],fontsize=8)
    plt.title(chart_title,fontsize=10)
    cbar = plt.colorbar()
    cbar.set_label('deviation [m]')
    plt.show()

def draw_ax_deviation(lm,chart,coastline_npy,ax,*,chart_title=' ',latlabel=True,lonlabel=True):
    d=np.array(lm.values)-np.array(chart.mean(axis=2))
    ax.contourf(lm.lonlabel,lm.latlabel,d,cmap='bwr',levels=np.arange(-5,5.1,0.5))
    dp2.draw_coastline(coastline_npy,fig=ax)
    ax.set_xlim(lm.lonrange)
    ax.set_ylim(lm.latrange)
    ax.set_yticks((35.62483,35.626))
    ax.set_xticks((139.78333,139.7852))
    if latlabel:
        ax.set_yticklabels(['35°37.49\'N','35°37.56\'N'],fontsize=6)
    else :
        ax.set_yticklabels([])
    if lonlabel:
        ax.set_xticklabels(['139°47.0\'E','139°47.11\'E'],fontsize=6)
    else :
        ax.set_xticklabels([])
    ax.set_title(chart_title,fontsize=10)


