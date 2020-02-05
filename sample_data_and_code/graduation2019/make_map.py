import numpy as np
import sys , os
import matplotlib.pyplot as plt
from helper import make_map_graduation2019, draw_ax, make_deviationmap,draw_ax_deviation  
sys.path.append('../../../')
import pdb
import deepy as dp2
pdb.set_trace()

coastline = dp2.cstline_xml2npy('FG-GML-533936-Cstline-20181001-0001.xml')
lm0 = make_map_graduation2019(['0806-1st0f60.csv','0806-1st25d2.csv'],coastline,show=True,chart_title='2019/8/6 1st')
lm1 = make_map_graduation2019(['0806-2nd_0f60.csv','0806-2nd_25D2.csv'],coastline,show=True,chart_title='2019/8/6 2nd')
lm2 = make_map_graduation2019(['0821map_25d2.csv'],coastline,show=True,chart_title='2019/8/21')
lm3 = make_map_graduation2019(['0906_barcheck1_map_25D2.csv'],coastline,show=True,chart_title='2019/9/6')
lm4 = make_map_graduation2019(['0910_1st25d2.csv'],coastline,show=True,chart_title='2019/9/10 1st')
lm5 = make_map_graduation2019(['0910_2nd25d2.csv'],coastline,show=True,chart_title='2019/9/10 2nd')
lm6 = make_map_graduation2019(['1106-25D2.csv'],coastline,show=True,chart_title='2019/11/6')

fig,ax = plt.subplots(4,2,figsize=(8,12))
draw_ax(lm0,coastline,ax[0,0],chart_title='2019/8/6 1st')
draw_ax(lm1,coastline,ax[0,1],chart_title='2019/8/6 2nd',latlabel=False)
draw_ax(lm2,coastline,ax[1,0],chart_title='2019/8/21')
draw_ax(lm3,coastline,ax[1,1],chart_title='2019/9/6',latlabel=False)
draw_ax(lm4,coastline,ax[2,0],chart_title='2019/9/10 1st')
draw_ax(lm5,coastline,ax[2,1],chart_title='2019/9/10 2nd',latlabel=False)
draw_ax(lm6,coastline,ax[3,0],chart_title='2019/11/6')

lm = dp2.LocalMesh()

fig.subplots_adjust(wspace=0.05,hspace=0.20,top=0.97,bottom=0.02,left=0.12,)
cbar = fig.colorbar(plt.contourf(lm0.values,cmap='viridis_r',levels=np.arange(0,14,0.5)),ax=ax.ravel().tolist(),pad=0.02,fraction=0.05)
print('ok')
cbar.set_label('depth [m]',fontsize=8)
ax[3,1].axis('off')
plt.show()

chart = np.dstack((lm0.values,lm1.values,lm2.values,lm3.values,lm4.values,lm5.values,lm6.values))
plt.contourf(lm.lonlabel,lm.latlabel,chart.mean(axis=2),cmap='viridis_r',levels=np.arange(0,14,0.5))
dp2.draw_coastline(coastline)
plt.xlim(lm.lonrange)
plt.ylim(lm.latrange)
plt.yticks(ticks=(35.62483,35.626),labels=['35°37.49\'N','35°37.56\'N'],fontsize=8)
plt.xticks(ticks=(139.78333,139.7852),labels=['139°47.0\'E','139°47.11\'E'],fontsize=8)
plt.title('mean',fontsize=10)
cbar = plt.colorbar()
cbar.set_label('depth [m]')
plt.show()

plt.contourf(lm.lonlabel,lm.latlabel,chart.var(axis=2,ddof=1),cmap='viridis_r',levels=np.arange(0,14,1))
dp2.draw_coastline(coastline)
plt.xlim(lm.lonrange)
plt.ylim(lm.latrange)
plt.yticks(ticks=(35.62483,35.626),labels=['35°37.49\'N','35°37.56\'N'],fontsize=8)
plt.xticks(ticks=(139.78333,139.7852),labels=['139°47.0\'E','139°47.11\'E'],fontsize=8)
plt.title('variance',fontsize=10)
cbar = plt.colorbar()
cbar.set_label('variance [m^2]')
plt.show()


make_deviationmap(lm0,chart,coastline,chart_title='Deviation  2019/8/6 1st')
make_deviationmap(lm1,chart,coastline,chart_title='Deviation  2019/8/6 2nd')
make_deviationmap(lm2,chart,coastline,chart_title='Deviation  2019/8/21')
make_deviationmap(lm3,chart,coastline,chart_title='Deviation  2019/9/6')
make_deviationmap(lm4,chart,coastline,chart_title='Deviation  2019/9/10 1st')
make_deviationmap(lm5,chart,coastline,chart_title='Deviation  2019/9/10 2nd')
make_deviationmap(lm6,chart,coastline,chart_title='Deviation  2019/11/6')
   
fig,ax = plt.subplots(4,2,figsize=(8,12))
draw_ax_deviation(lm0,chart,coastline,ax[0,0],chart_title='Deviation  2019/8/6 1st')
draw_ax_deviation(lm1,chart,coastline,ax[0,1],chart_title='Deviation  2019/8/6 2nd',latlabel=False)
draw_ax_deviation(lm2,chart,coastline,ax[1,0],chart_title='Deviation  2019/8/21')
draw_ax_deviation(lm3,chart,coastline,ax[1,1],chart_title='Deviation  2019/9/6',latlabel=False)
draw_ax_deviation(lm4,chart,coastline,ax[2,0],chart_title='Deviation  2019/9/10 1st')
draw_ax_deviation(lm5,chart,coastline,ax[2,1],chart_title='Deviation  2019/9/10 2nd',latlabel=False)
draw_ax_deviation(lm6,chart,coastline,ax[3,0],chart_title='Deviation  2019/11/6')

lm = dp2.LocalMesh()
d = np.array(lm0.values) - np.array(chart.mean(axis=2))

fig.subplots_adjust(wspace=0.05,hspace=0.20,top=0.97,bottom=0.02,left=0.12,)
cbar = fig.colorbar(plt.contourf(d,cmap='bwr',levels=np.arange(-5,5.1,0.5)),ax=ax.ravel().tolist(),pad=0.02,fraction=0.05)
print('ok')
cbar.set_label('deviation [m]',fontsize=8)
ax[3,1].axis('off')
plt.show()

