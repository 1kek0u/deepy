import matplotlib.pyplot as plt
from read_deeper import read_deeper, deeper_csv2localdata, deeper_csv2LocalMap, delete_landing
from coastline import cstline_xml2npy , draw_coastline
from local_map import LocalMap

def show_map(input_deeperfiles,coastline_xmlfile,*,localmap=LocalMap()):
    lm = localmap
    coast = cstline_xml2npy(coastline_xmlfile)
    odaiba =  deeper_csv2LocalMap(input_deeperfiles,localmap=lm)
    odaiba = delete_landing(odaiba,coast)
    plt.contourf(odaiba.lonlabel,odaiba.latlabel,odaiba.values)
    cbar = plt.colorbar()
    cbar.set_label(' depth [m]')
    draw_coastline(coast)
    plt.xlim((odaiba.lonrange[0],odaiba.lonrange[1]))
    plt.ylim((odaiba.latrange[0],odaiba.latrange[1]))
    plt.xlabel('longitude [deg]')
    plt.ylabel('latitude [deg]')
    plt.show()


