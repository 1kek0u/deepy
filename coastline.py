import numpy as np
import matplotlib.pyplot as plt
import codecs
import os

__all__ = ['cstline_xml2csv','cstline_xml2npy','draw_coastline']

def cstline_xml2csv(input_xml:str):
    """
    input : name of Coast Line xml file 
    output: None
    
    make Coast Line csvfiles devided by blocs
    
    """
    with codecs.open(input_xml,'r','utf-8','ignore')as f:
        s = f.read()
    dir_name = input_xml[:input_xml.find('.')-1]
    os.mkdir(dir_name)
    os.chdir(dir_name)
    P_beg = '<gml:posList>\r\n'  
    P_end = '</gml:posList>'
    coastlines = np.zeros([500,2,s.count(P_beg)])
    m = len(P_beg)
    i = 0
    n = 0
    while True:
        i = s.find(P_beg,i)
        if i == -1:
            break
        j = s.find(P_end,i)
        file_name = 'CoastLine'
        file_name = file_name + str(n) + '.csv' 
        fout = open(file_name,'w')
        fout.writelines(s[i+m : j])
        fout.close()
        n = n+1
        i = j + len(P_end) 


def cstline_xml2npy(input_xml:str):
    """
    input : name of CoastLine xmlfile
    output: numpy.ndarray[index,latlon,blocks]
    
    """
    with codecs.open(input_xml,'r','utf-8','ignore')as f:
        s = f.read()
    P_beg = '<gml:posList>\r\n'  
    P_end = '</gml:posList>'
    coastlines = np.zeros([500,2,s.count(P_beg)])
    m = len(P_beg)
    i = 0
    n = 0
    while True:
        i = s.find(P_beg,i)
        if i == -1:
            break
        j = s.find(P_end,i)
        pk = i + m
        cst = np.array([np.nan,np.nan])
        while True:
            k = s.find(' ',pk,j)
            if k == -1:
                cst = cst[1:,:]
                zero = np.empty((500 - len(cst),2))
                zero[:,:] = np.nan
                cst = np.vstack((cst,zero))
                coastlines[:,:,n] = cst
                n = n+1
                cst = np.array([np.nan,np.nan])
                i = pk 
                break
            lat = float(s[pk:k-1])
            pk = k + 1
            k = s.find('\n',pk,j)
            lon = float(s[pk:k-1])
            cst = np.vstack((cst,np.array([lat,lon])))
            pk = k + 1
    return coastlines 

def draw_coastline(coastlines,*,color='k'):
    c = color
    for i in range (len(coastlines[0,0,:])):
        lat = coastlines[:,0,i]
        lon = coastlines[:,1,i]
        lat = lat[lat !=0]
        lon = lon[lon !=0]
        plt.plot(lon,lat,'-',color=c)


