import numpy as np
import matplotlib.pyplot as plt
import codecs
import os

__all__ = ['cstline_xml2csv','cstline_xml2npy','draw_coastline']

def cstline_xml2csv(input_xml:str):
    """
    Convert Japanese coast-line imformation XML file into numerical CSV file devided by blocks.
    Japanese coast-line imformation XML file is avalable in here <http://>


    Parameters
    ----------
        input_xml:str
            Name of XML File that has coast-line imformations to convert into Numerical data that formatted CSV .
   
    """
    with codecs.open(input_xml,'r','utf-8','ignore')as f:
        s = f.read()
    dir_name = input_xml[:input_xml.find('.')]
    os.mkdir(dir_name)
    os.chdir(dir_name)
    P_beg = '<gml:posList>\r\n'  
    P_end = '</gml:posList>'
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
    
    Parameters
    ----------
        input_xml:str
            Name of CoastLine xmlfile.
    
    Return
    ------ 
        cstline_npy:numpy.ndarray 
            shape of (index,2,number of bloks) 
            you can acesses [['lat','lon]]' data in block 'n' , with typing 'cstline_npy[:,:,n] '
    """
    with codecs.open(input_xml,'r','utf-8','ignore')as f:
        s = f.read()
    P_beg = '<gml:posList>\r\n'  
    P_end = '</gml:posList>'
    cst = np.array([[np.nan,np.nan]])
    m = len(P_beg)
    i = 0
    n = 0
    coastlines = 0
    while True:
        i = s.find(P_beg,i)
        if i == -1:
            break
        j = s.find(P_end,i)
        pk = i + m
        while True:
            k = s.find(' ',pk,j)
            if k == -1:
                if type(coastlines) == int:
                    coastlines = cst[1:,:]
                    cst = np.array([np.nan,np.nan])
                    break
                cst = cst[1:,:]
                if len(cst) >= len(coastlines):
                    if np.ndim(coastlines) == 2:
                        zero = np.zeros([len(cst)-len(coastlines),2])
                        zero[:,:] = np.nan
                        coastlines = np.vstack((coastlines,zero))
                        coastlines = np.dstack((coastlines,cst))
                    else :
                        zero = np.zeros([len(cst)-len(coastlines),2,len(coastlines[0,0,:])])
                        zero[:,:,:] = np.nan
                        coastlines = np.vstack((coastlines,zero))
                        coastlines = np.dstack((coastlines,cst))
                if len(cst) < len(coastlines):
                    zero = np.zeros([len(coastlines)-len(cst),2])
                    zero[:,:] = np.nan
                    cst = np.vstack((cst,zero))
                    coastlines = np.dstack((coastlines,cst))
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
    """

def draw_coastline(coastlines,color='k',fig=np.nan):
    c = color
    for i in range (len(coastlines[0,0,:])):
        lat = coastlines[:,0,i]
        lon = coastlines[:,1,i]
        lat = lat[lat !=0]
        lon = lon[lon !=0]
        if type(fig) == float :
            plt.plot(lon,lat,'-',color=c)
        else :
            fig.plot(lon,lat,'-',color=c)
