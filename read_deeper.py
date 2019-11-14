import numpy as np
import pandas as pd
import codecs

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

def xmlcstline2npy(input_xml:str):
    with codecs.open(input_xml,'r','utf-8','ignore') as f:
        s = f.read()
    check = 'osList>\r\n'
    data = np.zeros([1000,2,s.count(check)])
    i=0
    ok = 0
    for i in range(len(s)):
        if s[i] == 'p':
            b = 0
            for n in range(len(check)):
                if s[i+1+n] == check[n]:
                    b = b + 1
                if b == len(check):
                    a = ''
                    ok = ok + 1
                    i = i + n + 1
                    while s[i] != '<':
                        if s[i] != '\r' and s[i] != '\n':
                            a = a + s[i]
                        if s[i] == '\r':
                            a = a + ' '
                            i = i +1
                    if len(a) > 10 :
                        fout = open("text.txt","w")
                        fout.writelines(a)
                        fout.close()
                        gml = pd.read_table('text.txt',delimiter=' ',dtype='f8',header=None )
                        gml = np.array(gml)
                        gml = gml[0,:-1]
                        gml = gml.reshape(int(len(gml)/2),2)
                        zero = np.zeros([1000-len(gml),2])
                        gml = np.vstack((gml,zero))
                        data[:,:,ok] = gml
    return data

