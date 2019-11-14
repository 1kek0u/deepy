import numpy as np
import pandas as pd

def kalmanfilter_2d(d:pd.core.frame.DataFrame,*,x:str='x',y:str='y', return_MSE:bool = False,
                    R = np.matrix([[0.420905,0,0,0],[0,1.24298,0,0],[0,0,0.420905,0],[0,0,0,1.24298]]),
                    Q = np.matrix([[0.21,0,0,0],[0,0.31,0,0],[0,0,0.21,0],[0,0,0,0.31]]),
                    P0 = np.matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])):
    G = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    I = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    H = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    xyMSE = 0
    P = P0
    f_d = pd.DataFrame(index = d.index , columns = d.columns)
    for j in range(len(d.columns)):
        if d.columns[j] != x and d.columns[j] != y:
            f_d.iloc[:,j] = d.iloc[:,j]
    f_d.at[0,x] = d.at[0,x]
    f_d.at[0,y] = d.at[0,y]
    Xhat = np.matrix([[f_d.at[0,x]],[f_d.at[0,y]],[0],[0]])
    for i in range(1,len(d)):
        dt =(d.at[i,'t'] - d.at[i-1,'t'])/1000
        if dt < 15 :
            F = np.matrix([[1,0,dt,0],[0,1,0,dt],[0,0,1,0],[0,0,0,1]])
            Xhat_k1 = np.dot(F,Xhat)
            P = np.dot(np.dot(F,P),F.T) + np.dot(np.dot(G,Q),G.T)
            Zvx = (d.at[i,x]-d.at[i-1,x])/dt
            Zvy = (d.at[i,y]-d.at[i-1,y])/dt
            Z = np.matrix([[d.at[i,x]],[d.at[i,y]],[Zvx],[Zvy]])
            e = Z - np.dot(H,Xhat_k1)
            S = R + np.dot(np.dot(H,P),H.T)
            K = np.dot(np.dot(P,H),np.linalg.inv(S))
            Xhat = Xhat_k1 + np.dot(K,e)
            P = np.dot((I-np.dot(K,H)),P)
            f_d.at[i,x] = Xhat[0,0]
            f_d.at[i,y] = Xhat[1,0]
            if return_MSE:
                xyMSE = xyMSE+((Xhat[0,0]-Xhat_k1[0,0])**2)+((Xhat[1,0]-Xhat_k1[1,0])**2)+((Xhat[0,0]-Z[0,0])**2)+((Xhat[1,0]-Z[1,0])**2)
        else:
            P = P0
            f_d.at[i,x] = d.at[i,x]
            f_d.at[i,y] = d.at[i,y]
            Xhat = np.matrix([[f_d.at[i,x]],[f_d.at[i,y]],[0],[0]])
    if return_MSE:
        xyMSE = xyMSE / len(d)
        return f_d,xyMSE
    else:
        return f_d


def kalmanfilter_3d(d:pd.core.frame.DataFrame,*,x:str='x',y:str='y',z:str='z', return_MSE:bool = False,
                    R = np.matrix([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]),
                    Q = np.matrix([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]),
                    P0 = np.matrix([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])):
    G = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
    I = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
    H = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
    xyzMSE = 0
    P = P0
    f_d = pd.DataFrame(index = d.index , columns = d.columns)
    for j in range(len(d.columns)):
        if d.columns[j] != x and d.columns[j] != y and d.columns[j] != z:
            f_d.iloc[:,j] = d.iloc[:,j]
    f_d.at[0,x] = d.at[0,x]
    f_d.at[0,y] = d.at[0,y]
    f_d.at[0,z] = d.at[0,z]
    Xhat = np.matrix([[f_d.at[0,x]],[f_d.at[0,y]],[f_d.at[0,z]],[0],[0],[0]])
    for i in range(1,len(d)):
        dt =(d.at[i,'t'] - d.at[i-1,'t'])/1000
        if dt < 15 :
            F = np.matrix([[1,0,0,dt,0,0],[0,1,0,0,dt,0],[0,0,1,0,0,dt],[0,0,0,0,1,0],[0,0,0,0,0,1]])
            Xhat_k1 = np.dot(F,Xhat)
            P = np.dot(np.dot(F,P),F.T) + np.dot(np.dot(G,Q),G.T)
            Zvx = (d.at[i,x]-d.at[i-1,x])/dt
            Zvy = (d.at[i,y]-d.at[i-1,y])/dt
            Zvz = (d.at[i,z]-d.at[i-1,z])/dt
            Z = np.matrix([[d.at[i,x]],[d.at[i,y]],[d.at[i,z]],[Zvx],[Zvy],[Zvz]])
            e = Z - np.dot(H,Xhat_k1)
            S = R + np.dot(np.dot(H,P),H.T)
            K = np.dot(np.dot(P,H),np.linalg.inv(S))
            Xhat = Xhat_k1 + np.dot(K,e)
            P = np.dot((I-np.dot(K,H)),P)
            f_d.at[i,x] = Xhat[0,0]
            f_d.at[i,y] = Xhat[1,0]
            f_d.at[i,z] = Xhat[2,0]
            if return_MSE:
                xyzMSE = xyzMSE+(((Xhat[0]-Xhat_k1[0])**2)+((Xhat[1]-Xhat_k1[1])**2)+((Xhat[2]-Xhat_k1[2])**2)+
                                 ((Xhat[0]-Z[0])**2)+((Xhat[1]-Z[1])**2)+((Xhat[2]-Z[2])**2))
        else:
            P = P0
            f_d.at[i,x] = d.at[i,x]
            f_d.at[i,y] = d.at[i,y]
            f_d.at[i,z] = d.at[i,z]
            Xhat = np.matrix([[f_d.at[i,x]],[f_d.at[i,y]],[f_d.at[i,z]],[0],[0],[0]])
    if return_MSE:
        xyMSE = xyzMSE / len(d)
        return f_d,xyzMSE
    else:
        return f_d
