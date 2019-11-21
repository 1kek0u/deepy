import numpy as np

__all__ = ['kalmanfilter_2d','judge_stability']

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


def kalmanfilter_3d(x,y,z,t,*,return_MSE:bool = False,
                    R = np.matrix([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]),
                    Q = np.matrix([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]),
                    P0 = np.matrix([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])):
    G = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
    I = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
    H = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
    xyzMSE = 0
    P = P0
    f_x = np.zeros(len(x))
    f_y = np.zeros(len(y))
    f_z = np.zeros(len(z))
    f_x[0] = x[0]
    f_y[0] = y[0]
    f_z[0] = z[0]
    Xhat = np.matrix([[f_x[0]],[f_y[0]],[f_z[0]],[0],[0],[0]])
    for i in range(1,len(d)):
        dt =(t[i] - t[i-1])/1000
        if dt < 15 :
            F = np.matrix([[1,0,0,dt,0,0],[0,1,0,0,dt,0],[0,0,1,0,0,dt],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
            Xhat_k1 = np.dot(F,Xhat)
            P = np.dot(np.dot(F,P),F.T) + np.dot(np.dot(G,Q),G.T)
            Zvx = (x[i]-x[i-1])/dt
            Zvy = (y[i]-y[i-1])/dt
            Zvz = (z[i]-z[i-1])/dt
            Z = np.matrix([[x[i]],[y[i]],[z[i]],[Zvx],[Zvy],[Zvz]])
            e = Z - np.dot(H,Xhat_k1)
            S = R + np.dot(np.dot(H,P),H.T)
            K = np.dot(np.dot(P,H),np.linalg.inv(S))
            Xhat = Xhat_k1 + np.dot(K,e)
            P = np.dot((I-np.dot(K,H)),P)
            f_x[i] = Xhat[0,0]
            f_x[i] = Xhat[1,0]
            f_z[i] = Xhat[2,0]
            if return_MSE:
                xyzMSE = xyzMSE+(((Xhat[0]-Xhat_k1[0])**2)+((Xhat[1]-Xhat_k1[1])**2)+((Xhat[2]-Xhat_k1[2])**2)+
                                 ((Xhat[0]-Z[0])**2)+((Xhat[1]-Z[1])**2)+((Xhat[2]-Z[2])**2))
        else:
            P = P0
            f_x[i] = x[i]
            f_y[i] = y[i]
            f_z[i] = z[i]
            Xhat = np.matrix([[f_x[i]],[f_y[i]],[f_z[i]],[0],[0],[0]])
    if return_MSE:
        xyzMSE = xyzMSE / len(d)
        return f_x,f_y,f_z,xyzMSE
    else:
        return f_x,f_y,f_z

def judge_stability(t,*,casting_interval=15,time2stable=5):
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

