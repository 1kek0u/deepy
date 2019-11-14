import numpy as np
import pandas as pd

class Ellipsoid:
    def __init__(self,*,model:str='wgs84',a:float=0 , f:float=0 ):
        if model =='wgs84':
            self.a = 6378137
            self.f = 1 / 298.257223563
            self.b = self.a * (1-self.f)
            self.e = np.sqrt(2*self.f -self.f**2)
        if model =='grs80':
            self.a = 6378137
            self.f = 1 / 298.257222101
            self.b = self.a * (1-self.f)
            self.e = np.sqrt(2*self.f -self.f**2)
        if model == 'other':
            self.a = a 
            self.f = f
            self.b = self.a * (1-self.f)
            self.e = np.sqrt(2*self.f -self.f**2)


def latlon2ecef(lat,lon,h,*,ell=Ellipsoid()):
    lat = np.radians(lat)
    lon = np.radians(lon) 
    N = ell.a/np.sqrt(1-((ell.e**2)*(np.sin(lat)**2)))
    x = (N + h) * np.cos(lat)*np.cos(lon)
    y = (N + h) * np.cos(lat)*np.sin(lon)
    z = (N*(1-ell.e**2)+h)*np.sin(lat)
    return x,y,z

def ecef2enu(x,y,z,lat0:float,lon0:float,h0:float,*,ell=Ellipsoid()):
    lat0 = np.radians(lat0)
    lon0 = np.radians(lon0)
    N = ell.a/np.sqrt(1-((ell.e**2)*(np.sin(lat0)**2)))
    x0 = (N + h0)*np.cos(lat0)*np.cos(lon0)
    y0 = (N + h0)*np.cos(lat0)*np.sin(lon0)
    z0 = (N*(1-ell.e**2)+h0)*np.sin(lat0)
    t = (x-x0)*np.cos(lon0) + (y-y0)*np.sin(lon0)
    e = (y-y0)*np.cos(lon0) - (x-x0)*np.sin(lon0)
    n = (z-z0)*np.cos(lat0) - t*np.sin(lat0) 
    u = (z-z0)*np.sin(lat0) + t*np.cos(lat0)
    return e,n,u

def enu2ecef(e,n,u,lat0,lon0,h0,*,ell=Ellipsoid()):
    lat0 = np.radians(lat0)
    lon0 = np.radians(lon0)
    N = ell.a/np.sqrt(1-((ell.e**2)*(np.sin(lat0)**2)))
    x0 = (N + h0)*np.cos(lat0)*np.cos(lon0)
    y0 = (N + h0)*np.cos(lat0)*np.sin(lon0)
    z0 = (N*(1-ell.e**2)+h0)*np.sin(lat0)
    t = u*np.cos(lat0)-n*np.sin(lat0)
    x = t*np.cos(lon0) - e*np.sin(lon0) + x0
    y = t*np.sin(lon0) + e*np.cos(lon0) + y0
    z = u*np.sin(lat0) + n*np.cos(lat0) + z0
    return x,y,z

def ecef2latlon(x,y,z,*,ell=Ellipsoid()):
    lat = np.arctan(z/((1-ell.e**2)*np.sqrt(x**2+y**2)))
    lon = np.arctan(y/x)
    lat = np.degrees(lat)
    lon = np.degrees(lon)
    lon[lon < 0] = lon[lon < 0] + 180 
    lon[lon > 180] = lon[lon > 180] -180 
    return lat , lon

def latlon2enu(lat,lon,h,lat0,lon0,h0,*,ell=Ellipsoid()):
    Ell = ell
    x,y,z = latlon2ecef(lat,lon,h,ell=Ell)
    e,n,u = ecef2enu(x,y,z,lat0,lon0,h0,ell=Ell)
    return e,n,u

def enu2latlon(e,n,u,lat0,lon0,h0,*,ell=Ellipsoid()):
    Ell = ell
    x,y,z = enu2ecef(e,n,u,lat0,lon0,h0,ell=Ell)
    lat,lon = ecef2latlon(x,y,z,ell=Ell)
    return lat,lon

