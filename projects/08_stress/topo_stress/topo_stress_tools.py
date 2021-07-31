import numpy as np
from scipy.fft import fft2,ifft2,fftshift,ifftshift,fftfreq
import pygmt
import netCDF4 as nc
from math import pi,exp

def read_GMT_netcdf(filename):
    f = nc.Dataset(filename,'r',format='NETCDF4')
    lon = f.variables['lon']
    lat = f.variables['lat']
    data = f.variables['z']
    my_lon = lon[:].copy()
    my_lat = lat[:].copy()
    my_data = data[:,:].copy()
    return my_lon, my_lat, my_data

def write_GMT_netcdf(filename,lons,lats,data):
    nlo = len(lons)
    nla = len(lats)
    f = nc.Dataset(filename,'w',format='NETCDF4')
    f.createDimension('lon',nlo)
    f.createDimension('lat',nla)
    lon = f.createVariable('lon','f8',('lon',))
    lon.long_name = 'longitude'
    lon.units = 'degrees_east'
    lon.standard_name = 'longitude'
    lon.axis = 'X'
    lon.actual_range = [lons.min(), lons.max()]
    lat = f.createVariable('lat','f8',('lat',))
    lat.long_name = 'latitude'
    lat.units = 'degrees_north'
    lat.standard_name = 'latitude'
    lat.axis = 'Y'
    lat.actual_range = [lats.min(), lats.max()]
    z = f.createVariable('z','f4',('lat','lon'))
    z.long_name = 'z'
    z.actual_range = [data.min(), data.max()]
    lat[:] = lats
    lon[:] = lons
    z[:,:] = data.astype('float32')
    f.close()
    return

def compute_window(grid):
    ni = grid.shape[0]
    nj = grid.shape[1]
    nsigy = ni/8
    nsigx = nj/8
    yidx = np.linspace(1,ni,ni)
    xidx = np.linspace(1,nj,nj)
    ywind = np.ones_like(yidx)
    xwind = np.ones_like(xidx)
    ywind[yidx < nsigy] = 0.5*(1-np.cos(pi*(yidx[yidx < nsigy] - 1)/nsigy))
    ywind[yidx > (ni-nsigy)] = 0.5*(1-np.cos(pi*(ni - yidx[yidx > (ni-nsigy)])/nsigy))
    xwind[xidx < nsigx] = 0.5*(1-np.cos(pi*(xidx[xidx < nsigx] - 1)/nsigx))
    xwind[xidx > (nj-nsigx)] = 0.5*(1-np.cos(pi*(nj - xidx[xidx > (nj-nsigx)])/nsigx))
    window = np.outer(ywind,xwind)
    return window

def fft_grd(to_fft, inv=False):
    if inv==True:
        out = ifft2(to_fft)
    else:
        out = fft2(to_fft)
    return out

def moho_topo(kx,ky,Te,rhoc,rhom,fk,young,rnu):
    grv = 9.81
    D = young*Te**3/12/(1-rnu**2)
    beta = 2*pi*np.sqrt(kx**2 + ky**2)
    transfer = 1/(1 + D*beta**4/grv/(rhom-rhoc)) 
    gk = fk*transfer
    return gk

def airy2load(kx,ky,zobs,H,fk,gk,rlam1,rmu1):
    
    beta = 2*pi*np.sqrt(kx**2 + ky**2)
    alpha = (rlam1+rmu1)/(rlam1+2*rmu1)
    eta = (rlam1+rmu1)/(3*rlam1 + 4*rmu1)
    bH = beta*H
    epH = np.exp(bH)
    ep2H = np.exp(2*bH)
    ep3H = np.exp(3*bH)
    ep4H = np.exp(4*bH)
    epH[bH < -50] = 0
    ep2H[bH < -50] = 0
    ep3H[bH < -50] = 0
    ep4H[bH < -50] = 0
    
    
    denom = beta**3*((eta+alpha*(-1+alpha+eta))**2 +
                    ep4H*(eta+alpha*(-1+alpha+eta))**2 -
                    2*ep2H*((eta+alpha*(-1+alpha+eta))**2 +
                    2*alpha**2*beta**2*(alpha-eta)**2*H**2))*rmu1
    
    A1 = gk*epH*eta*(eta-alpha*(1+eta+alpha*(-3+2*alpha+2*eta))
                    + alpha*beta*(alpha*(-1+3*alpha-5*eta)+eta)*H
                    + ep2H*(eta+alpha*(-1+alpha+eta))*(-1+alpha*(2+bH))) \
        + fk*eta*(-eta+alpha*(1+eta+alpha*(-3+2*alpha+2*eta))
                    + ep2H*(eta-alpha*(1+eta+2*beta*eta*H+2*(alpha+alpha*bH)**2
                    - alpha*(3+2*bH+2*eta*(-1+bH*(2+bH))))))
    
    B1 = -gk*alpha*epH*eta*(-eta-alpha*(-1+alpha+eta)
                           + ep2H*(eta+alpha*(-1+alpha+eta))+2*alpha*beta*(alpha-eta)*H) \
        + fk*(-(alpha*eta*(eta+alpha*(-1+alpha+eta))) + alpha*ep2H
             *eta*(eta+alpha*(-1+alpha+eta+2*beta*(alpha-eta)*H)))
    
    C1 = gk*(-(ep3H*eta*(-((-1+2*alpha)*(eta+alpha*(-1+alpha+eta)))
                        - alpha*beta*(alpha*(-1+3*alpha-5*eta)+eta)*H))
            + epH*eta*(eta+alpha*(-1+alpha+eta))*(1+alpha*(-2+bH))) \
        - fk*ep2H*eta*(eta+(-1+2*alpha)*ep2H*(eta+alpha*(-1+alpha+eta))
                      + alpha*(-1-2*alpha**2*(-1+bH)**2+eta*(-1+2*bH)
                              + alpha*(3-2*bH+2*eta*(-1+bH*(-2+bH)))))
    
    D1 = -fk*alpha*ep2H*eta*(-eta-alpha*(-1+alpha+eta)
                            + ep2H*(eta+alpha*(-1+alpha+eta))+2*alpha*beta*(alpha-eta)*H) \
        + gk*(-(alpha*epH*eta*(eta+alpha*(-1+alpha+eta))) + alpha
             *ep3H*eta*(eta+alpha*(-1+alpha+eta+2*beta*(alpha-eta)*H)))
    
    A1[beta != 0] = A1[beta != 0]/denom[beta != 0]
    B1[beta != 0] = B1[beta != 0]/denom[beta != 0]
    C1[beta != 0] = C1[beta != 0]/denom[beta != 0]
    D1[beta != 0] = D1[beta != 0]/denom[beta != 0]
    
    arg = beta*zobs
    ep = np.zeros_like(arg)
    en = np.exp(50*np.ones_like(arg))
    ep[arg >= -50] = np.exp(arg[arg >= -50])
    en[arg >= -50] = np.exp(-arg[arg >= -50])
    
    arg1 = beta*ep*(A1+B1+B1*beta*zobs)-beta*en*(C1-D1+beta*D1*zobs)
    arg2 = ep*(A1*alpha-2*B1+2*alpha*B1+alpha*B1*beta*zobs) \
            + en*(C1*alpha+2*D1-2*alpha*D1+alpha*D1*beta*zobs)
    arg3 = ep*(A1*alpha-2*B1+3*alpha*B1+alpha*B1*beta*zobs) \
            - en*(C1*alpha+2*D1-3*alpha*D1+alpha*D1*beta*zobs)
    arg4 = beta**2*ep*(A1+2*B1+B1*beta*zobs) \
            + beta**2*en*(C1-2*D1+beta*D1*zobs)
    
    cub = (2*pi*kx*alpha*arg1.imag + 1j*(-2*pi*kx*alpha*arg1.real))
    cvb = (2*pi*ky*alpha*arg1.imag + 1j*(-2*pi*ky*alpha*arg1.real))
    cwb = (-beta**2*arg2.real + 1j*(-beta**2*arg2.imag))
    cdwdz = (-beta**3*arg3.real + 1j*(-beta**3*arg3.imag))
    cdudz = (2*pi*kx*alpha*arg4.imag + 1j*(-2*pi*kx*alpha*arg4.real))
    cdvdz = (2*pi*ky*alpha*arg4.imag + 1j*(-2*pi*ky*alpha*arg4.real))
    cub[beta == 0] = (0 + 1j*0)
    cvb[beta == 0] = (0 + 1j*0)
    cwb[beta == 0] = (0 + 1j*0)
    cdwdz[beta == 0] = (0 + 1j*0)
    cdudz[beta == 0] = (0 + 1j*0)
    cdvdz[beta == 0] = (0 + 1j*0)
    
    return cub, cvb, cwb, cdwdz, cdudz, cdvdz

def disp2stress(kx,ky,cu,cv,cw,cdudz,cdvdz,cdwdz,rlam1,rmu1):
    # compute horizontal derivatives
    cux = (1j*2*pi*kx)*cu
    cvy = (1j*2*pi*ky)*cv
    cuy = (1j*2*pi*ky)*cu
    cvx = (1j*2*pi*kx)*cv
    cwx = (1j*2*pi*kx)*cw
    cwy = (1j*2*pi*ky)*cw
    cwz = cdwdz
    cuz = cdudz
    cvz = cdvdz
    
    # normal stresses
    cTxx = (rlam1+2*rmu1)*cux+rlam1*(cvy+cwz)
    cTyy = (rlam1+2*rmu1)*cvy+rlam1*(cux+cwz)
    cTzz = (rlam1+2*rmu1)*cwz+rlam1*(cux+cvy)
    
    # shear stresses
    cTxy = rmu1*(cuy+cvx)
    cTxz = rmu1*(cuz+cwx)
    cTyz = rmu1*(cvz+cwy)
    
    return cTxx, cTyy, cTzz, cTxy, cTxz, cTyz

def topo_stress(infile, zobs, H=7, Te=0, rhoc=2900):
    grv = 9.81
    rhom = 3300.0
    rhow = 1025.0
    young = 7.0e10
    rnu = 0.25
    rlam1 = young*rnu/(1+rnu)/(1-2*rnu)
    rmu1 = young/2/(1+rnu)
    bulk = young/3/(1-2*rnu)
    H = -abs(H*1000)
    Te = abs(Te*1000)
    zobs = -abs(zobs*1000)
    
    lon, lat, topo = read_GMT_netcdf(infile)
    
    ni = topo.shape[0]
    ni2 = int(ni/2+1)
    nj = topo.shape[1]
    nj2 = int(nj/2+1)
    
    window = compute_window(topo)
    load = topo*grv*rhoc*window/(ni*nj)
    
    load_k = fft_grd(load)
    dln = np.abs(lon[-1] - lon[0])/nj
    dlt = np.abs(lat[-1] - lat[0])/ni
    rlat2 = np.abs(lat[0]+ni*dlt/2)*pi/180
    xscl = np.cos(rlat2)
    dx = xscl*111000*dln
    dy = 111000*dlt
    width = nj*dx
    height = np.abs(ni*dy)
    
    ky_idx = np.linspace(1,ni2-1,ni2-1)
    kyP = (ni2 - 1 - np.flip(ky_idx))/height
    kyM = -1*np.linspace(1,ni2-1,ni2-1)/height
    ky = np.concatenate((kyP,np.flip(kyM)))
    
    kx_idx = np.linspace(1,nj2-1,nj2-1)
    kxP = (nj2 - 1 - np.flip(kx_idx))/width
    kxM = -1*np.linspace(1,nj2-1,nj2-1)/width
    kx = np.concatenate((kxP,np.flip(kxM)))
    kX, kY = np.meshgrid(kx,ky)
    
    g_k = moho_topo(kX,kY,Te,rhoc,rhom,load_k,young,rnu)
    cub, cvb, cwb, cdwdz, cdudz, cdvdz = airy2load(kX,kY,zobs,H,load_k,g_k,rlam1,rmu1)
    cTxx, cTyy, cTzz, cTxy, cTxz, cTyz = disp2stress(kX,kY,cub,cvb,cwb,cdudz,cdvdz,cdwdz,rlam1,rmu1)
    
    Txx = fft_grd(cTxx,inv=True)
    Tyy = fft_grd(cTyy,inv=True)
    Tzz = fft_grd(cTzz,inv=True)
    Txy = fft_grd(cTxy,inv=True)
    Txz = fft_grd(cTxz,inv=True)
    Tyz = fft_grd(cTyz,inv=True)
    
    return Txx, Tyy, Tzz, Txy, Txz, Tyz

if __name__ == '__main__':
    print('testing...')
    #Txx, Tyy, Tzz, Txy, Txz, Tyz = topo_stress('toymodel.nc',-0.5)