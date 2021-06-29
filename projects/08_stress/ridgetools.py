import numpy as np
import math
import pygmt
import pandas as pd
from scipy import integrate
import os

def load_one_ridge(filename):
    ridge = pd.read_table(filename,skiprows=1,names=['lon','lat'],header=0)
    return ridge

def sample_ridge_on_grid(points,grid):
    Ttype = os.path.basename(grid)[0:3] # this should be evident by the grid string
    dep = os.path.basename(grid)[5]
    gridtype = Ttype + dep
    z = pygmt.grdtrack(points,grid,newcolname=gridtype)
    return z

def get_avg_stress_tensor(stresses):
    cumStress = np.zeros((3,3))
    for ii in range(0,7):
        oneStress = np.array(
            [[stresses[ii*6], stresses[(ii*6)+3], stresses[(ii*6)+4]],
            [stresses[(ii*6)+3], stresses[(ii*6)+1], stresses[(ii*6)+5]],
            [stresses[(ii*6)+4], stresses[(ii*6)+5], stresses[(ii*6)+2]]]
        )
        cumStress = cumStress + oneStress*1000
    return cumStress

def get_ridge_azimuth(ridgeData):
    A = ridgeData.iloc[0,0:2]
    B = ridgeData.iloc[-1,0:2]
    d2r = math.pi/180
    dL = (B[0] - A[0])*d2r
    X = math.cos(B[1]*d2r)*math.sin(dL)
    Y = math.cos(A[1]*d2r)*math.sin(B[1]*d2r) - math.sin(A[1]*d2r)*math.cos(B[1]*d2r)*math.cos(dL)
    beta = math.atan2(X,Y)/d2r
    if beta < 0:
        beta = beta + 360
    return beta

def rotate_stress_tensor(azimuth, stress):
    ridgeNormal = azimuth + 90
    if ridgeNormal > 360:
        ridgeNormal = ridgeNormal - 360
    d2r = math.pi/180
    ridgeNormal=-1*ridgeNormal # remember, rotations are measured anti-clockwise
    rotMat = np.array([[math.cos(ridgeNormal*d2r), -1*math.sin(ridgeNormal*d2r), 0],
        [math.sin(ridgeNormal*d2r),math.cos(ridgeNormal*d2r),0],
        [0,0,1]]
    )
    rotStressTensor = rotMat@stress@rotMat.T
    return rotStressTensor

def get_loading_stress(ridgeDatum,azimuth):
    lon = ridgeDatum[0]
    lat = ridgeDatum[1]
    stresses = ridgeDatum[2:]
    avgStress = get_avg_stress_tensor(stresses)
    rotatedStress = rotate_stress_tensor(azimuth,avgStress)
    loadingStress = rotatedStress[1,1]
    return loadingStress

def lonlat2x(ridgeData):
    n = ridgeData.shape[0]
    r = 6371000
    d2r = math.pi/180
    lon = ridgeData.iloc[:,0]
    lat = ridgeData.iloc[:,1]
    x = np.zeros((n,))
    # could really speed this up by vectorizing
    for ii in range(0,n-1):
        dlat = abs(lat[ii+1] - lat[0])
        dlon = abs(lon[ii+1] - lon[0])
        dS = 2*math.asin(
            math.sqrt(
                math.sin(dlat*d2r/2)**2 + math.cos(lat[0])*math.cos(lat[ii+1])*math.sin(dlon*d2r/2)**2
            )
        )
        x[ii+1] = r*dS

    return x

def compute_K(x,p,H=7e3):
    a = (x[-1] - x[0]) + 1
    K_a = (2/H)*math.sqrt(a/math.pi)
    K_b = integrate.simps(p/np.sqrt(a**2 - x**2),x=x)
    K = K_a*K_b
    return K