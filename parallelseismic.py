from math import *
import numpy as np

"""
PS_forward_exact.py

calculates traveltimes (first arrival) for the parallel seismic method
for:
    R: pile radius in m
    D: distance pile edge-sensor in m
    zs: array of sensor depths in  m
    cp: wave velocity in concrete in m/s
    cs: wave velocity in soil in m/s
    L: pile length in m
    off: time offset in ms
    nsteps: steps for raytracing

using a simple ray tracing and interpolation approach
based on Niederleithinger 2010 (PhD thesis)

Niederleithinger BAM 8.2 2019
"""
def PSfe(R,D,zs,cp,cs,L,off,nsteps=500):
        n=zs.size
        traveltimes = np.zeros((n))
        zi = np.zeros((nsteps))
        ti = np.zeros((nsteps))
# calculation of traveltimes, upper part (z<=L)
        ni = int(round(nsteps/2))-1
        dz = L/ni
        for i in range(ni+1):
                z = (i)*dz
                ti[i] = sqrt(R*R+z*z)/cp
                alpha = atan(z/R)
                beta = asin(sin(alpha)*cs/cp)
                ti[i] = ti[i] + D/cos(beta)/cs
                zi[i] = z + tan(beta)*D
        z=zi[i-1];
# calculation of traveltimes, lower part (z<=L)        
        dz = (zs[n-1]-z)/ni
        t1 = sqrt(R*R+L*L)/cp
        for j in range(ni+1):
                z2 = z+j*dz
                ti[j+i-1] = t1 + sqrt(D*D+(z2-L)*(z2-L))/cs
                zi[j+i-1] = z2
        for i in range(n):
                j=1
                while (zi[j]<zs[i]):
                            j=j+1
                traveltimes[i] = ti[j]+(ti[j+1]-ti[j])/(zi[j+1]-zi[j])*(zs[i]-zi[j]);
        traveltimes = traveltimes + off
       
        return traveltimes