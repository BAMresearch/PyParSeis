"""
PyParSeis.py

PyParSeis contains routines to evaluate Parallel Seismic measurements
based on Ernst Niederleithinger's PhD thesis) 2010 and a subsequent 
publication:
Ernst Niederleithinger, 2012: "Improvement and extension of the parallel 
seismic method for foundation depth measurement", SOILS AND FOUNDATIONS 52 (6):
1093–101. https://doi.org/10.1016/j.sandf.2012.11.023

Parallel Seismic measurements are performed to determine foundation pile 
(or similar elements) lenght in the subsoil. A hammer blow is perfomed on the 
pile top-center to generate a downgoing stresswave and to trigger the recording. 
Arrriving signal are detected by geophyes od hydrophones in a nearby, parallel
borehole. For details see references.

While the original code was Matlab, it hase been recoded and extended
from 2019 to Python by the author. 

It is published at: https://github.com/BAMresearch/PyParSeis

The code is distributed under a GPL 3.0 license. Use at your own risk.
"""

"""
List of functions
1) PSfe(R,D,zs,cp,cs,L,off,nsteps=500)
calculates traveltimes (first arrival) for the parallel seismic method
for a borehole parallel to the pile.

2) PSfei(R,D,zs,cp,cs,L,off,inc, nsteps=500)
as PSfe, but for an inclined borehole with inclination inc in deg



PyParseis uses 
- math (for trigonometry)
- numpy (for array handling)

"""

from math import *
import numpy as np


"""
PSfe(R,D,zs,cp,cs,L,off,nsteps=500)
calculates traveltimes (first arrival) for the parallel seismic method. As no
direct anylytical solution exists, the travel paths (and thus the travel
times) are calculated for various angles (in nsteps steps, see parameters) 
from the vertical axes. As these paths might not "hit" the sensors exactly,
the travel time are interpolated by the nearest two travel paths. 
   
Parameters:
    R: pile radius in m
    D: distance pile edge-sensor in m
    zs: array of sensor depths in  m
    L: pile length in m
    cp: p-wave velocity in concrete in m/s (actually a stress wave, almost a bar wave)
    cs: p-wave velocity in soil in m/s
    nsteps: number od steps for raytracing
        
        
Return:
    traveltimes: first arrival travel time at sensor positions 

Note that z is not really the depth, the z-Axis is running parallel to the pile,
z=0 at pile top



PSfei(R,D,zs,cp,cs,L,off,inc, nsteps=500)
as PSfei, but for an inclined borehole with inclination inc in deg
        
        
Return:
    traveltimes: first arrival travel time at sensor positions 

Note that z is not really the depth, the z-Axis is running parallel to the pile,
z=0 at pile top

Based on Niederleithinger 2010 (PhD thesis)

Original Niederleithinger BAM 8.2 2019
16.08.2025: editing comments


"""



def PSfe(R,D,zs,cp,cs,L,off,nsteps=500):
        n=zs.size
        traveltimes = np.zeros((n))
        zi = np.zeros((nsteps))
        ti = np.zeros((nsteps))
    
# Calculation of traveltimes, upper part (zi<= L)
        ni = int(round(nsteps/2))-1
        dz = L/ni
        for i in range(ni+1): 
        	    z = (i)*dz
	            ti[i] = sqrt(R*R+z*z)/cp
	            alpha = atan(z/R)
	            beta = asin(sin(alpha)*cs/cp)
	            ti[i] = ti[i] + D/cos(beta)/cs
	            zi[i] = z + tan(beta)*D		
        z=zi[i-1]

# Calculation of traveltimes, lower part 
        dz = (zs[n-1]-z)/ni
        t1 = sqrt(R*R+L*L)/cp
        for j in range(ni+1):
	            z2 = z+j*dz
	            ti[j+i-1] = t1 + sqrt(D*D+(z2-L)*(z2-L))/cs
	            zi[j+i-1] = z2
		
# Interpolation für sensor depths

        for i in range(n):
	            j=1
	            while (zi[j]<zs[i]): 
                        j=j+1 
    
	            traveltimes[i] = ti[j]+(ti[j+1]-ti[j])/(zi[j+1]-zi[j])*(zs[i]-zi[j]);

        traveltimes = traveltimes + off

        return traveltimes
    
    
    
"""   
PSfei(R,D,zs,cp,cs,L,off,inc, nsteps=500)
as PSfe, but for an inclined borehole with inclination inc [deg]
        
        
Return:
    traveltimes: first arrival travel time at sensor positions 

Note that z is not really the depth, the z-Axis is running parallel to the pile,
z=0 at pile top

Note: the zs are sensors depth along the (inclined) borehole.

Note: in dissertation sin(inc) in calculation of zi[i] 
      instead of tan(inc) which has led to a small error

Note: PSfei does not consider boreholes inclined towards the pile in a way
      that some sensors are under the pile

Based on Niederleithinger 2010 (PhD thesis)

Original Niederleithinger BAM 8.2 2026
22.01.2026: first edition
30.01.2026: debugging


"""

def PSfei(R,D,zs,cp,cs,L,off,inc,nsteps=500):
        n=zs.size
        traveltimes = np.zeros((n))
        zi = np.zeros((nsteps))
        ti = np.zeros((nsteps))
    
# Calculation of traveltimes, upper part (zi<= L)
        ni = int(round(nsteps/2))-1
        dz = L/ni
        for i in range(ni+1): 
             za1 = i*dz
             ti[i] = sqrt(R*R+za1*za1)/cp
             alpha1 = atan(za1/R)
             theta1 = asin(sin(alpha1)*cs/cp)
             zi[i] = (za1 + D*tan(theta1))/(1-tan(theta1)*tan(radians(inc)))	               
             ti[i] = ti[i] + (D+zi[i]*sin(radians(inc)))/cos(theta1)/cs

        z=zi[i-1]

# Calculation of traveltimes, lower part 
        dz = (zs[n-1]-z)/ni
        t1 = sqrt(R*R+L*L)/cp
        for j in range(ni+1):
            z2 = z+j*dz
            d = D+z2*tan(radians(inc))
            ti[j+i-1] = t1 + sqrt(d*d+(z2-L)*(z2-L))/cs
            zi[j+i-1] = z2
		
# Interpolation für sensor depths

        for i in range(n):
            j=1
            while (zi[j]<zs[i]*cos(radians(inc))):
                j=j+1 
                traveltimes[i] = ti[j]+(ti[j+1]-ti[j])/(zi[j+1]-zi[j])*((zs[i]*cos(radians(inc)))-zi[j]);

        traveltimes = traveltimes + off

        return traveltimes


    
