# -*- coding: utf-8 -*-
#!/usr/bin/python

"""
PyParseisModel.py

Calculate synthetic parallel seismic first arrival data
Based on Niederleithinger 2010 (PhD thesis)

Niederleithinger BAM 8.2 2020
last change(s)
Niederleithinger 11.08.2025 Changing import filename
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.text as pltext

from PS_forward_exact import PSfe
        
# main script

#set parameters
R = 0.3 # pile radius in m
D = 0.5 # distance pile edge - sensor in m

z0 = 1.0 # shallowest sensor position in m (positive downwards)
dz = 0.5  # sensor interval in m 
z1 = 15 # deepest sensor position in m

z=np.arange(z0,z1*1.001,dz) # z1*1.001 to include z1 in array

cp = 4200.0 # concrete velocity in m/s
cs = 1500.0 # soil velocity  in m/s
L = 11.0 # pile length in m
off = 0.0005 # time offset in s
inc = 0.0 # inclination in deg

# evaluate start model
tcalc = PSfe(R,D,z,cp,cs,L,off)

#add noise 
noiseamp = 0.0/1000  # noise amplitude in s
iz = z.size
noise=np.random.random_sample((iz,))
noise= (noise-0.5)*noiseamp
tcalc =tcalc+noise

# plot
cm2inch=1/2.54
fig = plt.figure(figsize=(13*cm2inch,17*cm2inch))

plt.plot(tcalc*1000, -z, 'ok', label='data')
plt.xlabel("t in ms")
plt.ylabel("z below reference in m")
plt.legend()
xmin, xmax, ymin, ymax = plt.axis()
plt.text(xmin+(xmax-ymin)/100, ymin+(ymax-ymin)/100, 'cp=%5.0f, cs=%5.0f, L=%5.2f, off=%5.2f, inc=%5.2f' % (cp,cs,L,off*1000,inc))
plt.savefig('parseis-fit.png', dpi = 300)
plt.show()
plt.close()

# write first arrival data to ASCII text file

filename='PyParseisModel-R0.3-D0.5-cp4200-cs1500-L11-off0.5-inc0-noise0.dat' 
f = open(filename,'w')
for i in range (0,iz):
    f.write(str(z[i])+' '+str(tcalc[i]*1000)+'\n')
f.close()

 