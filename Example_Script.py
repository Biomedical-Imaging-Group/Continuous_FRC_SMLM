from FRCcont import *
import matplotlib.pyplot as plt
import time
#-----------------------------------------------------------
# Continuous Fourier Ring Correlation
#
# Exemple of use of the function FRCcont with and without the use of the
# lookup table. 
#
# Copyright (2018) Thanh-An Pham (thanh-an.pham@epfl.ch)
#                  Emmanuel Soubies (esoubies@gmail.com)
#-----------------------------------------------------------

## Parameters
fov=6400         # Field of view to generate points (nm)     
Nmol=1000        # Number of molecules
stdMol=8         # Standard deviation between the position of the molecules of the two sets
Nsamples = 200   # Number of points in the FRC curve
fmax = 0.1       # Maximal frequency. The FRC will be computed on [0,fmax]

## Load data
gt=np.random.rand(Nmol,2)*fov
loc=gt+np.random.randn(Nmol,2)*stdMol

## FRC computation
tic = time.time()
FRCcLook,freq = FRCcont(loc,gt,fov,Nsamples,fmax,True)
print('Computation with lookup table {:.2f} sec'.format(time.time() - tic))
FRCcNoLook,freq = FRCcont(loc,gt,fov,Nsamples,fmax,False)
print('Computation without lookup table {:.2f} sec'.format(time.time() - tic))

## Plots
plt.figure()
plt.plot(freq,FRCcLook,'r')
plt.plot(freq,FRCcNoLook,'g--')
plt.grid(True)
plt.xlabel('f')
plt.ylabel('FRC')
#plt.set(gca,'FontSize',14)
plt.legend(['Lookup Table','No Lookup Table'])
plt.xlim([0, 0.1])
plt.ylim([0, 1])
plt.show()

