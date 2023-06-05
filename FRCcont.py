#utils
import numpy as np
from scipy import special as sp

def FRCcont(S1,S2,distMax,Nsamples=1e3,fmax=0.1,useLookUp=False,step = 0.001):
    '''
    -----------------------------------------------------------
    Compute the FRC in continuous domain [1]

    Input:  S1        -> First set of points of size (numpy array, N1 x dim)
            S2        -> Second set of points of size (numpy array, N2 x dim)
            distMax   -> Maximal distance between two points (used if useLookUp=true)
            Nsamples  -> Number of samples for the FRC curve (default 1000)
            fmax      -> Maxiamal frequency for the FRC curve (default 0.1)
            useLookUp -> Boolean true to use a lookup table for sp.jv (default false)
            step      -> Step size for lookup table (default 0.001)

    Output: FRC  -> FRC curve
            freq -> frequencies where the FRC is evaluated 

    Reference:
    [1] Closed-Form Expression of the Fourier Ring-Correlation for Single-Molecule Localization Microscopy. 
        Proceeding of ISBI, 2019.
        Thanh-an Pham, Emmanuel Soubies, Daniel Sage, and Michael Unser.

    Copyright (2023) Thanh-An Pham (thanh-an.pham@epfl.ch)
                    Emmanuel Soubies (esoubies@gmail.com)
    -----------------------------------------------------------
    '''

    ## Compute Distances between each pair of points
    D = S1.shape[1]
    dist1 = 0.
    dist2 = 0.
    dist = 0.
    for d in range(0,D):
        dist1 = dist1 + (S1[:,d] -  np.atleast_2d(S1[:,d]).T)**2
        dist2 = dist2 + (S2[:,d] -  np.atleast_2d(S2[:,d]).T)**2
        dist = dist + (S1[:,d] -  np.atleast_2d(S2[:,d]).T)**2
    
    dist1 = np.sqrt(dist1[:,:,np.newaxis])
    dist2 = np.sqrt(dist2[:,:,np.newaxis])
    dist = np.sqrt(dist[:,:,np.newaxis])

    ## Compute FRC
    # - rho discretization
    k = 2*np.pi*fmax
    rho_set_orig = np.linspace(0,k, Nsamples)
    FRC = np.zeros(Nsamples)
    # - Lookup Table
    
    xquery = np.arange(0,distMax,step)[:,np.newaxis].T
    lookupB = sp.jv(0,xquery).squeeze()
    max_len = np.floor(2e9/(4*max(dist.size,max(dist1.size,dist2.size))))
    
    # - Main Loop
    for kk in range(0,np.ceil(Nsamples/max_len).astype(int)): #maybe not ceil
        cind = np.arange(kk*max_len,min((kk+1)*max_len,Nsamples),dtype=int)#.astype(int)
        rho_set = rho_set_orig[cind]
        if useLookUp:
            J0r = np.sum(np.sum(mybesselj(np.reshape(rho_set,(1,1,np.size(rho_set)))*dist,lookupB,step),0),0)
            S1N = np.sum(np.sum(mybesselj(np.reshape(rho_set,(1,1,np.size(rho_set)))*dist1,lookupB,step),0),0)
            S2N = np.sum(np.sum(mybesselj(np.reshape(rho_set,(1,1,np.size(rho_set)))*dist2,lookupB,step),0),0)
        else:
            J0r = np.sum(np.sum(sp.jv(0,np.reshape(rho_set,(1,1,np.size(rho_set)))*dist),0),0)
            S1N = np.sum(np.sum(sp.jv(0,np.reshape(rho_set,(1,1,np.size(rho_set)))*dist1),0),0)
            S2N = np.sum(np.sum(sp.jv(0,np.reshape(rho_set,(1,1,np.size(rho_set)))*dist2),0),0)
        
        FRCdenom = np.sqrt(S1N*S2N)
        FRC[cind] = J0r/FRCdenom
        
    freq=np.linspace(0,fmax,Nsamples)
    return FRC,freq


def mybesselj(x,lookupB,step):
    ## Evaluate Bessel with a lookup table
    sz_in = x.shape
    flx = np.floor(x.flatten()/step).astype(int)
    if np.any(flx>=lookupB.size):
        print('The limit of the lookup table is reached. Increase the maximal distance for the lookup table.')
    else:
        out = np.reshape(lookupB[flx],sz_in) # Faster but noisy
    return out
    