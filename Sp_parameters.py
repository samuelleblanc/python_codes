
# coding: utf-8

# Run the different import required at start

# In[ ]:

import numpy as np, h5py
import math
import os
import warnings
from scipy import interpolate
warnings.simplefilter('ignore', np.RankWarning)


# Prepare the different functions require to run the class Sp

# In[5]:

def closestindex(a,x):
    " Get the index from a of the closest value from x "
    return min(range(len(a)), key=lambda i: abs(a[i]-x))

def smooth(x, window):
    """Moving average of 'x' with window size 'window'."""
    return np.convolve(x, np.ones(window)/window, 'same')
        
def nanmasked(x):
    " Build an array with nans masked out and the mask output"
    mask = ~np.isnan(x)
    maskA = x[mask]
    return (maskA,mask)

def norm2max(x):
    " Returns a spectrum, x, that is normalized by its maximum value, ignores nans "
    return x/np.nanmax(x)    
    
def param(sp,wvlin):
    " Calculates the parameters from a spectrum."
    spc, mask = nanmasked(sp)
    wvl = wvlin[mask]
    norm = norm2max(spc)
    i1000 = closestindex(wvl,1000)
    norm2 = spc/spc[i1000]
    dsp = smooth(np.gradient(norm2,wvl/1000),2)
    npar = 16
    imaxwvl = np.argmax(spc)
    maxwvl = wvl[mask[imaxwvl]]
    # now calculate each parameter
    i1077 = closestindex(wvl,1077)
    fit0 = np.polyfit(wvl[i1000:i1077],norm2[i1000:i1077],1)
    fit0_fn = np.poly1d(fit0)
    i1493 = closestindex(wvl,1493)
    i1600 = closestindex(wvl,1600)
    fit7 = np.polyfit(wvl[i1493:i1600],norm2[i1493:i1600],1)
    fit7_fn = np.poly1d(fit7)
    fit8 = np.polyfit(wvl[i1000:i1077],dsp[i1000:i1077],1)
    i1200 = closestindex(wvl,1200)
    i1300 = closestindex(wvl,1300)
    fit9 = np.polyfit(wvl[i1200:i1300],dsp[i1200:i1300],1)
    i530 = closestindex(wvl,530)
    i610 = closestindex(wvl,610)
    fit10 = np.polyfit(wvl[i530:i610]/1000,norm[i530:i610],1)
    i1565 = closestindex(wvl,1565)
    i1634 = closestindex(wvl,1634)
    fit14 = np.polyfit(wvl[i1565:i1634],spc[i1565:i1634]/norm[i1565],1)
    par = [sum(norm2[i1000:i1077]-fit0_fn(wvl[i1000:i1077])),               # curvature of rad normed to 1000 nm for 1000 nm - 1077 nm
           dsp[closestindex(wvl,1193)],                                     # deriv of rad normed to 1000 nm at 1193 nm
           dsp[closestindex(wvl,1493)],                                     # deriv of rad normed to 1000 nm at 1493 nm
           norm[closestindex(wvl,1198)]/norm[closestindex(wvl,1236)],       # ratio of normalized rad of 1198 nm / 1236 nm
           np.nanmean(norm[closestindex(wvl,1248):closestindex(wvl,1270)]), # mean of normalized rad between 1248 nm - 1270 nm
           np.nanmean(norm[closestindex(wvl,1565):closestindex(wvl,1644)]), # mean of normalized rad between 1565 nm - 1644 nm
           np.nanmean(norm[closestindex(wvl,1000):closestindex(wvl,1050)]), # mean of normalized rad between 1000 nm - 1050 nm
           sum(norm2[i1493:i1600]-fit7_fn(wvl[i1493:i1600])),               # curvature of rad normed to 1000 nm for 1493 nm - 1600 nm
           fit8[1],                                                         # slope of deriv of rad normed to 1000 nm, 1000 nm - 1077 nm
           fit9[1],                                                         # slope of deriv of rad normed to 1000 nm, 1200 nm - 1300 nm
           fit10[1],                                                        # slope of normalized radiance between 530 nm - 610 nm
           norm[closestindex(wvl,1040)],                                    # normalized radiance at 1040 nm
           norm[closestindex(wvl,1000)]/norm[closestindex(wvl,1065)],       # ratio of normalized radiance at 1000 nm / 1065 nm
           norm[closestindex(wvl,600)]/norm[closestindex(wvl,870)],         # ratio of normalized radiance at 600 nm / 870 nm
           np.nanmax([0.003,fit14[1]]),                                     # slope of radiance / rad at 1565 between 1565 nm - 1634 nm
           spc[closestindex(wvl,515)]]                                      # radiance at 515 nm
    return par


# Create a new class for spectra that returns easy plotting, normalization, parameters, etc.

# In[10]:

class Sp:
    """ 
    Sp: spectrum class object that has all the tools for easy spectra analysis
    """    
    def params(self):
        " Runs through each spectrum in the sp array to calculate the parameters"
        sp = self.sp
        wvl = self.wvl
        self.npar = 16
        dims = [sp.shape[i] for i in range(len(sp.shape))]
        ui = [i for i in range(sp.ndim) if sp.shape[i] == len(wvl)]
        dims[ui[0]] = self.npar
        par = np.empty((dims))
        par[:] = np.NAN
        if sp.ndim == 5:
            for w in range(dims[0]):
                for t in range(len(self.tau)):
                    for r in range(len(self.ref)):
                        print w,t,r
                        par[w,:,0,r,t] = param(sp[w,:,0,r,t],wvl)
        elif sp.ndim == 2:
            for tt in range(len(self.utc)):
                par[tt,:] = param(sp[tt,:],wvl)
            
        else: raise LookupError
        self.par = par
        return par
    
    def sp_hires(self):
        " Interpolate the modeled sp to a finer resolution in tau and ref. Only works on modeled data with 5 dimensions"
        sp = self.sp
        wvl = self.wvl
        if sp.ndim !=5:
            warning('This method cannot be applied to this object, not enough dimensions')
            return
        tau_hires = np.arange(self.tau[0],0.9,0.1)+np.arange(1.0,3.5,0.5)+np.arange(4.0,self.tau[-1],1.0)
        ref_hires = np.arange(self.ref[0],self.ref[-1])
        ttau = lut.sp[0,0,0,:,:]*0.+lut.tau
        rref = np.transpose(lut.sp[0,0,0,:,:].T*0.+lut.ref.T)
        sp_hires = np.zeros([2,len(wvl),2,len(tau_hires),len(ref_hires)])
        for w in range(len(wvl)):
            for ph in [0,1]:
                 print w,ph
                 fx = interpolate.interp2d(ttau,rref,sp[ph,w,0,:,:])
                 sp_hires[ph,w,0,:,:] = fx(tau_hires,ref_hires)
        print('Overwriting the current sp, tau, and ref with the new high resolution values')
        self.sp = sp_hires
        self.tau = tau_hires
        self.res = ref_hires
        return       
    
    def norm_par(self,pcoef={}):
        """ 
        Normalize the parameters, if no keyword set, returns the normalized parameter values in self.parn
        if the keywords are not set, returns the normalized parameter values in self.parn and returns
        the pcoef dictionary containing the coefficients and additive values used for normalizing each parameter.
        
        """
        npar = self.npar
        par = self.par
        parn = zeros_like(self.par)
        if pceof is not None:
            # for defined coefficients
            if self.sp.dim == 5:
                for ph in [0,1]:
                    parn[ph,:,0,:,:] = np.transpose(self.par[ph,:,0,:,:]*pcoef['coef'].T-pcoef['add'].T)
            elif self.sp.dim == 2:
                for t in range(len(self.utc)):
                    parn[t,:] = self.par[t,:]*pcoef['coef']-pcoef['add']
            else:
                warning('There is a problem with the dimensions of the sp')
                return
        else:
            # for defining its own coefficients
            pco = np.empty(self.npar)
            padd = np.empty(self.npar)
            for p in range(self.npar):
                pco[p] = 1./(np.nanmax(self.par[:,p,0,:,:])-np.nanmin(self.par[:,p,0,:,:]))
                padd[p] = np.nanmin(self.par[:,p,0,:,:]*pco[p])
                for ph in [0,1]:
                    parn[ph,:,0,:,:] = np.transpose(self.par[ph,:,0,:,:]*pco.T-padd.T)
            
            pcoef = {'coef':pco,'add':padd}
        self.parn = parn
        return pcoef
    
    def __init__(self,s,**kwargs):
        if 'nm' in s:
            self.wvl = [item for sublist in s['nm'] for item in sublist]
        if 'zenlambda' in s:
            self.wvl = s['zenlambda']
        if 'wvl' in s:
            self.wvl = [item for sublist in s['wvl'] for item in sublist]
        self.iwvls = np.argsort(self.wvl)
        print len(self.iwvls), len(self.wvl)
        self.wvl = np.sort(self.wvl)
        self.nwvl = len(self.wvl)
        self.sp = self.wvlsort(s)
        self.norm = self.normsp(self.sp)
        print self.sp.shape
        if self.sp.ndim > 2:
            self.tau = s['tau']
            self.ref = s['ref']
       #     self.par = self.params(self.sp,self.wvl,tau=self.tau,ref=self.ref)
        else:
            self.utc = s['utc']
            self.good = s['good']
        #    self.par = self.params(self.sp,self.wvl,utc=self.utc)
        
    def wvlsort(self,s):
        "Function to sort spectra along the wavelength axis"
        iwvls = self.iwvls
        if 'sp' in s:
            print 'in sp'
            print s['sp'].shape
            ui = [i for i in range(s['sp'].ndim) if s['sp'].shape[i] == len(self.wvl)]
            if 1 in ui:
                sp = s['sp'][:,iwvls,:,:,:]
            else: raise LookupError
        if 'rad' in s:
            print 'in rad'
            print s['rad'].shape, s['rad'].ndim, len(iwvls)
            ui = [i for i in range(s['rad'].ndim) if s['rad'].shape[i] == len(self.wvl)]
            print ui
            if 1 in ui:
                print '1 in ui'
                sp = s['rad'][:,iwvls]
            else: 
                print 'not 1 in ui'
                sp = s['rad'][iwvls,:]
        return sp
    
    def normsp(self,sp):
        " Function to return the normalized spectra list"
        idxwvl = [i for i in range(sp.ndim) if sp.shape[i] == len(self.wvl)]
        if sp.ndim == 2:
            norm = sp/np.nanmax(sp,axis=idxwvl[0])[:,None]
        elif sp.ndim == 5:
            norm = sp/np.nanmax(sp,axis=idxwvl[0])[:,None,:,:,:]
        return norm
            


# In[ ]:



