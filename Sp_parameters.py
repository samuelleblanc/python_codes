# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Run the different import required at start

# <codecell>

from __future__ import division
import numpy as np
import math
import os
import warnings
import sys
warnings.simplefilter('ignore', np.RankWarning)

# <markdowncell>

# Prepare the different functions require to run the class Sp

# <codecell>

def closestindex(a,x):
    " Get the index from a of the closest value from x "
    return (np.abs(a - x)).argmin()

def find_closest(A, target):
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx

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
    from linfit import linfit
    spc, mask = nanmasked(sp)
    wvl = wvlin[mask]
    norm = norm2max(spc)
    [i1000,i1077,i1493,i1600,i1200,i1300,i530,i610,
     i1565,i1634,i1193,i1198,i1236,i1248,i1270,i1644,
     i1050,i1040,i1065,i600,i870,i515] = find_closest(wvl,np.array([1000,1077,1493,1600,1200,1300,530,
                                                     610,1565,1634,1193,1198,1236,1248,
                                                     1270,1644,1050,1040,1065,600,870,515]))
    norm2 = spc/spc[i1000]
    dsp = smooth(np.gradient(norm2,wvl/1000),2)
    npar = 16
    imaxwvl = np.argmax(spc)
    maxwvl = wvl[mask[imaxwvl]]
    # now calculate each parameter
    fit0 = np.polyfit(wvl[i1000:i1077],norm2[i1000:i1077],1)
    fit0_fn = np.poly1d(fit0)
    fit7 = np.polyfit(wvl[i1493:i1600],norm2[i1493:i1600],1)
    fit7_fn = np.poly1d(fit7)
    fit8,z = linfit(wvl[i1000:i1077],dsp[i1000:i1077])
    fit9,z = linfit(wvl[i1200:i1300],dsp[i1200:i1300])
    fit10,z = linfit(wvl[i530:i610]/1000,norm[i530:i610])
    fit14,z = linfit(wvl[i1565:i1634],spc[i1565:i1634]/norm[i1565])
    par = [sum(norm2[i1000:i1077]-fit0_fn(wvl[i1000:i1077])),   # curvature of rad normed to 1000 nm for 1000 nm - 1077 nm
           dsp[i1193],                                          # deriv of rad normed to 1000 nm at 1193 nm
           dsp[i1493],                                          # deriv of rad normed to 1000 nm at 1493 nm
           norm[i1198]/norm[i1236],                             # ratio of normalized rad of 1198 nm / 1236 nm
           np.nanmean(norm[i1248:i1270]),                       # mean of normalized rad between 1248 nm - 1270 nm
           np.nanmean(norm[i1565:i1644]),                       # mean of normalized rad between 1565 nm - 1644 nm
           np.nanmean(norm[i1000:i1050]),                       # mean of normalized rad between 1000 nm - 1050 nm
           sum(norm2[i1493:i1600]-fit7_fn(wvl[i1493:i1600])),   # curvature of rad normed to 1000 nm for 1493 nm - 1600 nm
           fit8[1],                                             # slope of deriv of rad normed to 1000 nm, 1000 nm - 1077 nm
           fit9[1],                                             # slope of deriv of rad normed to 1000 nm, 1200 nm - 1300 nm
           fit10[1],                                            # slope of normalized radiance between 530 nm - 610 nm
           norm[i1040],                                         # normalized radiance at 1040 nm
           norm[i1000]/norm[i1065],                             # ratio of normalized radiance at 1000 nm / 1065 nm
           norm[i600]/norm[i870],                               # ratio of normalized radiance at 600 nm / 870 nm
           np.nanmax([0.003,fit14[1]]),                         # slope of radiance / rad at 1565 between 1565 nm - 1634 nm
           spc[i515]]                                           # radiance at 515 nm
    return par

# <markdowncell>

# For fancy ouputting of progress bars

# <codecell>

def startprogress(title):
    """
    Creates a progress bar 40 chars long on the console
    and moves cursor back to beginning with BS character
    """
    global progress_x, title_global
    title_global = title
    sys.stdout.write(title + ": [" + "-" * 40 + "]")
    sys.stdout.flush()
    progress_x = 0


def progress(x):
    """
    Sets progress bar to a certain percentage x.
    Progress is given as whole percentage, i.e. 50% done
    is given by x = 50
    """
    global progress_x, title_global
    x = int(x * 40 // 100)                      
    sys.stdout.write("\r" + title_global + ": [" + "#" * x + "-" * (40 - x) + "]")
    sys.stdout.flush()
    progress_x = x


def endprogress():
    """
    End of progress bar;
    Write full bar, then move to next line
    """
    global title_global
    sys.stdout.write("\r" + title_global + ": [" +"#" * 40 + "]\n")
    sys.stdout.flush()

# <markdowncell>

# Create a class for handling the luts

# <codecell>

class lut:
    """
    lut: look up table class object that is used for specific parameter tables
    Has sub objects similar to Sp class
    """
    def __init__(self,par,tau,ref,phase=None):
        self.par = par
        self.tau = tau
        self.ref = ref
        self.npar = list(set(par.shape) - set([tau.size,ref.size,2]))[0] # 2 for the possible 2 phases
        self.phase = phase
        
    def norm_par(self,pcoef=None):
        """ 
        Normalize the parameters, if no keyword set, returns the normalized parameter values in self.parn
        if the keywords are not set, returns the normalized parameter values in self.parn and returns
        the pcoef dictionary containing the coefficients and additive values used for normalizing each parameter.
        
        """
        import numpy as np
        npar = self.npar
        par = self.par
        parn = np.zeros_like(self.par)
        if pcoef is not None:
            # for defined coefficients
            if self.phase is 'liq' or 'ice':
                parn = np.transpose(self.par*pcoef['coef'].T-pcoef['add'].T)
            elif self.phase is None:
                if self.par.ndim == 4:
                    for ph in [0,1]:
                        parn[ph,:,:,:] = np.transpose(self.par[ph,:,:,:]*pcoef['coef'].T-pcoef['add'].T)
                elif self.par.ndim == 5:
                    for ph in [0,1]:
                        parn[ph,:,0,:,:] = np.transpose(self.par[ph,:,0,:,:]*pcoef['coef'].T-pcoef['add'].T)
            else:
                warning('There is a problem with the phase of parameters')
                return
        else:
            # for defining its own coefficients
            pco = np.empty(self.npar)
            padd = np.empty(self.npar)
            if self.par.ndim == 3:
                for p in xrange(self.npar):
                    pco[p] = 1./(np.nanmax(self.par[:,:,p])-np.nanmin(self.par[:,:,p]))
                    padd[p] = np.nanmin(self.par[:,:,p]*pco[p])
                    parn[:,:,p] = self.par[:,:,p]*pco[p]-padd[p]
            elif self.par.ndim == 4:
                for p in xrange(self.npar):
                    pco[p] = 1./(np.nanmax(self.par[:,:,:,p])-np.nanmin(self.par[:,:,:,p]))
                    padd[p] = np.nanmin(self.par[:,:,:,p]*pco[p])
                    for ph in [0,1]:
                        parn[ph,:,:,p] = self.par[ph,:,:,p]*pco[p]-padd[p]
            elif self.par.ndim == 5:
                for p in xrange(self.npar):
                    pco[p] = 1./(np.nanmax(self.par[:,p,0,:,:])-np.nanmin(self.par[:,p,0,:,:]))
                    padd[p] = np.nanmin(self.par[:,p,0,:,:]*pco[p])
                    for ph in [0,1]:
                        parn[ph,p,0,:,:] = self.par[ph,p,0,:,:]*pco[p]-padd[p]
            pcoef = {'coef':pco,'add':padd}
        self.parn = parn
        return pcoef
        

# <markdowncell>

# Create a new class for spectra that returns easy plotting, normalization, parameters, etc.

# <codecell>

class Sp:
    """ 
    Sp: spectrum class object that has all the tools for easy spectra analysis.
    Includes tools for building look up tables
    """    
    def params(self):
        " Runs through each spectrum in the sp array to calculate the parameters"
        from Sp_parameters import param
        import itertools
        sp = self.sp
        wvl = self.wvl
        if sp.ndim == 5:
            w,r,t = np.mgrid[0:2,0:len(self.ref),0:len(self.tau)]
            partemp = map(lambda w,r,t:param(sp[w,:,0,r,t],wvl),w.ravel(),r.ravel(),t.ravel())
            par = np.reshape(partemp,[2,len(self.ref),len(self.tau),-1])
            self.npar = par.shape[3]
        elif sp.ndim == 2:
            par = np.array(map(lambda tt:param(sp[tt,:],wvl),xrange(len(self.utc))))
            self.npar = par.shape[1]
        else: raise LookupError
        self.par = par
        return par
    
    def sp_hires(self):
        " Interpolate the modeled sp to a finer resolution in tau and ref. Only works on modeled data with 5 dimensions"
        #from Sp_parameters import startprogress, progress, endprogress
        from scipy import interpolate
        import numpy as np
        sp = self.sp
        wvl = self.wvl
        tau = self.tau
        ref = self.ref
        if sp.ndim !=5:
            warning('This method cannot be applied to this object, not enough dimensions')
            return
        tau_hires = np.concatenate((np.arange(tau[0],1.0,0.1),np.arange(1.0,4.0,0.5),np.arange(4.0,tau[-1]+1.0,1.0)))
        ref_hires = np.arange(ref[0],ref[-1]+1.0)
        print tau_hires.shape
        print ref_hires.shape
        #ttau = sp[0,0,0,:,:]*0.+tau
        #rref = np.transpose(sp[0,0,0,:,:].T*0.+ref.T)
        ttau = tau
        rref = ref
        sp_hires = np.zeros([2,len(wvl),2,len(ref_hires),len(tau_hires)])
        startprogress('Running interpolation')
        for w in xrange(len(wvl)):
            for ph in [0,1]:
                fx = interpolate.RectBivariateSpline(rref,ttau,sp[ph,w,0,:,:],kx=1,ky=1)
                sp_hires[ph,w,0,:,:] = fx(ref_hires,tau_hires)
            progress(w/len(wvl)*100.0)
        endprogress()
        print('Overwriting the current sp, tau, and ref with the new high resolution values')
        self.sp = sp_hires
        self.tau = tau_hires
        self.ref = ref_hires
        return       
    
    def norm_par(self,pcoef=None):
        """ 
        Normalize the parameters, if no keyword set, returns the normalized parameter values in self.parn
        if the keywords are not set, returns the normalized parameter values in self.parn and returns
        the pcoef dictionary containing the coefficients and additive values used for normalizing each parameter.
        
        """
        import numpy as np
        npar = self.npar
        par = self.par
        parn = np.zeros_like(self.par)
        if pcoef is not None:
            # for defined coefficients
            if self.par.ndim == 4:
                for ph in [0,1]:
                    parn[ph,:,:,:] = np.transpose(self.par[ph,:,:,:]*pcoef['coef'].T-pcoef['add'].T)
            elif self.par.ndim == 2:
                for t in xrange(len(self.utc)):
                    parn[t,:] = self.par[t,:]*pcoef['coef']-pcoef['add']
            elif self.par.ndim == 5:
                for ph in [0,1]:
                    parn[ph,:,0,:,:] = np.transpose(self.par[ph,:,0,:,:]*pcoef['coef'].T-pcoef['add'].T)
            else:
                warning('There is a problem with the dimensions of the sp')
                return
        else:
            # for defining its own coefficients
            pco = np.empty(self.npar)
            padd = np.empty(self.npar)
            if self.par.ndim == 4:
                for p in xrange(self.npar):
                    pco[p] = 1./(np.nanmax(self.par[:,:,:,p])-np.nanmin(self.par[:,:,:,p]))
                    padd[p] = np.nanmin(self.par[:,:,:,p]*pco[p])
                    for ph in [0,1]:
                        parn[ph,:,:,p] = self.par[ph,:,:,p]*pco[p]-padd[p]
            elif self.par.ndim == 5:
                for p in xrange(self.npar):
                    pco[p] = 1./(np.nanmax(self.par[:,p,0,:,:])-np.nanmin(self.par[:,p,0,:,:]))
                    padd[p] = np.nanmin(self.par[:,p,0,:,:]*pco[p])
                    for ph in [0,1]:
                        parn[ph,p,0,:,:] = self.par[ph,p,0,:,:]*pco[p]-padd[p]
            pcoef = {'coef':pco,'add':padd}
        self.parn = parn
        return pcoef
    
    def __init__(self,s,**kwargs):
        import numpy as np
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
        import numpy as np
        idxwvl = [i for i in range(sp.ndim) if sp.shape[i] == len(self.wvl)]
        if sp.ndim == 2:
            norm = sp/np.nanmax(sp,axis=idxwvl[0])[:,None]
        elif sp.ndim == 5:
            norm = sp/np.nanmax(sp,axis=idxwvl[0])[:,None,:,:,:]
        return norm
    
    def reshape_lut(self,phase=None):
        " Function that reformats the lut such that only liquid clouds exist "
        if phase is None:
            warning('No phase selected, returning nothing')
            return
        elif phase is 'liq':
            ref = self.ref[self.ref <= 30]
            if self.par.ndim == 4:
                par = self.par[0,self.ref<=30,:,:]
            elif self.par.ndim == 5:
                par = np.rollaxis(self.par[0,:,0,self.ref<=30,:],0,3)
            else:
                warning('Problem with par dimensions')
            self.liq = lut(par,tau,ref,phase=phase)
        elif phase is 'ice':
            ref = self.ref[self.ref >= 10]
            if self.par.ndim == 4:
                par = self.par[1,self.ref>=10,:,:]
            elif self.par.ndim == 5:
                par = np.rollaxis(self.par[1,:,0,self.ref>=10,:],0,3)
            else:
                warning('Problem with par dimensions')
            self.ice = lut(par,tau,ref,phase=phase)
        return
            

# <codecell>


