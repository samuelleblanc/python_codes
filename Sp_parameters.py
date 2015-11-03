
# coding: utf-8

# Run the different import required at start

# In[1]:

from __future__ import division
import numpy as np
import math
import os
import warnings
import sys
warnings.simplefilter('ignore', np.RankWarning)


# Prepare the different functions require to run the class Sp

# In[1]:

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

def smooth(x, window,nan=True):
    """
    Moving average of 'x' with window size 'window'.
    If the nan keyword is set to True (default), the array is returned with the nans removed and substituted by an interpolated value
    """
    if nan:
        from Sp_parameters import nanmasked
        from scipy import interpolate
        ix = np.arange(len(x))
        xmasked, mask = nanmasked(x)
        fx = interpolate.interp1d(ix[mask],xmasked,bounds_error=False,fill_value=0)
        xinterp = fx(ix)
        xout = np.convolve(xinterp, np.ones(window)/window, 'same')
    else:
        xout = np.convolve(x, np.ones(window)/window, 'same')
    return xout

def deriv(y,x):
    """
    Numerical derivative based on IDL's deriv function: returns the derivative of a three-point (quadratic) Lagrangian interpolation:
    For evenly space: Y'[0] = (–3*Y[0] + 4*Y[1] – Y[2]) / 2
                      Y'[i] = (Y[i+1] – Y[i–1]) / 2 ; i = 1...N–2
                      Y'[N–1] = (3*Y[N–1] – 4*Y[N–2] + Y[N–3]) / 2
    For uneven: for the middle with x = x1:
              y' = y0x12/(x01x02) + y1(1/x12 – 1/x01) – y2x01/(x02x12)
            The first point is computed at x = x0:
              y' = y0(x01 + x02)/(x01x02) – y1x02/(x01x12) + y2x01/(x02x12)
            The last point is computed at x = x2:
              y' = –y0x12/(x01x02) + y1x02/(x01x12) – y2(x02 + x12)/(x02x12)
            with x12 = x1-x2 and x=[x0,x1,x2] (for the 3 neighbors)
    """
    import numpy as np
    d = np.zeros_like(x)*np.nan
    if len(x) != len(y):
        print '*** warning x is not the same length as y nothing done ***'
        return
    if len(x) < 3:
        print '*** Not enough points ***'
        return
    d[0] = y[0]*((x[0]-x[1])+(x[0]-x[2]))/((x[0]-x[1])*(x[0]-x[2])) -            y[1]*(x[0]-x[2])/((x[0]-x[1])*(x[1]-x[2])) +            y[2]*(x[0]-x[1])/((x[0]-x[2])*(x[1]-x[2]))
    d[-1] = - y[-3]*(x[-2]-x[-1])/((x[-3]-x[-2])*(x[-3]-x[-1])) +             y[-2]*(x[-3]-x[-1])/((x[-3]-x[-2])*(x[-2]-x[-1])) -             y[-1]*((x[-3]-x[-1])+(x[-2]-x[-1]))/((x[-3]-x[-1])*(x[-2]-x[-1]))
    for i in xrange(1,len(x)-1):
        d[i] = y[i-1]*(x[i]-x[i+1])/((x[i-1]-x[i])*(x[i-1]-x[i+1])) +                y[i]*(1.0/(x[i]-x[i+1]) - 1.0/(x[i-1]-x[i])) -                y[i+1]*(x[i-1]-x[i])/((x[i-1]-x[i+1])*(x[i]-x[i+1]))
    return d
        
def nanmasked(x):
    """
    Build an array with nans masked out and the mask output
    Input: x=array to be masked
    Output: maskedarray,maskindices
    """
    mask = ~np.isnan(x)
    maskA = x[mask]
    return (maskA,mask)

def doublenanmask(x,y,return_mask=False):
    """
    Build two arrays that have no nans, in either. Returns smaller array
    if return_mask is set to True (default is False), the mask is also returned
    """
    if len(x) != len(y):
        print "The two arrays don't match sizes, returning"
        return
    mask = ~np.isnan(x) & ~np.isnan(y)
    if return_mask:
        return x[mask],y[mask],mask
    else:
        return x[mask],y[mask]

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.
    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """
    return np.isnan(y), lambda z: z.nonzero()[0]

def norm2max(x,iws=None):
    " Returns a spectrum, x, that is normalized by its maximum value, ignores nans "
    import warnings
    warnings.simplefilter("ignore",RuntimeWarning)
    try:
        if not iws:
            import numpy as np
            iws = np.where(x)[0]
    except ValueError:
        pass
    return x/np.nanmax(x[iws])    
    
def param(sp,wvlin,iws=None):
    " Calculates the parameters from a spectrum."
    from linfit import linfit
    from Sp_parameters import nanmasked, norm2max, smooth, deriv, find_closest
    npar = 16
    spc, mask = nanmasked(sp)
    wvl = wvlin[mask]
    try:
        norm = norm2max(spc,iws=iws)
    except ValueError:
        par = np.zeros(npar)*np.nan
        return par
    [i1000,i1077,i1493,i1600,i1200,i1300,i530,i610,
     i1565,i1634,i1193,i1198,i1236,i1248,i1270,i1644,
     i1050,i1040,i1065,i600,i870,i515] = find_closest(wvl,np.array([1000,1077,1493,1600,1200,1300,530,
                                                     610,1565,1634,1193,1198,1236,1248,
                                                     1270,1644,1050,1040,1065,600,870,515]))
    if np.isnan(spc[i1000]) or not spc[i1000]:
        par = np.zeros(npar)*np.nan
        return par
    norm2 = spc/spc[i1000]
    dsp = smooth(deriv(norm2,wvl/1000),2,nan=False)
    imaxwvl = np.argmax(spc)
    maxwvl = wvl[mask[imaxwvl]]
    # now calculate each parameter
    fit0 = np.polyfit(np.array([wvl[i1000],wvl[i1077]]),np.array([norm2[i1000],norm2[i1077]]),1)
    fit0_fn = np.poly1d(fit0)
    fit7 = np.polyfit(np.array([wvl[i1493],wvl[i1600]]),np.array([norm2[i1493],norm2[i1600]]),1)
    fit7_fn = np.poly1d(fit7)
    fit8,z = linfit(wvl[i1000:i1077],dsp[i1000:i1077])
    fit9,z = linfit(wvl[i1200:i1300],dsp[i1200:i1300])
    fit10,z = linfit(wvl[i530:i610]/1000,norm[i530:i610])
    fit14,z = linfit(wvl[i1565:i1634],spc[i1565:i1634]/norm[i1565])
    par = [sum(norm2[i1000:i1077]-fit0_fn(wvl[i1000:i1077])),   # 1 curvature of rad normed to 1000 nm for 1000 nm - 1077 nm
           dsp[i1198],                                          # 2 deriv of rad normed to 1000 nm at 1198 nm (!=IDL version)
           dsp[i1493],                                          # 3 deriv of rad normed to 1000 nm at 1493 nm
           norm[i1198]/norm[i1236],                             # 4 ratio of normalized rad of 1198 nm / 1236 nm
           np.nanmean(norm[i1248:i1270]),                       # 5 mean of normalized rad between 1248 nm - 1270 nm
           np.nanmean(norm[i1565:i1644]),                       # 6 mean of normalized rad between 1565 nm - 1644 nm
           np.nanmean(norm[i1000:i1050]),                       # 7 mean of normalized rad between 1000 nm - 1050 nm
           sum(norm2[i1493:i1600]-fit7_fn(wvl[i1493:i1600])),   # 8 curvature of rad normed to 1000 nm for 1493 nm - 1600 nm
           fit8[0],                                             # 9 slope of deriv of rad normed to 1000 nm, 1000 nm - 1077 nm
           fit9[0],                                             # 10 slope of deriv of rad normed to 1000 nm, 1200 nm - 1300 nm
           fit10[0],                                            # 11 slope of normalized radiance between 530 nm - 610 nm
           norm[i1040],                                         # 12 normalized radiance at 1040 nm
           norm[i1000]/norm[i1065],                             # 13 ratio of normalized radiance at 1000 nm / 1065 nm
           norm[i600]/norm[i870],                               # 14 ratio of normalized radiance at 600 nm / 870 nm
           np.nanmin([0.003,fit14[0]]),                         # 15 slope of radiance / rad at 1565 between 1565 nm - 1634 nm
           spc[i515]]                                           # 16 radiance at 515 nm
    # do a check for bad points
    if np.all(np.isnan(par[0:13])): 
        par[14] = np.nan
        par[15] = np.nan
    return par


# For fancy ouputting of progress bars

# In[3]:

def startprogress(title):
    """
    Creates a progress bar 40 chars long on the console
    and moves cursor back to beginning with BS character
    """
    global progress_x, title_global
    title_global = title
    sys.stdout.write(title + ": [" + "-" * 40 + "] 00% ")
    sys.stdout.flush()
    progress_x = 0


def progress(x):
    """
    Sets progress bar to a certain percentage x.
    Progress is given as whole percentage, i.e. 50% done
    is given by x = 50
    """
    global progress_x, title_global
    xp = int(x * 40 // 100)                      
    sys.stdout.write("\r" + title_global + ": [" + "#" * xp + "-" * (40 - xp) + "] %02d%%" %x)
    sys.stdout.flush()
    progress_x = x


def endprogress():
    """
    End of progress bar;
    Write full bar, then move to next line
    """
    global title_global
    sys.stdout.write("\r" + title_global + ": [" +"#" * 40 + "]100% -- Done! \n")
    sys.stdout.flush()


# Create a class for handling the luts

# In[4]:

class lut:
    """
    lut: look up table class object that is used for specific parameter tables
    Has sub objects similar to Sp class
    """
    def norm_par(self,pcoef=None):
        """ 
        Normalize the parameters, if no keyword set, returns the normalized parameter values in self.parn
        if the keywords are not set, returns the normalized parameter values in self.parn and returns
        the pcoef dictionary containing the coefficients and additive values used for normalizing each parameter.
        
        """
        import numpy as np
        npar = self.npar
        par = self.par
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
            parn = np.zeros_like(par)*np.nan
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
    
    def __init__(self,par,tau,ref,phase=None):
        self.par = par
        self.tau = tau
        self.ref = ref
        self.npar = list(set(par.shape) - set([tau.size,ref.size,2]))[0] # 2 for the possible 2 phases
        self.phase = phase
        self.pcoef = self.norm_par()
        


# Create a new class for spectra that returns easy plotting, normalization, parameters, etc.

# In[1]:

class Sp:
    """ 
    Sp: spectrum class object that has all the tools for easy spectra analysis.
    Includes tools for building look up tables
    
    Modified (v1.1): Samuel LeBlanc, NASA Ames, 2015-05-14
                    - added in _init_ for checking if using radiance subset instead of all radiances (rad vs. rads), in matlab loaded values
    Modified (v1.2): Samuel LeBlanc, NASA Ames, 2015-10-01
                    - added subset of wavelengths for normalization
    Modified (v1.3): Samuel LeBlanc, NASA Ames, 2015-11-02
                    - added __getitem__ method for calling like a dict
                    - added method for obtaining the valid ref ranges for ice or liquid, making Sp more general.
    """    
    import numpy as np
    def __init__(self,s,irrad=False):
        import numpy as np
        if 'nm' in s:
            self.wvl = np.array([item for sublist in s['nm'] for item in sublist])
        if 'zenlambda' in s:
            self.wvl = s['zenlambda']
        if 'wvl' in s:
            self.wvl = np.array([item for sublist in s['wvl'] for item in sublist])
        if 'w' in s:
            self.wvl = np.array([item for sublist in s['w'] for item in sublist])
        if self.wvl[0] < 100.0:
            self.wvl = self.wvl*1000.0
        self.iwvls = np.argsort(self.wvl)
        print len(self.iwvls), len(self.wvl)
        self.wvl = np.sort(self.wvl)
        self.wvlsort(s,irrad)
        self.isubwvl = self.wvl_for_norm(self.wvl,wrange=[315.0,940.0])
        self.norm = self.normsp(self.sp,iws=self.isubwvl)
        print self.sp.shape
        if self.sp.ndim > 2:
            self.tau = s['tau']
            self.ref = s['ref']
        else:
            self.utc = s['utc'][self.iset]
            if 'good' in s.keys():
                if isinstance(s['good'],tuple):
                    self.good = s['good'][0]
                else:
                    self.good = s['good']
            else:
                print 'No indexed good values, choosing all times that are greater than 0'
                self.good = np.where(self.utc>0.0)[0]
            if 'Alt' in s:
                self.alt = s['Alt'][self.iset]
            if 'Lat' in s:
                self.lat = s['Lat'][self.iset]
            if 'Lon' in s:
                self.lon = s['Lon'][self.iset]
        if 'sza' in s:
            self.sza = s['sza']
        # initiate with NANs the values that need to be populated
        self.npar = np.nan
        self.par = np.zeros_like(self.sp)*np.nan
        self.parn = np.zeros_like(self.par)*np.nan
        self.meansp = np.zeros_like(self.wvl)*np.nan
        self.stdsp = np.zeros_like(self.wvl)*np.nan
        self.stdpar = np.zeros((16))*np.nan
        self.stdparn = np.zeros((16))*np.nan
        self.sphires = False
        self.parhires = False
        self.hiresirrad = False
        self.irrad = irrad
        self.reflect = np.zeros_like(self.sp)*np.nan
    
    def __getitem__(self,i):
        'Method to call only the variables in the SP class like a dict'
        return self.__dict__get(i)
    
    def params(self):
        " Runs through each spectrum in the sp array to calculate the parameters"
        useezmap = False
        from Sp_parameters import param
        import warnings
        print('Running Parameters')
        if useezmap:
            import ezmap
            from ezmap import map as zmap
        warnings.filterwarnings("ignore",category=RuntimeWarning,module="param")
        warnings.filterwarnings("ignore",category=RuntimeWarning,module="Sp_parameters.normsp")
        sp = self.sp
        wvl = self.wvl
        if sp.ndim == 5:
            w,r,t = np.mgrid[0:2,0:len(self.ref),0:len(self.tau)]
            if useezmap:
                gx = lambda a:param(sp[a[0],:,0,a[1],a[2]],wvl)
                args = ((aw,ar,at) for aw in w.ravel() for ar in r.ravel() for at in t.ravel())
                partemp = zmap(gx,args,progress=True,ncpu=2)
            else:
                applypar = lambda w,r,t:param(sp[w,:,0,r,t],wvl)
                partemp = map(applypar,w.ravel(),r.ravel(),t.ravel())
            par = np.reshape(partemp,[2,len(self.ref),len(self.tau),-1])
            self.npar = par.shape[3]
        elif sp.ndim == 2:
            applypartime = lambda tt:param(sp[tt,:],wvl)
            if useezmap:
                par = np.array(zmap(applypartime,xrange(len(self.utc)),ncpu=2,progress=True))
            else:
                part = map(applypartime,xrange(len(self.utc)))
                par = np.array(part)
            self.npar = par.shape[1]
        else: raise LookupError
        self.par = par
        return
    
    def param_hires(self):
        " Runs through the parameter space and interpolates to create a hires res version, should be run instead of sp_hires "
        from Sp_parameters import param
        from scipy import interpolate
        print('Running parameter hires')
        if np.all(np.isnan(self.par)):
            print 'Run params() before'
            return
        tau = self.tau
        ref = self.ref
        par = self.par
        if self.sp.ndim == 5:
            tau_hires = np.concatenate((np.arange(self.tau[0],1.0,0.1),np.arange(1.0,4.0,0.5),np.arange(4.0,self.tau[-1]+1.0,1.0)))
            ref_hires = np.arange(self.ref[0],self.ref[-1]+1.0)
            print tau_hires.shape
            print ref_hires.shape
            import gc; gc.collect()
            par_hires = np.zeros([2,len(ref_hires),len(tau_hires),self.npar])*np.nan
            startprogress('Running interpolation on params')
            refranges = self.get_refrange()
            for ph in [0,1]:
                for tt in xrange(1,self.tau.size-1):
                    for rr in refranges[ph]:
                        #look for NaNs and remove them by interpolating over the neighbors
                        if np.any(np.isnan(par[ph,rr,tt,:])) or not np.any(par[ph,rr,tt,:]):
                            fs = interpolate.interp1d([tau[tt-1],tau[tt+1]],[par[ph,rr,tt-1,:],par[ph,rr,tt+1,:]],axis=0)
                            par[ph,rr,tt,:] = fs(tau[tt])
                for pp in xrange(self.npar):
                    fx = interpolate.RectBivariateSpline(ref[refranges[ph]],tau,par[ph,refranges[ph],:,pp],kx=1,ky=1)
                    par_hires[ph,:,:,pp] = fx(ref_hires,tau_hires)
                    progress(float(pp+self.npar*ph)/(self.npar*2)*100.0)
            endprogress()
        self.par = par_hires
        self.tausp = tau
        self.refsp = ref
        self.tau = tau_hires
        self.ref = ref_hires
        self.parhires = True
        return
    
    def get_refrange(self):
        """
        Simple program that returns the ref ranges with valid first parameter for two phases
        """
        if np.all(np.isnan(self.par)):
            print 'Run params() before'
            return
        ice_r = [r for r in xrange(len(self.ref)) if ~ np.isnan(self.par[1,r,10,0])]
        liq_r = [r for r in xrange(len(self.ref)) if ~ np.isnan(self.par[0,r,10,0])]
        return (liq_r,ice_r)
    
    def sp_hires(self,doirrad=False):
        """
        Interpolate the modeled sp to a finer resolution in tau and ref. 
        Only works on modeled data with 5 dimensions
        When doirrad is set to True, run the interpolation of reflectance at z=1
        """
        #from Sp_parameters import startprogress, progress, endprogress
        from scipy import interpolate
        from Sp_parameters import nanmasked
        import numpy as np
        import warnings
        if doirrad:
            if not self.irrad:
                warnings.warn("Irradiance was not set when first intializing the Sp object")
                return
            sp = self.sp_irrup/self.sp_irrdn
            z = 1
        else:
            sp = self.sp
            z = 0
        wvl = self.wvl
        if self.parhires:
            tau = self.tausp
            ref = self.refsp
        elif self.irradhires and not self.sphires:
            tau = self.tausp
            ref = self.refsp
        else:
            tau = self.tau
            ref = self.ref
        if doirrad and self.sphires:
            tau = self.tausp
            ref = self.refsp
        if sp.ndim !=5:
            warnings.warn('This method cannot be applied to this object, not enough dimensions')
            return
        if self.sphires and not doirrad:
            print('sp is already hires')
            return
        tau_hires = np.concatenate((np.arange(tau[0],1.0,0.1),np.arange(1.0,4.0,0.5),np.arange(4.0,tau[-1]+1.0,1.0)))
        ref_hires = np.arange(ref[0],ref[-1]+1.0)
        print tau_hires.shape 
        print ref_hires.shape
        import gc; gc.collect()
        sp_hires = np.zeros([2,len(wvl),2,len(ref_hires),len(tau_hires)])*np.nan
        startprogress('Running interpolation')
        refranges = (range(0,23),range(3,35))
        for ph in [0,1]:
            for tt in xrange(1,len(tau)-1):
                for rr in refranges[ph]:
                    #look for NaNs and remove them by interpolating over the neighbors
                    if np.any(np.isnan(sp[ph,:,z,rr,tt])) or not np.any(sp[ph,:,z,rr,tt]):
                        fs = interpolate.interp1d([tau[tt-1],tau[tt+1]],[sp[ph,:,z,rr,tt-1],sp[ph,:,z,rr,tt+1]],axis=0)
                        sp[ph,:,z,rr,tt] = fs(tau[tt])
            for w in xrange(len(wvl)):
                fx = interpolate.RectBivariateSpline(ref[refranges[ph]],tau,sp[ph,w,z,refranges[ph],:],kx=1,ky=1)
                sp_hires[ph,w,z,:,:] = fx(ref_hires,tau_hires)
                progress((w+len(wvl)*ph)/(len(wvl)*2)*100.0)
        endprogress()
        if doirrad:
            print('Creating reflect spectra with tau, and ref of the new high resolution values')
            self.reflect = sp_hires
            self.hiresirrad = True
            self.tausp = tau
            self.refsp = ref
        else:
            print('Overwriting the current sp, tau, and ref with the new high resolution values')
            self.sp = sp_hires
            self.sphires = True
            self.tausp = tau
            self.refsp = ref
        self.tau = tau_hires
        self.ref = ref_hires  
    
    def norm_par(self,pcoef=None,std=False,vartitle=None):
        """ 
        Normalize the parameters, if no keyword set, returns the normalized parameter values in self.parn
        if the keywords are not set, returns the normalized parameter values in self.parn and returns
        the pcoef dictionary containing the coefficients and additive values used for normalizing each parameter.
        Applies the normalization to std if set to true, and requires the pceof to be set and creates stdparn, overrides vartitle
        Saves to vartitle if vartitle is set.
        
        """
        import numpy as np
        npar = self.npar
        par = self.par
        parn = np.zeros_like(self.par)
        if pcoef is not None:
            # for defined coefficients
            if not std:
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
                    import warnings; warnings.warn('There is a problem with the dimensions of the sp')
                    return
                if vartitle is not None:
                    setattr(self,vartitle,parn)
                    return
                else:
                    self.parn = parn
                    return
            else:
                self.stdparn = self.stdpar*pcoef['coef']-pcoef['add']
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
        
    def wvlsort(self,s,irrad):
        """
        Function to sort spectra along the wavelength axis
        select processing depending on the keys in the original dictionary input
        
        Creates the sp element in self, and iset for measurements
        """
        iwvls = self.iwvls
        if sorted(iwvls) is iwvls:
            print '*** wvls are already sorted, there may be a problem! ***'
        if 'sp' in s:
            print 'in sp'
            print s['sp'].shape
            ui = [i for i in range(s['sp'].ndim) if s['sp'].shape[i] == len(self.wvl)]
            if 1 in ui:
                sp = s['sp'][:,iwvls,:,:,:]
            else: 
                raise LookupError
        if 'rads' in s:
            print 'in rads'
            print s['rads'].shape, s['rads'].ndim, len(iwvls)
            ui = [i for i in range(s['rads'].ndim) if s['rads'].shape[i] == len(self.wvl)]
            #print ui
            if 1 in ui:
                print '1 in ui'
                sp = s['rads'][:,iwvls]
            else: 
                print 'not 1 in ui'
                sp = s['rads'][iwvls,:]
            if 'iset' in s:
                if s['iset'].ndim>1:
                    self.iset = s['iset'][:,0]
                else:
                    self.iset = s['iset']
            else:
                print '** Problem, rads present (radiance subset), but not the subset integers **'
        elif 'rad' in s:
            print 'in rad'
            print s['rad'].shape, s['rad'].ndim, len(iwvls)
            ui = [i for i in range(s['rad'].ndim) if s['rad'].shape[i] == len(self.wvl)]
            #print ui
            if 1 in ui:
                print '1 in ui'
                sp = s['rad'][:,iwvls]
            else: 
                print 'not 1 in ui'
                sp = s['rad'][iwvls,:]
            self.iset = np.where(s['rad'][:,0])[0]
        if irrad:
            print 'in irrad'
            ui = [i for i in range(s['sp_irrdn'].ndim) if s['sp_irrdn'].shape[i] == len(self.wvl)]
            if 1 in ui:
                self.sp_irrdn = s['sp_irrdn'][:,iwvls,:,:,:]
                self.sp_irrup = s['sp_irrup'][:,iwvls,:,:,:]
            else: 
                raise LookupError
        self.sp = sp
        
    def wvl_for_norm(self,wvl,wrange=[315.0,940.0]):
        " Function that gets the indices for the wavelengths to be used in normalization "
        import numpy as np
        return np.where((wvl>=wrange[0])&(wvl<=wrange[1]))[0]
    
    def normsp(self,sp,iws=None):
        " Function to return the normalized spectra list"
        import numpy as np
        import warnings
        warnings.simplefilter("ignore",RuntimeWarning)
        try:
            if not iws:
                iws = np.where(self.wvl)[0]
        except ValueError:
            pass
        idxwvl = [i for i in range(sp.ndim) if sp.shape[i] == len(self.wvl)]
        if sp.ndim == 2:
            norm = sp/np.nanmax(sp[:,iws],axis=idxwvl[0])[:,None]
        elif sp.ndim == 5:
            norm = sp/np.nanmax(sp[:,iws,:,:,:],axis=idxwvl[0])[:,None,:,:,:]
        return norm
    
    def reshape_lut(self,phase=None):
        " Function that reformats the lut such that only liquid or ice clouds values exist "
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
            self.liq = lut(par,self.tau,ref,phase=phase)
        elif phase is 'ice':
            ref = self.ref[self.ref >= 10]
            if self.par.ndim == 4:
                par = self.par[1,self.ref>=10,:,:]
            elif self.par.ndim == 5:
                par = np.rollaxis(self.par[1,:,0,self.ref>=10,:],0,3)
            else:
                warning('Problem with par dimensions')
            self.ice = lut(par,self.tau,ref,phase=phase)
        return
    
    def mean(self):
        " function that returns the mean spectra."
        sp = self.sp
        notaxis = tuple(np.where(np.array(self.sp.shape) != self.wvl.shape[0])[0])
        meansp = np.nanmean(sp,axis=notaxis)
        self.meansp = meansp
        return meansp
    
    def std(self,sp=None):
        " function that returns the standard deviation spectra."
        if sp is None:
            sp = self.sp
            savetoself = True
        else:
            savetoself = False
        notaxis = tuple(np.where(np.array(self.sp.shape) != self.wvl.shape[0])[0])
        stdsp = np.nanstd(sp,axis=notaxis)
        if savetoself:
            self.stdsp = stdsp
        return stdsp 
    
    def build_std(self):
        """
        Function that creates a set of uncertainty for each parameter.
        Currently just uses white noise (gaussian)
        """
        from Sp_parameters import param
        meansp = self.mean()
        stdsp = self.std()
        num_noise = 200
        noise = np.random.normal(1,0.005,(num_noise,self.wvl.size)) # add 0.5% variance to signal at all wavelengths
        # should be at every sp in utc, but for now, use mean sp
        sp_arr = meansp*noise
        #import code; code.interact(local=locals())
        par_noisy = np.array(map(lambda tt:param(sp_arr[tt,:],self.wvl),xrange(num_noise)))
        notaxis = tuple(np.where(par_noisy.shape != self.npar)[0])
        stdpar = np.nanstd(par_noisy,axis=notaxis)
        self.stdpar = stdpar
        return stdpar
            


# ## Create some plotting functions for use with the Sp class

# In[1]:

def plt_zenrad(meas):
    """
    Purpose:
        Plot each zzenith radiance spectra on the same scale. 
    Input:
        meas : Sp_parameters.Sp object with params already run
    Output: 
        matplotlib fig value
    Keywords:
        None
    Dependencies:
        - matplotlib
        - mpltools for color
    Modification History:
        Writtten: Samuel LeBlanc, 2015-10-27, NASA Ames, CA
    """
    import matploltib.pyplot as plt
    from mpltools import color
    fig = plt.figure()
    color.cycle_cmap(len(meas.utc),cmap=plt.cm.gist_ncar,ax=plt.gca())
    for i in range(len(meas.utc)):
        plt.plot(meas.wvl,meas.sp[i,:]/1000.0)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Radiance [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')
    plt.title('All radiance spectra')
    scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar)
    scalarmap.set_array(meas.utc)
    cba = plt.colorbar(scalarmap)
    cba.set_label('UTC [h]')
    return fig


# In[2]:

def plt_norm_zenrad(meas):
    """
    Purpose:
        Plot each normalized zenith radiance spectra on the same scale. 
    Input:
        meas : Sp_parameters.Sp object with params already run
    Output: 
        matplotlib fig value
    Keywords:
        None
    Dependencies:
        - matplotlib
        - mpltools for color
    Modification History:
        Writtten: Samuel LeBlanc, 2015-10-27, NASA Ames, CA
    """
    import matploltib.pyplot as plt
    from mpltools import color
    fig = plt.figure()
    color.cycle_cmap(len(meas.utc),cmap=plt.cm.gist_ncar,ax=plt.gca())
    for i in range(len(meas.utc)):
        plt.plot(meas.wvl,meas.norm[i,:])
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Normalized Radiance')
    plt.ylim([0,1])
    plt.title('All normalized radiance spectra')
    scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar)
    scalarmap.set_array(meas.utc)
    cba = plt.colorbar(scalarmap)
    cba.set_label('UTC [h]')
    return fig


# In[3]:

def curtain_zenrad(meas,utc=True):
    """
    Purpose:
     Create a figure of a curtain containing each zenith radiance spectra on the same color scale along the time axis in y. 
    Input:
        meas : Sp_parameters.Sp object with params already run
    Output: 
        matplotlib fig value
    Keywords:
        utc: (default true) if set plots the curtain as a function of UTC time. 
             If false, plots as a function of measurements number
    Dependencies:
        - matplotlib
        - mpltools for color
    Modification History:
        Writtten: Samuel LeBlanc, 2015-10-27, NASA Ames, CA
    """
    import matploltib.pyplot as plt
    from mpltools import color
    fig,ax = plt.subplots(1,1,figsize=(8,14))
    if utc:
        pco = ax.pcolorfast(meas.wvl,meas.utc,meas.sp[:-1,:-1]/1000.0,cmap='gist_ncar',vmin=0,vmax=0.8)
        ax.set_ylabel('UTC [h]')
    else:
        pco = ax.pcolorfast(meas.wvl,np.where(meas.utc)[0],meas.sp[:-1,:-1]/1000.0,cmap='gist_ncar',vmin=0,vmax=0.8)
        ax.set_ylabel('Measurement number')
    cba = plt.colorbar(pco)
    cba.set_label('Radiance [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')
    ax.set_xlabel('Wavelength [nm]')
    ax.set_title('All radiance spectra [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')
    return fig


# In[5]:

def curtain_norm_zenrad_num(meas, utc=True):
    """
    Purpose:
     Create a figure of a curtain containing each *Normalized* zenith radiance spectra on the same color scale along the time axis in y. 
    Input:
        meas : Sp_parameters.Sp object with params already run
    Output: 
        matplotlib fig value
    Keywords:
         utc: (default true) if set plots the curtain as a function of UTC time. 
             If false, plots as a function of measurements number
    Dependencies:
        - matplotlib
        - mpltools for color
    Modification History:
        Writtten: Samuel LeBlanc, 2015-10-27, NASA Ames, CA
    """
    import matploltib.pyplot as plt
    from mpltools import color
    fig,ax = plt.subplots(1,1,figsize=(8,14))
    if utc:
        pco = ax.pcolorfast(meas.wvl,meas.utc,meas.norm[:-1,:-1],cmap='gist_ncar',vmin=0,vmax=1.0)
        ax.set_ylabel('UTC [h]')
    else:
        pco = ax.pcolorfast(meas.wvl,np.where(meas.utc)[0],meas.norm[:-1,:-1],cmap='gist_ncar',vmin=0,vmax=1.0)
        ax.set_ylabel('Measurement number')
    cba = plt.colorbar(pco)
    cba.set_label('Normalized Radiance')
    ax.set_xlabel('Wavelength [nm]')
    ax.set_title('All radiance spectra [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')
    return fig

