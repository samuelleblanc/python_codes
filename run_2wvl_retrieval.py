# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Module that runs the 2-wavelength reflectance retrieval of cloud properties

# <codecell>

def run_2wvl_retrieval(meas,model,wvls=[500.0,1600.0],sza_norm=True):
    """ 
    Name:

        run_2wvl_retrieval
    
    Purpose:

        Retrieves cloud optical properties from reflectance 2 wavelength method (Nakajima & King, 1990)
        A function that uses the Sp class to run through each utc point in meas, 
        and finds the closest model values.
        If sza is present then the reflectance is calculated with sza_factor (1/cos(sza))
    
    Calling Sequence:

        (ta,re,ki) = run_2wvl_retrieval(meas,model,wvls=[500.0,1600.0],sza_norm=True)
    
    Input: 
  
        meas : object with keys of utc (time in fractional hours), Rvis (reflectance in vis), Rnir (reflectance in nir)
        model: Sp object with irradiance values set
     
    Output:

        ta: 2d array of cloud optical thickness, for [:,0] = liquid, for [:,1] = ice
        re: 2d array of cloud particle effective radius, for [:,0] = liquid, for [:,1] = ice
        ki: 2d array of minimal ki^2 values, for [:,0] = liquid, for [:,1] = ice
    
    Keywords: 

        wvls: (optional) array of two values with the vis wavlength in first, and nir wavelength in second [nm]
        sza_norm: (optional) if set to True (default) then searches for sza in meas and model. 
                 If there, normalizes meas and model reflectances by 1/cos(sza) 
    
    Dependencies:

        Sp_parameters
        numpy
        warnings
        gc: for clearing the garbage
        run_2wvl_retrieval (this file)
        pdb: for debugging when there is an error
        math
    
    Required files:
   
        none
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2014-12-08, NASA Ames

    """
    import Sp_parameters as Sp
    import numpy as np
    import run_2wvl_retrieval as rw
    import warnings
    import math
    if not model.irrad:
        warnings.warn('model lut does not have irradiances set! please rerun with irradiances')
        return
    if sza_norm & hasattr(meas,'sza'):
        meas.sza_factor = 1.0/np.cos(np.radians(meas.sza))
    else:
        meas.sza_factor = 1.0
    meas.Rvis = meas.Rvis*meas.sza_factor
    meas.Rnir = meas.Rnir*meas.sza_factor
    model.lut_2wvl = rw.build_lut_2wvl(model,wvls,sza_norm)    
    #start loop over measurement
    Sp.startprogress('Retrieval progress over times')
    ki = np.ndarray((len(meas.utc),2))*np.nan
    ta = np.ndarray((len(meas.utc),2))*np.nan
    re = np.ndarray((len(meas.utc),2))*np.nan
    #make good filter
    meas.good = np.where((np.isfinite(meas.Rvis)) & (meas.Rvis > 0) & (np.isfinite(meas.Rnir)) & (meas.Rnir > 0))[0]
    for tt in meas.good:
        for ph in [0,1]:
            ki_ref_tau = ((meas.Rvis[tt]-model.lut_2wvl[ph,0,:,:])/meas.Rvis[tt])**2+((meas.Rnir[tt]-model.lut_2wvl[ph,1,:,:])/meas.Rnir[tt])**2
            try:
                ki_minin = np.unravel_index(np.nanargmin(ki_ref_tau),ki_ref_tau.shape)
            except:
                import pdb; pdb.set_trace()
            ki[tt,ph] = np.nanmin(ki_ref_tau)
            (ta[tt,ph],re[tt,ph]) = (model.tau[ki_minin[1]],model.ref[ki_minin[0]])
         #   if (ph == 1) & (meas.utc[tt] > 18.4):
         #       import pdb; pdb.set_trace()
        Sp.progress(float(tt)/len(meas.utc)*100.0)
    Sp.endprogress()
    return (ta,re,ki)

# <codecell>

def build_lut_2wvl(model,wvls,sza_norm):
    " function that builds the 2 wavelength Reflectance lut from the model Sp object (lut) "
    import numpy as np
    import math
    from Sp_parameters import find_closest
    iw = find_closest(model.wvl,np.array(wvls))
    if hasattr(model,'sza') & sza_norm:
        sza_factor = 1.0/math.cos(math.radians(model.sza))
    else:
        sza_factor = 1.0
    return np.reshape(model.reflect[:,iw,1,:,:],(2,2,model.ref.size,model.tau.size))*sza_factor #ph, wvl, z, ref, tau

