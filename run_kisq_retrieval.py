# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Functions to run the ki^2 retrieval based on the parameters calculated in the Sp class, from modeled and measured source

# <codecell>

def assure_input(sp):
    " A function that checks the input and runs all the required functions to generate parameters"
    import numpy as np
    if not(sp.hires) and hasattr(sp,'tau'):
        print('Running the interpolation to hi-res')
        sp.sp_hires()
    if not(hasattr(sp,'par')):
        print('Running the params on current class')
        sp.params()        
    elif (np.isnan(sp.par)).all():
        print('Running the params on current class')
        sp.params() 

# <codecell>

def phase(parn,model,stdparn):
    """
    Function to determine the phase of the measurement by checking if the 
    measured parameters are within the possible model values.
    Lists the default parameter boundaries. 
    **** to be updated ****
    """
    import numpy as np
    ph = 2 # set default to uncertain
    try:
        if parn[1]-stdparn[1] > np.nanmax(model.parn[1,:,:,1]):
            ph = 0
    except IndexError:
        import pdb; pdb.set_trace()
    if parn[1]+stdparn[1] < np.nanmin(model.parn[0,:,:,1]):
        ph = 1
    if parn[0]+stdparn[0] < np.nanmin(model.parn[0,:,:,0]):
        ph = 1
    if parn[8]-stdparn[8] > np.nanmax(model.parn[0,:,:,8]):
        ph = 1
    if parn[9]-stdparn[9] > np.nanmax(model.parn[0,:,:,9]):
        ph = 1
    return ph

# <codecell>

def run_retrieval(meas,model,subp=range(15)):
    """ 
    A function that uses the Sp class to run through each utc point in meas, 
    and finds the closest model values, 
    uses only the subset of parameters defined by subp, which defaults to first 15 parameters
    """
    import Sp_parameters as Sp
    import numpy as np
    import run_kisq_retrieval as rk
    rk.assure_input(meas)
    rk.assure_input(model)
    # normalize the parameters
    pcoef = model.norm_par()
    meas.norm_par(pcoef=pcoef)
    
    # build the measurement error
    meas.build_std() #in future must reference a set of measurements to establishe the uncertainty, right now only white noise at 0.5%
    meas.stdparn = meas.stdpar/pcoef['coef'] #for creating the normalized version of the stdpar
    
    # build the reshaped model for ice and liquid
    model.reshape_lut(phase='liq')
    model.reshape_lut(phase='ice')
    
    # build the normalized measurement parameters based on liquid or ice values
    meas.norm_par(pcoef=model.ice.pcoef,vartitle='parn_ice')
    meas.norm_par(pcoef=model.liq.pcoef,vartitle='parn_liq')
    
    # build the weights from the stdparn of measurements
    wg = np.sqrt(meas.stdpar)/pcoef['coef']
    
    #start loop over measurement
    Sp.startprogress('Retrieval progress over times')
    ph = np.zeros_like(meas.utc)*np.nan
    ki = np.zeros_like(meas.utc)*np.nan
    ta = np.zeros_like(meas.utc)*np.nan
    re = np.zeros_like(meas.utc)*np.nan
    #ki_2ph = np.zeros_like(meas.utc) #kisq with 2 phase
    for tt in meas.good:
        if np.all(np.isnan(meas.parn[tt,subp])):
            continue
        #first get the phase in first method
        #import pdb; pdb.set_trace()
        ph[tt] = rk.phase(meas.parn[tt,:].ravel(),model,meas.stdparn)
        if ph[tt] == 2: # undecided, must do the kisq
            ki_2ph = np.nansum(wg[subp]*(meas.parn[tt,subp]-model.parn[:,:,:,subp])**2,axis=3)
            ki_minin = np.unravel_index(np.nanargmin(ki_2ph),ki_2ph.shape)
            ph[tt] = ki_minin[0]
        if ph[tt] == 0: #liquid
            goodpi = [i for i in subp if (meas.parn_liq[tt,i]+meas.stdparn[i]>=0) and (meas.parn_liq[tt,i]-meas.stdparn[i]<=1)]
            if len(goodpi) < 4:
                continue
            ki_arr = np.nansum(wg[goodpi]*(meas.parn_liq[tt,goodpi]-model.liq.parn[:,:,goodpi])**2,axis=2)
            ki[tt] = np.nanmin(ki_arr)
            #print ki[tt]
            ki_minin = np.unravel_index(np.nanargmin(ki_arr),ki_arr.shape)
            (ta[tt],re[tt]) = (model.liq.tau[ki_minin[1]],model.liq.ref[ki_minin[0]])
        elif ph[tt] == 1: #ice
            goodpi = [i for i in subp if (meas.parn_ice[tt,i]+meas.stdparn[i]>=0) and (meas.parn_ice[tt,i]-meas.stdparn[i]<=1)]
     #       print goodpi
            if len(goodpi) < 4:
                continue
            ki_arr = np.nansum(wg[goodpi]*(meas.parn_ice[tt,goodpi]-model.ice.parn[:,:,goodpi])**2,axis=2)
            ki[tt] = np.nanmin(ki_arr)
            ki_minin = np.unravel_index(np.nanargmin(ki_arr),ki_arr.shape)
            (ta[tt],re[tt]) = (model.ice.tau[ki_minin[1]],model.ice.ref[ki_minin[0]])
            #print ki[tt]
        else:
            warning('Problem with phase!')
        Sp.progress(float(tt)/len(meas.utc)*100.0)
    Sp.endprogress()
    
    #save the file
    return (ta,re,ph,ki)
    

# <codecell>


