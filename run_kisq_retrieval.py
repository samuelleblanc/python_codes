# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Functions to run the ki^2 retrieval based on the parameters calculated in the Sp class, from modeled and measured source

# <codecell>

def run_retrieval(meas,model):
    """ A function that uses the Sp class to run through each utc point in meas, and finds the closest model values """
    import Sp_parameters as Sp
    assure_input(meas)
    assure_input(model)
    # normalize the parameters
    pcoef = model.par_norm()
    pcoef = meas.par_norm(pcoef=pcoef)
    
    # build the model and measurement error
    # model.std()
    # model.std.par_norm(pcoef=pcoef)
    
    # build the reshaped model for ice and liquid
    model.reshape_lut(phase='liq')
    model.reshape_lut(phase='ice')
    pcoef_liq = model.liq.par_norm()
    pcoef_ice = model.ice.par_norm()
    # meas.std()
    # meas.std.par_norm(pcoef=pcoef)
    
    # build the weights
    wg,wg_liq,wg_ice = par_weights(meas,model)
    #start loop over measurement
    Sp.startprogress('Progress over times')
    ph = np.zeros_like(meas.utc)
    ki = np.zeros_like(meas.utc)
    ta = np.zeros_like(meas.utc)
    re = np.zeros_like(meas.utc)
    for tt in xrange(len(meas.utc)):
        #first get the phase
        ph[tt] = kisq_phase(meas.parn[tt,:],model,wg[tt,:])
        if ph[tt] == 0: #liquid
            goodpi = reject_outofbounds(meas.parn[tt,:],model.liq.lut)
            ki[tt],ta[tt],re[tt] = kisq(meas.parn[tt,:],model.liq,wg,goodpi)
        elif ph[tt] == 1: #ice
            goodpi = reject_outofbounds(meas.parn[tt,:],model.ice.lut)
            ki[tt],ta[tt],re[tt] = kisq(meas.parn[tt,:],model.ice,wg,goodpi)
        else:
            warning('Problem with phase!')
        Sp.progress(tt/len(meas.utc))
        
    Sp.endprogress()
    
    #save the file
    
    

# <codecell>

def assure_input(sp):
    " A function that checks the input and runs all the required functions to generate parameters"
    if sp.sp.ndim == 5:
        print('Running the interpolation to hi-res')
        sp.sp_hires()
    if not(hasattr(sp,'par')):
        print('Running the params on current class')
        sp.params()        

# <codecell>

def par_weights(meas,model):
    """
    A function to combine the standard deviation uncertainty in measurement
    and model to build a set of weigths for the ki^2 regression.  
    """
    pstd = np.sqrt(meas.std.parn**2.0+model.std.parn**2.0)
    return 1./pstd

# <codecell>

def kisq_phase(parn,model,weights,goodpi):
    """
    Function that first gets the phase from the unique possibilities, 
    and then calculates the ki^2 value with phase.
    """
    ki = -999.9
    return ki

# <codecell>

def reject_outofbounds(parn,lut):
    """
    functions that builds a mask of bad parameters that are deemed out of bounds of the model luts
    """
    maskgood = [1,0]
    return maskgood

# <codecell>

def kisq(parn,model,wg):
    """ 
    function to calculate the ki squared value for a single set of 
    measured parameter compared to model lut. Returns the tau min and ref min
    """
    taumin = -999
    refmin = -999
    return taumin,refmin

