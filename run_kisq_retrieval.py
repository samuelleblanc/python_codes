
# coding: utf-8

## Functions to run the ki^2 retrieval based on the parameters calculated in the Sp class, from modeled and measured source

# In[ ]:

def run_retrieval(meas,model):
    """ A function that uses the Sp class to run through each utc point in meas, and finds the closest model values """
    assure_input(meas)
    assure_input(model)
    # normalize the parameters
    pcoef = model.par_norm()
    pcoef = meas.par_norm(pcoef=pcoef)
    
    #start loop over measurement
    


# In[ ]:

def assure_input(sp):
    " A function that checks the input and runs all the required functions to generate parameters"
    if not(hasattr(sp,'par')):
        print('Running the params on current class')
        sp.params()
        

