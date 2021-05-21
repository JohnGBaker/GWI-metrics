#This is a place for background/foreground models
import numpy as np

#functions for PSD
def compute_background_PSD(f,model):
    '''
    This function computes a set of backgrounds/foregrounds as noises to add to the instrument PSD.
    For a nontrival result, the 'background' key should be specified in the model dict
    The allowed formats for its value are:
     'background':func
     'background':[func,parsdict]
     'background':[[func1,parsdict1],func2,...]
    BG/FG model functions are expected to meet the folowing interface:
    The interface function arguments are:
                 f  -  a frequency (or array) 
       instrParams  -  the general instrument params dict
      nativeParams  -  dict with any additional params specific to this function 
    The return value is the PSD in fractional frequency units.
    '''
    S_bg=0
    if "background" in model:
        bgspec=model['background']
        if not isinstance(bgspec,list) or ( len(bgspec)==2 and isinstance(bgspec[1],dict) ): bgspec=[bgspec] #allow simpler specification of just one bg func
        for bg in bgspec:
            if isinstance(bg,list):
                S_bg+=bg[0](f,model,bg[1])
            else:
                S_bg+=bg(f,model)
    return S_bg

def GBF_model_TN(f,instrParams,nativeParams=None):
    '''
    Model for the unresolved galactic binary foreground from the LISA Sensitivity Tech Note.
    This conforms to a general interface for background/foreground models
    The interface argumentss are:
                 f  -  a frequency (or array) 
       instrParams  -  the general instrument params dict
      nativeParams  -  dict with any additional params specific to this function 
    For the TechNote function, there are no additional params.
    The function returns the PSD of the foreground (in fractional frequency units).
    '''
    #set constants
    A=1.14e-44
    alpha=1.8
    f2=0.31 
    a1=-0.25
    b1=-2.7
    ak=-0.27
    bk=-2.47
    #import
    Tobs=instrParams.get("SciDuration",4.0) #default to 4yrs observation 
    #derive
    f1=Tobs**a1*10**b1
    fk=Tobs**ak*10**bk
    SGal=  A * f**(-7/3) * np.exp(-(f/f1)**alpha) * (1+np.tanh(-(f-fk)/f2))/2
    return SGal

#reference bg defintion example
LISA_Galaxy={'label':"galaxy",'spec':GBF_model_TN}

#convenience function to add background info to an instrument model based on definition (see example above)
def add2model(model,bg=LISA_Galaxy):
    result=model.copy()
    result['label']+=' + '+bg['label']
    result['background']=bg['spec']
    return result


