#These are functions used for GW mission sensitivity calculations and other figures of merit related to GW iimaging or sky localization

import numpy as np
import constants
import PhenomWaveform_nonspinning as chirp
import subsystems
import background

#Probably adapt more for GW Imager concepts 
def PSD_noise_components(fr, model):
    '''
    Make PSD Noise Components
    We follow LISA in computing the noise. In that case the noise is more-or-less directly derived from 
    - a measure of the acceleration noise in $m/s^2/\sqrt(Hz)$ with some additional reddening factors
    - a measure of optical measurement system noise in $m/\sqrt(Hz)$
    The noise PSD is reported in fractional frequency units
    '''
        
    c=constants.c
    
    ### Acceleration noise
    if 'sqSacc_ASD' in model:
        sqSacc_ASD = model.get('sqSacc_ASD') 
        Sa_a = subsystems.F_Noise_PSD(fr,sqSacc_ASD,True)**2
    else:
        Sa_a = model.get('sqSacc_func')(fr,model) #Can provide a func here instesad of a value
    Sa_d = Sa_a*(2.*np.pi*fr)**(-4.)
    Sa_nu = Sa_d*(2.0*np.pi*fr/c)**2

    ### Optical Metrology System
    if 'sqSoms_ASD' in model:
        sqSoms_ASD = model.get('sqSoms_ASD')
        Soms_d = subsystems.F_Noise_PSD(fr,sqSoms_ASD,True)**2
    else:
        Soms_d = model.get('sqSoms_func')(fr, model) #Can provide a func based Jeff's calculations
    Soms_nu = Soms_d*(2.0*np.pi*fr/c)**2
    
    return [Sa_nu, Soms_nu]


#Component for the LISA TN computation:
def AvFXp2_approx(fr,L):
    return 16.*(3./20.)*(1./(1.+0.6*(2.*np.pi*fr*L)**2))

#Components for the Larson computation
gamma = 60*180/np.pi
sg = np.sin(gamma)
cg = np.cos(gamma)
def fIte(eps,th1,u,cu,su,sg,cg):
    cth1 = np.cos(th1)
    sth1 = np.sin(th1)
    cth2 = cg*cth1+sg*sth1*np.cos(eps)
    sa = sg*np.sin(eps)/np.sqrt(1-cth2*cth2)
    eta = (cu-np.cos(u*cth1))*(cu-np.cos(u*cth2))*cth1*cth2+(su-cth1*np.sin(u*cth1))*(su-cth2*np.sin(u*cth2))
    #print(eps,th1,cth1,sth1,sa,eta,sa)
    return np.sin(th1)*(1-2*sa*sa)*eta

def Tarm(f,L):
    from scipy.integrate import dblquad
    w = 2*np.pi*f
    u = w*L
    cu = np.cos(u)
    su = np.sin(u)
    
    
    Iet = dblquad(fIte, 0, 2*np.pi, lambda x: 0, lambda x: np.pi, args=(u,cu,su,sg,cg))
    #print(Iet)
    return ((1+cu*cu)*(1./3.-2./(u*u))+su*su+4.*su*cu/(u*u*u))/(u*u) - Iet[0]/(4.*np.pi)

    
#Compute sensitivty curve
def makeSensitivity(fr, model,style='TN'):
    '''
    Using the semi-analytical average response, the semi-analytical sensitivity for TDI X 4 links is:
    $$
    S_{h,X} =  \frac{ S_{OMS} + \left( 3 + \cos \left( \frac{2 \omega L}{c} \right)  \right)  S_{acc} }
    { \left( {\omega L \over c} \right)^2 \ R_{\Sigma}^2(f, L) }
    $$
    $$
    S_{h,X} 
    =  \frac{S_{n,X_{2.0}}}
    {<R_{L, X_{2.0}}(f)>} 
    = \frac{ 64 \sin^2 \left( \omega L \right) \sin^2 \left(2 \omega L \right) (S_{OMS} + \left( 3 + \cos \left( \frac{2 \omega L}{c} \right)  \right)  S_{acc}) }
    { (4\omega L)^2 \sin^2{(\omega L)} (2 \sin{(2\omega L)})^2 <(F^{+}_{X})^2> }
    = \frac{ S_{OMS} + \left( 3 + \cos \left( \frac{2 \omega L}{c} \right)  \right)  S_{acc} }
    { \left( \omega L \right)^2  <(F^{+}_{X})^2>  }
    '''
    [Sa_nu,Soms_nu] = PSD_noise_components(fr, model)
    L=model.get('Lconst')
    N=model.get('Nindep')
    c=constants.c
    phiL = 2*np.pi*fr*L/c
    if style=='TN':
        AvFXp2 = AvFXp2_approx(fr,L/c)
        S_hX = (Soms_nu + Sa_nu*(3.+np.cos(2*phiL)) ) / (phiL**2 * AvFXp2/4**2)#LISA TN
    elif style=='Larson': #This is very slow!
        yTarm = np.zeros(len(fr))
        for i in range(len(fr)):
            yTarm[i] = Tarm(fr[i],L)
        #renormalize as in TN
        yTarmN = yTarm/yTarm[0]
        yTarmRN = yTarmN * (np.sin(np.pi/3.)**2/5)
        S_hX = (Soms_nu + Sa_nu*(3+np.cos(2*phiL)) ) / (phiL**2 * yTarmRN)
    S_h = S_hX / N
    #add optional background
    S_h += background.compute_background_PSD(fr,model)
    return S_h

# Get the baseline for a particualr model
def getBaseline(model,t=0,tstart=0):
    '''
    Compute the effective baseline for a concept depending on the parameters in the model. 
    If a multi-constellation model, the baseline is always defined by the DSep parameter (e.g Dsep ~= 0)
    For single constellations, the baseline starts as the detector size and increases with time according to the orbital parameters
    The parameter t defines the time(s) at which to evaluate the orbit.
    The parameter tstart is the time at which to start including orbital effects (e.g. after some SNR threshold)
    '''
    
    # minimum baseline is the detector size
    B = np.zeros_like(t) + model.get('Lconst')/constants.c
    
    # if you have multiple constellations, use the Dsep parameter (which is in what units?)
    if 'Dsep' in model:
        Dsep = model.get('Dsep')
        if Dsep > 0:
            return Dsep*constants.AU + B
        
    # otherwise start using your orbit from tstart
    istart = np.argmin(np.abs(t-tstart))
    theta = 2*np.pi*np.clip((t[istart:]-t[istart])/(model.get('Torbit')*constants.year),0,0.5)
    Borbit = model.get('Rorbit')*constants.AU*np.sqrt(2*(1-np.cos(theta)))
    B[istart:]= B[istart:]+Borbit
    return B
    
        
#Make SNR for continuous-wave source
def getCWsnr(f0,h0,T,model,style='TN'):
    '''
    Compute the SNR for a monochromatic GW source based on the source frequency f0, source ampltiude h0, obsertvation time T, and an instrument model.
    
    The calculation follows (83) from the LISA-LCST-SGS-TN-001 to compute the inclinaiton, polarization, and sky-position averaged SNR for a monochomatic source:
    $$
    \left<SNR^2\right>_{\iota,\psi,sky} = 10 \frac{\left(\frac{2}{5}h_0\right)^2T}{S_h\left(f_0\right)}
    $$
    
    where the instrument sensitivity is computed from the provided model and the makeSensitivity method
    '''
    
    
    # compute sensitivity at the GW frequency from the model
    S_h_f0 = makeSensitivity(f0,model,style)
    
    # apply equation (83)
    rho2 = 10*(((2./5.)*h0)**2)*T/S_h_f0 
    
    # return SNR
    return np.sqrt(rho2)

#Make SNR for continuous-wave source
def getChirpSNR(mtot,eta,dl,model,tstart=-constants.year,Npts = 1000,style='TN'):
    '''
    Compute the SNR for a chirping source
    '''
    # set up time vector
    # stop time is when we get to merger frequency / 3 to avoid PN blow-up
    tstop = chirp.tFromF(0.3*chirp.getFmerge(mtot,eta),mtot,eta)
    tvals = -np.flip(np.logspace(np.log10(-tstop),np.log10(-tstart),Npts))
    
    # get the corresponding frequency vector
    fvals = chirp.fFromT(tvals,mtot,eta)
    
    # add in the final frequency and time
    tvals = np.append(tvals,0)
    fvals = np.append(fvals,chirp.getFcut(mtot,eta))
    
    # get the corresponding amplitude 
    hvals = chirp.binaryAmp(fvals,mtot,eta,dl)
    Sh = makeSensitivity(fvals, model)
    snri = 4*np.real(hvals*np.conjugate(hvals)/Sh)
    snr = np.sqrt(np.cumsum(np.diff(fvals)*snri[1:]))
    tvals = tvals[1:]
    
    return snr, tvals, fvals, hvals
    
#Make snr from a source
def getSourceSnr(source,model,T = 4*constants.year, Npts = 1000,style='TN'):
    '''
    Compute the SNR given a GW source description and a GW model, both as dictionaries. 
    
    The calculaiton will use the SNR computation method that is relevant for that source
    '''
    
    stype = source.get('type')
    # continuous-wave source
    if stype == 'CW':
        
        # just do the frequency and amplitude
        if 'h0' in source:
            h0 = source.get('h0')
            f0 = source.get('f0')
        
        else :

            # get the chirp mass, first we try chirp mass directlycomponent masses
            if 'mchirp' in source:
                mchirp = source.get('mchirp')
                mtot = source.get('mtot')
            # failing that, we try the component masses
            else:
                m1 = source.get('m1')
                m2 = source.get('m2')
                mchirp = ((m1*m2)**(2./5.))/((m1+m2)**(1./5.))
                mtot = m1+m2

            # convert to seconds
            mtot = mtot * constants.MSun2s
            mchirp = mchirp * constants.MSun2s

            # get the frequency, compute the semi-major axis
            if 'f0' in source:
                f0 = source.get('f0')
                a = (mtot / (np.pi*f0)**2)**(1./3.)
            # get the semi-major axis, compute the frequency
            else:
                a = source.get('a')*constants.AU
                f0 = (1./np.pi)*(mtot/(a**3.))**(1./2.)

            # get the luminosity distance
            dl = source.get('dl')*constants.kpc2s

            # compute the ampltiude
            h0 = (2./dl)*(mchirp**(5./3.))*((np.pi*f0)**(2./3.))


        # compute the SNR
        rho = getCWsnr(f0,h0,T,model,style)
        
        # copy source to return with computed parameters
        sourceOut = source.copy()
        sourceOut['snr'] = rho
        sourceOut['T'] = T
        if not('f0' in sourceOut):
            sourceOut['f0']=f0
        
        if not('h0' in sourceOut):
            sourceOut['h0']=h0
        
        return rho, sourceOut
        
    # chirping source
    elif stype == 'chirp':
            # get the total mass
            if 'mtot' in source:
                mtot = source.get('mtot')*constants.MSun2s
            else:
                mtot = (source.get('m1') + source.get('m2'))*constants.MSun2s
                
            if 'eta' in source:
                eta = source.get('eta')
            else:
                eta = (source.get('m1')*source.get('m2'))/((source.get('m1')+source.get('m2'))**2)

            ds = source.get('dl')*constants.kpc2s
            
    
            print('mtot = %3.2g, eta = %3.2g, ds = %3.2g, T = %3.2g' % (mtot,eta,ds,T))
            snrt, tvals, fvals, hvals = getChirpSNR(mtot,eta,ds,model,T,Npts,style)
        
            i10 = np.argmin(np.abs(snrt-10))
            t10 = tvals[i10]

            observation = {
                'source' : source.copy(),
                'model' : model.copy(),
                't' : tvals,
                'f' : fvals,
                'h' : hvals,
                'SNR of t' : snrt,
                'SNR' : snrt[-1],
                'observation time' : t10
            }
            
            return observation
            
        
    # unsupported source, maybe need to throw an error/warning
    else: 
        print('Unsupported source type')
        return -1.0
    
    

### Imaging

def dResRange(fr,model):
    '''
    Here we construct two elementary figures of merit relevant for imaging, relevant for our imaging incoherent and our astrometric notions of imaging. 

    The first is basically diffraction limited resolution for short or long duration sources concentrated at some rerference frequency. The estimate is:
    $$
    \Delta \theta_\mathrm{diff} \approx F \frac{\lambda}{D}
    $$
    Where we suppose $F\approx 1$ and
    $$
    \max(L_\mathrm{constellation},D_\mathrm{sep})\leq D \leq \max(2R_{orbit},L_\mathrm{constellation},D_\mathrm{sep})
    $$
    depending on how long the source lasts compared to $T_\mathrm{orbit}$.

    The other is astrometric localization which is scaled by the SNR:
    $$
    \Delta \theta_\mathrm{am} \approx \Delta \theta_\mathrm{diff}/\rho
    $$

    There are a number of different ways we can think about making plots using these, including horizon distances for reference classes of obervations, etc. More thought is needed on what makes sense..
    '''
    
    D=model['Lconst']
    c=constants.c
    au=constants.AU*c
    if 'Dsep' in model: D=max([D,model['Dsep']*au])
    Dshort=D
    if 'Rorbit' in model: D=max([D,2*model['Rorbit']*au])
    Dlong=D
    dtheta_long=c/fr/Dlong
    dtheta_short=c/fr/Dshort
    return dtheta_long,dtheta_short
