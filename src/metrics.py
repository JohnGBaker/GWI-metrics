#These are functions used for GW mission sensitivity calculations and other figures of merit related to GW iimaging or sky localization

import numpy as np
import constants
import PhenomWaveform_nonspinning as Phenom
import subsystems
import background
import sources

#Probably adapt more for GW Imager concepts 
def PSD_noise_components(fr, model):
    '''
    Make PSD Noise Components
    We follow LISA in computing the noise. In that case the noise is more-or-less directly derived from 
    - a measure of the acceleration noise in m/s^2/sqrt(Hz) with some additional reddening factors
    - a measure of optical measurement system noise in m/sqrt(Hz)
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

#Components for the Larson computation from LISA TN
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

    
#Compute sensitivity curve
def makeSensitivity(fr, model,style='TN'):
    #'''
    #Using the semi-analytical average response, the semi-analytical sensitivity for TDI X 4 links is:
    #$$
    #S_{h,X} =  \frac{ S_{OMS} + \left( 3 + \cos \left( \frac{2 \omega L}{c} \right)  \right)  S_{acc} }
    #{ \left( {\omega L \over c} \right)^2 \ R_{\Sigma}^2(f, L) }
    #$$
    #$$
    #S_{h,X} 
    #=  \frac{S_{n,X_{2.0}}}
    #{<R_{L, X_{2.0}}(f)>} 
    #= \frac{ 64 \sin^2 \left( \omega L \right) \sin^2 \left(2 \omega L \right) (S_{OMS} + \left( 3 + \cos \left( \frac{2 \omega L}{c} \right)  \right)  S_{acc}) }
    #{ (4\omega L)^2 \sin^2{(\omega L)} (2 \sin{(2\omega L)})^2 <(F^{+}_{X})^2> }
    #= \frac{ S_{OMS} + \left( 3 + \cos \left( \frac{2 \omega L}{c} \right)  \right)  S_{acc} }
    #{ \left( \omega L \right)^2  <(F^{+}_{X})^2>  }
    #'''
    #print('makeSens:',fr[0],'< f <',fr[-1],'style=',style, 'model:')
    #display(model)
    
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

    # initialize baseline array
    B = np.zeros_like(t)
    
    # minimum baseline is the detector size
    Bc = model.get('Lconst')/constants.c

    haveOrbit=False
    if 'Torbit' in model and 'Rorbit' in model and not np.isnan(tstart):
        haveOrbit=True
        # start using your orbit from tstart
        istart = np.argmin(np.abs(t-tstart))
        theta = 2*np.pi*np.clip((t[istart:]-t[istart])/(model.get('Torbit')*constants.year),0,0.5)
        Borbit = model.get('Rorbit')*constants.AU*np.sqrt(2*(1-np.cos(theta)))
        B[istart:]= B[istart:] + Borbit
        if np.isnan(sum(B)):
            print('Got ',np.count_nonzero(np.isnan(B)),' NaNs (out of',len(B),' for orbit baseline. istart',istart,' tstart',tstart)

    # if you have multiple constellations, use the Dsep parameter (which is in what units?)
    if 'Dsep' in model:
        Dsep = model.get('Dsep')
        if Dsep*constants.AU > Bc: # essentially a zero check
            Bmin = Dsep*constants.AU
        else:
            Bmin = Bc
    else:
        Bmin = Bc

    if haveOrbit:
        # use the greater of the available baselines for each timestep
        igreater = np.where(B < Bmin)
        B[igreater] = Bmin
    else:
        B=Bmin

    return B

#Make SNR for continuous-wave source
def getCWsnr(f0,h0,T,model,style='TN'):
    #'''
    #Compute the SNR for a monochromatic GW source based on the source frequency f0, source ampltiude h0, obsertvation time T, and an instrument model. 
    #
    #The calculation follows (83) from the LISA-LCST-SGS-TN-001 (https://arxiv.org/abs/2108.01167) to compute the inclinaiton, polarization, and sky-position averaged SNR for a monochomatic source:
    #$$
    #\left<SNR^2\right>_{\iota,\psi,sky} = 10 \frac{\left(\frac{2}{5}h_0\right)^2T}{S_h\left(f_0\right)}
    #$$
    #
    #where the instrument sensitivity is computed from the provided model and the makeSensitivity method
    #'''
    
    
    # compute sensitivity at the GW frequency from the model
    S_h_f0 = makeSensitivity(f0,model,style)
    #print('T,f0,S_h_f0',T,f0,S_h_f0)
    
    # apply equation (83)
    rho2 = 10*(((2./5.)*h0)**2)*T/S_h_f0 
    
    # return SNR
    return np.sqrt(rho2)

#Make SNR for continuous-wave source
def getChirpSNR(mtot,eta,dl,model,tstart=-constants.year,Npts = 1000,style='TN',tstop=None):
    '''
    Compute the SNR for a chirping source
    '''
    # set up time vector
    # stop time is when we get to merger frequency / 3 to avoid PN blow-up
    if tstop is None:
        tstop = Phenom.tFromF(0.3*Phenom.getFmerge(mtot,eta),mtot,eta)
    tvals = -np.flip(np.logspace(np.log10(-tstop),np.log10(-tstart),Npts))
    
    # get the corresponding frequency vector
    fvals = Phenom.fFromT(tvals,mtot,eta)
    
    # get the corresponding amplitude 
    hvals = Phenom.binaryAmp(fvals,mtot,eta,dl)
    Sh = makeSensitivity(fvals, model)
    snri = 4*np.real(hvals*np.conjugate(hvals)/Sh)
    snrt = np.sqrt(np.cumsum(np.diff(fvals)*snri[1:]))
    tvals = tvals[1:]
    
    # for the total SNR, we integrate over all the frequencies to better sample the waveform
    fsnr = np.logspace(np.log10(fvals[0]),np.log10(Phenom.getFcut(mtot,eta)),Npts)
    hsnr = Phenom.binaryAmp(fsnr,mtot,eta,dl)
    Shsnr = makeSensitivity(fsnr,model)
    snri = 4*np.real(hsnr*np.conjugate(hsnr)/Shsnr)
    snr = np.sqrt(np.sum(np.diff(fsnr)*snri[1:]))
    
    # add in the final frequency and time
    tvals = np.append(tvals,0)
    snrt = np.append(snrt,snr)

    
    return snrt, tvals, fsnr, hsnr
    
#Make snr from a source
def getSourceSnr(source,model,T = 4*constants.year, Npts = 1000,style='TN'):  
    '''
    Compute the SNR given a GW source description and a GW model, both as dictionaries. 
    
    The calculaiton will use the SNR computation method that is relevant for that source
    '''
    
    stype = source.get('type')
    # continuous-wave source
    if stype == 'CW':

        #T needs to be positive
        T=abs(T)
        
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
                f0 = np.array(source.get('f0'))
                a = (mtot / (np.pi*f0)**2)**(1./3.)
            # get the semi-major axis, compute the frequency
            else:
                a = source.get('a')*constants.AU
                f0 = np.array((1./np.pi)*(mtot/(a**3.))**(1./2.))

            # get the luminosity distance
            dl = source.get('dl')*constants.kpc2s

            # compute the ampltiude
            h0 = (2./dl)*(mchirp**(5./3.))*((np.pi*f0)**(2./3.))

        if np.size(T) == 1:
            T = np.linspace(1,T,Npts)        #JGB REVIEW: What is going on here? Seem to get an array of mostly NANs

        # compute the SNR
        snrt = getCWsnr(f0,h0,T,model,style)
        #print('snrt: len,min,max=',len(snrt),min(snrt),max(snrt))
        
        i10 = np.argmin(np.abs(snrt-10))
        t10 = T[i10]                         #JGB REVIEW: Does this work for low SNR cases?
        #print('snrt[-1]=',snrt[-1])
        #print(snrt)
        
        observation = {
                'source' : source.copy(),
                'model' : model.copy(),
                't' : T,
                'f' : f0,
                'h' : h0,
                'SNR of t' : snrt,
                'SNR' : snrt[-1],
                'observation time' : t10
            }
        
    # chirping source
    elif stype == 'chirp':
        
            #T needs to be negative
            T=-abs(T)
        
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
            
    
            #print('mtot = %3.2g, eta = %3.2g, ds = %3.2g, T = %3.2g' % (mtot,eta,ds,T))
            tstop=source.get('timecut',None)
            if tstop is not None:tstop*=counstants.year
            snrt, tvals, fvals, hvals = getChirpSNR(mtot,eta,ds,model,T,Npts,style,tstop=tstop)
        
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
            

            
        
    # unsupported source, maybe need to throw an error/warning
    else: 
        print('Unsupported source type')
        observation = {}
    
    
    return observation

    
    
    

### Imaging

def getResolution(obsIn):
    '''
    Compute the angular resolution as a function of time for an observation.
    
    '''
    
    obsOut = obsIn.copy()
    t = obsOut.get('t')
    snr = obsOut.get('SNR of t')
    f = obsOut.get('f')
    isnr2 = np.clip(np.argmin(np.abs(snr-0.5*snr[-1])),0,len(t)-1)
    tsnr2 = t[isnr2]
    snr2 = snr[isnr2]
    if(np.isnan(tsnr2)):
        print('tsnr2 is nan. snr=',snr[-1],'isnr2=',isnr2)
    
    if np.size(f)==1:
        fsnr2 = f
    else:
        fsnr2 = f[isnr2]
    
    obsOut['t half SNR'] = tsnr2
    obsOut['f half SNR'] = fsnr2
    
    # estimate the diffraciton limit
    lamGW = constants.c / fsnr2
    B = getBaseline(obsOut.get('model'),t,tsnr2)
    deltaThetaDiff = (lamGW/constants.c)/B
    obsOut['Baseline'] = B
    obsOut['Diffraction Limit'] = deltaThetaDiff
    
    # estimate the angular resolution
    deltaTheta = deltaThetaDiff/snr
    obsOut['Angular Resolution of t'] = deltaTheta
    obsOut['Angular Resolution'] = deltaTheta[-1]
    return obsOut
    

def dResRange(fr,model):
    #'''
    #Here we construct two elementary figures of merit relevant for imaging, relevant for our imaging incoherent and our astrometric notions of imaging. 
    #
    #The first is basically diffraction limited resolution for short or long duration sources concentrated at some rerference frequency. The estimate is:
    #$$
    #\Delta \theta_\mathrm{diff} \approx F \frac{\lambda}{D}
    #$$
    #Where we suppose $F\approx 1$ and
    #$$
    #\max(L_\mathrm{constellation},D_\mathrm{sep})\leq D \leq \max(2R_{orbit},L_\mathrm{constellation},D_\mathrm{sep})
    #$$
    #depending on how long the source lasts compared to $T_\mathrm{orbit}$.

    #The other is astrometric localization which is scaled by the SNR:
    #$$
    #\Delta \theta_\mathrm{am} \approx \Delta \theta_\mathrm{diff}/\rho
    #$$
    #
    #There are a number of different ways we can think about making plots using these, including horizon distances for reference classes of obervations, etc. More thought is needed on what makes sense..
    #'''
    
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

######################################
## Revised calculations Summer 2022 ##
######################################

# This is new variation on what was in getChirpSNR above. The part focuses on
# providing what is needed for the new localization calculations (also for SNR).
def computeChirp(source,model,tstart=-1,Npts = 1000,fstop=None,tstop=None, resp_style='TN',fstart=None):
    '''
    source -A dict containing:
      eta    Unitless reduced mass
      dl     Luminosity distance in kpc
    model  -Observatory model
    tstart -Initial time before coalescence in years (tstart>tstop or tstart<0)
    Npts   -Number of time samples in output
    fstop  -highest freq (also reference for t=0) [default=fRingdown, but decreasing for EMRI's]
    tstop  -a time before 0 at fmerge (a version of merger) to cut the signal
    fstart If present, this overrides tstart

    By comparison with getChirpSNR, this version works a little differently.
    Some changes are to support using the full signal through merger. 
    Previously, the fFromT function was a PN expression which didn't support 
    merger. Instead we detive t_of_f (a new function in Phenom) directly from
    the model phase. (See FoundationCalcs notebook for some tests of that.
    This is only defined in the f->t direction, so this version differs in 
    that the grid is explicitly defined in f, with t derived. By default the 
    frequency terminates at ringdown (for comparable masses).  

    For more extreme mass ratios, the models nere fail at high frequencies 
    (which also are not usually relevant with insufficient signal power.)
    It would be natural to terminate EMRI inspirals at ISCO.  Here we just pick
    0.3Fmerge (which was used in the previous code.  There is a blending
    to transition between these terminal frequencies based on mass ratio.

    A consequence of the grid being implicit in time is that we cannot direct
    impose a stop time.  Instead we have a search for the latest f,t sample
    which does not exceed th desired cutoff.  We have implemented a search 
    scheme, but the result is approximate, depending on the grid density.
    '''

    mtot,eta=sources.get_mtot_eta(source)
    dl = source.get('dl')*constants.kpc2s
    tstart = tstart*constants.year
    fmerge=Phenom.getFmerge(mtot,eta)
    if fstop is None:
        fring=1.0*Phenom.getFring(mtot,eta)
        fcut_emri=Phenom.getFisco(mtot,eta)
        blend=eta/(eta+.05)
        fstop=fring-(1-blend)*(fring-fcut_emri)
        #print('fstop=',fstop)
    else:
        fstop=fstop

    fstartin=fstart
    if fstart is None:
        fstart=Phenom.fFromT(tstart,mtot,eta)*0.95 #coarse estimate
        if np.isnan(fstart): fstart=1e-6
        print('rough tstart =',tstart,fstart,'< f <',fstop)
        fvals=np.logspace(np.log10(fstart),np.log10(fstop),Npts)
        #print('fstop,fmerge,fring',fstop,fmerge,fring)
        tvals=Phenom.t_of_f(fvals,mtot,eta,zero_at_f=fmerge)
        #Now trim to fine-tune start
        ncut=find_ncut(tvals,tstart,-1)
        tvals=tvals[ncut:]
        fvals=fvals[ncut:]
        fstart=fvals[0]
        print('ncut=',ncut, 'n=',len(tvals))
    fvals=np.logspace(np.log10(fstart),np.log10(fstop),Npts)
    tvals=Phenom.t_of_f(fvals,mtot,eta,zero_at_f=fmerge)
    print('tstart =',tstart,'t[0]=',tvals[0],fstart,'< f <',fstop)
    
    
    #implement a stop time
    if tstop is not None:
        tstop=tstop*constants.year
        #we search for the closest cut point in the time series
        ncut=find_ncut(tvals,tstop)
        tvals=tvals[:-ncut]
        fvals=fvals[:-ncut]
    
    # get the corresponding amplitude 
    hvals = Phenom.binaryAmpOnly(fvals,mtot,eta,dl)
    Sh = makeSensitivity(fvals, model, style=resp_style)
    #snri = 4*np.real(hvals*np.conjugate(hvals)/Sh)
    snri = 4*(hvals**2/Sh)
    #snrt = np.sqrt(np.cumsum(np.diff(fvals)*snri[1:]))
    tvals = tvals[1:]
    #snri=np.concatenate(([0],snri))
    return [tvals,fvals,snri,Sh]

def find_ncut(times,tstop,sign=1):
    '''
    Find the number of elements to cut from the end so that the remaining values have time<=tstop.  Assumes monotonicity in the relevant range.
    '''
    verbose=False
    #we search for the closest cut point in the time series
    if verbose: print("looking for t=",tstop,"in [",times[0],",",times[-1],"]")
    ncut=0
    stepfac=4
    step=stepfac**4
    count=0
    while step>0:
        if verbose:
            print('step=',step)
            if (ncut+step)<(len(times)-1): print('ttest:',times[-sign*(ncut+step)],'tlimit',tstop,'-->',times[-sign*(ncut+step)]*sign>=tstop*sign)
        while (ncut+step)<(len(times)-1) and times[-sign*(ncut+step)]*sign>=tstop*sign: 
            ncut+=step
            if verbose: print('ncut=',ncut)
            count+=1
            if verbose:
                if (ncut+step)<(len(times)-1): print('ttest:',times[-sign*(ncut+step)],'tlimit',tstop,'-->',times[-sign*(ncut+step)]*sign>=tstop*sign)

        step=step//stepfac
    if verbose: print('found ncut=',ncut,'in',count,'  ',times[-(ncut+1)],'<=',tstop,times[-ncut])
    return ncut

def get_rho2om2(chirp):
    tvals,fvals,snri,Sh=chirp
    f=fvals[1:]
    return np.sum(np.diff(fvals)*snri[1:]*f**2)*4*np.pi**2

#Compute the sky-angle resolution for a chirping signal
def get_chirp_snr_sky_sigma2_orbit(chirp,model):
    tvals,fvals,snri,Sh=chirp
    f=fvals[1:]
    df=np.diff(fvals)
    t=tvals
    wt=snri[1:]*f**2
    if len(wt)<1:
        print('get_chirp_snr_sky_sigma2_orbit: Empty chirp, returning NaNs')
        return {'sig2':float('nan'),'rho2':float('nan'),'rho2om2':float('nan')}
    i0=np.argmax(wt)
    fac=sum(df*wt)
    wt*=df/fac
    fac=fac*4*np.pi**2
    t0=t[i0]
    #print('get_...orbit:t0=',t0,'fac=',fac)
    haveOrbit=False
    if 'Rorbit' in model and 'Torbit' in model:
        R=model['Rorbit']*constants.AU
        Om=2*np.pi/(model['Torbit']*constants.year)
        tau=t+R*np.sin(Om*t)
        gam=+R*np.cos(Om*t)
        tauav=sum(wt*tau)
        gamav=sum(wt*gam)
        tau-=tauav
        gam-=gamav
        #print('get_...orbit:tauav,gamav',tauav,gamav)
        I11=sum(wt*tau**2)
        I12=sum(wt*tau*gam)
        I22=sum(wt*gam**2)
        sig2=1/fac/(I22-I12**2/I11)
    else:
        sig2=float('inf')
    #print('get_...orbit:I11,I12,I22,sig2',I11,I12,I22,sig2)
    rho2=np.sum(df*snri[1:])
    #rho2om2=np.sum(df*snri[1:]*f**2)*4*np.pi**2
    res={'sig2':sig2,'rho2':rho2,'rho2om2':fac}
    return res

def get_CW_sky_res_orbit(f0,h0,rho, model):
    '''
    Compute CW sky angle variance via orbial motion using simplified 
    Fisher-based formulation.
    
    See FomulationCalcs notebook for derivation.  Note that the SNR rho may be
    provided as a value or an nd.array.
    '''
    if 'Rorbit' in model and 'Torbit' in model:
        R=model['Rorbit']*constants.AU
        Om=2*np.pi/(model['Torbit']*constants.year)
        rhoom=2*np.pi*f0*rho
        sig2=2/(R*rhoom)**2
    else:
        sig2=float('inf')
    return sig2

def get_sky_sigma2_Lconst(rho2om2,model):
    '''
    Compute sky angle variance via constellation size using simplified 
    Fisher-based formulation.
    
    See FomulationCalcs notebook for derivation.  The same formula applies
    for CW and chirping sources. Note that the rho2om2 may be provided as 
    a value or an nd.array.
    '''
    L=model['Lconst']/constants.c
    return 4/rho2om2/L**2     

def get_sky_sigma2_Dsep(rho2om2,model):
    '''
    Compute sky angle variance via constellation size using simplified 
    Fisher-based formulation.
    
    See FomulationCalcs notebook for derivation.  The same formula applies
    for CW and chirping sources. Note that the rho2om2 may be provided as 
    a value or an nd.array and that 'rho' here is for the total, not for 
    each constellation.  This formula assumes that each constellation has
    the same sensitivity.
    '''
    Dsep=model.get('Dsep',0)*constants.AU
    if Dsep==0:Dsep=1e-20
    return 4/rho2om2/Dsep**2   

def getSNRandSkyResolution(source,model, Nsamp=0, Nres = 1000, Tmax = None, resp_style='TN',SNRdetect=10,SNRcut=None):  
    '''
    Compute the SNR and angular sky resolution for an observation, applying 
    the method appropriate for the source type.

    Args:
      source       Source descriptor dict
      model        Concept descriptor dict
      Nsamp        Number of time/freq samples to include in output grid.
      Nres         Number of samples to use in resolving chirp signals
      Tmax         Maximum observation duration (yr). For chirps, this 
                   specifies the approximate start time before merger. 
      resp_style   Style option for response calculation.
      SNRdetect    Starting SNR from which to estimate detection
      SNRcut       Force sig->2pi as SNR->SNRcut, small SNR cases
    
    If Nsamp>1, then result includes nd.array with a grid of computations for 
    multiple durations of observation. The grid is chosen appropriately for 
    the source type.
    '''
    
    stype = source.get('type')
    
    observation = {
        'source' : source.copy(),
        'model' : model.copy()
    }
    
    if stype == 'CW':
        # continuous-wave source
        
        h0, f0 = sources.get_CW_h0_f0(source)

        if Tmax is None: Tmax=100
        Tdur=min([Tmax,model.get('SciDuration',Tmax)])

        if Nsamp>0:
            #Create a time grid uniformly spaced in sqrt(t), not including 0
            t = np.linspace(0,np.sqrt(Tdur),Npts+1)[1:]**2
        else:
            t=Tdur
        f=f0
        snr = getCWsnr(f0,h0,t,model,resp_style)
        rho2om2= snr**2*(2*np.pi*f0)**2 # added square on second factor [2025-04-25]
        #rho2om2= snr**2*(2*np.pi*f0)
        sig2orb=get_CW_sky_res_orbit(f0,h0,snr, model)

        observation['f']=f0
        observation['h']=h0

        
    elif stype == 'chirp':
        # chirping source

        #First specify the relevant time period
        Tdur=Tmax
        Tdur=min((Tdur,model.get('SciDuration',Tdur)))
        tstart=-Tdur
        #print('Tdur',Tdur,'tstart',tstart)
        daysperyear=constants.year/(24*3600)
        tstop_source=source.get('timecut',None)
        if tstop_source is not None: tstop_source*=-1/daysperyear
        #print('tstop:',tstop_source)
        #Make first pass on chirp, based on source
        tvals,fvals,snri,Sh=computeChirp(source,model,tstart=tstart,Npts=Nres,fstop=None,tstop=tstop_source,resp_style=resp_style)

        #Now refine for model observation period
        if len(tvals) < Nres/2:
            print('Recomputing chirp because net resolution is low.',len(tvals),'<',Nres/2)
            #print('initial chirp:',tvals,fvals,snri,Sh)
            #Under-resolved so we will redo the chirp
            #first find the cutoff freq
            fstop=None
            if len(fvals)>0: fstop=fvals[-1]*1.2 #add a buffer
            #print('fstop',fstop)
            #compute the net stop time
            #redo the chirp only up to cutoff freq
            tvals,fvals,snri,Sh=computeChirp(source,model,tstart=tstart,Npts=Nres,fstop=fstop,tstop=tstop_source,resp_style=resp_style)
            #print('recomputed chirp:',tvals,fvals,snri,Sh)
        
        chirp=tvals,fvals,snri,Sh
        #print(chirp)
        
        if Nsamp>0:
            #set up arrays for output
            snr=np.zeros(Nsamp)
            t=np.zeros(Nsamp)
            f=np.zeros(Nsamp)
            rho2om2=np.zeros(Nsamp)
            sig2orb=np.zeros(Nsamp)
            #set the final values to correspond to the full signal
            res=get_chirp_snr_sky_sigma2_orbit(chirp,model)
            t[-1]=tvals[-1]
            f[-1]=fvals[-1]
            snr[-1]=np.sqrt(res['rho2'])
            rho2om2[-1]=res['rho2om2']
            sig2orb[-1]=res['sig2']
            
            #Create a log-uniform output time grid
            tcuts = -np.flip(np.logspace(np.log10(max([100,-tvals[-1]])),np.log10(-tvals[0]),Nsamp))
            tcuts=tcuts[1:]
            #Fill in grid
            for i in range(Nsamp-1):
                #first trim the chirp
                tcut=tcuts[i]
                ncut=find_ncut(tvals,tcut)
                if ncut>=len(tvals):
                    t[i]=float('nan')
                    f[i]=float('nan')              
                    snr[i]=float('nan')
                    rho2om2[i]=float('nan')
                    sig2orb[i]=float('nan')
                    continue
                if ncut>0:
                    chirp=tvals[:-ncut],fvals[:-ncut],snri[:-ncut],Sh[:-ncut]
                else:
                    chirp=tvals,fvals,snri,Sh
                #print(len(tvals),ncut,len(chirp[0]))
                res=get_chirp_snr_sky_sigma2_orbit(chirp,model)
                t[i]=chirp[0][-1]
                f[i]=chirp[1][-1]
                snr[i]=np.sqrt(res['rho2'])
                rho2om2[i]=res['rho2om2']
                sig2orb[i]=res['sig2']
        else:
            #single sample chirp
            #print('chirp:',chirp[0][0],'< t <',chirp[0][-1],' range',(chirp[0][-1]-chirp[0][0])/constants.year,'yr  ',chirp[1][0],'< f <',chirp[1][-1])
            res=get_chirp_snr_sky_sigma2_orbit(chirp,model)
            #print('res')
            #display(res)
            if len(tvals)==0:
                t=float('nan')
                f=float('nan')              
                snr=float('nan')
                rho2om2=float('nan')
                sig2orb=float('nan')
            else:
                t=tvals[-1]
                f=fvals[-1]
                snr=np.sqrt(res['rho2'])
                rho2om2=res['rho2om2']
                sig2orb=res['sig2']
            
        #completion of chirp info
        observation['f']=f
    else: 
        print('Unsupported source type')

    #Now, for both CW and chirping,
    #compute the constellation and separation contirbutions to resolution
    sig2con=get_sky_sigma2_Lconst(rho2om2,model)
    sig2sep=get_sky_sigma2_Dsep(rho2om2,model)
    #Assemble the total adding in quadrature
    sig2=(sig2orb**-1+np.minimum(sig2con,sig2sep)**-1)**-1
    if SNRcut is not None:
        sig2+=4*np.pi**2/(1+(snr/SNRcut)**4)
    sig=np.sqrt(sig2)

    if Nsamp>0:
        #print('t',t)
        #print('SNR',snr)
        idet = np.argmin(np.abs(snr-SNRdetect))
        if snr[idet]>=SNRdetect:
            tdet = t[idet]
            observation['detection time']=tdet                            
        observation['SNR of t']=snr
        observation['Angular Resolution of t'] = sig
        snr=snr[-1]
        sig=sig[-1]
        observation['Angle Variance (orbit) of t'] = sig2orb
        observation['Angle Variance (const) of t'] = sig2con
        observation['Angle Variance (sep) of t'] = sig2sep
        observation['Effective omega of t'] = np.sqrt(rho2om2)/snr
        observation['Effective D of t'] = 4/rho2om2/sig**2
        sig2orb=sig2orb[-1]
        sig2con=sig2con[-1]
        sig2sep=sig2sep[-1]
        rho2om2=rho2om2[-1]
    observation['t']=t
    observation['SNR']=snr
    observation['Angular Resolution'] = sig
    observation['Angle Variance (orbit)'] = sig2orb
    observation['Angle Variance (const)'] = sig2con
    observation['Angle Variance (sep)'] = sig2sep
    observation['Effective omega'] = np.sqrt(rho2om2)/snr
    observation['Effective D'] = 4/rho2om2/sig**2
    
    return observation

    
    
    

