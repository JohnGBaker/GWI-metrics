# This file defines fiducial sources for GW imaging metric calculaitons
import constants


# Parameters for AmCV taken from CLHT demo
AmCV = {
    'type' : 'CW',
    'label' : 'AmCV',
    'm1' : 0.68,               # mass of object 1 in solar masses
    'm2' : 0.125,              # mass of object 2 in solar masses
    'dl' : 0.3,              # luminosity distance (kpc)
    'a'  : 0.5*9.5e-4,         # binary separation
}

NSNS = {
    'type' : 'chirp',
    'label' : 'NS-NS',
    'm1' : 1.5,
    'm2' : 1.5,
    'dl' : 500*1e3   # distance is 500Mpc
    
}

NSNS1_1a = {
    'type' : 'chirp',
    'label' : 'NS-NS POC 1.1.a',
    'm1' : 1.5,
    'm2' : 1.5,
    'dl' : 230*1e3   # distance is 230Mpc
    
}

NSBH = {
    'type' : 'chirp',
    'label' : r'NS-BH @ 1Gpc',
    'm1' : 10,
    'm2' : 1.5,
    'dl' : 1e6   # distance is 1Gpc
    
}

NSBH1_2a = {
    'type' : 'chirp',
    'label' : r'NS-BH @ 1Gpc (POC 1.2.a)',
    'm1' : 10,
    'm2' : 1.5,
    'dl' : 1e6   # distance is 1Gpc
    
}

SOBH = {
    'type' : 'chirp',
    'label' : r'$30\,M_{\odot}+30\,M_{\odot}\:@\:1\:Gpc$',
    'm1' : 30,
    'm2' : 30,
    'dl' : 1e6   # distance is 1Gpc
    
}

SOBH1_3a = {
    'type' : 'chirp',
    'label' : r'$30\,M_{\odot}+30\,M_{\odot}\:@\:1\:Gpc$ (POC 1.3.a)',
    'm1' : 30,
    'm2' : 30,
    'dl' : 1e6   # distance is 1Gpc
    
}

MBH2_1a = {
    'type' : 'chirp',
    'label': r'$10^{6}\,M_{\odot}+10^{6}\,M_{\odot}\:@\:z=1$ (POC 2.1.a)',
    'm1' : 1e6,
    'm2' : 1e6,
    'dl' : constants.z1kpc # z = 1
}

MBH2_2a = {
    'type' : 'chirp',
    'label': r'$5 \times 10^{5}\,M_{\odot}+5 \times 10^{5}\,M_{\odot}\:@\:z=2$ (POC 2.2.a)',
    'm1' : 5e5,
    'm2' : 5e5,
    'dl' : constants.z2kpc # z = 2
}

MBHB61 = {
    'type' : 'chirp',
    'label': r'$10^{6}\,M_{\odot}+10^{6}\,M_{\odot}\:@\:z=1$',
    'm1' : 1e6,
    'm2' : 1e6,
    'dl' : constants.z1kpc
}

MBHB33 = {
    'type' : 'chirp',
    'label': r'$10^{3}\,M_{\odot}+10^{3}\,M_{\odot}\:@\:z=3$',
    'm1' : 1e3,
    'm2' : 1e3,
    'dl' : constants.z3kpc
}

MBHB31 = {
    'type' : 'chirp',
    'label': r'$10^{3}\,M_{\odot}+10^{3}\,M_{\odot}\:@\:z=1$',
    'm1' : 1e3,
    'm2' : 1e3,
    'dl' : constants.z1kpc
}

MBHB43 = {
    'type' : 'chirp',
    'label': r'$10^{4}\,M_{\odot}+10^{4}\,M_{\odot}\:@\:z=3$',
    'm1' : 1e4,
    'm2' : 1e4,
    'dl' : constants.z3kpc
}

MBHB41 = {
    'type' : 'chirp',
    'label': r'$10^{4}\,M_{\odot}+10^{4}\,M_{\odot}\:@\:z=1$',
    'm1' : 1e4,
    'm2' : 1e4,
    'dl' : constants.z1kpc
}

IMBH313 = {
    'type' : 'chirp',
    'label': r'$10^{3}\,M_{\odot}+10\,M_{\odot}\:@\:z=3$',
    'm1' : 1e3,
    'm2' : 10,
    'dl' : constants.z3kpc
}

IMBH413 = {
    'type' : 'chirp',
    'label': r'$10^{4}\,M_{\odot}+10\,M_{\odot}\:@\:z=3$',
    'm1' : 1e4,
    'm2' : 10,
    'dl' : constants.z3kpc
}

EMRI3_1a = {
    'type' : 'chirp',
    'label': r'$10^{8}\,M_{\odot}+100\,M_{\odot}\:@\:z=0.5$ (POC 3.1.a)',
    'm1' : 1e8,
    'm2' : 100,
    'dl' : 2.92e6 # distance is 2.92 Gpc
}

EMRI3_2a = {
    'type' : 'chirp',
    'label': r'$10^{7}\,M_{\odot}+1\,M_{\odot}\:@\:z=0.5$ (POC 3.2.a)',
    'm1' : 1e7,
    'm2' : 1,
    'dl' : 2.92e6 # distance is 2.92 Gpc
}

EMRI3_3a = {
    'type' : 'chirp',
    'label': r'$5 \times 10^{5}\,M_{\odot}+10\,M_{\odot}\:@\:z=1.5$ (POC 3.3.a)',
    'm1' : 5e5,
    'm2' : 10,
    'dl' : 1.22e7 # distance is 12.2 Gpc
}

MBHB73 = {
    'type' : 'chirp',
    'label': r'$10^{7}\,M_{\odot}+10^{7}\,M_{\odot}\:@\:z=3$',
    'm1' : 1e7,
    'm2' : 1e7,
    'dl' : constants.z3kpc
}

MBHB81 = {
    'type' : 'chirp',
    'label': r'$10^{8}\,M_{\odot}+10^{8}\,M_{\odot}\:@\:z=1$',
    'm1' : 1e8,
    'm2' : 1e8,
    'dl' : constants.z1kpc
}

GW150914 = {
    'type' : 'chirp',
    'label': r'GW150914',
    'm1' : 36,
    'm2' : 29,
    'dl' : 410*1e3
}

BHB31G = {
    'type' : 'chirp',
    'label': r'BHB31G',
    'm1' : 1e3,
    'm2' : 1e3,
    'dl' : 1e6
}

PSB4_1a = {
    'type' : 'persistent',
    'label': r'WDB $0.6\,M_{\odot}+0.6\,M_{\odot}\:@\:z=0.2$ (POC 4.1.a)',
    'm1' : 0.6,
    'm2' : 0.6,
    'dl' : 1e6 # distance is 1 Gpc
}

PSB4_1b = {
    'type' : 'persistent',
    'label': r'WDB $0.6\,M_{\odot}+0.6\,M_{\odot}\:@\:z=0.2$ (POC 4.1.a)',
    'm1' : 0.6,
    'm2' : 0.6,
    'dl' : 1e3 # distance is 1 Mpc
}

PSB4_2a = {
    'type' : 'persistent',
    'label': r'WDB $0.6\,M_{\odot}+0.6\,M_{\odot}\:@\:z=0.2$ (POC 4.1.a)',
    'm1' : 0.6,
    'm2' : 0.6,
    'dl' : 1e3 # distance is 1 Mpc
}


def get_mtot_eta(source):
    '''
    A utility for revtieving mtot and eta from available info in the source. 

    Generally intended for chirping sources.
    '''

    if 'mtot' in source:
        mtot = source.get('mtot')*constants.MSun2s
    else:
        mtot = (source.get('m1') + source.get('m2'))*constants.MSun2s                    
    if 'eta' in source:
        eta = source.get('eta')
    else:
        eta = (source.get('m1')*source.get('m2'))/((source.get('m1')+source.get('m2'))**2)
    return mtot,eta


def get_CW_h0_f0(source):
    '''
    A small utility for revtieving h0 and f0 from available CW source data.
    '''
    
    if 'h0' in source:
        # if available just take the frequency and amplitude
        h0 = source.get('h0')
        f0 = source.get('f0')
        
    else :
        #compute from source physical params

        # convert to seconds
        mtot = mtot * constants.MSun2s
        mchirp = mchirp * constants.MSun2s

        if 'f0' in source:
            # get the frequency, compute the semi-major axis
            f0 = np.array(source.get('f0'))
            a = (mtot / (np.pi*f0)**2)**(1./3.)
        else:
            # get the semi-major axis, compute the frequency
            a = source.get('a')*constants.AU
            f0 = np.array((1./np.pi)*(mtot/(a**3.))**(1./2.))

        # get the luminosity distance
        dl = source.get('dl')*constants.kpc2s
        
        # compute the ampltiude
        h0 = (2./dl)*(mchirp**(5./3.))*((np.pi*f0)**(2./3.))

    return h0,f0
