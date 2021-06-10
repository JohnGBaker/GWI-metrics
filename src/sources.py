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

NSBH = {
    'type' : 'chirp',
    'label' : r'NS-BH @ 1Gpc',
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