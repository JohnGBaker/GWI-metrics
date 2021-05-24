# This file defines fiducial sources for GW imaging metric calculaitons
import constants


# Parameters for AmCV taken from CLHT demo
AmCV = {
    'type' : 'CW',
    'label' : 'AmCV',
    'm1' : 0.68,               # mass of object 1 in solar masses
    'm2' : 0.125,              # mass of object 2 in solar masses
    'dkpc' : 0.3,              # luminosity distance (kpc)
    'a'  : 0.5*9.5e-4,         # binary separation
}


MBHB61 = {
    'type' : 'chirp',
    'label': r'$10^{6}\,M_{\odot}+10^{6}\,M_{\odot}\:@\:z=1$',
    'm1' : 1e6,
    'm2' : 1e6,
    'dl' : constants.z1kpc
}

MBHB73 = {
    'type' : 'chirp',
    'label': r'$10^{7}\,M_{\odot}+10^{7}\,M_{\odot}\:@\:z=3$',
    'm1' : 1e7,
    'm2' : 1e7,
    'dl' : constants.z3kpc
}