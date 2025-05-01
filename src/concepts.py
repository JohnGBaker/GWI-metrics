# This file provides descriptors for reference GW imager mission concepts.
# Input parameters initially based on the Gravitational Wave Imager/Mission Concepts Design/Mission Architecture Trades 2.xlsx spreadsheet 
# Compared and updates with the parameters of concepts published and available in the literature as of February 2025
# last updated Feb 2025

import subsystems

##################################################
## LISA Concepts
##################################################

# Baseline LISA concept
# SciRD document, reference LISA-LCST-SGS-TN-001, https://arxiv.org/abs/2108.01167
# hardcoded ACC and OMS noises
LISASciRDv1 = {
    'label' : 'LISA(SciRDv1)',
    'sqSacc_ASD' : [[3e-15,.4e-3*3e-15],[0,-1]],
    'sqSoms_ASD' : [[15e-12,15e-12*4.e-6],[0,-2]],
    'Lconst' : 2.5e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4
}

# LISA SciRD with low-level noise
# Replicates SciRDv1 but ACC and OMS noise are evaluated through the appropriate functions and not hardcoded 
LISASciRDLowLev = {
    'label' : 'LISA(SciRDLowLev)',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 2.0,
    'lambdaOMS' : 1064,
    'D_Tx' : 0.3,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : [[15e-12],[0]], 
    'Lconst' : 2.5e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4,
    'TMxOmega2' : -8e-7, 
    'OBGRSOmega2' : -7e-7, 
    'TMsize' : 0.046,
    'TMmat' : 'AuPt',
    'VacuumPressure' : 1e-6, 
    'ACCEL_other_ASD' : [[1e-18, 2e-15],[-1, 0]], # accounting for acceleration noise terms that we did not include in the model to match LISA SciRDv1
}

LISAold = {
    'label' : 'LISA(old 5Mkm arms)',
    #'sqSacc_ASD' : [[3e-15,.4e-3*3e-15],[0,-1]],
    'sqSoms_ASD' : [[15e-12,15e-12*4.e-6],[0,-2]],
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    #'sqSoms_func' : subsystems.OMS_Noise_PSD,
     'P_Tx' : 2.0,
    'lambdaOMS' : 1064,
    'D_Tx' : 0.3,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'Lconst' : 5e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4
}

# LISA CBE
LISACBE = {
    'label' : 'LISA(CBE)',
    'sqSacc_ASD' : [[3e-15,.4e-3*3e-15],[0,-1]],
    #'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 2.0,
    'lambdaOMS' : 1064,
    'D_Tx' : 0.3,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'Lconst' : 2.5e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4
}

##################################################
## Concepts under discussion in our study
##################################################

# LISA Grande
LISAGrande = {
    'label' : 'LISA Grande',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 4, # doubling LISA SciRD P_Tx (was 3 when we assumed LISA P_Tx = 1.5)
    'lambdaOMS' : 1064,
    'D_Tx' : 0.5,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'ACCEL_other_ASD' : [[1e-18, 2e-15],[-1, 0]], # accounting for acceleration noise terms that we did not include in the model to match LISA SciRDv1
    'Lconst' : 25e9, # 10 times bigger than LISA
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'Dsep' : 0,
    'SciDuration' : 4
}

# LISA-like mission in Folkner-like orbit
LISAU = {
    'label' : 'LISA-AU',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 4, # doubling LISA SciRD P_Tx (was 3 when we assumed LISA P_Tx = 1.5)
    'lambdaOMS' : 1550,
    'D_Tx' : 1.0,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'ACCEL_other_ASD' : [[1e-18, 2e-15],[-1, 0]], # accounting for acceleration noise terms that we did not include in the model to match LISA SciRDv1
    'Lconst' : 2.55e11,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'Dsep' : 0,
    'SciDuration' : 4
}

# Twin LISA
TwinLISA = {
    'label' : 'Twin LISA',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 2.0,
    'lambdaOMS' : 1064,
    'D_Tx' : 0.3,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'ACCEL_other_ASD' : [[1e-18, 2e-15],[-1, 0]], # accounting for acceleration noise terms that we did not include in the model to match LISA SciRDv1
    'Lconst' : 2.5e9,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 4, # adding an additional constellation doubles number of independent links
    'Dsep' : 1, # adding an additional constellation at 1 AU
    'SciDuration' : 4
}

# Baseline LISA concept GoBIGLowF
GoBIGLISA = {
    'label' : 'GoBIG LISA',
    'description' : 'Two LISA-like constellations on near-radial trajectories (in near opposing directions) to larger distance in solar system. This version uses a LISA-sized constellation',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 2.0,
    'lambdaOMS' : 1064,
    'D_Tx' : 0.3,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'ACCEL_other_ASD' : [[1e-18, 2e-15],[-1, 0]], # accounting for acceleration noise terms that we did not include in the model to match LISA SciRDv1
    # 'Lconst' : 5e10, # according to the spreadsheet
    'Lconst' : 2.5e9, # according to the LISA-like constellation
    'Dsep' : 30,
    'Rorbit' : 0.0,
    'Torbit' : 0.0,
    'Nindep' : 4,
    'SciDuration' : 4
}

# Baseline LISA concept GoBIGLowF
GoBIGLowF = {
    'label' : 'GoBig(LowF)',
    'description' : 'Two LISA-like constellations on near-radial trajectories (in near opposing directions) to larger distance in solar system. This version optimized for the low-ff side of LISA\'s band.',
    'sqSacc_ASD' : [[3e-15,.4e-3*3e-15],[0,-1]],
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 5,
    'lambdaOMS' : 1064,
    'D_Tx' : 1.0,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'Lconst' : 50e9,
    'Dsep' : 30,
    'Rorbit' : 0,
    'Torbit' : 0,
    'Nindep' : 2,
    'SciDuration' : 4
}

# Baseline LISA concept GoBIGLowF2
GoBIGLowF2 = {
    'label' : 'GoBig(LowF2)',
    'description' : 'Two LISA-like constellations on near-radial trajectories (in near opposing directions) to larger distance in solar system. This version optimized for the low-ff side of LISA\'s band.',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 5,
    'lambdaOMS' : 1064,
    'D_Tx' : 1.0,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'Lconst' : 50e9,
    'Dsep' : 30,
    'Rorbit' : 0,
    'Torbit' : 0,
    'Nindep' : 2,
    'TMsize' : .46,
    'TMmat'  : 'W',
    'SciDuration' : 4
}

# Baseline ALIA concept
ALIA = {
    'label' : 'ALIA',
    'sqSacc_ASD' : [[6e-16,.4e-3*6e-16],[0,-1]],
    'sqSoms_ASD' : [[5e-13,5e-13*4.e-6],[0,-2]],
    'Lconst' : 0.5e9,
    'lambdaOMS' : 1064, 
    'Responsivity' : 0.7,
    'P_Tx' : 10,
    'D_Tx' : 1.0,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4
}

# ALIA low-level concept
ALIAlowL = {
    'label' : 'ALIA low level',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 10,
    'lambdaOMS' : 1064,
    'D_Tx' : 1.0,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 0.05*10e-12, # LISA: 10e-12,
    'TMxOmega2' : .5*-8e-7, # LISA: -8e-7,
    'OBGRSOmega2' : .5*-7e-7, # LISA: -7e-7,
    'TMsize' : 2*.046, # LISA: .046
    'Lconst' : 0.5e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4,
}

# Baseline ALIA concept
ALIAbender = {
    'label' : 'ALIA Bender',
    'sqSacc_ASD' : [[3e-16,.4e-3*3e-16],[0,-1]],
    # 'sqSoms_ASD' : [[5e-13,5e-13*4.e-6],[0,-2]],
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'Lconst' : 0.5e9,
    'lambdaOMS' : 1064, 
    'Responsivity' : 0.7,
    # 'OMS_other_ASD': 1e-13,
    'P_Tx' : 30,
    'D_Tx' : 1.0,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4
}

# Baseline ALIA concept
ALIAcornish = {
    'label' : 'ALIA Cornish',
    'sqSacc_ASD' : [[3e-16,.4e-3*3e-16],[0,-1]],
    'sqSoms_ASD' : [[1e-13,1e-13*4.e-6],[0,-2]],
    # 'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'Lconst' : 0.5e9,
    'lambdaOMS' : 1064, 
    'Responsivity' : 0.7,
    'P_Tx' : 30,
    'D_Tx' : 1.0,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4
}


# Twin ALIA concept
ALIAtwin = {
    'label' : 'Twin ALIA',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 10,
    'lambdaOMS' : 1064,
    'D_Tx' : 1.0,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 0.05*10e-12, # LISA: 10e-12,
    'TMxOmega2' : .5*-8e-7, # LISA: -8e-7,
    'OBGRSOmega2' : .5*-7e-7, # LISA: -7e-7,
    'TMsize' : 2*.046, # LISA: .046
    'Lconst' : 0.5e9,
    'Dsep' : 1.0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 4,
    'SciDuration' : 4,
}

# Baseline LISA concept GoBIGLowF
GoBIGALIA = {
    'label' : 'GoBIG ALIA',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 10,
    'lambdaOMS' : 1064,
    'D_Tx' : 1.0,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 0.05*10e-12, # LISA: 10e-12,
    'TMxOmega2' : .5*-8e-7, # LISA: -8e-7,
    'OBGRSOmega2' : .5*-7e-7, # LISA: -7e-7,
    'TMsize' : 2*.046, # LISA: .046
    'Lconst' : 0.5e9,
    'Dsep' : 30,
    'Rorbit' : 0.0,
    'Torbit' : 0.0,
    'Nindep' : 4,
    'SciDuration' : 4
}

## other concepts not under discussion in our study, but useful as reference and comparison

# AMIGO
AMIGO = {
    'label' : 'AMIGO',
    'sqSacc_ASD' : [[3e-16,.4e-3*3e-16],[0,-1]],
    'sqSoms_ASD' : [[15e-13,15e-13*4.e-6],[0,-2]],
    'Lconst' : 2.5e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4
}

# Baseline DECIGO
DECIGO = {
    'label' : 'DECIGO',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 10,
    'lambdaOMS' : 515,
    'D_Tx' : 1.0,
    'Responsivity' : 0.7,
    # 'OMS_other_ASD' : .5e-12, #10e-12,
    'Lconst' : 1e6,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Dsep' : 1.73,
    'Nindep' : 6, # not sure how to specify the 4th constellation which is co-located with one of the other three, so skipping it.
    'SciDuration' : 4
}

# LISAMax
LISAMax = {
    'label' : 'LISAMax',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 1.0,
    'lambdaOMS' : 1064,
    'D_Tx' : 1.0,
    # 'Responsivity' : 0.7, # "optical_efficiency": 0.3,
    'QE' : 0.3,
    'Lconst' : 259e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4,
    'TMsize' : 0.046 * 2 ** (1 / 3), # Assuming a test mass twice as massive as LISA
    'TMmass': 1.928 * 2,
    'VacuumPressure' : 2e-6 / 10, # Assuming a pressure 10 times lower than the LISA allocation
    'ACCEL_other_ASD' : [[0],[0]],
}

# mu-Ares
muAres = {
    'label' : 'muAres',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 10.0,
    'lambdaOMS' : 1064,
    'D_Tx' : 1.0,
    # 'Responsivity' : 0.7, # "optical_efficiency": 0.3,
    'QE' : 0.5,
    'VacuumPressure' : 2e-6 / 10, # Assuming a pressure 10 times lower than the LISA allocation
    'Lconst' : 430e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4,
}

# DO-conservative
DOcons = {
    'label' : 'DOcons',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 10.0,
    'lambdaOMS' : 532,
    'D_Tx' : 1.0,
    # 'Responsivity' : 0.7, # "optical_efficiency": 0.3,
    'QE' : 0.5,
    'TMsize' : 0.046 * 2 ** (1 / 3), # Assuming a test mass twice as massive as LISA
    'TMmass': 1.928 * 2,
    'VacuumPressure' : 2e-6 / 10, # Assuming a pressure 10 times lower than the LISA allocation
    'Lconst' : 1e8,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4,
}

menuNames='LISASciRDv1,LISACBE,TwinLISA,LISASciRDLowLev,LISAGrande,LISAU,GoBIGLISA,GoBIGLowF,GoBIGLowF2,ALIA,ALIAbender,ALIAcornish,ALIAtwin,ALIAlowL,GoBIGALIA,DECIGO,AMIGO'.split(',')

menu={name:globals()[name] for name in menuNames}    
