#This file provides descriptors for reference GW imager mission concepts.
import subsystems


#Baseline LISA concept
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

#AMIGO
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

#LISA CBE
LISACBE = {
    'label' : 'LISA(CBE)',
    'sqSacc_ASD' : [[3e-15,.4e-3*3e-15],[0,-1]],
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 1.5,
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

#Twin LISA
TwinLISA = {
    'label' : 'Twin LISA',
    'sqSacc_ASD' : [[3e-15,.4e-3*3e-15],[0,-1]],
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 1.5,
    'lambdaOMS' : 1064,
    'D_Tx' : 0.3,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'Lconst' : 2.5e9,
    'Dsep' : 1,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 4,
    'SciDuration' : 4
}

#LISA SciRD with low-level noise
LISASciRDLowLev = {
    'label' : 'LISA(SciRDLowLev)',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 1.5,
    'lambdaOMS' : 1064,
    'D_Tx' : 0.3,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : [[15e-12],[0]],
    'Lconst' : 2.5e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'TMxOmega2' : -8e-7,
    'OBGRSOmega2' : -7e-7,
    'TMsize' : 0.046,
    'TMmat' : 'AuPt',
    'VacuumPressure' : 1e-6,
    'ACCEL_other_ASD' : [[0],[0]],
    'SciDuration' : 4
}

LISAGrande = {
    'label' : 'LISA Grande',
    'sqSacc_ASD' : [[3e-15,.4e-3*3e-15],[0,-1]],
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 3,
    'lambdaOMS' : 1064,
    'D_Tx' : 0.5,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'Lconst' : 25e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4
}

# LISA-like mission in Folkner-like orbit
LISAU = {
    'label' : 'LIS-AU',
    'sqSacc_ASD' : [[3e-15,.4e-3*3e-15],[0,-1]],
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 3,
    'lambdaOMS' : 1550,
    'D_Tx' : 1.0,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'Lconst' : 2.6e11,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4
}

#Baseline LISA concept GoBIGLowF
GoBIGLISA = {
    'label' : 'GoBig(LISA)',
    'description' : 'Two LISA-like constellations on near-radial trajectories (in near opposing directions) to larger distance in solar system. This version uses a LISA-sized constellation',
    'label' : 'GOBIG',
    'sqSacc_ASD' : [[3e-15,.4e-3*3e-15],[0,-1]],
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 1.5,
    'lambdaOMS' : 1550,
    'D_Tx' : 0.3,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : 10e-12,
    'Lconst' : 5e9,
    'Dsep' : 30,
    'Nindep' : 4,
    'SciDuration' : 4
}
    


#Baseline LISA concept GoBIGLowF
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

#Baseline LISA concept GoBIGLowF2
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

#Baseline ALIA concept
ALIA = {
    'label' : 'ALIA',
    'sqSacc_ASD' : [[6e-16,.4e-3*6e-16],[0,-1]],
    'sqSoms_ASD' : [[5e-13,5e-13*4.e-6],[0,-2]],
    'Lconst' : 0.5e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4
}

#Twin ALIA concept
ALIAtwin = {
    'label' : 'ALIAtwin',
    'sqSacc_ASD' : [[6e-16,.4e-3*6e-16],[0,-1]],
    'sqSoms_ASD' : [[5e-13,5e-13*4.e-6],[0,-2]],
    'Lconst' : 0.5e9,
    'Dsep' : 1,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 4,
    'SciDuration' : 4
}

#ALIA low-level concept
ALIAlowL = {
    'label' : 'ALIA low level',
    'Lconst' : 0.5e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2,
    'SciDuration' : 4,
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 10,
    'lambdaOMS' : 1064,
    'TMxOmega2' : .4*-8e-7, # -8e-7,
    'OBGRSOmega2' : .4*-7e-7, # -7e-7,
    'D_Tx' : 1.0,
    'Responsivity' : 0.7,
    'OMS_other_ASD' : .5e-12, #10e-12,
    'VacuumPressure' : 1e-6,
    'TMsize' : .046*2.5 # .046
}

#Baseline LISA concept GoBIGLowF
GoBIGALIA = {
    'label' : 'GoBig(ALIA)',
    'description' : 'Two ALIA-like constellations on near-radial trajectories (in near opposing directions) to larger distance in solar system. This version uses a LISA-sized constellation',
    'sqSacc_ASD' : [[6e-16,.4e-3*6e-16],[0,-1]],
    'sqSoms_ASD' : [[5e-13,5e-13*4.e-6],[0,-2]],
    'Lconst' : 0.5e9,
    'Dsep' : 30,
    'Nindep' : 4,
    'SciDuration' : 4
}

#Baseline DECIGO
DECIGO = {
    'label' : 'DECIGO',
    'sqSacc_ASD' : [[2e-18,.4e-3*2e-18],[0,-1]],
    'sqSoms_ASD' : [[15e-18,15e-18*4.e-6],[0,-2]],
    'Lconst' : 10e6,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Dsep' : 1.73,
    'Nindep' : 6, #not sure how to specify the 4th constellation which is co-located with one of the other three, so skipping it.
    'SciDuration' : 4
}

menuNames='LISASciRDv1,LISACBE,TwinLISA,LISASciRDLowLev,LISAGrande,LISAU,GoBIGLISA,GoBIGLowF,GoBIGLowF2,ALIA,ALIAtwin,ALIAlowL,GoBIGALIA,DECIGO,AMIGO'.split(',')
menu={name:globals()[name] for name in menuNames}    
