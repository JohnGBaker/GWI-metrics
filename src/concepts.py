#This file provides descriptors for reference GW imager mission concepts.
import subsystems

#Baseline LISA concept
LISASciRDv1 = {
    'label' : 'LISA(SciRDv1)',
    'sqSacc_level' : 3e-15,
    'sqSoms_level' : 15e-12,
    'Lconst' : 2.5e9,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2
}
#Baseline LISA concept
LISASciRDLowLev = {
    'label' : 'LISA(SciRDLowLev)',
    'sqSacc_func' : subsystems.ACC_Noise_PSD,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 1.5,
    'lambdaOMS' : 1064,
    'D_Tx' : 0.3,
    'Responsivity' : 0.7,
    'OMS_other_PSD' : 15e-12,
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
    'ACCEL_other_PSD' : 0
}

#LISA Grande
LISAGrande = {
    'label' : 'LISAGrande',
    'sqSacc_level' : 1e-15,
    'sqSoms_level' : 18e-12,
    'Lconst' : 5e10,
    'Dsep' : 0,
    'Rorbit' : 1.0,
    'Torbit' : 1.0,
    'Nindep' : 2
}

#Baseline LISA concept GoBIGLowF
GoBIGLowF = {
    'label' : 'GoBig(LowF)',
    'description' : 'Two LISA-like constellations on near-radial trajectories (in near opposing directions) to larger distance in solar system. This version optimized for the low-ff side of LISA\'s band.',
    'sqSacc_level' : 3e-15,
    'sqSoms_func' : subsystems.OMS_Noise_PSD,
    'P_Tx' : 5,
    'lambdaOMS' : 1064,
    'D_Tx' : 1.0,
    'Responsivity' : 0.7,
    'OMS_other_PSD' : 10e-12,
    'Lconst' : 50e9,
    'Dsep' : 30,
    'Rorbit' : 0,
    'Torbit' : 0,
    'Nindep' : 2
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
    'OMS_other_PSD' : 10e-12,
    'Lconst' : 50e9,
    'Dsep' : 30,
    'Rorbit' : 0,
    'Torbit' : 0,
    'Nindep' : 2,
    'TMsize' : .46,
    'TMmat'  : 'W'
}
