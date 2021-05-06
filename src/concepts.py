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
    'sqSacc_level' : 3e-15,
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
    'Nindep' : 2
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
