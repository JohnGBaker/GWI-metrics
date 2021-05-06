#These are functions for defining GW mission subsystems. They are used in sensitivity calculations and possibly other places.
import constants
import numpy as np

#OMS Optical Metrology System
#The initial version of this captures the contents of a spreadsheet that Jeff put together

#OMS parameters
'''
P_Tx          | W           | Transmitted power
lambdaOMS     | nm          | Laser wavelength
D_Tx          | m           | Transmitter mirror diameter
D_Rx          | m           | Receiver mirror diameter (optional, otherwise=D_Tx)
Responsivity  | A/W         | Photodetector responsivity [need Resp or QE]
QE            | -           | Ouantum Efficiency [need Resp or QE]
OMS_other     | pm/sqrt(Hz) | Other sources of OMS noise to add in quadrature
'''

def OMS_Noise_PSD(fr, model):
    '''
    Return OMS subsystem noise power in m^/Hz on freqs in fr
    Parameters in model dictionary should include:
    P_Tx          | W           | Transmitted power
    lambdaOMS     | nm          | Laser wavelength
    D_Tx          | m           | Transmitter mirror diameter
    D_Rx          | m           | Receiver mirror diameter (optional, otherwise=D_Tx)
    Responsivity  | A/W         | Photodetector responsivity [need Resp or QE]
    QE            | -           | Ouantum Efficiency [need Resp or QE]
    OMS_other_PSD | m/sqrt(Hz)  | Other sources of OMS noise to add in quadrature  
    '''
    P_Tx=model['P_Tx']
    lambdaOMS=model['lambdaOMS']*1e-9
    D_Tx=model['D_Tx']
    if 'D_Rx' in model: D_Rx=model['D_Rx']
    else: D_Rx=D_Tx
    L_arm=model['Lconst']

    if ( 'QE' in model and 'Responsivity' in model ) or not ( 'QE' in model or 'Responsivity' in model ):
        raise ValueError('Need QE or Responsivity, not both')
    if 'Responsivity' in model:
        Responsivity=model['Responsivity']
        QE=Responsivity*constants.h_planck*constants.c/constants.e_charge/lambdaOMS
    if 'QE' in model:
        QE=model['QE']
    
    Sn_shot=constants.h_planck*constants.c*lambdaOMS**3*L_arm**2/(2*np.pi**2*QE*P_Tx*D_Tx**2*D_Rx**2) #Jeff's formula
    #sqSn_shot=constants.h_planck*constants.c*lambdaOMS**3*L_arm**2/(QE*P_Tx*D_Tx**2*D_Rx**2) #Differs by 2pi to agree with LISA value
    print('QE:',QE)
    
    print('P_Rx:',OMS_received_power(model))
    print('sqSn_shot:',np.sqrt(Sn_shot))
    
    OMS_other=0
    if 'OMS_other_PSD' in model: OMS_other=model['OMS_other_PSD']
    
    Sn = Sn_shot + OMS_other**2
    return Sn

def OMS_received_power(model):
    P_Tx=model['P_Tx']
    lambdaOMS=model['lambdaOMS']*1e-9
    D_Tx=model['D_Tx']
    if 'D_Rx' in model: D_Rx=model['D_Rx']
    else: D_Rx=D_Tx
    L_arm=model['Lconst']
    
    radius_over_omega = 2.24
    scaling_const = 0.5*np.pi**2/2/radius_over_omega**2

    P_Rx=scaling_const*P_Tx*D_Tx**2*D_Rx**2/(L_arm*lambdaOMS)**2
    return P_Rx
    
           
