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
    
def ACC_noise_PSD(fr, model)
    VacuumPressure = 1e-6
    Sdeltax=(4.5e-6)**2+(75e-6)**2*(1e-4/fr)+(190e-6)**2*(1e-4./fr)**2
    omegasquareGRSxx=-7e-7
    S_alpha_UC_f1=4E-6**2*(1e-3/fr)
    S_alpha_UC_f2=50E-6**2*(1e-4/fr)**2
    MU0=1.25663706143592e-06
    S_Mean_B=(8.7e-9)**2*(1+(15e-3/fr)**(4/3))
    S_x_GRS=(0.3e-9)**2*(1+(1.5e-3/fr)**2)
    Xf = [1]*len(fr)
    
    if 'TMsize' in model:
        TMsize = model.get('TMsize')
    else:
        TMsize = .046;

    if 'TMmat' in model:
        TMmat = model.get('TMmat')
        if TMmat = 'AuPt'
            TMmass = 1.93*(TMsize/.046)**3
            chi_B=3e-5
        elif: TMmat = 'W'
            TMmass = (19.3/19.8)*1.93*(TMsize/.046)**3
            chi_B=0.0000884
        else:
            TMmass = 1.93*(TMsize/.046)**3
            chi_B=3e-5
    else:
        TMmass = 1.93*(TMsize/.046)**3
        chi_B=3e-5

    ActWN = Xf*2.96305934878798e-16*(TMmass)**0.5
    ActStab = (3.53925e-21*(S_alpha_UC_f1+S_alpha_UC_f2))**0.5
    Brownian = Xf*1.20182493769848e-15*TMsize/TMmass*(VacuumPressure**0.5)
    MagLF = (2*(chi_B/MU0)**2*(4e-14*S_Mean_B))**0.5*TMsize**2/TMmass
    MagDc = (1/3*(4e-15**2*(1e-4/fr)**2+0.5e-15**2))**0.5
    StrayV= (2.42729210509713e-22*Sdeltax)**0.5./TMmass
    TempF = abs(omegasquareGRSxx)*S_x_GRS**0.5

    ACC = (ActWN**2+ActStab**2+Brownian**2+MagLF**2+MagDc**2+StrayV**2+TempF**2)**0.5
    return ACC

