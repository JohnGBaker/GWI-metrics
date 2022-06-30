################################################################################
# These are functions for defining GW mission subsystems.
# They are used in sensitivity calculations and possibly other places.
# Written by: J Slutsky
# Some edits by: E Castelli 06/2022
################################################################################
import constants
import numpy as np


#### OMS Optical Metrology System
# The initial version of this captures the contents of a spreadsheet that Jeff put together
# [Ref. to spreadsheet?]

# OMS parameters
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
    OMS_other_ASD | m/sqrt(Hz)  | Other sources of OMS noise to add in quadrature
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
    if 'OMS_other_ASD' in model:OMS_other=model['OMS_other_ASD']
    OMS_other_ASD = F_Noise_PSD(fr,OMS_other)

    Sn = Sn_shot + OMS_other_ASD**2
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

#### TM Acceleration noise subsystem
# The initial version of this is derived from the LISA Performance model documents
# Ref. LISA-LCST-INST-TN-003 and LISA-LCST-INST-RP-003 and related references

# Acceleration parameters
# do we also want to include actuation parameters or other stuff here?
'''
Return ACCEL subsystem noise power in m^2/s^4/Hz on freqs in fr
Parameters in model dictionary should include:
TMxOmega2          | s^-2           | Test Mass x-axis stiffness
OBGRSOmega2        | s^-2           | OB/GRS TMx stiffness (?)
TMsize             | m              | TM linear dimension
TMmat              |                | TM material composition
VacuumPressure     | Pa             | Residual gas vacuum presure
ACCEL_other_ASD    | m/s^2          | Catch-all unknown acceleration noise input
'''

def ACC_Noise_PSD(fr, model):

    # make unit array of same dimensions as frequency array
    Xf = F_Noise_PSD(fr,[[1],[0]])

    # INITIALIZE internal functions and values
    MU0 = constants.MU0
    KB = constants.KB
    H2Omo = constants.H20mo

    # READ PARAMETERS FROM INPUT
    # Test Mass x-axis stiffness
    omegasquarexx = -8e-7
    if 'TMxOmega2' in model:
        omegasquarexx = model.get('TMxOmega2');
    # OB/GRS TMx stiffness (?)
    omegasquareGRSxx=-7e-7
    if 'OBGRSOmega2' in model:
        omegasquareGRSxx = model.get('OBGRSOmega2');
    # TM cube linear dimension
    if 'TMsize' in model:
        TMsize = model.get('TMsize')
    else:
        TMsize = .046;
    # TM material
    TMmassAu = 1.93*(TMsize/.046)**3
    if 'TMmat' in model:
        TMmat = model.get('TMmat')
        if TMmat == 'AuPt':
            TMmass = TMmassAu
            chi_B=3e-5
        elif TMmat == 'W':
            TMmass = (19.3/19.8)*TMmassAu
            chi_B=0.0000884
        else:
            TMmass = TMmassAu
            chi_B=3e-5
    else:
        TMmass = TMmassAu
        chi_B=3e-5
    # Vacuum Pressure
    VacuumPressure = 1e-6
    if 'VacuumPressure' in model:
        VacuumPressure=model.get('VacuumPressure')
    # Margin/misc
    ACCEL_other=0
    if 'ACCEL_other_ASD' in model:ACCEL_other=model['ACCEL_other_ASD']

    ### Parameters from the Perf. model and NOT read from input
    # + Noise spectra of known quantities used in the evaluation of the single noise terms
    # Ref. LISA-LCST-INST-RP-003

    ## Actuation noise
    # electrodes effective armlength between x-face electrodes
    R_star = 0.033 # ref. ESA-L3-EST-INST-DD-001 v2.2, LISA Payload Definition Document and upcoming LPF actuation paper
    # partial derivative of C_x TM w.r.t. x
    DCX_DX_GRS1 = 2.91139e-10 # ref. ESA-L3-EST-INST-DD-001 and J.Phys.Conf.Ser.154:012008
    # partial derivative of C_x EH w.r.t. x
    DCXH_DX_GRS1 = -6.9675e-11 # ref. ESA-L3-EST-INST-DD-001 and J.Phys.Conf.Ser.154:012008
    # maximum torque acting on phi
    gamma_c = 1e-9 # Ref. LISA-UTN-INST-TN-006 v1 grav specs analysis and upcoming LPF actuation paper
    # authority margin on phi
    authority_margin = 1.5 # Ref. LISA-UTN-INST-TN-006 v1 grav specs analysis and upcoming LPF actuation paper
    # Uncorrelated gain fluctuations (1/f)
    # Ref. upcoming LPF actuation paper +  J.Phys.Conf.Ser. 840:012006 (2017) + CQG 33:235015 (2016)
    S_alpha_UC_f1 = 4E-6**2*(1e-3/fr)
    # Uncorrelated gain fluctuations (1/f**2)
    # Ref. upcoming LPF actuation paper +  J.Phys.Conf.Ser. 840:012006 (2017) + CQG 33:235015 (2016)
    S_alpha_UC_f2 = 50E-6**2*(1e-4/fr)**2
    # additive actuation voltage noise
    S_V = (2e-6)**2 # ref. ESA-L3-EST-INST-DD-001 and J.Phys.Conf.Ser.154:012008

    ## Brownian noise
    # Brownian noise amplification factor due to the proximity of the EH surfaces to the TM
    alphanoise = 13 # Ref. PRL 108:140601 + PRD 84:063007
    # GRS Temperature
    T_GRS = 293.15

    ## Stray electrostatics noise
    # Total TM capacitance
    CTOT = 3.42e-11 # Ref. PRL 108:181101
    # Number of elementary Test Mass charges
    Nq = 15e6 # Ref. PRL 108:181101
    # effective single electrode noise which is a sum of three terms depending on [0 -1 -2] powers of the frequency
    # Sdeltax = A*VTM_meansq*f^0 + B*VTM_meansq*f^-1 + C*VTM_meansq*f^-2
    # VTM_meansq is the mean squared value of TM voltage evaluated from the RMS charge measured on the spacecraft
    # A, B, C come from a charged TM voltage noise test (need to recover the pest values by running a Matlab script)
    # Ref. PRL 108:181101
    Sdeltax = (10e-6)**2 + (75e-6)**2*(1e-4/fr) + (190e-6)**2*(1e-4/fr)**2

    ## Magnetic quantities
    # Square modulus of the difference of static (average over one day) magnetic field between the X faces of the TM, averaged over x TM face
    delta_B2_DC = 4e-14 # Ref. lisa-utn-inst-tn-015 RevMagReqs + OHB technical reports (don't understand if there is smtg published or in preparation)
    # Mean of the fluctuating magnetic field of the two X faces of the TM [T^2 Hz^-1]
    S_Mean_B = (8.7e-9)**2*((15e-3/fr)**(4/3)) # Ref. lisa-utn-inst-tn-015 RevMagReqs + OHB technical reports (don't understand if there is smtg published or in preparation)
    # Maximum absolute value of each component of the static magnetic field at any point on the x faces of the TM
    Bmax = 2.5e-13 # Ref. lisa-utn-inst-tn-015 RevMagReqs
    # Maximum noise of any component of the magnetic field generated by sources on board the SC
    S_dB_max = (1e-9)**2*(1 + (15e-3/fr)**(4/3)) # Ref. lisa-utn-inst-tn-015 RevMagReqs

    ## Temperature fluctuation noise
    # GRS-OB baseline deformation contribution to x jitter
    S_x_GRS = (0.3e-9)**2*(1+(1.5e-3/fr)**2)
    # TestMass Jitter along x with respect to MOSA
    Sx_tm = (0.95e-9)**2*(1+(2e-4/fr)**2)/(1+(fr/8e-3/(10**0.5))**4);

    #######################
    # CALCULATE NOISE TERMS
    #######################
    ## Actuation noises
    # Actuation white noise
    #ActWN = sqrt(2./TMmass.*(DCX_DX_GRS1+DCXH_DX_GRS1).*(R_star.*authority_margin.*gamma_c).*S_V)
    ActWN = Xf*(2/TMmass.*(DCX_DX_GRS1+DCXH_DX_GRS1)*(R_star*authority_margin*gamma_c)*S_V)**0.5
    # ActWN = Xf*2.96305934878798e-16/(TMmass)**0.5                           #M**-0.5

    # Actuation stability noise
    # = sqrt(R_star.^2.*(gamma_c.^2+(authority_margin.*gamma_c).^2).*(S_alpha_UC_f1+S_alpha_UC_f2));
    ActStab = (R_star**2*(gamma_c**2+(authority_margin*gamma_c)**2)*(S_alpha_UC_f1+S_alpha_UC_f2))**0.5
    # ActStab = (3.53925e-21*(S_alpha_UC_f1+S_alpha_UC_f2))**0.5              #[]gain fluctuations, should this scale like Act white noise? (it's a force, right?)

    ## Brownian noise
    # = sqrt(alphanoise.*(1+pi./8).*(TMsize./TMmass).^2.*VacuumPressure.*(512.*H2Omo.*KB.*T_GRS/pi).^0.5)
    Brownian = (alphanoise*(1+np.pi/8)*(TMsize/TMmass)**2*VacuumPressure*(512*H2Omo*KB*T_GRS/np.pi)**0.5)**0.5
    # Brownian = Xf*5.0424e-11*TMsize/TMmass*(VacuumPressure**0.5)            #M**-1 S

    ## Stray electrostatics noise
    # Low frequency fluctuations in average stray electrostatic fields arise from the interaction of TM charge with the residual
    # effective DC potential bias applied on a single GRS electrode, which explains the residual sensitivity of the electrodes to a change in the TM charge.
    # Therefore, the presence of a non-zero residual charge on the TM would couple to the effective potential on the surface of single x-axis electrode,
    # giving rise to a noisy force acting on the TM. The surface effect is companioned by fluctuations in the electrode potential generated
    # by the Front End Electronics actuation, which will add to the stray voltages.
    # = sqrt(1./TMmass.^2 .* ((DCX_DX_GRS1 + DCXH_DX_GRS1)./CTOT .* Nq .* E).^2.* Sdeltax)
    StrayV = (1/TMmass^2 * ((DCX_DX_GRS1 + DCXH_DX_GRS1)/CTOT * Nq * E)**2 * Sdeltax)
    # StrayV = (2.42729210509713e-22*Sdeltax)**0.5/TMmass                      #M**-1 maybe need TMsize scaling?

    ## Magnetic noises
    # Low frequency interplanetary B
    MagLF = (2/TMmass**2*(chi_B*TMsize**2/MU0)**2*(delta_B2_DC*S_Mean_B))**0.5
    # MagLF = (2*(chi_B/MU0)**2*(4e-14*S_Mean_B))**0.5*TMsize**2/TMmass       #S M**-1
    # noise from local source of B field fluctuations (still low freq)
    # (just added this because it's comparable to MagLF)
    # Ref. lisa-utn-inst-tn-015 RevMagReqs
    MAGfluct = 3./TMmass**2*(2*chi B*TMsize**2/MU0)**2 *(Bmax * S_dB_max))**0.5
    #  Down-conversion of audio-frequency magnetic fields into a low frequency force, via their coupling with the corresponding induced audio-frequency eddy currents
    # Ref. lisa-utn-inst-tn-015 RevMagReqs
    # = sqrt(1./3.*(4e-15.^2.*(1e-4./f).^2+0.5e-15.^2))
    MagDc = (1/3*(4e-15**2*(1e-4/fr)**2+0.5e-15**2))**0.5                   #[] ??

    ## Temperature force noise
    # OB/GRS baseline deformation
    # = sqrt(S_x_GRS).*abs(omegasquareGRSxx)
    TempF = abs(omegasquareGRSxx)*S_x_GRS**0.5                              #[] no obvious scaling, MSS material need to be included?
    # TM x stiffness
    # = sqrt(Sx_tm.*abs(omegasquarexx).^2)
    Xstiff = (Sx_tm*abs(omegasquarexx)**2)**0.5                             #[]TM jitter, should this be mass dependent?

    # any other noise
    AO = F_Noise_PSD(fr,ACCEL_other)

    # SUM SQUARES
    ACC = ActWN**2 + ActStab**2 + Brownian**2 + MagLF**2 + Magfluct**2 + MagDc**2 + StrayV**2 + TempF**2 + Xstiff**2 + AO**2
    return ACC

def F_Noise_PSD(fr, pLaws, QUAD=False):
    # Generate an output array with length equal to that of input frequency axis 'fr'. Output array values are determined by one
    # or more powerlaws.  If input pLaws is an int or float, then a constant amplitude output at that value is generated.  If
    # input is a two-element list, then the first element contains one or more amplitudes, and the second element contains the
    # same number of respective power laws. If the two lists are of different lengths, '0's are padded onto the end, for example
    # [[1] [0 2]] is the same as [[1 0] [0 2]], in effect [[1] [0]] which is a an array of where every element has the value '1'.

    if (type(pLaws) is int) or (type(pLaws) is float):
        freqPower = [0]
        freqAmp = [pLaws]
    elif type(pLaws) is list:
        freqPower = pLaws[:][1]
        freqAmp = pLaws[:][0]

    fMatch = np.size(freqPower) - np.size(freqAmp)
    if fMatch > 0: freqAmp += [0]*fMatch
    if fMatch < 0: freqPower += [0]*(-fMatch)

    outPSD = np.zeros(np.size(fr))
    if QUAD:
        for ii in range(np.size(freqPower)):
            outPSD += (freqAmp[ii]*fr**freqPower[ii])**2
            ii += 1
        outPSD = outPSD**0.5
    else:
        for ii in range(np.size(freqPower)):
            outPSD += freqAmp[ii]*fr**freqPower[ii]
            ii += 1

    return outPSD
