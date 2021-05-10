import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# simple, early phenom frequency-space waveform according to Ajith et al. (2008)
# https://arxiv.org/abs/0710.2335
# Just calculating the amplitude and phase here (latter is relative, so you have to choose an offset)

def AmpLow(f,C,f1):

    return C*(f/f1)**(-7.0/6.0)


def AmpMid(f,C,f1):

    return C*(f/f1)**(-2.0/3.0)


def AmpHigh(f,C,f1,f2,sigma):

    w = (np.pi*sigma/2.0)*(f2/f1)**(-2.0/3.0)
    return C*w*(1.0/2.0/np.pi)*sigma/( (f-f2)**2 + sigma**2/4.0)


##### Main
def getIntrinsic(f,M,eta,d=1.0):


    a1 = 2.9740e-1
    a2 = 5.9411e-1
    a3 = 5.0801e-1
    a4 = 8.4845e-1

    b1 = 4.4810e-2
    b2 = 8.9794e-2
    b3 = 7.7515e-2
    b4 = 1.2848e-1

    c1 = 9.5560e-2
    c2 = 1.9111e-1
    c3 = 2.2369e-2
    c4 = 2.7299e-1

    #alpha1
    fmerg = (a1*eta**2 +b1*eta + c1)/np.pi/M

    fring = (a2*eta**2 +b2*eta + c2)/np.pi/M

    sigma = (a3*eta**2 +b3*eta + c3)/np.pi/M

    fcut  = (a4*eta**2 +b4*eta + c4)/np.pi/M

    # Eq. 4.17
    C = M**(5.0/6.0)*fmerg**(-7.0/6.0)*np.sqrt(5.0*eta/24.0)/d/np.pi**(2.0/3.0)

    freqs = f

    amp1_vals = AmpLow(freqs,C,fmerg)
    amp2_vals = AmpMid(freqs,C,fmerg)
    amp3_vals = AmpHigh(freqs,C,fmerg,fring,sigma)

    f_in_bin1 = freqs < fmerg
    f_in_bin2 = (freqs >= fmerg) & (freqs < fring)
    f_in_bin3 = freqs >= fring

    amp_vals = f_in_bin1*amp1_vals + f_in_bin2*amp2_vals + f_in_bin3*amp3_vals

    # phase calculation
    x0 =  1.7516e-1
    x2 = -5.1571e1
    x3 =  6.5866e2
    x4 = -3.9031e3
    x6 = -2.4874e4
    x7 =  2.5196e4

    y0 =  7.9483e-2
    y2 = -1.7595e1
    y3 =  1.7803e2
    y4 = -7.7493e2
    y6 = -1.4892e3
    y7 =  3.3970e2

    z0 = -7.2390e-2
    z2 =  1.3253e1
    z3 = -1.5972e2
    z4 =  8.8195e2
    z6 =  4.4588e3
    z7 = -3.9573e3

    psi0 = (x0*eta**2 + y0*eta + z0)/eta/(np.pi*M)**( (5.0-0)/3.0)
    psi2 = (x2*eta**2 + y2*eta + z2)/eta/(np.pi*M)**( (5.0-2)/3.0)
    psi3 = (x3*eta**2 + y3*eta + z3)/eta/(np.pi*M)**( (5.0-3)/3.0)
    psi4 = (x4*eta**2 + y4*eta + z4)/eta/(np.pi*M)**( (5.0-4)/3.0)
    psi6 = (x6*eta**2 + y6*eta + z6)/eta/(np.pi*M)**( (5.0-6)/3.0)
    psi7 = (x7*eta**2 + y7*eta + z7)/eta/(np.pi*M)**( (5.0-7)/3.0)

    t0 = 0.0
    phi0 = -5400.0

    phase_vals = 2.0*np.pi*freqs*t0 + phi0
    phase_vals = phase_vals + psi0*freqs**( (0-5.0)/3.0)
    phase_vals = phase_vals + psi2*freqs**( (2-5.0)/3.0)
    phase_vals = phase_vals + psi3*freqs**( (3-5.0)/3.0)
    phase_vals = phase_vals + psi4*freqs**( (4-5.0)/3.0)
    phase_vals = phase_vals + psi6*freqs**( (6-5.0)/3.0)
    phase_vals = phase_vals + psi7*freqs**( (7-5.0)/3.0)

    h = amp_vals*np.exp(1j*phase_vals)

    return h

