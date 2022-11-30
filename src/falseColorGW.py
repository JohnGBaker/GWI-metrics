# These are functions used to make flase-color GW images from a set of individual GW sources


# modules
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from astropy.visualization import make_lupton_rgb
#import astropy.constants as constants

# from lisacattools.catalog import GWCatalogs, GWCatalogType
# from lisacattools import convert_ecliptic_to_galactic, HPhist
# from lisacattools import OFF # could be CRITICAL, FATAL, ERROR, WARNING, TRACE, INFO, DEBUG, NOTSET
# logger = logging.getLogger("lisacattools")
# logger.setLevel(OFF) 



# Define function for making super-ensemble of "GW photons" from catalog of sources
def makeSuperEnsemble(catalogIn, cdlColName = 'Diffraction Limit',snrColName = 'SNR', Nsamp = 1e3):
    # load in the catalog.  Check that the right fields exist (amplitude, freuqency, and sky positions)
    # loop through catalog entries, for each catalog draw from a Gaussian in sky position
    # concatenate to make superensemble
    
    from scipy.special import expit
    
    
    source_chain = list()
    # pick out the sources and start looping
    sources = list(catalogIn.index)

    for source in sources:
        N = int(np.ceil(Nsamp*expit(catalogIn.loc[source][snrColName]-8)))
        
        source_lats = np.random.normal(catalogIn.loc[source]['Galactic Latitude'], catalogIn.loc[source][cdlColName],N)
        source_lons = np.random.normal(catalogIn.loc[source]['Galactic Longitude'], catalogIn.loc[source][cdlColName],N)
        source_amp = np.zeros_like(source_lats)+catalogIn.loc[source]['Amplitude']
        source_freq = np.zeros_like(source_lats)+catalogIn.loc[source]['Frequency']
        source_df = pd.DataFrame.from_dict({'Amplitude' : source_amp, 'Frequency' : source_freq, 'Galactic Latitude' : source_lats, 'Galactic Longitude' : source_lons})
        source_chain.append(source_df)

    ensembleOut = pd.concat(source_chain)
    return(ensembleOut)

# define function for making frequency v intensity 
def makeGWluminosity(ensemble, fmax = 12e-3):
    # Map the frequencies to a color, take the rgb values of that color:
    ensemble['Power'] = ensemble['Amplitude']**2*ensemble['Frequency']**2
    cmap_name = 'gist_rainbow'
    cmap = plt.get_cmap(cmap_name)
    norm = mpl.colors.Normalize(vmax=fmax)
    scalarMap = cm.ScalarMappable(norm=norm, cmap=cmap)
    fmap = scalarMap.to_rgba(ensemble['Frequency'])
    f_bins = np.linspace(0,fmax,255)
    
    # return ensemble with frequency
    return ensemble, f_bins, fmap

    
    

# Define function for making Image from SuperEnsemble
def makeGWimage(ensemble, fmap, dynamicRange=8, skyRange = [-180,180,-90,90], skyBins = [512,1024]):
    
    # Do sky area binning
    Nlat = skyBins[1]
    coscolat_range = np.cos((90-np.array([skyRange[2],skyRange[3]]))*np.pi/180)
    colat_edges = np.arccos(np.linspace(coscolat_range[0],coscolat_range[1],Nlat))*180/np.pi
    lat_edges = 90 - colat_edges
    lat_centers = lat_edges[:-1]+0.5*(lat_edges[1:]-lat_edges[:-1])
    
    Nlon = skyBins[0]
    lon_edges = np.linspace(skyRange[0],skyRange[1],Nlon)
    lon_centers = lon_edges[:-1]+0.5*(lon_edges[1:]-lon_edges[:-1])

    
    # Generate an rgbi array for the image: 
    fcol = np.empty((4,(len(lon_edges)-1)*(len(lat_edges)-1)))
    for id in range(0,4): 
        fnew = fmap[:,id]*ensemble['Power']
        ftmp, xedges, yedges = np.histogram2d(ensemble['Galactic Longitude'],ensemble['Galactic Latitude'], bins=(lon_edges,lat_edges), normed=False, weights=fnew)
        fcol[id,:] = ftmp[:].ravel()

    # Scale to max value and reshape
    fcol /= np.max(fcol[3,:])
    fcol[3,:] = 1.
    fcol = fcol.reshape((4,len(lon_edges)-1,len(lat_edges)-1))
    fcol = np.transpose(fcol,(0,2,1))
    
    ftmp = np.transpose(fcol,(1,2,0))
    fmin = np.min(ftmp[ftmp.nonzero()])
    ftmp += fmin * 1e-10
    fmax = np.max(ftmp)
    ftmp /= fmax

    flog = np.log10(ftmp)
    fout =  (flog/dynamicRange)+1
    fout[fout < 0.] = -1. 
    fout[fout > 1.] = 1. 
    return(fout)

