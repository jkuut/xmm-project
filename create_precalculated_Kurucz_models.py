import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import integrate, interpolate
from scipy.optimize import curve_fit


filters = ['gaia_BPmag', 'gaia_Gmag', 'gaia_RPmag',
                '2mass_Jmag', '2mass_Hmag', '2mass_Kmag', 
                'W1mag', 'W2mag', 'W3mag', 'W4mag']

filterList = ['gaia_BPmag', 'Gaia_e_BPmag', 'gaia_Gmag', 'Gaia_e_Gmag', 'gaia_RPmag', 'Gaia_e_RPmag', 
                '2mass_Jmag', '2mass_e_Jmag', '2mass_Hmag', '2mass_e_Hmag', '2mass_Kmag', '2mass_e_Kmag', 
                'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag', 'W3mag', 'e_W3mag', 'W4mag', 'e_W4mag']
                
centralWavelengths = [5035.75, 5822.39, 7619.96, 
                        12350.00, 16620.00, 21590.00, 
                        33526.00, 46028.00, 115608.00, 220883.00]
                        
                        
def magToFlux(datapoints):

    # F = F0 * 10^(-mag/2.5), for erg/cm2/s/A F0 is as below
    zeroPoints = np.array([4.08e-9, 2.50e-9, 1.27e-9, 3.13e-10, 1.13e-10, 4.28e-11, 
                            8.18e-12, 2.42e-12, 6.52e-14, 5.09e-15])

    fluxes = zeroPoints * np.power(10, -np.array(datapoints)/2.5)
    return fluxes


def get_fnus(wavel, flux, response):

    filter_wavel, filter_response = np.loadtxt("responses/" + response + ".dat", unpack=True)
    fint = interpolate.interp1d(wavel, flux)
    fluxInFilterRange = fint(filter_wavel)
    response2 = filter_response / integrate.simps(filter_response, filter_wavel)    
    fnu = integrate.simps(fluxInFilterRange*response2, filter_wavel)
    return fnu
    
    

def get_Kurucz_data(temperature, metallicity, logg):

    metaldict = {0.5: "p05", 0.2: "p02", 0.0: "p00", -0.5: "m05", -1.0: "m10", -1.5: "m15", -2.0: "m20", -2.5:"m25"}
    met = metaldict[metallicity]
    datafile = "/home/kuuttila/3ML/Kurucz/ck" + met + "_" + str(int(temperature)) + ".fits"   
    if logg == 0.0:
        logg_index = "g00"
    elif logg == 0.5:
        logg_index = "g05"
    else:
        logg_index = "g" + str(int(logg*10))
        
    with fits.open(datafile) as hdul:
        alldata = hdul[1].data
        wavel = alldata['WAVELENGTH']
        modeldata = alldata[logg_index]
    fnus = []
    for fil in filters:
        fnus.append(get_fnus(wavel, modeldata, fil))
    return np.array(fnus)


temperature = 4000
metallicity = -1.0
logg = 1.0

alldata = np.loadtxt("KuruczConvolvedFilter.dat")
temps = alldata[:,0]
mets = alldata[:, 1]
loggs = alldata[:, 2]

iind = np.where((temps == temperature) & (mets == metallicity) & (loggs == logg))[0][0]
print(alldata[iind,:])


exit()

temp_grid = np.concatenate((np.arange(3500, 13000, 250), np.arange(13000, 51000, 1000)))
metal_grid = np.array([0.5, 0.2, 0.0, -0.5, -1.0, -1.5, -2.0, -2.5])
logg_grid = np.arange(0.0, 5.5, 0.5)


output = open("KuruczConvolvedFilter.dat", 'w')

for a in temp_grid:
    for b in metal_grid:
        for c in logg_grid:
            modelData = get_Kurucz_data(a, b, c)
            output.write(str(a) + " " + str(b) + " " + str(c) + " ")
            for d in modelData:
                output.write(str(d) + " ")
            output.write("\n")
            
output.close()


















