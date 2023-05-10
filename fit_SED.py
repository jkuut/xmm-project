import os
import time
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from astropy.io import fits
from scipy import integrate, interpolate
from scipy.optimize import curve_fit
from astropy.modeling import models
from astropy import units as u
from multiprocessing import Pool



# Fit extinction and blackbody components?
fit_Av = True
fit_BB = True
verbose = False     # Print parameters during fitting?
createPlots = True  # Create plots?

# Catalog file
datafile = "/home/kuuttila/Topcat/XMM_DR11/WISE_distance_simbad_GoodSources.fits"




# List of filter response file names
filters = ['gaia_BPmag', 'gaia_Gmag', 'gaia_RPmag',
                '2mass_Jmag', '2mass_Hmag', '2mass_Kmag', 
                'W1mag', 'W2mag', 'W3mag', 'W4mag']
# List of photometric band names and their errors, as appears in our catalog
filterList = ['gaia_BPmag', 'Gaia_e_BPmag', 'gaia_Gmag', 'Gaia_e_Gmag', 'gaia_RPmag', 'Gaia_e_RPmag', 
                '2mass_Jmag', '2mass_e_Jmag', '2mass_Hmag', '2mass_e_Hmag', '2mass_Kmag', '2mass_e_Kmag', 
                'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag', 'W3mag', 'e_W3mag', 'W4mag', 'e_W4mag']
                
# Central wavelengths for the different photometric filters
#   the order needs to be the same as in the list "filters"    
centralWavelengths = [5035.75, 5822.39, 7619.96, 
                        12350.00, 16620.00, 21590.00, 
                        33526.00, 46028.00, 115608.00, 220883.00]
                        
                        
# Castelli & Kurucz model parameter grids
temp_grid = np.concatenate((np.arange(3500, 13000, 250), np.arange(13000, 21000, 1000)))
metal_grid = np.array([0.5, 0.2, 0.0, -0.5, -1.0, -1.5, -2.0, -2.5])
logg_grid = np.arange(0.0, 5.5, 0.5)

# Read the extinction model
# The extinction at each wavelength is calculated as: A_λ = A_V * k_λ/k_V, 
#  where k_λ is the opacity for a given λ and k_V=211.4. A_V is the fit parameter.
extinctionWavel, extinction = np.loadtxt("responses/ExtinctionCurve.dat", unpack=True)
extinction = extinction/211.4 
f = interpolate.interp1d(extinctionWavel, extinction)
extinctionAtCentralWavel = f(centralWavelengths)                       
                       
                        
print(extinctionAtCentralWavel)
exit()
                 
                        
def magToFlux(datapoints):
    # Convert Vega class magnitudes to fluxes
    # F = F0 * 10^(-mag/2.5), for erg/cm2/s/A the F0 is as below
    zeroPoints = np.array([4.08e-9, 2.50e-9, 1.27e-9, 3.13e-10, 1.13e-10, 4.28e-11, 
                            8.18e-12, 2.42e-12, 6.52e-14, 5.09e-15])
    fluxes = zeroPoints * np.power(10, -0.4*np.array(datapoints))
    return fluxes


def get_fnus(wavel, flux, response):
    # Convolve the model spectrum with a filter to produce the flux in a certain photometric filter
    filter_wavel, filter_response = np.loadtxt("responses/" + response + ".dat", unpack=True)
    fint = interpolate.interp1d(wavel, flux)
    fluxInFilterRange = fint(filter_wavel)
    response2 = filter_response / integrate.simps(filter_response, filter_wavel)    
    fnu = integrate.simps(fluxInFilterRange*response2, filter_wavel)
    return fnu
    
    

def calculate_Kurucz_data(temperature, metallicity, logg):
    # Read the Castelli & Kurucz stellar model data from files
    # Datafile name and column number depend on the parameters 
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
    # Apply the photometric filters to the models
    for fil in filters:
        fnus.append(get_fnus(wavel, modeldata, fil))
    return np.array(fnus)
    

# Instead of calculating the models and filters each time, read pre-calculated numbers from a file
# The I/O is slow, it's much faster to read the data into memory and the function accesses only this memory, not the file
alldata = np.loadtxt("KuruczConvolvedFilter.dat")
temps = alldata[:, 0]
mets = alldata[:, 1]
loggs = alldata[:, 2]
def read_Kurucz_data(temperature, metallicity, logg):
    iind = np.where((temps == temperature) & (mets == metallicity) & (loggs == logg))[0][0]
    return alldata[iind, 3:]



                        
def blackbody(temp, scale):
    # Blackbody in the data range with the photometric filters applied
    bb = models.BlackBody(temperature=temp*u.K, scale = scale *u.erg / (u.cm**2 * u.s * u.AA * u.sr))
    wav = np.arange(3100, 280600) * u.AA    # BB wavelength limits, based on the Gaia G and Wise W4 filter ranges
    bbfluxes = bb(wav)*4.0*3.141592653*u.sr
    fnus = []                 
    for fil in filters:
        fnus.append(get_fnus(wav, bbfluxes.value, fil))
    return np.array(fnus)
    
    
def getSourceName(srcid):
    # Get the source name based on the numerical source id for plotting
    with fits.open(datafile) as hdul:
        fullcatalog = hdul[1].data
        sources = fullcatalog['xmm_IAUNAME'] 
        ids = fullcatalog['XMM_SRCID'] 
        sourceName = sources[np.where(ids==srcid)][0].replace(" ", "")
    return sourceName


def fit_Kurucz_models(sourceData):
    
    # This is the main fitting function. Takes as an input the sourceData, where first column is source id, 
    #  second column is the distance to the source, and the rest are magnitudes and associated errors in various filters. 
    #  The observations and errors are in the same order as in the pre-defined list "filterList"

    obsData = magToFlux(sourceData[2::2])
    obsErrors = obsData * sourceData[3::2]
    
    #Check for NaNs in errors and replace them with 100% errors
    for qwe in range(len(obsErrors)):
        if math.isnan(obsErrors[qwe]):
            obsErrors[qwe] = obsData[qwe]
        # If any observations are NaN then skip these sources and return a zero array
        if math.isnan(obsData[qwe]):
            return np.array([sourceData[0], sourceData[1], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
           
    #Second error array for "visual" goodness of fit, similar to VOSA
    # if original error less than 10%, set the error to 10% 
    obsErrors2 = []      
    for qw in range(len(obsErrors)):
        if obsErrors[qw] < 0.1 * obsData[qw]:
            obsErrors2.append(0.1 * obsData[qw])
        else:
            obsErrors2.append(obsErrors[qw])            
    obsErrors2 = np.array(obsErrors2)
    

    params = []
    K_array = []
    Av_array = []
    gofs = []
    Vgofs = []
    # Loop through the parameter grids, and for each parameter combination fit the K and Av 
    # Then calculate gof for each combination and at the end choose the best fittin model
    for a in temp_grid:
        for b in metal_grid:
            for c in logg_grid:
                #Read pre-calculated model data for set of parameters. It's faster to read the data rather than always calculate it
                modelData = read_Kurucz_data(a, b, c)
                
                # Some parameter combinations are not included in the models, skip these
                if all(modelData <= 0.0):
                    continue

                if fit_Av:  
                    #Fit K and Av for each parameter combination
                    fitfunc = lambda x, *p: p[0] * x * np.power(10, -0.4*p[1]*extinctionAtCentralWavel)
                    popt, pcov = curve_fit(fitfunc, modelData, obsData, p0=[1e-20, 0.0], bounds=((1e-30, 0.0),(1e-10, 3.0)), sigma=obsErrors)
                    Av_array.append(popt[1])
                    #Calculate actual goodness of fit, and "visual" gof with min 10% errors
                    gofs.append(np.sum((fitfunc(modelData, *popt) - obsData)**2/obsErrors**2)/5.0)
                    Vgofs.append(np.sum((fitfunc(modelData, *popt) - obsData)**2/obsErrors2**2)/5.0)     
                else:
                    fitfunc = lambda x, *p: p[0] * x 
                    popt, pcov = curve_fit(fitfunc, modelData, obsData, p0=[1e-20], sigma=obsErrors)
                    gofs.append(np.sum((fitfunc(modelData, *popt) - obsData)**2/obsErrors**2)/6.0)
                    Vgofs.append(np.sum((fitfunc(modelData, *popt) - obsData)**2/obsErrors2**2)/6.0)                    
                    
                params.append([a, b, c])
                K_array.append(popt[0])
                
    # The best fit model is the combination of parameters with smallest goodness of fit value     
    min_gof = min(gofs)
    minIndex = gofs.index(min_gof)    
    min_Vgof = Vgofs[minIndex]
    bestfit = params[minIndex]
    # best fitting K and Av for the best parameter combination
    bestK = K_array[minIndex]
    if fit_Av:  
        bestAv = Av_array[minIndex]
    else:
        bestAv = 0    
        
    if verbose:
        print("XMM source ID = ", str(getSourceName(sourceData[0])))
        print("Best K = ", str(bestK))
        print("Best Av = ", str(bestAv))
        print("Best fit = ", str(bestfit))
        print("Min gof = ", str(min_gof))
        print(" Vgof = " + str(Vgofs[minIndex]))
    
    # If fitting the BB, and the original fit is bad, then add a cold BB component and fit again
    if fit_BB:
        if min_gof > 3.0 and min_Vgof > 3.0: 
        
            # Get the best fitting Kurucz models (with or without extinction) and add a BB component to that       
            if fit_Av:
                bestModel = bestK * read_Kurucz_data(*bestfit) * np.power(10, -0.4*bestAv*extinctionAtCentralWavel)
            else:
                bestModel = bestK * read_Kurucz_data(*bestfit)

            fitfunc_BB = lambda x, *bbp: x + blackbody(*bbp)
            BB_popt, pcovBB = curve_fit(fitfunc_BB, bestModel, obsData, p0=[100, bestK], bounds=((30, bestK*1.e-6),(1000, bestK*1.e6)), sigma=obsErrors)
            BB_gof = (np.sum((fitfunc_BB(bestModel, *BB_popt) - obsData)**2/obsErrors**2)/5.0)
            BB_Vgof = (np.sum((fitfunc_BB(bestModel, *BB_popt) - obsData)**2/obsErrors2**2)/4.0)
            BB_temp, BB_scale = BB_popt[0], BB_popt[1]
            
            # Did the BB improve the fit? 1 for yes, 0 for no 
            if BB_gof < min_gof and BB_Vgof < min_Vgof: 
                BB_fitImprov = 1
            else:
                BB_fitImprov = 0
        else:
            # If not doing BB fit, set the BB values to 0 anyways
            BB_temp, BB_scale, BB_gof, BB_Vgof = 0, 0, 0, 0
            BB_fitImprov = 0
    else:
        # If not doing BB fit, set the BB values to 0 anyways
        BB_temp, BB_scale, BB_gof, BB_Vgof = 0, 0, 0, 0
        BB_fitImprov = 0
    
    if verbose:
        print("BB temp, scale, gof, vgof = ")
        print(BB_temp, BB_scale, BB_gof, BB_Vgof)
   
    
    if createPlots:
        # Plot the data and models and save in a file
        sourceName = getSourceName(sourceData[0])
        fig = plt.figure()
        ax1 = fig.add_subplot(111, xscale='log', yscale='log')
        if BB_fitImprov:
            ax1.set_title(sourceName + ": T, M, Logg=" + str(bestfit) + f"\nK = {bestK:.2e}" + ", Av=" + str(round(bestAv,2)) + ", " + \
                           "BB temp = " + str(int(BB_temp)) + f", scale = {BB_scale:.2e}")
        else:
            ax1.set_title(sourceName + ": T, M, Logg=" + str(bestfit) + f"\nK = {bestK:.2e}" + ", Av=" + str(round(bestAv,2)) )
            
        ax1.errorbar(centralWavelengths, obsData*centralWavelengths, yerr=obsErrors*centralWavelengths, fmt='bo')
        legend_elements = [Line2D([0], [0], color='w', marker='o', markerfacecolor='b', label='Observations'),
                       Line2D([0], [0], marker='*', color='r', label='Kurucz model')]
        if fit_Av:
            plotdata = bestK * read_Kurucz_data(*bestfit) * np.power(10, -0.4*bestAv*extinctionAtCentralWavel) * centralWavelengths
            ax1.plot(centralWavelengths, plotdata, 'r*')
            ax1.plot(centralWavelengths, plotdata, 'r--')
        else:
            plotdata = bestK * read_Kurucz_data(*bestfit) * centralWavelengths
            ax1.plot(centralWavelengths, plotdata, 'r*')
            ax1.plot(centralWavelengths, plotdata, 'r--')        
        if fit_BB and BB_fitImprov:
            BB_plotdata = blackbody(*BB_popt) * centralWavelengths
            ax1.plot(centralWavelengths, BB_plotdata, 'ms')
            ax1.plot(centralWavelengths, BB_plotdata, 'm--')      
            legend_elements.append(Line2D([0], [0], color='m', marker='s', label='Blackbody T=' +str(int(BB_temp))))
            
            ax1.plot(centralWavelengths, plotdata + BB_plotdata, 'k-.')
            legend_elements.append(Line2D([0], [0], color='k', linestyle='-.', label='Kurucz + BB'))
            
        ax1.legend(handles=legend_elements)  
        ax1.set_ylim(min(obsData*centralWavelengths)*0.1, max(obsData*centralWavelengths)*10.0)
        ax1.set_xlabel(r'Wavelength ($\AA$)')
        ax1.set_ylabel(r"$\nu$F$_{\nu}$ (erg s$^{-1}$ cms$^{-2}$)")  
          
        fig.savefig("figs/" + sourceName + ".png", dpi=200)
        plt.close()

    # This is the output array with all the fit parameters
    outputAr = np.array([sourceData[0], sourceData[1], bestfit[0], bestfit[1], bestfit[2], \
        bestK, bestAv, min_gof, min_Vgof, BB_temp, BB_scale, BB_gof, BB_Vgof, BB_fitImprov])

    return outputAr
    
    
    
    
    
if __name__ == '__main__':    

    start = time.time()

    # Read the data from the catalog and select the photometric data for fitting
    with fits.open(datafile) as hdul:
        fullcatalog = hdul[1].data
        dataTable = np.zeros((len(fullcatalog), len(filterList)+2))
        sources = fullcatalog['xmm_IAUNAME'] 
        dataTable[:, 0] = fullcatalog['XMM_SRCID']  # XMM ID for source identification
        dataTable[:, 1] = fullcatalog['rpgeo']      # Distance for converting normalization K to radius
        for i in range(len(filterList)):            # Get the photometric data for pre-defined filters
            dataTable[:, i+2] = fullcatalog[filterList[i]]
            
            
    # Create figs folder if it doesn't exist
    if not os.path.isdir('figs'):
        os.mkdir('figs')

    # Initialize output file
    with open("fit_output.txt", 'w') as outf:
        outf.write("# Source, Distance, Temperature, metallicity, log g, K, Av, gof, visual gof, BB temp, BB scale, BB gof, BB visual gof, Fit improved by BB?\n")
        
        # Simple multi-processing method for faster calculation
        pool = Pool(4)
        
        # Run the fitting function parallel in arbitrary order for all the sources in the catalog
        for outputArray in pool.imap(fit_Kurucz_models, (dataTable[ij, :] for ij in range(len(sources)))):
            # Save the results in a file
            sourceName = getSourceName(outputArray[0])
            for jk in range(len(outputArray)):
                if jk == 0:
                    outf.write(sourceName + " ")
                else:
                    outf.write(str(outputArray[jk]) + " ")
            outf.write("\n")
            
            #np.savetxt(outf, np.c_[sourceName, outputArray[1:]], fmt="%s")
            #print(outputArray)
    
        end = time.time()
        print("Total run time = ", str(end - start))
    


    #Convert output txt file to fits file
    
    names = ["xmm_IAUNAME", "rpgeo", "Temperature", "Metallicity", "Log g", "K Normalisation", "Extinction Av", "Goodness of fit", "Visual goodness of fit", 
        "Blackbody temperature", "BB normalisation", "BB goodness of fit", "BB Visual goodness of fit", "BB improves fit?"]
    formats = ['A21', 'F', 'F', 'F', 'F', 'E', 'F', 'F', 'F', 'F', 'E', 'F', 'F', 'A3']

    outputdataArray = np.loadtxt("fit_output.txt", dtype='str')

    assert len(names) == len(formats)
    assert len(names) == len(outputdataArray[0,:])

    allColumns = []
    for i in range(len(names)):
        if i == 0:
            arr = np.array([aa.replace("J", " J") for aa in outputdataArray[:, i]])
        elif i == len(names)-1:
            arr1 = np.array([aa.replace("1.0", "True") for aa in outputdataArray[:, i]])
            arr = np.array([a.replace("0.0", "False") for a in arr1])
        else:
            arr = np.array(outputdataArray[:, i])

        col = fits.Column(name=names[i], array=arr, format=formats[i])
        allColumns.append(col)

    t = fits.BinTableHDU.from_columns(allColumns)
    t.writeto('fit_results.fits', overwrite=True)












