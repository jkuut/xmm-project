import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def distance_formula(lon1, lat1, lon2, lat2, deg_or_rad='rad'):
    # Angular distance with a approximate formula
    # good enough here for small distances
    if deg_or_rad == 'deg':
        lon1 = np.radians(lon1)
        lat1 = np.radians(lat1)
        lon2 = np.radians(lon2)
        lat2 = np.radians(lat2)
    x = (lon2 - lon1) * np.cos( 0.5*(lat2+lat1) )
    y = lat2 - lat1
    d = np.sqrt( x*x + y*y )
    return d

datafile = '/home/kuuttila/Topcat/XMM_DR11/Catalogs/clean_catalog_slim.fits'


#Splice the catalog based on previously created coordinate file
RAs, DECs, radius = np.loadtxt('obsidlist_for_xmatch.txt', unpack=True)
RAs = np.radians(RAs)
DECs = np.radians(DECs)
radius = np.radians(radius)

with open("fitsfiles.txt", 'w') as fasfs:
    fasfs.write("# OBSID RA DEC RADIUS\n")

with fits.open(datafile) as hdul:

    alldata = hdul[1].data
    slim_RAs = np.radians(alldata['SC_RA'])
    slim_DECs = np.radians(alldata['SC_DEC'])

    used_RAs = []
    used_DECs = []
    total_size = 0

    for i in range(len(RAs)):
        
        temp_RAs = slim_RAs[np.where(( np.abs( slim_RAs - RAs[i] ) < 0.035) & ( np.abs( slim_DECs - DECs[i] ) < 0.035))]
        temp_DECs = slim_DECs[np.where(( np.abs( slim_RAs - RAs[i] ) < 0.035) & ( np.abs( slim_DECs - DECs[i] ) < 0.035))]
    
        distances = distance_formula(RAs[i], DECs[i], temp_RAs, temp_DECs)

        save_RAs = temp_RAs[np.where((distances <= radius[i]))]
        save_DECs = temp_DECs[np.where((distances <= radius[i]))]

        mask = np.in1d(slim_RAs, save_RAs)

        mask2 = np.in1d(slim_RAs, used_RAs, invert=True)
        mask3 = np.in1d(slim_DECs, used_DECs, invert=True)

        savedata = alldata[mask & mask2 & mask3]


        if len(savedata) > 0:
            hdu = fits.BinTableHDU(data=savedata)
            hdu.writeto(str(i) + '.fits', overwrite=True)
            total_size = total_size + len(savedata)

            with open("fitsfiles.txt", 'a') as fasfs:
                fasfs.write(str(i) + " " + str(np.degrees(RAs[i])) + " " + str(np.degrees(DECs[i])) + " " + str(np.degrees(radius[i])) + "\n")


        for j in range(len(save_RAs)):
            used_RAs.append(save_RAs[j])
            used_DECs.append(save_DECs[j])


    #np.savetxt(str(len(alldata)) + "_" + str(total_size) + ".txt", np.c_[used_RAs, used_DECs])



