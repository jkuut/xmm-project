import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from astropy.io import fits
from astropy.table import Table, Column

MASS2_filters = {'2MASS/2MASS.H' : ['2mass_Hmag', '2mass_e_Hmag'], 
                '2MASS/2MASS.J' : ['2mass_Jmag', '2mass_e_Jmag'], 
                '2MASS/2MASS.Ks' : ['2mass_Kmag', '2mass_e_Kmag']}
GAIA_filters =  {'GAIA/GAIA3.G' : ['gaia_Gmag', 'Gaia_e_Gmag'], 
                'GAIA/GAIA3.Gbp' : ['gaia_BPmag', 'Gaia_e_BPmag'], 
                'GAIA/GAIA3.Grp' : ['gaia_RPmag', 'Gaia_e_RPmag']}
WISE_filters =  {'WISE/WISE.W1' : ['W1mag', 'e_W1mag'], 
                'WISE/WISE.W2' : ['W2mag', 'e_W2mag'], 
                'WISE/WISE.W3' : ['W3mag', 'e_W3mag'], 
                'WISE/WISE.W4' : ['W4mag', 'e_W4mag']}


def create_data_line(alldata, index, dfilter):

    if dfilter.startswith('GAIA'):
        mag = GAIA_filters[dfilter][0]
        err = GAIA_filters[dfilter][1]
    elif dfilter.startswith('2MASS'):
        mag = MASS2_filters[dfilter][0]
        err = MASS2_filters[dfilter][1]
    elif dfilter.startswith('WISE'):
        mag = WISE_filters[dfilter][0]
        err = WISE_filters[dfilter][1]
    else:
        print("Wrong filter")
        exit()

    objName = alldata['main_id'][index].replace('*', '_AST_').replace(' ', '_')
    outputString = objName + " " + str(alldata['gaia__RAJ2000'][index])  + " " + str(alldata['gaia__DEJ2000'][index])

    if alldata['rpgeo'][index] > 0:
        outputString = outputString + " " + str(alldata['rpgeo'][index])
    else: 
        outputString = outputString + " ---"

    outputString = outputString + " ---" # extinction set to zero for now

    outputString = outputString + " " + dfilter

    if alldata[mag][index] > 0:
        outputString = outputString + " " + str(alldata[mag][index])
    else: 
        outputString = outputString + " ---"

    if alldata[err][index] > 0:
        outputString = outputString + " " + str(alldata[err][index])
    else: 
        outputString = outputString + " ---"

    outputString = outputString + " mag Av:0.0/5.0"
    return outputString



#sources = []
#fbad = open("./VOSA_sources.txt", 'r')
#for line in fbad:
#    sources.append(str(line.rstrip('\n')))
#fbad.close()

sources = []
with fits.open("../BEST_WISE_SPTYPE.fits") as hdul2:
    alldata2 = hdul2[1].data
    main_ids = alldata2['Numerical_Sp_type']
    allsources = alldata2['main_id']
    for aa in range(len(main_ids)):
        if len(main_ids[aa]) > 0:
            sources.append(allsources[aa])

print(len(sources))

sources = ["2MASX J00032769+0207014"]

with fits.open("../WISE_distance_simbad.fits") as hdul:
    alldata = hdul[1].data
    main_ids = alldata['main_id']
    mask_id = np.in1d(main_ids, sources, invert=False)
    newdata = alldata[mask_id]

newsources = newdata['main_id']
print(len(newsources))

filters = list(GAIA_filters.keys()) + list(MASS2_filters.keys()) + list(WISE_filters.keys())

output = open("test.txt", 'w') 
for i in range(len(newsources)):
    for j in range(len(filters)):
        output.write(create_data_line(newdata, i, filters[j]) + "\n")
output.close()





























