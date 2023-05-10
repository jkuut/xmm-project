import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from astropy.io import fits
from astropy.table import Table, Column

tempType = { 'O' : 1, 'B' : 2, 'A' : 3, 'F' : 4, 'G' : 5, 'K' : 6, 'M' : 7 }

with fits.open("fullTable_fitResults.fits") as hdul:
    alldata = hdul[1].data
    spectral_type = alldata['sp_type']
    fitTemp = alldata['Temperature']

fig = plt.figure()
plt.rcParams.update({'font.size': 20})
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['ytick.major.width'] = 2
lsize=20
plt.rcParams['xtick.labelsize'] = lsize
plt.rcParams['ytick.labelsize'] = lsize
plt.rcParams['xtick.minor.size'] = 8
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.minor.size'] = 8
plt.rcParams['ytick.minor.width'] = 1
ax1 = fig.add_subplot(111)
ax1.set_xlabel(r"VOSA best fit temperature (K)", size=lsize)
ax1.set_ylabel("Numerical spectral type", size=lsize)
#ax1.set_title("Nearby main sequence stars")



temperatures = []
fit_temps = []



for i in range(len(spectral_type)):
    sp = spectral_type[i]
    temp = fitTemp[i]
    
    if len(sp) > 0:

        if len(sp) == 1:
            if sp in tempType.keys():
                temperatures.append(tempType[sp]+0.5)
                fit_temps.append(temp)
        else:
            if sp[0] in tempType.keys():
                maintype = tempType[sp[0]]
                if sp[1] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                    maintype = maintype + int(sp[1])/10.0
                else:
                    maintype = maintype + 0.5
                temperatures.append(maintype)
                fit_temps.append(temp)

print(len(temperatures))


np.savetxt("aaaaa.txt", np.c_[fit_temps, temperatures])

#ax2 = ax1.twinx()
labels = ['O', 'B', 'A', 'F', 'G', 'K', 'M']
ax1.set_ylim(8.0, 1.0)
ax1.set_ylim(8.0, 1.0)
ax1.set_xlim(3500, 20000)
ax1.set_yticklabels(labels)
ax1.set_ylabel("Simbad spectral type")

ax1.plot(fit_temps, temperatures, 'ko')

fig.set_size_inches(16, 9)
fig.savefig('temperature.png', dpi=200)

plt.show()
    























