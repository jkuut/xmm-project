import os
import numpy as np
import time
import shutil
from astropy.io import fits
from astropy.table import Table, vstack
import subprocess


XMATCH_template = """
get FileLoader file={OBSID}.fits
set pos ra=SC_RA dec=SC_DEC
set poserr type=CIRCLE param1=SC_POSERR/sqrt(2)
set cols *
prefix xmm_

get VizieRLoader tabname=I/355/gaiadr3 mode=cone center="{RA} {DEC}" radius={RAD}arcmin allcolumns
set pos ra=RAJ2000 dec=DEJ2000
set poserr type=COR_ELLIPSE param1=e_RAJ2000/1000.0 param2=e_DEJ2000/1000.0 param3=RADEcorJ2000
set cols *
prefix gaia_

get VizieRLoader tabname=II/246/out mode=cone center="{RA} {DEC}" radius={RAD}arcmin allcolumns
set pos ra=RAJ2000 dec=DEJ2000
set poserr type=ELLIPSE param1=errMaj param2=errMin param3=errPA
set cols *
prefix 2mass_

xmatch probaN_v1 joins=II completeness=0.9973 area={AREA}
merge dist mec
save {OBSID}_xmm_gaia_2mass.fits fits
"""



def create_XMATCH_scripts(OBSIDs, RAs, DECs, radius):
    RADII = radius*60.0 # arc minutes
    rAREA = 3.141592 * RADII**2 * 8.46138e-8 # radians^2
    files = []
    for i in range(len(OBSIDs)):
        OBSID = int(OBSIDs[i])
        RA = RAs[i]
        DEC = DECs[i]
        RAD = RADII[i]
        AREA = "{:.8f}".format(rAREA[i]) #correct format for xmatch script
        XMATCHscript = {
            "OBSID" : OBSID,
            "RA"    : RA,
            "DEC"   : DEC,
            "RAD"   : RAD,
            "AREA"  : AREA
        }
        with  open("XMATCH_" + str(OBSID) + ".txt", 'w') as myfile:
            myfile.write(XMATCH_template.format(**XMATCHscript))
        files.append("XMATCH_" + str(OBSID) + ".txt")
    return files
    
    

def run_XMATCH(OBSID, login=False):
    if login:
        subprocess.call('printf "anonymous\n anonymous" | ./arches.bash i', shell=True)
    print('./arches.bash put {}.fits'.format(OBSID))
    subprocess.call('./arches.bash put {}.fits'.format(OBSID), shell=True)
    print('./arches.bash x XMATCH_{}.txt'.format(OBSID))
    subprocess.call('./arches.bash x XMATCH_{}.txt'.format(OBSID), shell=True)   
    print('./arches.bash g {}_xmm_gaia_2mass.fits'.format(OBSID))
    subprocess.call('./arches.bash g {}_xmm_gaia_2mass.fits'.format(OBSID), shell=True)
    print('./arches.bash r {}_xmm_gaia_2mass.fits'.format(OBSID))
    subprocess.call('./arches.bash r {}_xmm_gaia_2mass.fits'.format(OBSID), shell=True)
    print('./arches.bash r {}.fits'.format(OBSID))
    subprocess.call('./arches.bash r {}.fits'.format(OBSID), shell=True)
    #shutil.move("{}.log".format(OBSID), "logs/{}.log".format(OBSID))
    return



if not os.path.isdir("logs"):
    os.mkdir("logs")
if not os.path.isdir("XmatchFiles"):
    os.mkdir("XmatchFiles")
if not os.path.isdir("XmatchScripts"):
    os.mkdir("XmatchScripts")
if not os.path.isdir("skyFitsFiles"):
    os.mkdir("skyFitsFiles")
    


OBSIDs, RAs, DECs, radius = np.loadtxt("fitsfiles.txt", unpack=True, dtype='float')

#for iij in range(len(radius)):
#    if radius[iij] > 0.49 and radius[iij] < 0.51:
#        radius[iij] = 0.37
    
files = create_XMATCH_scripts(OBSIDs, RAs, DECs, radius)

starttime = time.time()
failedFiles = []
for i in range(len(OBSIDs)):
    obsid = str(int(OBSIDs[i]))
    if os.path.exists("XMATCH_" + obsid + ".txt") and os.path.exists(obsid + ".fits"):
        if i < 0:
            run_XMATCH(obsid, False)
            cwd = os.path.dirname(os.path.realpath(__file__))
            shutil.copy2(cwd + '/' + obsid + '_xmm_gaia_2mass.fits', cwd + '/joined_results.fits')
            shutil.move(obsid + '_xmm_gaia_2mass.fits', 'XmatchFiles/' + str(obsid) + '_xmm_gaia_2mass.fits')
            shutil.move("XMATCH_{}.txt".format(obsid), "XmatchScripts/XMATCH_{}.txt".format(obsid))
            shutil.move("{}.fits".format(obsid), "skyFitsFiles/{}.fits".format(obsid))
        else:
            run_XMATCH(obsid, False)
            try:
                newTable = Table.read(obsid + '_xmm_gaia_2mass.fits', format='fits')
                fullTable = Table.read('joined_results.fits', format='fits')
                joinedTable = vstack([fullTable, newTable])
                joinedTable.write('joined_results.fits', format='fits', overwrite=True)
                shutil.move(obsid + '_xmm_gaia_2mass.fits', 'XmatchFiles/' + str(obsid) + '_xmm_gaia_2mass.fits')
                shutil.move("XMATCH_{}.txt".format(obsid), "XmatchScripts/XMATCH_{}.txt".format(obsid))
                shutil.move("{}.fits".format(obsid), "skyFitsFiles/{}.fits".format(obsid))
            except:
                failedFiles.append(obsid)


endtime = time.time()
print(endtime - starttime)
#print("FAILED: " + failedFiles)










