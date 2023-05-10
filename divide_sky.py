import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import random

#### OUTPUT FILE
# Create output file
outputFile = "obsidlist_for_xmatch.txt"
with open(outputFile, 'w') as outf:
    outf.write("#ID, RA, DEC, radius in radians\n")


def distance_formula(lon1, lat1, lon2, lat2, deg_or_rad='rad'):
    # Angular distance based on haversine formula
    if deg_or_rad == 'deg':
        lon1 = np.radians(lon1)
        lat1 = np.radians(lat1)
        lon2 = np.radians(lon2)
        lat2 = np.radians(lat2)
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))     
    return c


def quicker_distance(lon1, lat1, lon2, lat2, deg_or_rad='rad'):
    # Angular distance with a quicker approximate formula
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

def combine_nearby_old(RAs, DECs, radius=5):
    rad = radius/60.0/180.0 * np.pi 
    new_RAs = []
    new_DECs = []
    for i in range(len(RAs)):
        temp_RAs = []
        temp_DECs = []
        distances = np.sqrt( np.power((RAs - RAs[i]), 2) + np.power((DECs - DECs[i]), 2))
        for j in range(len(distances)):
            if distances[j] < rad:
                temp_RAs.append(RAs[j])
                temp_DECs.append(DECs[j])  
        new_RAs.append(np.round(np.mean(temp_RAs), 6))
        new_DECs.append(np.round(np.mean(temp_DECs), 6))   
    nnew_RAs = []
    nnew_DECs = []
    return np.array(new_RAs), np.array(new_DECs)
    

def combine_nearby(RAs, DECs, radius=5):
    rad = radius/60.0/180.0*np.pi
    new_RAs = []
    new_DECs = []
    used_RAs = np.array([])
    used_DECs = np.array([])
    for a, b in zip(RAs, DECs):
        if a in used_RAs and b in used_DECs:
            continue      
        distances = quicker_distance(a, b, RAs, DECs)
        close_RAs = RAs[distances <= rad]
        close_DECs = DECs[distances <= rad] 
        new_RAs.append(np.round(np.mean(close_RAs), 6))
        new_DECs.append(np.round(np.mean(close_DECs), 6))
        used_RAs = np.concatenate((used_RAs, close_RAs))
        used_DECs = np.concatenate((used_DECs, close_DECs))
    return np.array(new_RAs), np.array(new_DECs)


def separate_singles(RAs, DECs, radius=30):
    # separate single pointings from crowded areas and save singles into the outputfile
    rad = radius/60.0/180.0 * np.pi 
    single_RA = []
    single_DEC = []
    other_RA = []
    other_DEC = []
    for a, b in zip(RAs, DECs):
        flag = 1
        for x, y in zip(RAs, DECs):
            if np.abs(a-x) < 0.0174533 and np.abs(b-y) < 0.0174533 and a != x and b != y:
                dist = quicker_distance(a, b, x, y)
                if dist < rad:
                    other_RA.append(a)
                    other_DEC.append(b)
                    flag = 0
                    break
        if flag:
            single_RA.append(a)
            single_DEC.append(b)    
    radii = np.zeros(len(single_RA)) + rad
    with open(outputFile, 'a') as outf:
        np.savetxt(outf, np.c_[np.degrees(single_RA), np.degrees(single_DEC), np.degrees(radii)], fmt='%.6f')
    
    return np.array(other_RA), np.array(other_DEC)



def separate_small_groups(RAs, DECs, radius=55):
    rad = radius/60.0/180.0*np.pi
    largeGroup_RAs = []
    largeGroup_DECs = []
    false_RAs = []
    false_DECs = []
    mean_RAs = []
    mean_DECs = []
    used_RAs = np.array([])
    used_DECs = np.array([])
    for a, b in zip(RAs, DECs):
        if a in used_RAs and b in used_DECs:
            continue
        distances = quicker_distance(a, b, RAs, DECs)
        close_RAs = RAs[distances < rad]
        close_DECs = DECs[distances < rad]    
        flag = 1
        for c, d in zip(close_RAs, close_DECs):
            newdist = quicker_distance(c, d, RAs, DECs)
            candRAs = RAs[newdist < rad]
            candDECs = DECs[newdist < rad]
            for e, f in zip(candRAs, candDECs):
                if e not in close_RAs and f not in close_DECs:
                    flag = 0
                    break                      
        if flag:
            mean_RA = np.mean(close_RAs)
            mean_DEC = np.mean(close_DECs)
            mean_RAs.append(mean_RA)
            mean_DECs.append(mean_DEC)
            used_RAs = np.concatenate((used_RAs, close_RAs))
            used_DECs = np.concatenate((used_DECs, close_DECs)) 
            for g, h in zip(largeGroup_RAs, largeGroup_DECs):
                lrdist = quicker_distance(g, h, mean_RA, mean_DEC)
                if lrdist < rad:
                    false_RAs.append(g)
                    false_DECs.append(h)
        else:
            largeGroup_RAs.append(a)
            largeGroup_DECs.append(b)
    radii = np.zeros(len(mean_RAs)) + rad
    with open(outputFile, 'a') as outf:
        np.savetxt(outf, np.c_[np.degrees(mean_RAs), np.degrees(mean_DECs), np.degrees(radii)], fmt='%.6f')       
    returnRAs = []
    returnDECs = []     
    for n, m in zip(largeGroup_RAs, largeGroup_DECs):
        if n not in false_RAs and m not in false_DECs:
            returnRAs.append(n)
            returnDECs.append(m)
    return np.array(returnRAs), np.array(returnDECs)



def splice_large_groups(RAs, DECs, radius=55):
    rad = radius/60.0/180.0*np.pi
    mean_RAs = np.array([])
    mean_DECs = np.array([])
    used_RAs = np.array([])
    used_DECs = np.array([])
    for a, b in zip(RAs, DECs):
        if a in used_RAs and b in used_DECs:
            continue
        distances = quicker_distance(a, b, RAs, DECs)
        close_RAs = RAs[distances < rad]
        close_DECs = DECs[distances < rad]               
        groupRAs = list(close_RAs)
        groupDECs = list(close_DECs)
        flag = True
        while flag:
            flag = False
            for c, d in zip(groupRAs, groupDECs):
                dists = quicker_distance(c, d, RAs, DECs)
                close_RAs = RAs[distances < rad]
                close_DECs = DECs[distances < rad]                   
                for e, f in zip(close_RAs, close_DECs):
                    if e not in groupRAs and f not in groupDECs:
                        groupRAs.append(e)
                        groupDECs.append(f)
                        flag = True
        minLength = 1e5
        bestMeanRAs = []
        bestMeanDECs = []
        for rndi in range(len(groupRAs)*1000):
            qw = list(zip(groupRAs, groupDECs))
            random.shuffle(qw)
            groupRAs, groupDECs = zip(*qw)
            groupRAs = np.array(groupRAs)
            groupDECs = np.array(groupDECs)
            used_groupRAs = np.array([])
            used_groupDECs = np.array([])
            mean_groupRAs = []
            mean_groupDECs = []
            for g, h in zip(groupRAs, groupDECs):
                if g not in used_groupRAs and h not in used_groupDECs:
                    groupdists = quicker_distance(g, h, groupRAs, groupDECs)
                    closeRAs = groupRAs[groupdists < rad]
                    closeDECs = groupDECs[groupdists < rad] 
                    meanRAs = []
                    meanDECs = []
                    for n, m in zip(closeRAs, closeDECs):
                        if n not in used_groupRAs and m not in used_groupDECs:
                            meanRAs.append(n)
                            meanDECs.append(m)
                    mean_groupRAs.append(np.mean(meanRAs))
                    mean_groupDECs.append(np.mean(meanDECs))
                    used_groupRAs = np.concatenate((used_groupRAs, closeRAs))
                    used_groupDECs = np.concatenate((used_groupDECs, closeDECs))
            if len(mean_groupRAs) < minLength:
                minLength = len(mean_groupRAs)
                bestMeanRAs = mean_groupRAs
                bestMeanDECs = mean_groupDECs
        used_RAs = np.concatenate((used_RAs, np.array(groupRAs)))
        used_DECs = np.concatenate((used_DECs, np.array(groupDECs)))
        mean_RAs = np.concatenate((mean_RAs, np.array(bestMeanRAs)))
        mean_DECs = np.concatenate((mean_DECs, np.array(bestMeanDECs)))
    radii = np.zeros(len(mean_RAs)) + rad
    with open(outputFile, 'a') as outf:
        np.savetxt(outf, np.c_[np.degrees(mean_RAs), np.degrees(mean_DECs), np.degrees(radii)], fmt='%.6f')
    aa = zip(list(used_RAs), list(used_DECs))
    res = sorted(list(set(aa)))
    usRAs, usDECs = list(map(list, zip(*res)))
    return len(RAs) - len(usRAs)
        
        
# Read the list of XMM OBSIDs
# This list includes all unique OBSIDs with REFCAT > 0.0 (i.e. good observations)
OBSIDs, allRAs, allDECs = np.loadtxt('/home/kuuttila/Topcat/XMM_DR11/Catalogs/clean_obslist.txt', unpack=True)
allRAs = np.radians(allRAs)
allDECs = np.radians(allDECs)

# First separate single pointings from the grouped ones
# Single ones meaning distance to other pointings > 30 arcmin
# These are saved into file to be xmatched later and grouped ones processed further
RAs, DECs = separate_singles(allRAs, allDECs)  

# Then combine nearby ( < 5 arcmin) observations into one mean location 
RAs, DECs = combine_nearby(RAs, DECs)
RAs, DECs = combine_nearby(RAs, DECs)

# Repeat single separation
RAs, DECs = separate_singles(RAs, DECs)  

# Next separate the small groups that can be combined into one single fields (within 1 degree)
# and leave the large groups for last
RAs, DECs = separate_small_groups(RAs, DECs)  


#np.savetxt("OTHERS.txt", np.c_[np.degrees(RAs), np.degrees(DECs)], fmt='%.6f') 
#exit()
#RAs, DECs = np.loadtxt('OTHERS.txt', unpack=True)
#RAs = np.radians(RAs)
#DECs = np.radians(DECs)

# Finally splice the remaining large groups in 1 degree circles
exitValue = splice_large_groups(RAs, DECs)
if exitValue:
    print("Spliced unsuccessfully: ", exitValue)


### Check radii between large and small pointings
allRAs, allDECs, allRadii = np.loadtxt(outputFile, unpack=True, dtype='float')
allRAs = np.radians(allRAs)
allDECs = np.radians(allDECs)
largerRadii = allRadii[allRadii > 0.6]
largerRAs = allRAs[allRadii > 0.6]
largerDECs = allDECs[allRadii > 0.6]
smallerRadii = allRadii[allRadii < 0.6]
smallerRAs = allRAs[allRadii < 0.6]
smallerDECs = allDECs[allRadii < 0.6]
allRadii = np.radians(allRadii)
largerRadii = np.radians(largerRadii)
smallerRadii = np.radians(smallerRadii)
newRadii = list(largerRadii)
newRAs = list(largerRAs)
newDECs = list(largerDECs)
maxdist = 45.0/60.0/180.0*np.pi
for a, b, c in zip(smallerRAs, smallerDECs, smallerRadii):
    distances = quicker_distance(a, b, largerRAs, largerDECs)
    if all(distances > maxdist):
        newRAs.append(a)
        newDECs.append(b)
        newRadii.append(c)      
newRAs = np.array(newRAs)
newDECs = np.array(newDECs)
newRadii = np.array(newRadii)
#with open(outputFile, 'w') as outf:
#    outf.write("#ID, RA, DEC, radius in radians\n")
#np.savetxt(outputFile, np.c_[newRAs, newDECs, newRadii], fmt='%.6f')



### Check data and pointing distances

datafile = '/home/kuuttila/Topcat/XMM_DR11/Catalogs/clean_catalog_slim.fits'
with fits.open(datafile) as hdul:
    alldata = hdul[1].data
    slim_RAs = np.radians(alldata['SC_RA'])
    slim_DECs = np.radians(alldata['SC_DEC'])

for a, b in zip(slim_RAs, slim_DECs):
    distances = quicker_distance(a, b, newRAs, newDECs)
    dist2 = sorted(set(distances))
    minRad = newRadii[distances == dist2[0]]
    minRad2 = newRadii[distances == dist2[1]]
    minRad3 = newRadii[distances == dist2[2]]   
    if dist2[0] > minRad and dist2[1] > minRad2 and dist2[2] > minRad3:
        minRA = newRAs[distances == dist2[0]]
        minDEC = newDECs[distances == dist2[0]]
        minindex = np.where((newRAs == minRA) & (newDECs == minDEC))
        newRadii[minindex] = newRadii[minindex] + 15.0/60.0/180.0*np.pi
            
with open(outputFile, 'w') as outf:
    outf.write("#RA, DEC, radius in radians\n")          
np.savetxt("obsidlist_for_xmatch.txt", np.c_[np.degrees(newRAs), np.degrees(newDECs), np.degrees(newRadii)], fmt='%.6f')


