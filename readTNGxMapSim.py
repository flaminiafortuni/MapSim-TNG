#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 13:21:00 2019

@author: flaminia.fortuni@inaf.it
"""

import illustris_python as il
import numpy as np
import h5py as h5
import sys
from bisect import bisect_left
if 'threading' in sys.modules:
    del sys.modules['threading']


def take_closest(myList, myNumber):
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0
    if pos == len(myList):
        return len(myList)-1
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return pos
    else:
       return (pos-1)



def FreadTNG(i):
    
    h=0.6774
    
    #PATHS#
    basePath = '/YOUR_PATH/MapSim-TNG/output'
    AGEPath="/YOUR_PATH/MapSim-TNG/age_bc03.txt" #age bc03 - not filtered
    ZLISTPath="/YOUR_PATH/MapSim-TNG/illustris_python/redshift_list_TNG_pyth.txt" #TNG z,age,lbt list - not filtered
    OUTPath='/YOUR_PATH/MapSim-TNG/sn_%d.txt' #output snapshot
    HEADERPath='/YOUR_PATH/MapSim-TNG/header_%d.txt'
    f_1dig='/YOUR_PATH/MapSim-TNG/output/snapdir_00%d/snap_00%d.0.hdf5' #snapshot with 1 digit
    f_2dig='/YOUR_PATH/MapSim-TNG/output/snapdir_0%d/snap_0%d.0.hdf5' #snapshot with 2 digits    
    fields_s=['Masses','Coordinates', 'GFM_InitialMass', 'GFM_StellarFormationTime', 'GFM_Metallicity']#,'GFM_StellarPhotometrics']

################(a snaps=i)###########################################

    stars=il.snapshot.loadSubset(basePath, i, 'stars', fields=fields_s)

    #STAR PARTICLE   
    coords=stars['Coordinates'] #Spatial comoving coordinates of this star particle [ckpc/h]
    xs=coords[:,0].tolist()
    ys=coords[:,1].tolist() 
    zs=coords[:,2].tolist()
    GFM_InMass=stars['GFM_InitialMass']* 1.e10 #Mass of this star particle when it was formed [Msun/h]
    GFM_InTime=stars['GFM_StellarFormationTime'] #The exact time (given as the scale factor) when this star was formed
    #GFM_Phot=stars['GFM_StellarPhotometrics'] #Stellar magnitudes of this star particle in eight bands: U,B,V,K,g,r,i,z (rest frame) [mag] (available only for "big snapshots")
    StellarMass=stars['Masses']*1.e10 #Mass of this star particle [Msun/h]
    
    ZF=[] #Formation redshift of the star particle, computed as (1/GFM_InTime)-1

    GFM_sMet=stars['GFM_Metallicity'] #Ratio Mz/Mtot of this star particle (with Mz total mass of all metal elements above He). Not in solar metallicity. To convert in solar met, divide by 0.0127
    met_bc03=[0.0001,0.0004,0.004,0.008, 0.02,0.05] #bc03 available metallicities
    met_app=[] #Metallicity of this star particle - GFM_sMet approximated to the closest value in met_bc03
    
    age_bc03=(np.loadtxt(AGEPath,dtype="double",usecols=(0),unpack=True)*1.e-9).tolist() #bc03 available ages in [Gyr]
    age=[] #Age of this star particle computed from DT
    age_app=[] #Age of this star particle - age approximated to the closest value in age_bc03

    
    #Sorting z_list from the end to the beginning, otherwise take_closest() doesn't work   
    rlist=np.loadtxt(ZLISTPath,dtype="double",usecols=(0,1,2,3),skiprows=1, unpack=True)
    snap_list=np.transpose(rlist[0,::-1]).tolist() #TNG snapshot list
    z_list=np.transpose(rlist[1,::-1]).tolist() #TNG redshift list
    ageuni_list=np.transpose(rlist[2,::-1]).tolist() #TNG age of the Uni corresponding to z_list [Gyr]
    lbt_list=np.transpose(rlist[3,::-1]).tolist() #TNG lookback time of the Uni corresponding to z_list [Gyr]

    
    #Computing the age of a star particle [Gyr] - in TNG it is given as a scale factor
    for l in range(len(GFM_InTime)):
            ZF.append(1./GFM_InTime[l]-1)
            
    for l in range(len(ZF)):
        if ZF[l]<0:
            age.append(1000)
        elif ZF[l]>0:
            idx1=take_closest(z_list,ZF[l])
            idx2=snap_list.index(i)
            age.append(ageuni_list[idx2]-ageuni_list[idx1]) #Age of ssp (calculated as DT between zf and zobs) 
    
    for l in range(len(GFM_sMet)):
        idx3=take_closest(met_bc03, GFM_sMet[l])
        met_app.append(met_bc03[idx3]) #Approximated metallicity closer to met_bc03
        
    for l in range(len(age)):
        idx4=take_closest(age_bc03, age[l])
        age_app.append(age_bc03[idx4]) #Approximated age closer to age_bc03
        
    
    STARS=np.row_stack((xs,ys,zs,GFM_InMass,StellarMass, met_app, age_app, ZF))
    

    ###Creation of STAR PARTICLES file (at snapshot i)###
    #values printed in the sn_*.txt file = (xs, ys, zs, GFM_InMass, StellarMass, met_app, age_app, ZF)#
    outfile_i=open(OUTPath %(i), 'w')
    for l in range(len(STARS[0])):
        	outfile_i.write('%.6f \t  %.6f \t %.6f \t %.6f \t %.6f \t %.6f \t %.6f \t %.6f \n' %(STARS[0][l],STARS[1][l],STARS[2][l],STARS[3][l],STARS[4][l],STARS[5][l], STARS[6][l], STARS[7][l]) )
    outfile_i.close()
    
    
            
    ###Creation of HEADER file (at snapshot i)###
    #values printed in the header_*.txt file = (box size, mass_gas, mass_dm, mass_tracer, mass_star, mass_bh, n_part_gas, n_part_dm, n_part_tracer, n_part_star, n_part_bh, Om0, OmL, time, redshift)# 
    if i<10:
        ff=h5.File(f_1dig %(i,i), 'r')
    elif i>=10:
        ff=h5.File(f_2dig %(i,i), 'r')
    dset=ff[u'Header']
    bs=dset.attrs['BoxSize']
    masst=dset.attrs['MassTable']*1.e10 #Masses of particle (of each type) which have a constant mass (only DM - other types are empty) [Msun/h]
    npart=dset.attrs['NumPart_Total'] #Total number of particles (of each type) in this snapshot
    Om0=dset.attrs['Omega0']
    OmL=dset.attrs['OmegaLambda']
    time=dset.attrs['Time'] #The scale factor a = 1/(1 + z) corresponding to the current snapshot
    rs=dset.attrs['Redshift']
    
    
    outfile_j=open(HEADERPath %(i), 'w')
    outfile_j.write('%f \t  %f \t %f \t %f \t %f \t %f \t %f \t %d \t %d \t %d \t %d \t %d \t %d \t %f \t %f \t %f \t %f \n' %(bs, masst[0], masst[1], masst[2], masst[3], masst[4], masst[5], npart[0], npart[1], npart[2], npart[3], npart[4], npart[5],Om0,OmL,time,rs))
    outfile_j.close()
        

