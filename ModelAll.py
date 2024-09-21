''' 
**** ADJUSTED FOR SUBMM DATA **** 


Note: 
ms values are not source unlensed AB mags as per original code, but source unlensed log S from Cai data;
msrc values not source lensed AB mags as per original code, but source lensed flux.
 - CFW  10/3/2019

Upgraded to Python 3.8
  - CJS Aug 2021
  
Submm, survey loops removed (not used - this is to simplify code) 
  - CJS 16 Aug 2021
'''


################################
# imports

import os
import sys
import time
import shutil
import pickle
import matplotlib.pyplot as plt
from astropy.io import fits
from __init__ import *


#################################
# set parameters


if len(sys.argv)==3:                # input for criterion, frac added CJS 20/9/2024
    criterion = float(sys.argv[1]) 
    frac = float(sys.argv[2])   

else:
    criterion = input('Enter criterion e.g. 10 or 80: ')
    criterion = float(criterion)
    frac = input('Enter frac (fraction of sky) e.g. 0.25: ')
    frac = float(frac)
    print()

survey = "submm"  # anachronism - no survey needed for submm (but used in file structure)
sourcepop = "submm"

sigfloor=200
a=20  # SN threshold - irrelevant for submm
b=3   # Magnification threshold
c=1000
d=1000
firstod=1
nsources=1
t0=time.perf_counter()
chunk=0
Si=0
SSPL={}
foundcount={}
foundcount[survey]=0    

nall=3700000       # cw check this (should be >idealisedlenses)      
nall=int(nall*frac)  
print('fraction of whole sky: ',frac)
print('criterion >:           ',criterion, ' mJy')

# re images
n_image = 0   # counter for no. of images saved (lines ~235+ below)
numim = 20    # sets limit on number of images (all images will be a lot of data!)

image_folder="images_submm"
if os.path.exists(image_folder):
    shutil.rmtree(image_folder)
os.mkdir(image_folder)
    

#################################
# create L and S[survey]

print('generating lens sample ...')
L=LensSample(reset=False,sigfloor=sigfloor) 

S={}
n={}

print('reading in source catalogue ...')
S[survey]=FastLensSim(survey, fractionofseeing=1)
S[survey].bfac=float(2)
S[survey].rfac=float(2)

    
#################################
# loop for lens i

print("about to load first lens ...")

for i in range(nall):       
    if i%10000==0:

        L.LoadLensPop(i,sourcepop)
        print('Count: ', i, ' of ',nall)

    if i!=0:
        if i%10000==0 or i==300 or i==1000 or i==3000:
            t1=time.perf_counter()
            ti=(t1-t0)/float(i)
            tl=(nall-i)*ti
            tl/=60  #mins
            hl=numpy.floor(tl/(60))
            ml=tl-(hl*60)
            print('time left: '+str(int(hl))+'h '+str(int(ml))+'m ... ',end='')
            print('completed: '+str(i))

    lenspars=L.lens[i]
    if lenspars["lens?"]==False:
        del L.lens[i]
        continue

    lenspars["rl"]["VIS"]=(lenspars["rl"]["r_SDSS"]+\
                           lenspars["rl"]["i_SDSS"]+lenspars["rl"]["z_SDSS"])/3
    for mi in [lenspars["ml"],lenspars["ms"][1]]:
        mi["VIS"]=(mi["r_SDSS"]+mi["i_SDSS"]+mi["z_SDSS"])/3  
    

    lenspars["mag"]={}
    lenspars["msrc"]={}
    lenspars["mag"]={}
    lenspars["msrc"]={}
    lenspars["SN"]={}
    lenspars["bestband"]={}
    lenspars["pf"]={}
    lenspars["resolved"]={}
    lenspars["poptag"]={}
    lenspars["seeing"]={}
    lenspars["rfpf"]={}
    lenspars["rfsn"]={}

    lastsurvey="non"

    S[survey].setLensPars(lenspars["ml"],lenspars["rl"],\
                              lenspars["ql"],reset=True)
        
    #  ms below are read in directly as log of unlensedflux (mJy) from Cai data
    
    for j in range(nsources): 
        S[survey].setSourcePars(lenspars["b"][j+1],lenspars["ms"][j+1],\
                                lenspars["xs"][j+1],lenspars["ys"][j+1],\
                                lenspars["qs"][j+1],lenspars["ps"][j+1],\
                                lenspars["rs"][j+1],sourcenumber=j+1    )

    if survey[:3]+str(i)!=lastsurvey:
        model=S[survey].makeLens(stochasticmode="MP")
        SOdraw=numpy.array(S[survey].SOdraw)
        if type(model)!=type(None):
            lastsurvey=survey[:3]+str(i)
        if S[survey].seeingtest=="Fail":
            lenspars["pf"][survey]={}
            lenspars["rfpf"][survey]={}
            for src in S[survey].sourcenumbers:
                lenspars["pf"][survey][src]=False
                lenspars["rfpf"][survey][src]=False
            continue  #try next survey
    else: 
        S[survey].loadModel(model)
        S[survey].stochasticObserving(mode="MP",SOdraw=SOdraw)
        if S[survey].seeingtest=="Fail":
            lenspars["pf"][survey]={}
            for src in S[survey].sourcenumbers:
                lenspars["pf"][survey][src]=False
            continue  #try next survey
        S[survey].ObserveLens()

    mag,msrc,SN,bestband,pf=S[survey].SourceMetaData(criterion, SNcutA=a,magcut=b,SNcutB=[c,d])
    lenspars["SN"][survey]={}
    lenspars["bestband"][survey]={}
    lenspars["pf"][survey]={}
    lenspars["resolved"][survey]={}
    lenspars["poptag"][survey]=i
    lenspars["seeing"][survey]=S[survey].seeing
    rfpf={}
    rfsn={}
    for src in S[survey].sourcenumbers:
        rfpf[src]=False      
        rfsn[src]=[0]
        lenspars["mag"][src]=mag[src]     # magnification
        lenspars["msrc"][src]=msrc[src]   # lensed flux S500 
        lenspars["SN"][survey][src]=SN[src]
        lenspars["bestband"][survey][src]=bestband[src]
        lenspars["pf"][survey][src]=pf[src]
        lenspars["resolved"][survey][src]=S[survey].resolved[src]
    lenspars["rfpf"][survey]=rfpf
    lenspars["rfsn"][survey]=rfsn


    #############################################
    # generate and save images - CJS 27/9/2020
            
    if lenspars["pf"][survey][1] and n_image < numim:
        image_folder="images_submm"
        image_folder = "%s/%i"%(image_folder, i)
        os.mkdir(image_folder)
        for band in S[survey].bands:
            try:
                img=S[survey].image[band]
                sig=S[survey].sigma[band]
                psf=S[survey].psf[band]
                resid=S[survey].fakeResidual[0][band] # the lens subtracted

        #resid contains the lensed source, with the lens subtracted
        #assuming the subtraction is poisson noise limited (i.e. ideal)

                fits.PrimaryHDU(img).writeto("%s/image_%s.fits"%(image_folder,band),\
                                         overwrite=True)
                fits.PrimaryHDU(sig).writeto("%s/sigma_%s.fits"%(image_folder,band),\
                                         overwrite=True)
                fits.PrimaryHDU(psf).writeto("%s/psf_%s.fits"%(image_folder,band),\
                                         overwrite=True)
                fits.PrimaryHDU(resid).writeto("%s/galsub_%s.fits"%(image_folder,band),\
                                           overwrite=True)
                n_image+=1
                print("saved image",i)
            
            except:
                print("image error, lens",i)


    #############################################
        

    L.lens[i]=None      #delete used data for memory saving
            
    accept=False
    if lenspars["pf"][survey][1]:
        accept=True

    if accept:
        Si+=1
        SSPL[Si]=lenspars.copy() 
        if (Si+1)%1000==0:
            f=open("LensStats/%s_%s_%s_Lens_stats_%i.pkl"%(sourcepop,frac,criterion,chunk),"wb")
            pickle.dump([frac,SSPL],f,2)
            f.close()
            SSPL={} # reset SSPL or memory fills up
            print('chunk ',chunk, 'saved as pkl file')
            chunk+=1 
            
    del L.lens[i]

##########################################

if not os.path.exists("LensStats"):
    os.mkdir("LensStats")
                                          
f=open("LensStats/%s_%s_%s_Lens_stats_%i.pkl"%(sourcepop,frac,criterion,chunk),"wb")
pickle.dump([frac,SSPL],f,2)
f.close()
print('residual chunk ',chunk, 'saved as pkl file')
    

##########################################
    
print('Predicted lenses for', frac, 'of whole sky: ',  Si)

bl=[]
for j in SSPL.keys():
    try: 
        if SSPL[j]["rfpf"][survey][1]:
            bl.append(SSPL[j]["b"][1])
    except KeyError:pass
    
######################    
print("Routine ended")  
