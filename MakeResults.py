''' 
ADJUSTED FOR SUBMM DATA


Note: ms values are not source unlensed AB mags as per original code, but source unlensed logS from Cai data. Also, msrc values not source lensed AB mags, but source lensed flux.
  - amendment by CW, 2019

Updated to Python 3.8, removed source loop etc. 
   - CJS, Aug 2021
   
'''


import pickle
import sys,os
import matplotlib.pyplot as plt
import glob
from __init__ import *


#############################
# Set parameters

sourcepop="submm"   
frac=0.25            # cjs 
criterion = 100      # cjs
survey = "Euclid"   # anachronistic ... to remove


###############################
# Set up file names and parameter structure  

#filename="%s_%s_%s_lists.pkl"%(sourcepop, frac, criterion)   # not used later
filename="%s_%s_%s_results.pkl"%(sourcepop, frac, criterion)  # not used later

lensparsfile="lenses_%s_%s_%s.txt"%(sourcepop, frac, criterion)
f=open(lensparsfile,"w")

    
#initialise parameters as dictionaries    
bl={}
zs={}
zl={}
sigl={}
ql={}
rs={}
ms={}
mag={}
weights={}
    
#create keywords for each parameter dictionary    
for key in ["resolved","rfpf"]:
    bl[key]=[]
    zs[key]=[]
    rs[key]=[]
    ms[key]=[]
    zl[key]=[]
    sigl[key]=[]
    ql[key]=[]
    mag[key]=[]
    rs[key]=[]
    weights[key]=[]

bands=["VIS"]    # note: only one band: VIS (from Euclid) is anachronism


# cw - amended this for Herschel sky survey area for submm (not used now)
# frac = 42000.*1./600.   # or  frac=42000.*1./15000.
       

    
########################################
# read in pickle files, write txt files


# find pathnames for all relevant pickle files of lens stats
filelist=glob.glob("LensStats/%s_%s_%s_Lens_stats_*.pkl"%(sourcepop, frac, criterion))

# start loop for filelist    
chunki=0
ilist=[]

for chunk in filelist:
    print("chunk", chunki)
    chunki+=1        
        
    # read in data from pickle file (chunk)  
    # cw - fracsky is the fraction of idealised lenses examined in MAll
    #   (default = 0.1); frac (here) is 1/(fraction of sky covered by survey)        
    f2=open(chunk,"rb")
    fracsky,sspl=pickle.load(f2)
    fract=frac*fracsky   
    f2.close()
        
        
    I=0
    for i in sspl.keys():
        if i in ilist:  
            continue 
        else:
            try:
                sspl[i]["seeing"][survey]
            except KeyError:
                continue
            f.write("%.2f "%sspl[i]["zl"])
            f.write("%.2f "%sspl[i]["zs"][1])
            f.write("%.2f "%sspl[i]["b"][1])
            f.write("%.2f "%sspl[i]["sigl"])
            f.write("%.2f "%sspl[i]["ql"])
            f.write("%.2f "%sspl[i]["rl"]["g_SDSS"])
            for band in bands:       #cw - only VIS band used for Euclid
                f.write("%.2f "%sspl[i]["ml"][band])
            f.write("%.2f "%sspl[i]["xs"][1])
            f.write("%.2f "%sspl[i]["ys"][1])
            f.write("%.2f "%sspl[i]["qs"][1])
            f.write("%.2f "%sspl[i]["ps"][1])
            f.write("%.2f "%sspl[i]["rs"][1])
            for band in bands:                  
                f.write("%.12f "%sspl[i]["ms"][1][band]) # unlensed source logS values
            f.write("%.2f "%sspl[i]["mag"][1])   
            f.write("%.2f "%sspl[i]["msrc"][1]['VIS'])   # lensed source total flux
            for band in bands:
                f.write("%.2f "%sspl[i]["seeing"][survey][band]) 
            if survey!="Euclid":
                f.write("%.2f "%sspl[i]["rfsn"][survey][1][0])
            f.write("\n")
                
            ilist.append(str(i))
            if sspl[i]["pf"][survey][1]==False:continue

                
            try:
                if sspl[i]["resolved"][survey][1][sspl[i]["bestband"][survey][1]]:
                    bb=sspl[i]["bestband"][survey][1]
                    if sspl[i]["mag"][1]<3:
                        continue
                    if sspl[i]["SN"][survey][1][bb][0]<20:
                        continue 
                    bl["resolved"].append(sspl[i]["b"][1])
                    weights["resolved"].append(1./fract)
                    zs["resolved"].append(sspl[i]["zs"][1])
                    rs["resolved"].append(sspl[i]["rs"][1])
                    zl["resolved"].append(sspl[i]["zl"])
                    sigl["resolved"].append(sspl[i]["sigl"])
                    ql["resolved"].append(sspl[i]["ql"])
                    mag["resolved"].append(sspl[i]["mag"][1])
                    ms["resolved"].append(sspl[i]["ms"][1]["g_SDSS"])

                if sspl[i]["rfpf"][survey][1]:
                    if sspl[i]["rfsn"][survey][1][0]<20:continue
                    if sspl[i]["resolved"][survey][1]["RF"]==False:continue
                    bl["rfpf"].append(sspl[i]["b"][1])
                    weights["rfpf"].append(1./fract)
                    zs["rfpf"].append(sspl[i]["zs"][1])
                    rs["rfpf"].append(sspl[i]["rs"][1])
                    zl["rfpf"].append(sspl[i]["zl"])
                    sigl["rfpf"].append(sspl[i]["sigl"])
                    ql["rfpf"].append(sspl[i]["ql"])
                    mag["rfpf"].append(sspl[i]["mag"][1])
                    ms["rfpf"].append(sspl[i]["ms"][1]["g_SDSS"])

            except KeyError:
                pass

f.close()        
   
    
# write pickle file    
f=open(filename,"wb")
pickle.dump([weights,bl,zs,rs,ms,zl,sigl,ql,mag],f,2)
f.close()
         

##################################################


print('txt file saved: ', lensparsfile)
print('pkl file saved: ', filename)
print("Routine ended!")

