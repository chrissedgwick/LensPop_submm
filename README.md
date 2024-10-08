# LensPop_submm
(LensPop model adapted for simulating lensed submillimetre galaxies)
​

The original LensPop model by Tom Collett (https://github.com/tcollett/LensPop) has been adapted for simulating lensed submillimetre galaxies.
​
For details of changes to the original model, please see: 

(a) Charles Weiner, 2019, Strong Gravitational Lens Modelling, Master Thesis, The Open University; and

(b) Chris Sedgwick, Stephen Serjeant, and Charles Weiner, 2024, The detection of strongly-lensed submillimetre galaxies (in preparation)


​
The main changes to the original model have been:
    
(1) A submillimetre source population (generated from data in Cai et al., 2013, ApJ, 768, 21, as described in the above two sources) has replaced the simulated LSST population in the original code. This is supplied in the file submmdataCai.txt (note size 3Gb: still to upload).

​
(2) The criteria for accepting a lens is detectable have changed; in addition to the Einstein radius criterion, a much simpler criterion using the 500 micron flux level is applied. In the original paper proposing this (Negrello et al., 2010, Science, 330, 800) the flux cutoff was set at >100mJy; later work has used >80mJy; some of the results in the Sedgwick et al. paper cited above use >10mJy to estimate the total number of submm lenses.


(3) The code has been updated to Python 3.8, although the resulting changes were fairly minor.

Other minor amendments to the code are indicated in the individual routines. The structure of the model and the names of the various routines are unchanged.

A Jupyter Notebook is provided which includes the run sequence, some example results, and flowcharts showing the relationship between components of the model, and the flow of data input and output between routines. These flowcharts are equally applicable to the original code.

Chris Sedgwick

22 September 2024

===============================================================================

LICENSE

​
The original code is open access, but please email thomas.collett@port.ac.uk to tell him that you are using it. Please cite Collett, T., 2015, ApJ, 811, 20, if you make use of these codes in your work.

The amended code for submillimatre lenses is also open access. Please email chris@sedgwick.uk.net and cite Sedgwick et al. 2024 (in prep.) if you make use of this in your work. 

===============================================================================

INSTALLATION NOTE FROM ORIGINAL CODE


​
The code is mostly python, so should work out of the box (if you have standard astrophysical libraries installed already) except the deflection angles code which must be compiled using:
​
cd pylens
f2py -c -m powerlaw powerlaw.f

===============================================================================    
   
HOW TO REPRODUCE RESULTS in Sedgwick et al. 2024 


python   MakeLensPop.py (takes ~7 hours, makes all the lenses on the sky) to generate idealized foreground lens population

python ModelAll.py (takes ~24 hrs for submm 10mJy and frac=0.25, fraction of sky) to generate and observe the idealized lens population

python MakeResults.py (a few minutes) to generate a .txt file with all parameters of results


These are shown, with earlier results for specified parameters, in the Jupyter Notebook provided. Timimgs shown are for a reasonably recent laptop.
