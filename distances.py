"""

*** Astrodistances.py ***

REPLACES THE DISTANCES MODULE WITH FlatLambdaCDM ASTROPY COSMOLOGY

A module to compute cosmological distances, including:
    comoving_distance (Dc)
    angular_diameter_distance (Da)
    luminosity_distance (Dl)
    comoving_volume (volume)
 - CW

"""


import astropy.cosmology  
from astropy.cosmology import FlatLambdaCDM  

'''
this is the class needed for Astropy distance functions; we could alternatively import
FlatwCDM, rather than wCDM, and not bother with Ode parameter (as FLAT => Ode = 1-Om).
'''

import warnings   # original Collett code and comment
#warnings.warn("Default cosmology is Om=0.3,Ol=0.7,h=0.7,w=-1 and distance units are Mpc!",ImportWarning)

cosmo=(0.324,0.667)  # NOTE COSMOLOGY (FLAT with LAMBDA)

class Distance(FlatLambdaCDM):   # INTRODUCTION OF FlatLambdaCDM COSMOLOGY

# Distance is the existing class name used in Collett code; to avoid having to change the name throughout
# the code I have simply kept the classname Distance but have it inherit the FlatLambda/w CDM properties

    def __init__(self):
        FlatLambdaCDM.__init__(self,H0=(cosmo[1]*100),Om0=cosmo[0])
        self.OMEGA_M = cosmo[0]
        print('Om0 = ', self.OMEGA_M)
        print('w = -1.0')
        print('h = ',cosmo[1])
        
       # self.w = -1.
       # self.wpars = None          used in Collett code to accommodate variable w
       # self.w_analytic = False    used in Collett code to accommodate variable w
        self.Dc = self.co_distance
        self.Dt = self.co_t_distance
        self.Dm = self.co_t_distance
        self.Da = self.ang_diam_distance
        self.Dl = self.lum_distance
        self.dm = self.dist_modulus
        self.volume = self.co_volume
        
     # original Collett code has h defined as an attribute, but this is autmatically defined within Astropy

    def co_distance(self,z):
        return self.comoving_distance(z).value
        
    def co_t_distance(self,z):
        return self.comoving_distance(z).value    
        
    def ang_diam_distance(self,z1,z2=0):
        if z2<z1:
            z1,z2 = z2,z1
        return self.angular_diameter_distance_z1z2(z1,z2).value   
    
    def lum_distance(self,z):
        return self.luminosity_distance(z).value 
        
    def dist_modulus(self,z):
        if z>0:                          # needed to avoid runtime 'divide by zero' error
            return self.distmod(z).value    
        else:
            return 0 
        
    def co_volume(self,z):
        return self.comoving_volume(z).value 