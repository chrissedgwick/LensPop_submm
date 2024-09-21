'''
**** ADJUSTED FOR SUBMM DATA **** 

CFW  13/3/2019
'''


import pickle,numpy
class Survey():
    def  __init__(self,Name):
        self.zeroexposuretime=1
        self.strategy="resolve"
        self.strategyx=1

        if Name=="submm":
            self.pixelsize=0.1    
            self.side=200
            self.bands=['VIS']
            self.zeropoints=[25.5]
            self.zeroexposuretime=1.
            self.skybrightnesses=[22.2]
            self.exposuretimes=[1610]
            self.gains=[1]
            self.seeing=[.2]
            self.nexposures=4
            self.degrees_of_survey=600   # amended to reflect Herschel
            self.readnoise=(4.5)
            twodVIS=numpy.array([[0.17,22.2],[0.17,22.2]])
            self.stochasticobservingdata=[twodVIS]

        elif Name=="ideal":
            self.pixelsize=0.05
            self.side=400
            self.bands=['g','r','i']
            self.zeropoints=[24.5,24.5,24.5]
            self.zeroexposuretime=4.
            self.skybrightnesses=[220,220,220]
            self.exposuretimes=[22400000000,22400000000,22400000000]
            self.gains=[1,1,1]
            self.seeing=[.05,0.05,0.05]
            self.nexposures=1
            self.readnoise=(.005)
            twodr=numpy.array([[0.1,220],[0.1,220]])
            twodg=numpy.array([[0.1,220],[0.1,220]])
            twodi=numpy.array([[0.1,220],[0.1,220]])
            self.stochasticobservingdata=[twodg,twodr,twodi]
            self.degrees_of_survey=41253
        else:
            print("I don't know that survey")
            exit()



        #convert bandnames into the required formats
        for i in range(len(self.bands)):
            bandname=self.bands[i]
            if bandname=="g":
                self.bands[i]="g_SDSS"
            if bandname=="r":
                self.bands[i]="r_SDSS"
            if bandname=="i":
                self.bands[i]="i_SDSS"
            if bandname=="z":
                self.bands[i]="z_SDSS"
            if bandname=="F814":
                self.bands[i]="F814W_ACS"

        degrees_of_whole_sky=41253.
        self.f_sky=float(self.degrees_of_survey)/degrees_of_whole_sky

