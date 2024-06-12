# a script to write out pythia cards
# performs a scan of lifetimes and masses for the higgs and vector portal models

from dark_shower import dark_shower


# make a dark shower object, specifying the decay portal using the mandatory "portal" flag. 
# The value of this flag must be "photon", "vector", "gluon", "higgs" or "darkphoton"
vector_portal=dark_shower("vector")
higgs_portal=dark_shower("higgs")

# writing out a pythia card

xio = 1.0
xil = 1.0
masses = [1,2,3,5,10]
ctaus = [3, 10, 30, 100] # mm
production = "wh" # ggf or zh

def to_cm(mm):
    return mm/10.

for mass in masses:
    for ctau in ctaus:
    	vector_portal.pythia_card("cards/{}_vectorportal_m-{}_xio-{}_xil-{}_ctau-{}.cmnd".format(production,mass,xio,xil,ctau),
    	                         m=mass,
    	                         xi_omega=1.0,
    	                         xi_Lambda=1.0,
    	                         ctau=to_cm(ctau), #takes in cm
    	                         Nc=3,
    	                         Nf=1,
    	                         mH=125,
    	                         mayDecay=True,
    	                         production=production,
    	                         userName="Karri DiPetrillo")

    	higgs_portal.pythia_card("cards/{}_higgsportal_m-{}_xio-{}_xil-{}_ctau-{}.cmnd".format(production,mass,xio,xil,ctau),
    	                         m=mass,
    	                         xi_omega=1.0,
    	                         xi_Lambda=1.0,
    	                         ctau=to_cm(ctau), #takes in cm
    	                         Nc=3,
    	                         Nf=1,
    	                         mH=125,
    	                         mayDecay=True,
    	                         production=production,
    	                         userName="Karri DiPetrillo")