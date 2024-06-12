#################################################################################################
## Package to generate self-consistent configuration cards for the pythia 8 hidden valley module.
## Based on arXiv xxxxxxx. Please cite this paper when using this package.
##
## Writen by Simon Knapen (smknapen@lbl.gov)
## on 02/24/2023 
#################################################################################################
import os 
import sys
import sys
import glob
import math
import numpy as np
from scipy.interpolate import interp1d
from scipy.special import exp1
from scipy import integrate
from datetime import date

data_dir = os.path.dirname(__file__)+"/data/"

GeVtimescm=1.98e-14
TeV_in_GeV=1e3
eV_in_GeV=1e-9


def round_sig(x, sig=3):
  if float(x)==0.0:
    return 0.0
  else:
    return round(float(x), sig-int(math.floor(math.log10(abs(float(x)))))-1)


class dark_shower_two_flavor():
    def __init__(self,m_pi2,m_eta,Lambda,sin_theta,mA,g,epsilon): 
        
      #######################
      # Universal variables
      #######################
      
      # constants
      self.m_h=125.0 # GeV
      
      self.m_e=0.511e-3 # GeV
      self.m_mu=0.105   # GeV
      self.m_tau=1.778  # GeV
      self.m_s=0.095  # GeV
      self.m_c=1.19  # GeV
      self.m_b=4.19  # GeV
      self.m_t=174.0  # GeV
      
      self.m_pi=0.139  # GeV
      self.m_K=0.494  # GeV
      self.m_D=1.870  # GeV
      self.m_B=5.279 # GeV
      
      self.m_W=80.4 # GeV
      self.m_Z=91.2 # GeV
      
      self.G_F=1.166e-5 # GeV^-2
      self.Nc=3.0 # SM colors
      self.alpha_EM=1.0/137. # fine structure constant
      self.vev=246.0 # higgs vev
      
      # number of flavors as function of mass scale
      self.Nf_L=(lambda m: 3.0+float(m>self.m_c)+float(m>self.m_b))
      self.Nf_H=(lambda m: 6.0-self.Nf_L(m))
      
      # running of alpha_s
      [alphasx,alphasy]=np.loadtxt(data_dir+"gluon_portal/alphas_running.txt").T
      self.alphas=interp1d(alphasx,alphasy)
      
      # load and check model parameters
      self.load_params(m_pi2,m_eta,Lambda,sin_theta,mA,g,epsilon)
      # load dark photon
      self.load_dark_photon()
      # compute dark sector dependent parameters and branching ratios
      self.compute_dark_sector()  
    
    
    def load_params(self,m_pi2=0,m_eta=0,Lambda=0,sin_theta=0,mA=0,g=0,epsilon=0):
      #######################
      # Input model variables
      #######################
      if m_pi2>0:
        self.m_pi2=m_pi2              # pi_2 mass
      if m_eta>0:
        self.m_eta=m_eta              # eta mass
      if Lambda>0:
        self.Lambda=Lambda            # confinement scale
      if sin_theta>0:
        self.sin_theta=sin_theta      # pi_3 - pi_0 mixing angle, parametrizes isospin violation
        self.cos_theta=np.sqrt(1.0-self.sin_theta**2)
      if mA>0:
        self.mA=mA                    # dark photon mass
      if g>0:
        self.g=g                      # dark photon gauge coupling
      if epsilon>0:
        self.epsilon=epsilon          # dark photon mixing
      
      self.errors=[]
      self.warnings=[]
      if(self.m_pi2>self.m_eta):
        self.errors.append("Error: please choose m_pi2 < m_eta")
      if(self.m_pi2>self.Lambda):
        self.errors.append("Error: please choose m_pi2 < Lambda")  
      if(self.sin_theta<0 or self.sin_theta>np.sqrt(2.0)/2.0):
        self.errors.append("Error: please choose sin_theta between 0 and \sqrt(2)/2")  
      if(self.epsilon<0 or self.epsilon>1):
        self.errors.append("Error: please choose epsilon between 0 and 1")
      if(self.g<0 or self.g>4*np.pi):
        self.errors.append("Error: please choose g between 0 and 4\pi")
      if(self.g>2):
        self.warnings.append("Warning: dark U(1) may be non-perturbative")
        
    # routine to change model parameters
    def update_params(self,m_pi2=0,m_eta=0,Lambda=0,sin_theta=0,mA=0,g=0,epsilon=0): 
      # load and check new model parameters
      self.load_params(m_pi2=m_pi2,
                       m_eta=m_eta,
                       Lambda=Lambda,
                       sin_theta=sin_theta,
                       mA=mA,
                       g=g,
                       epsilon=epsilon)
      # recompute dark sector dependent parameters and branching ratios
      self.compute_dark_sector()
      
    def load_dark_photon(self):
      #########################
      # A' Branching ratios
      #########################
      # PDG codes of decay products
      self.decay_PDG_codes={"d dbar":"1 -1",\
                              "u ubar":"2 -2",\
                              "s sbar":"3 -3",\
                              "c cbar":"4 -4",\
                              "b bbar":"5 -5",\
                              "e+ e-":"11 -11",\
                              "mu+ mu-":"13 -13",\
                              "tau+ tau-":"15 -15",\
                              "pi+ pi-":"211 -211",\
                              "K+ K-":"321 -321",\
                              "pi+ pi- pi0":"211 -211 111",\
                              "pi+ pi- pi+ pi-":"211 -211 211 -211",\
                              "pi+ pi- pi0 pi0":"211 -211 111 111"
                             }
      
      self.eta_decay_PDG_codes={"pi1 pi1 pi3":"4900211 4900211 4900111",\
                              "pi2 pi2 pi3":"-4900211 -4900211 4900111",\
                              "A' A'":"999999 999999",\
                              "pi2 A'":"4900211 999999",
                              "pi2 d dbar":"4900211 1 -1",\
                              "pi2 u ubar":"4900211 2 -2",\
                              "pi2 s sbar":"4900211 3 -3",\
                              "pi2 c cbar":"4900211 4 -4",\
                              "pi2 b bbar":"4900211 5 -5",\
                              "pi2 e+ e-":"4900211 11 -11",\
                              "pi2 mu+ mu-":"4900211 13 -13",\
                              "pi2 tau+ tau-":"4900211 15 -15",\
                              "pi2 pi+ pi-":"4900211 211 -211",\
                              "pi2 K+ K-":"4900211 321 -321",\
                              "pi2 pi+ pi- pi0":"4900211 211 -211 111",\
                              "pi2 pi+ pi- pi+ pi-":"4900211 211 -211 211 -211",\
                              "pi2 pi+ pi- pi0 pi0":"4900211 211 -211 111 111"  
                               }
      self.eta_decay_PDG_codes.update(self.decay_PDG_codes)
        
      # scale at which we switch to quarks instead of exclusive channels
      self.pert_QCD=1.65 
        
      # loop over all files in the folder and interpolate
      Brfilelist=glob.glob(data_dir+"darkphoton_portal/*")
      self.branching_interpol={}
      for filename in Brfilelist:
        with open(filename, "r") as f:
          f.readline() # skip first line
          PDGstring=f.readline().replace("PDG codes: ","").replace("\n","") 
        # verify if the PDG codes in the file match those in self.decay_PDG_codes, throw exception if not
        try:  
          channel=list(self.decay_PDG_codes.keys())[list(self.decay_PDG_codes.values()).index(PDGstring)]
          datatemp=np.loadtxt(filename,skiprows=2).T
          self.branching_interpol[channel]=interp1d(datatemp[0],datatemp[1],fill_value=0.0,bounds_error=False,copy=True)
        except:
          print(f"PDG codes in decay file {filename} do not match known decay channel.")
        
      # formula for perturbative partial widths, setting epsilon=1
      self.pert_width=(lambda m,m_f,Nc,q: 0.0 if m<2*m_f else Nc*q**2*self.alpha_EM/3.0*m*np.sqrt(1.0-4.0*m_f**2/m**2)*(1+2*m_f**2/m**2))
        
      # perturbative partial widths, setting epsilon=1
      self.partial_width={}
      self.partial_width["e+ e-"]=(lambda m: self.pert_width(m,self.m_e,1.0,1.0))
      self.partial_width["mu+ mu-"]=(lambda m: self.pert_width(m,self.m_mu,1.0,1.0))
      self.partial_width["tau+ tau-"]=(lambda m: self.pert_width(m,self.m_tau,1.0,1.0))
      self.partial_width["u ubar"]=(lambda m: 0.0 if m<=self.pert_QCD else self.pert_width(m,0.0,3.0,2.0/3.0))
      self.partial_width["c cbar"]=(lambda m: 0.0 if m<=self.pert_QCD else self.pert_width(m,self.m_c,3.0,2.0/3.0))
      self.partial_width["d dbar"]=(lambda m: 0.0 if m<=self.pert_QCD else self.pert_width(m,0.0,3.0,1.0/3.0))
      self.partial_width["s sbar"]=(lambda m: 0.0 if m<=self.pert_QCD else self.pert_width(m,self.m_s,3.0,1.0/3.0))
      self.partial_width["b bbar"]=(lambda m: 0.0 if m<=self.pert_QCD else self.pert_width(m,self.m_b,3.0,1.0/3.0))
        
      # total width, reinstate epsilon
      self.A_total_width=(lambda m, eps: eps**2*self.partial_width["e+ e-"](m)/self.branching_interpol["e+ e-"](m) if m<=self.pert_QCD else eps**2*sum([self.partial_width[channel](m) for channel in self.partial_width.keys()]))
      self.A_ctau=GeVtimescm/self.A_total_width(self.mA,self.epsilon)
      
      # branching ratios
      self.branching_ratios={}
      for channel in self.partial_width.keys():
        self.branching_ratios[channel]=(lambda m, channel=channel: self.partial_width[channel](m)/self.A_total_width(m,1.0))
        
        self.branching_ratios["pi+ pi-"]=(lambda m: 0.0 if (m>self.pert_QCD or m<2.0*self.m_pi) else self.branching_interpol["pi+ pi-"](m))
        self.branching_ratios["K+ K-"]=(lambda m: 0.0 if (m>self.pert_QCD or m<2.0*self.m_K)else self.branching_interpol["K+ K-"](m))
        self.branching_ratios["pi+ pi- pi0"]=(lambda m: 0.0 if (m>self.pert_QCD or m<3.0*self.m_pi)else self.branching_interpol["pi+ pi- pi0"](m))
        self.branching_ratios["pi+ pi- pi+ pi-"]=(lambda m: 0.0 if (m>self.pert_QCD or m<4.0*self.m_pi) else self.branching_interpol["pi+ pi- pi+ pi-"](m))
        self.branching_ratios["pi+ pi- pi0 pi0"]=(lambda m: 0.0 if (m>self.pert_QCD or  m<2.0*self.m_pi) else self.branching_interpol["pi+ pi- pi0 pi0"](m))
              
    def compute_dark_sector(self): # calculate spectrum and branching ratios
      ###########################################################
      # Derived hidden sector variables and consistency checks
      ###########################################################
      self.alpha_D=self.g**2/(4.0*np.pi)        # g^2/4pi
      self.m0=self.Lambda                       # mass parameter for the eta. This choice is not valid for large Nc.
      self.f=self.Lambda/(4.0*np.pi)            # estimate, really f is an additional free parameter.
      self.mbar=self.m_pi2**2/self.Lambda       # sum quark masss
      self.delta_m=self.m_pi2**2/(2.0*self.Lambda)*np.tan(2.0*np.arcsin(self.sin_theta)) # difference quark masss
      self.mquark1=(self.mbar-self.delta_m)/2.0
      self.mquark2=(self.mbar+self.delta_m)/2.0
      
      if(self.mA<self.g*self.f):
        print("Fatal error: mA too low. Please increase mA or decreasing g.")
        sys.exit()
      else:  
        self.v_dark=np.sqrt(self.mA**2/self.g**2-self.f**2)   # vev of dark sector higgs
      
      self.y11=self.delta_m/(np.sqrt(2)*self.v_dark)
      self.y12=self.mbar/(np.sqrt(2)*self.v_dark)
      
      if(self.y11>4.0*np.pi or self.y12>4.0*np.pi):
        self.errors.append("Error: yukawas are too large. Please increase mA or decreasing g.")
      if(self.y11>2.0 or self.y12>2.0):
        self.warnings.append("Warning: yukawas are almost non-perturbative, proceed with caution. Can be addressed by increasing mA or decreasing g.")  
      
      self.m_pi1=self.m_pi2*(1.0+self.f**2/self.v_dark**2)
      self.m_pi3=self.m_pi2*(1.0-np.tan(np.arcsin(self.sin_theta))*self.delta_m/self.mbar)
      if(np.abs(self.m_pi1-self.m_pi2)/self.m_pi2>0.1):
        self.warnings.append("Warning: pi1 and pi2 masses less degenerate than 10\%. pi1 and pi2 are taken to be degenerate in the pythia cards.")
      
      self.pi1_mayDecay=False
      self.pi2_mayDecay=False
      if(2.0*self.mA<self.m_pi3):
        self.pi3_mayDecay=True
      else:
        self.pi3_mayDecay=False
      
      ########################
      # dark photon decays
      ########################
      self.A_ctau=GeVtimescm/self.A_total_width(self.mA,self.epsilon)
      
      #########################
      # eta decays
      #########################
      x_pi=self.m_pi2/self.m_eta
      x_A=self.mA/self.m_eta
      x_e=self.m_e/self.m_eta
      # specify kinematic endpoints for phase space integrals in eta -> pi f f
      kin_endpoint={'e+ e-': 2.0*self.m_e, 
              'pi+ pi-': 2.0*self.m_pi,
              'K+ K-': 2.0*self.m_K,
              'pi+ pi- pi0': 3.0*self.m_pi,
              'pi+ pi- pi+ pi-': 4.0*self.m_pi,
              'pi+ pi- pi0 pi0': 4.0*self.m_pi,
              'mu+ mu-': 2.0*self.m_mu,
              'tau+ tau-': 2.0*self.m_tau,
              'u ubar': 0.0,
              'c cbar': 2.0*self.m_c, 
              'd dbar': 0.0, 
              's sbar': 2.0*self.m_s,
              'b bbar': 2.0*self.m_b}
      
        
      if(self.m_pi3+2.0*self.m_pi2<self.m_eta): # eta -> 3pi is open and is assumed to always dominate
        self.eta_ctau=0.0 # always prompt
        self.eta_branching_ratios={"pi1 pi1 pi3":0.5,"pi2 pi2 pi3":0.5}
          
      elif(self.m_pi2+self.mA<self.m_eta): # eta -> A' pi2 is open
        Gamma_eta_to_Api=self.g**2*self.sin_theta**2*self.m_eta**3/(16*np.pi*self.mA**2)\
        *((1-x_pi**2-x_A**2)**2-4.0*x_pi**2*x_A**2)**(1.5)
          
        if(2.0*self.mA<self.m_eta): # eta -> A' A' is also open
          Gamma_eta_to_AA=self.alpha_D**2*self.Nc**2/(64.0*np.pi**3)*self.cos_theta**2*self.m_eta**3/(self.f**2)
          Gamma_eta=Gamma_eta_to_Api+Gamma_eta_to_AA
          self.eta_ctau=GeVtimescm/Gamma_eta
          self.eta_branching_ratios={"A' A'":Gamma_eta_to_AA/Gamma_eta,"pi2 A'":Gamma_eta_to_Api/Gamma_eta}
        else:
          self.eta_ctau=GeVtimescm/Gamma_eta_to_Api
          self.eta_branching_ratios={"pi2 A'":1.0}
        
      elif(2.0*self.mA<self.m_eta): # only eta -> A' A' is open
          Gamma_eta_to_AA=self.alpha_D**2*self.Nc**2/(64.0*np.pi**3)*self.cos_theta**2*self.m_eta**3/(self.f**2)
          self.eta_ctau=GeVtimescm/Gamma_eta_to_AA
          self.eta_branching_ratios={"A' A'":1}
          
      else: # only eta -> pi f f is open. All Br computed relative to eta -> pi e e
          self.eta_branching_ratios={}
          self.eta_partial_widths={}
          Gamma_eta=0.0
          prefactor=self.epsilon**2*self.alpha_D*self.alpha_EM*self.sin_theta**2/(12*np.pi)
          dfdx=(lambda x: (x+2*x_e**2)/(x-x_A**2)**2*np.sqrt(x-4*x_e**2)*(((x-1)**2-2*(x+1)*x_pi**2+x_pi**4)/x)**(1.5))
          ### loop over all final states
          # use A' branching ratios to compute eta branching ratios
          for Ap_channel in self.branching_ratios.keys(): # Ap branching ratio
            channel="pi2 "+Ap_channel
            xstart=kin_endpoint[Ap_channel]**2/self.m_eta**2
            xend=(1-x_pi)**2
            if(xstart<xend):
              pwidth=prefactor*integrate.quad((lambda x: dfdx(x)*self.branching_ratios[Ap_channel](np.sqrt(self.m_eta**2*x))/self.branching_ratios["e+ e-"](np.sqrt(self.m_eta**2*x))), xstart, xend)[0]
              self.eta_partial_widths.update({channel: pwidth})
            else:
              self.eta_partial_widths.update({channel: 0.0})
            Gamma_eta=Gamma_eta+self.eta_partial_widths[channel]
          self.eta_ctau=GeVtimescm/Gamma_eta  
          for channel in self.eta_partial_widths.keys():
            self.eta_branching_ratios.update({channel:self.eta_partial_widths[channel]/Gamma_eta})
    
      #########################
      # pi decays
      #########################
    
      if(2.0*self.mA<self.m_pi3): # pi3 -> A' A' is open
        Gamma_pi3_to_AA=self.alpha_D**2*self.Nc**2/(64*np.pi**3)*self.sin_theta**2*self.m_pi3**3/(self.f**2)
        self.pi3_ctau=GeVtimescm/Gamma_pi3_to_AA #temp
      
      #########################
      # warnings and errors
      #########################
      
      for error in self.errors:
        print(error)
      for warning in self.warnings:
        print(warning)  

    #########################
    # output routines
    #########################
   
    def get_spectrum(self):
      return f"""
##############################
# meson/dark photon masses
# m_pi1:    {round_sig(self.m_pi1)}  GeV
# m_pi2:    {round_sig(self.m_pi2)}  GeV
# m_pi3:    {round_sig(self.m_pi3)}  GeV
# m_eta:    {round_sig(self.m_eta)}  GeV
# m_A:      {round_sig(self.mA)}    GeV
##############################"""
    def print_spectrum(self):
      print(self.get_spectrum())
        
    def get_param(self):
      return f"""
##############################
# meson/dark photon masses
# m_pi1:     {round_sig(self.m_pi1)}  GeV
# m_pi2:     {round_sig(self.m_pi2)}  GeV
# m_pi3:     {round_sig(self.m_pi3)}  GeV
# m_eta:     {round_sig(self.m_eta)}  GeV
# m_A:       {round_sig(self.mA)}    GeV
# other parameters
# Nc:        3
# Nf:        2
# m quark 1: {round_sig(self.mquark1)}   GeV
# m quark 2: {round_sig(self.mquark2)}   GeV
# Lambda:    {round_sig(self.Lambda)}   GeV
# m0:        {round_sig(self.m0)}   GeV
# f:         {round_sig(self.f)}   GeV
# v_dark:    {round_sig(self.v_dark)}   GeV
# g:         {round_sig(self.g)}
# epsilon:   {round_sig(self.epsilon)}
# y11:       {round_sig(self.y11)}
# y12:       {round_sig(self.y12)}
##############################"""
    def print_param(self):
      print(self.get_param())
      
      
    def get_decay_table(self):
      string=f"""
##############################
# decay table
# pi1:    stable
# pi2:    stable\n"""
      if(self.pi3_mayDecay):
        string=string+f"""# pi3:    A' A'   1.0
#         ctau = {round_sig(self.pi3_ctau)} cm\n"""
      else:
         string=string+f"""# pi3:    stable\n"""
      for i, mode in enumerate(self.eta_branching_ratios.keys()):
        if i==0:
          string=string+("{:<35}".format(f"""# eta:    {mode}"""))
          string=string+f"""{round_sig(self.eta_branching_ratios[mode])}\n"""
        else:
          string=string+("{:<35}".format(f"""#         {mode}"""))
          string=string+f"""{round_sig(self.eta_branching_ratios[mode])}\n"""
      string=string+f"""#         ctau = {round_sig(self.eta_ctau)} cm\n"""
      for i, mode in enumerate(self.branching_ratios.keys()):
        if i==0:
          string=string+("{:<35}".format(f"""# A':     {mode}"""))
          string=string+f"""{round_sig(self.branching_ratios[mode](self.mA))}\n"""
        else:
          string=string+("{:<35}".format(f"""#         {mode}"""))
          string=string+f"""{round_sig(self.branching_ratios[mode](self.mA))}\n"""
          #string=string+f"""{self.branching_ratios[mode](self.mA)}\n"""
      string=string+f"""#         ctau = {round_sig(self.A_ctau)} cm
############################## """
      return string
      
    def print_decay_table(self):
      print(self.get_decay_table())
    
    #########################
    # output routines
    #########################  
      
    def pythia_card(self,filename,userName="Simon"):
      file= open(filename,"w")
      today = date.today()
      file.write(
        f"""
##############################################################################################
# Pythia 8 card for dark shower models with flavor violation, as defined in arXiv xxxx.xxxxx""")
      file.write(self.get_param())
      file.write(self.get_decay_table())
      file.write("""\n""")
      for error in self.errors:
        file.write(f"""# """+error+"""\n""")
      for warning in self.warnings:
        file.write(f"""# """+warning+"""\n""")  
      file.write(f"""
# Generated by {userName} on {today.strftime("%B %d, %Y")}
# lines starting with ! or # are commented out
##############################################################################################
        """)
      
      file.write(f"""  
! Settings used in the main program.
Main:numberOfEvents = 1000                     ! number of events to generate
Beams:eCM = 13000.                              ! CM energy of collision

! Settings related to output in init(), next() and stat().
Init:showChangedSettings = on                   ! list changed settings
Init:showChangedParticleData = on              ! list changed particle data
Next:numberCount = 1000                         ! print message every n events
Next:numberShowInfo = 0                         ! print event information n times
Next:numberShowProcess = 0                      ! print process record n times
Next:numberShowEvent = 2                        ! print event record n times

! For debugging purposes only
! PartonLevel:ISR = off
! PartonLevel:FSR = off
! PartonLevel:MPI = off
! HadronLevel:all = off

! Production settings
! decay the Higgs to two dark quark, turn off all SM branching ratios
HiggsSM:gg2H = on
25:m0 =125
25:addChannel = 1 0.5 102 4900101 -4900101
25:addChannel = 1 0.5 102 4900102 -4900102
25:0:onMode=0
25:1:onMode=0
25:2:onMode=0
25:3:onMode=0
25:4:onMode=0
25:5:onMode=0
25:6:onMode=0
25:7:onMode=0
25:8:onMode=0
25:9:onMode=0
25:10:onMode=0
25:11:onMode=0
25:12:onMode=0
25:13:onMode=0

! HiddenValley Settings
HiddenValley:Ngauge = 3                          ! number of colors
HiddenValley:nFlav = 2                           ! number of flavors
HiddenValley:fragment = on
HiddenValley:FSR = on
HiddenValley:alphaOrder = 1                      ! use running coupling
HiddenValley:setLambda = on                      ! only for pythia 8.309 and higher
HiddenValley:Lambda = {self.Lambda}              ! dark confinement scale
HiddenValley:pTminFSR = {1.1*self.Lambda}        ! pT cut off on dark shower (IR regulator)
HiddenValley:spinFv=0                            ! spin of bifundamentals: not used, but set for consistency
HiddenValley:probVector=0.75                     ! fraction of hadronization to spin 1 
HiddenValley:separateFlav = on                   ! allow for non-degenerate mesons
HiddenValley:probKeepEta1 = 1.0                  ! suppression factor for eta hadronization

! dark sector meson mass spectrum
4900101:m0 = {self.Lambda+self.mquark1}          ! Dark Quark Mass, following arXiv:2203.09503
4900102:m0 = {self.Lambda+self.mquark2}          ! Dark Quark Mass, following arXiv:2203.09503
4900111:m0 = {self.m_pi2}                        ! Setting pion Mass
4900211:m0 = {self.m_pi2}                        ! Setting pion Mass
4900221:m0 = {self.m_eta}                        ! Setting eta Mass
4900113:m0 = {self.Lambda}                       ! Setting rho Mass
4900213:m0 = {self.Lambda}                       ! Setting rho Mass

! vector meson decay chains
""")
      if(self.Lambda > 2.0*self.m_pi2): ## This is a little approximate, as pythia assumes pi^+ and pi^- are anti-particles
        file.write(f"""4900113:addChannel = 1 1.00 91 4900211 -4900211 \n""")
        file.write(f"""4900213:addChannel = 1 1.00 91 4900211 4900111 \n""")
      else:
        file.write("""4900113:onMode = 0\n""")
        file.write("""4900213:onMode = 0\n""")

      if(3.0*self.m_pi2>self.m_eta>self.mA+self.m_pi2 or self.mA<0.5*self.m_pi3):
        file.write(f"""
! define the dark photon A'""")
        file.write(f"""
999999:all = GeneralResonance void 1 0 0 {self.mA} 0.001 0. 0. 0.        ! dark photon A'
""")
        file.write("""
! A' decay modes""")
        for channel in self.decay_PDG_codes:
          if round_sig(self.branching_ratios[channel](self.mA))>0.0:
            file.write(f"""
999999:addChannel = 1 {round_sig(self.branching_ratios[channel](self.mA))} 91 {self.decay_PDG_codes[channel]}              ! A' -> {channel}""")
        file.write(f"""
999999:tau0 = {round_sig(self.A_ctau)*10.0}                            ! proper lifetime, in mm
""")
      
      file.write("""
! eta decay chains""")  
      for channel in self.eta_branching_ratios.keys():
        if round_sig(self.eta_branching_ratios[channel])>0.0:
          file.write(f"""
4900221:addChannel = 1 {round_sig(self.eta_branching_ratios[channel])} 91 {self.eta_decay_PDG_codes[channel]}              ! eta -> {channel}""")
      file.write(f"""
4900221:tau0 = {round_sig(self.eta_ctau)*10.0}                            ! proper lifetime, in mm
""")
      
      file.write("""
! pion decay chains""")  
      if(self.m_pi3>2.0*self.mA):
        file.write(f"""
4900111:addChannel = 1 1.0 91 999999 999999              ! pi3 -> A'A'""")
        file.write(f"""
4900111:tau0 = {round_sig(self.pi3_ctau)*10.0}                            ! proper lifetime, in mm
""")
      else:
        file.write("""
4900111:onMode = 0""")
      file.write("""
4900211:onMode = 0""")  
      
      file.close()
    
      
      