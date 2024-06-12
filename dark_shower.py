#################################################################################################
## Package to generate self-consistent configuration cards for the pythia 8 hidden valley module.
## Based on arXiv 2103.01238. Please cite this paper when using this package, in addition to the theory
## papers included in the "citation" string of each model.
##
## Writen by Simon Knapen (smknapen@lbl.gov)
## on 02/02/2021 
#################################################################################################
import os 
import sys
import glob
import math
import numpy as np
from scipy.interpolate import interp1d
from scipy.special import exp1
from datetime import date

data_dir = os.path.dirname(__file__)+"/data/"

GeVtimescm=1.98e-14
TeV_in_GeV=1e3
eV_in_GeV=1e-9

# the parameter "m" is the mass of the dark meson which decays to the SM: \tilde omega for the vector portal, \tilde \eta otherwise.

class dark_shower():
    def __init__(self,portal=""): # if no portal is specified, use photon portal
      self.portal=portal # which portal is considered
      
      if(not (self.portal=="photon" or self.portal=="vector"  or self.portal=="gluon"  or self.portal=="higgs"  or self.portal=="darkphoton" )):
        print("Error: portal flag should be one of the following: photon, vector, gluon, higgs or darkphoton")
      
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
      
      # fragmentation probabilities, see appendix A
      self.spin1_fraction=(lambda xi_Lambda,xi_omega: 3.0*exp1(xi_omega**2/xi_Lambda**2)/(3.0*exp1(xi_omega**2/xi_Lambda**2)+exp1(1.0/xi_Lambda**2)))
      
      #########################################
      # Model specific variables and functions
      #########################################
      if self.portal=="photon":
        # Variable and functions related to IR physics
        # ---------------------------------------------
        self.Gamma=(lambda m,f: self.alpha_EM**2/(4.0*np.pi)*m**3/f**2)
        # ctau(m,f) in cm, with m and f in GeV
        self.ctau=(lambda m,f: GeVtimescm/self.Gamma(m,f))
       
        # Variables related to UV completion
        # -----------------------------------
        # mass of heavy charged particles
        self.M_chi=700. 
        # yukawas
        self.y_chi=3.0
        self.y_psi=3.0
        self.NQ=1.0
        # minimum decay constant 
        self.f_min=(lambda m,xi_Lambda=1.0,m_a=-1.0: self.M_chi*np.max([10.0,2.0*m])**2/(2.0*self.y_chi*self.y_psi*self.NQ*m**2*xi_Lambda**2) if m_a<0.0 else self.M_chi*m_a**2/(2.0*self.y_chi*self.y_psi*self.NQ*m**2*xi_Lambda**2))
        # minimum ctau for UV complete model
        self.ctau_min=(lambda m,xi_Lambda=1.0,m_a=-1.0: self.ctau(m,self.f_min(m,xi_Lambda=xi_Lambda,m_a=m_a)))
        
        # Branching ratios
        # ------------------
        # PDG codes of decay products
        self.decay_PDG_codes={"gamma gamma":"22 22"}
        # Branching ratios of decays channels
        self.branching_ratios={"gamma gamma":(lambda m:1)}
        
        
        # Citation String
        #-----------------
        self.citation_string=""
        
      if self.portal=="gluon":
        
        # Variables and functions related to IR physics
        # ---------------------------------------------
        # scale at which we switch to a->gg instead of exclusive channels
        self.pert_QCD=1.84
        # taken from Fig S1 of 1811.03474
        [gluon_width_x,gluon_width_y]=np.loadtxt(data_dir+"gluon_portal/alp_lifetime_soreq.txt").T
        # convert normalization from 1811.03474
        gluon_ctau_y=np.array(list(map(lambda Ga: GeVtimescm/(Ga*eV_in_GeV)*(32.0*np.pi**2/TeV_in_GeV)**2,gluon_width_y))) 
        gluon_ctau_interpol=interp1d(gluon_width_x,np.log10(gluon_ctau_y), kind='linear')
        # ctau(m,f) in cm, with m and f in GeV
        self.ctau=(lambda m,f: 10**gluon_ctau_interpol(m)*f**2 if m<self.pert_QCD else GeVtimescm/(self.alphas(m)**2/(32.0*np.pi**3)*m**3/(f**2)*(1+(97.0/4.0-7.0/6.0*self.Nf_L(m))*self.alphas(m)/(np.pi))))
                   
        # Variables related to UV completion
        # -----------------------------------
        # UV parameters
        self.M_Q=2000. 
        self.y_Q=3.0
        self.y_psi=3.0
        self.N_Q=1.0
        # minimum decay constant
        self.f_min=(lambda m,xi_Lambda=1.0,m_a=400.0: m_a**2*self.M_Q/(self.N_Q*self.y_Q*self.y_psi*xi_Lambda**2*m**2)) 
        # minimum ctau for UV complete model
        self.ctau_min=(lambda m,xi_Lambda=1.0,m_a=400.0: self.ctau(m,self.f_min(m,xi_Lambda=xi_Lambda,m_a=m_a)))
      
      
        # Branching ratios
        # ------------------
        # PDG codes of decay products
        self.decay_PDG_codes={"gamma gamma":"22 22",\
                              "pi+ pi- pi0":"211 -211 111",\
                              "pi0 pi0 pi0":"111 111 111",\
                              "pi+ pi- gamma":"211 -211 22",\
                              "eta pi+ pi-":"221 211 -211",\
                              "eta pi0 pi0":"221 111 111",\
                              "etap pi+ pi-":"331 211 -211",\
                              "etap pi0 pi0":"331 111 111",\
                              "K+ K- pi0":"321 -321 111",\
                              "K+ KL pi-":"321 130 -211",\
                              "K- KL pi+":"-321 130 211",\
                              "K+ KS pi-":"321 310 -211",\
                              "K- KS pi+":"-321 310 211",\
                              "KS KL pi0":"310 130 111",\
                              "rho+ rho-":"213 -213",\
                              "rho0 rho0":"113 113",\
                              "omega omega":"223 223",\
                              "gluon gluon":"21 21"
                             }
        # Branching ratios of decays channels
        self.branching_ratios={}
        self.branching_ratios["gluon gluon"]=(lambda m: float(m>self.pert_QCD))
        
        datatemp=np.array(list(filter(lambda x: x[0]<self.pert_QCD*1.01,
                                      np.loadtxt(data_dir+"gluon_portal/gluonportal_2gamma.txt",skiprows=2)))).T
        self.branching_ratios["gamma gamma"]= interp1d(datatemp[0],datatemp[1],fill_value=0.0,bounds_error=False)
        
        datatemp=np.array(list(filter(lambda x: x[0]<self.pert_QCD*1.01,
                                      np.loadtxt(data_dir+"gluon_portal/gluonportal_3pi.txt",skiprows=2)))).T
        self.branching_ratios["pi+ pi- pi0"]= interp1d(datatemp[0],0.4*datatemp[1],fill_value=0.0,bounds_error=False)
        self.branching_ratios["pi0 pi0 pi0"]= interp1d(datatemp[0],0.6*datatemp[1],fill_value=0.0,bounds_error=False)
        
        datatemp=np.array(list(filter(lambda x: x[0]<self.pert_QCD*1.01,
                                      np.loadtxt(data_dir+"gluon_portal/gluonportal_2pi_gamma.txt",skiprows=2)))).T
        self.branching_ratios["pi+ pi- gamma"]= interp1d(datatemp[0],datatemp[1],fill_value=0.0,bounds_error=False)
        
        datatemp=np.array(list(filter(lambda x: x[0]<self.pert_QCD*1.01,
                                      np.loadtxt(data_dir+"gluon_portal/gluonportal_2pi_eta.txt",skiprows=2)))).T
        self.branching_ratios["eta pi+ pi-"]= interp1d(datatemp[0],0.66*datatemp[1],fill_value=0.0,bounds_error=False)
        self.branching_ratios["eta pi0 pi0"]= interp1d(datatemp[0],0.33*datatemp[1],fill_value=0.0,bounds_error=False)
        
        datatemp=np.array(list(filter(lambda x: x[0]<self.pert_QCD*1.01,
                                      np.loadtxt(data_dir+"gluon_portal/gluonportal_2pi_etap.txt",skiprows=2)))).T
        self.branching_ratios["etap pi+ pi-"]= interp1d(datatemp[0],0.66*datatemp[1],fill_value=0.0,bounds_error=False)
        self.branching_ratios["etap pi0 pi0"]= interp1d(datatemp[0],0.33*datatemp[1],fill_value=0.0,bounds_error=False)
        
        datatemp=np.array(list(filter(lambda x: x[0]<self.pert_QCD*1.01,
                                      np.loadtxt(data_dir+"gluon_portal/gluonportal_2K_pi.txt",skiprows=2)))).T
        self.branching_ratios["K+ K- pi0"]= interp1d(datatemp[0],0.16*datatemp[1],fill_value=0.0,bounds_error=False)
        self.branching_ratios["K+ KL pi-"]= interp1d(datatemp[0],0.16*datatemp[1],fill_value=0.0,bounds_error=False)
        self.branching_ratios["K- KL pi+"]= interp1d(datatemp[0],0.16*datatemp[1],fill_value=0.0,bounds_error=False)
        self.branching_ratios["K+ KS pi-"]= interp1d(datatemp[0],0.16*datatemp[1],fill_value=0.0,bounds_error=False)
        self.branching_ratios["K- KS pi+"]= interp1d(datatemp[0],0.16*datatemp[1],fill_value=0.0,bounds_error=False)
        self.branching_ratios["KS KL pi0"]= interp1d(datatemp[0],0.16*datatemp[1],fill_value=0.0,bounds_error=False)

        datatemp=np.array(list(filter(lambda x: x[0]<self.pert_QCD*1.01,
                                      np.loadtxt(data_dir+"gluon_portal/gluonportal_2rho.txt",skiprows=2)))).T
        self.branching_ratios["rho+ rho-"]= interp1d(datatemp[0],0.66*datatemp[1],fill_value=0.0,bounds_error=False)
        self.branching_ratios["rho0 rho0"]= interp1d(datatemp[0],0.33*datatemp[1],fill_value=0.0,bounds_error=False)
        
        datatemp=np.array(list(filter(lambda x: x[0]<self.pert_QCD*1.01,
                                      np.loadtxt(data_dir+"gluon_portal/gluonportal_2omega.txt",skiprows=2)))).T
        self.branching_ratios["omega omega"]= interp1d(datatemp[0],datatemp[1],fill_value=0.0,bounds_error=False)

        
        # Citation String
        #-----------------
        self.citation_string="""
        @article{Aloni:2018vki, 
                author = \"Aloni, Daniel and Soreq, Yotam and Williams, Mike\",
                title = \"{Coupling QCD-Scale Axionlike Particles to Gluons}\",
                eprint = \"1811.03474\",
                archivePrefix = \"arXiv\",
                primaryClass = \"hep-ph\",
                reportNumber = \"CERN-TH-2018-237, MIT-CTP/5080, MIT-CTP-5080\",
                doi = \"10.1103/PhysRevLett.123.031803\",
                journal = \"Phys. Rev. Lett.\",
                volume = \"123\",
                number = \"3\"
                pages = \"031803\",
                year = \"2019\"}
        @article{Chetyrkin:1998mw,
                author = \"Chetyrkin, K. G. and Kniehl, Bernd A. and Steinhauser, M. and Bardeen, William A.\",
                title = \"{Effective QCD interactions of CP odd Higgs bosons at three loops}\",
                eprint = \"hep-ph/9807241\",
                archivePrefix = \"arXiv\",
                reportNumber = \"TTP-98-21, NYU-TH-98-04-02, MPI-PHT-98-32, FERMILAB-PUB-98-126-T\",
                doi = \"10.1016/S0550-3213(98)00594-X\",
                journal = \"Nucl. Phys. B\",
                volume = \"535\",
                pages = \"3--18\",
                year = \"1998\"
                }        
                """
        
        
        
      if self.portal=="darkphoton" or self.portal=="vector":
  
        print("darkphoton/vector portal Warning: only most important branching ratios were included. For full list, see 1505.07459")
        
        # Branching ratios
        # ------------------
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
        self.total_width=(lambda m, eps: eps**2*self.partial_width["e+ e-"](m)/self.branching_interpol["e+ e-"](m) if m<=self.pert_QCD else eps**2*sum([self.partial_width[channel](m) for channel in self.partial_width.keys()]))
        
        # branching ratios
        self.branching_ratios={}
        for channel in self.partial_width.keys():
          self.branching_ratios[channel]=(lambda m, channel=channel: self.partial_width[channel](m)/self.total_width(m,1.0))
        
        self.branching_ratios["pi+ pi-"]=(lambda m: 0.0 if (m>self.pert_QCD or m<2.0*self.m_pi) else self.branching_interpol["pi+ pi-"](m))
        self.branching_ratios["K+ K-"]=(lambda m: 0.0 if (m>self.pert_QCD or m<2.0*self.m_K)else self.branching_interpol["K+ K-"](m))
        self.branching_ratios["pi+ pi- pi0"]=(lambda m: 0.0 if (m>self.pert_QCD or m<3.0*self.m_pi)else self.branching_interpol["pi+ pi- pi0"](m))
        self.branching_ratios["pi+ pi- pi+ pi-"]=(lambda m: 0.0 if (m>self.pert_QCD or m<4.0*self.m_pi) else self.branching_interpol["pi+ pi- pi+ pi-"](m))
        self.branching_ratios["pi+ pi- pi0 pi0"]=(lambda m: 0.0 if (m>self.pert_QCD or  m<2.0*self.m_pi) else self.branching_interpol["pi+ pi- pi0 pi0"](m))
              
        # Variable and functions related to IR physics
        # ---------------------------------------------
        # ctau(m,f) in cm, with m and f in GeV
        self.ctau=(lambda m,eps: GeVtimescm/self.total_width(m,eps))
        
        
        # Citation String
        #-----------------
        self.citation_string="""@article{Buschmann:2015awa,
        author         = \"Buschmann, Malte and Kopp, Joachim and Liu, Jia and
                        Machado, Pedro A. N.\",
        title          = \"{Lepton Jets from Radiating Dark Matter}\",
        journal        = \"JHEP\",
        volume         = \"07\",
        year           = \"2015\",
        pages          = \"045\",
        doi            = \"10.1007/JHEP07(2015)045\",
        eprint         = \"1505.07459\",
        archivePrefix  = \"arXiv\",
        primaryClass   = \"hep-ph\",
        reportNumber   = \"MITP-15-036, IFT-UAM-CSIC-15-047, FTUAM-15-13\",
        SLACcitation   = \"%%CITATION = ARXIV:1505.07459;%%\"
        }"""


        
        if self.portal=="darkphoton":
          # Variables related IR constraints
          # -----------------------------------
          self.eps_max=3.0e-4
          # minimum ctau subject to IR constraints
          self.ctau_min=(lambda m,xi_Lambda=1.0,xi_Ap=0.4: self.ctau(xi_Ap*m,self.eps_max) if xi_Ap <0.5 else "m_{A'} must be < 1/2 dark meson mass")
        
        
        if self.portal=="vector":
          EWPTdata=np.loadtxt(data_dir+"vector_portal/EWPT.txt",skiprows=1).T
          EWPT_interpol=interp1d(EWPTdata[0],EWPTdata[1],fill_value="extrapolate",bounds_error=False,copy=True)
          self.eps_EWPT=(lambda m: np.min([EWPT_interpol(m),1.0]))
          
          # Variables related to UV completion
          # -----------------------------------
          self.g=1.0
          # m_a refers to the mass of the heavy mediator, in this case the dark photon 
          self.eps_max=(lambda m,m_a=20.0,eps=-1.0: self.g*0.3*self.eps_EWPT(m_a)*m**2/m_a**2 if eps<0.0 else self.g*0.3*eps*(xi_Lambda*m)**2/m_a**2) 
          # minimum ctau for UV complete model. xi_Lambda not used, just included for uniformity
          self.ctau_min=(lambda m,xi_Lambda=1.0,m_a=20.0,eps=-1.0: self.ctau(m,self.eps_max(m,m_a=m_a,eps=eps)))
          
          self.citation_string_EWPT="""@article{Curtin:2014cca,
    author = \"Curtin, David and Essig, Rouven and Gori, Stefania and Shelton, Jessie\",
    title = \"{Illuminating Dark Photons with High-Energy Colliders}\",
    eprint = \"1412.0018\",
    archivePrefix = \"arXiv\",
    primaryClass = \"hep-ph\",
    reportNumber = \"YITP-SB-14-49\",
    doi = \"10.1007/JHEP02(2015)157\",
    journal = \"JHEP\",
    volume = \"02\",
    pages = \"157\",
    year = \"2015\"
}"""
          
          
      if self.portal=="higgs":
        print("Higgs portal warning: 4pi, 2eta and 2rho branching ratios are unknown and were not included. As a result, branching ratio to Kaons in the 1 GeV to 2 GeV mass regime is likely slightly inflated. See 1809.01876")
        
        # Branching ratios
        # -----------------
        
        # PDG codes of decay products
        self.decay_PDG_codes={
                              "pi+ pi-":"211 -211",\
                              "pi0 pi0":"111 111",\
                              "K+ K-":"321 -321",\
                              "K0 K0bar":"311 -311",\
                              "gluon gluon":"21 21",\
                              "s sbar":"3 -3",\
                              "c cbar":"4 -4",\
                              "b bbar":"5 -5",\
                              "e+ e-":"11 -11",\
                              "mu+ mu-":"13 -13",\
                              "tau+ tau-":"15 -15",\
                              "gamma gamma":"22 22"
                             }
        
        # scale at which we switch to quarks instead of exclusive channels 
        self.pert_QCD=2.0
        
        self.loop_f=(lambda x: np.arcsin(np.sqrt(1.0/np.real(x)))**2 if x>=1 else -0.25*(np.log((1.0+np.sqrt(1.0-np.real(x)))/(1.0-np.sqrt(1.0-np.real(x))))-1.0j*np.pi)**2)
        self.loop_g=(lambda x:x*(1.0+(1.0-x)*self.loop_f(x)))
        
        self.partial_width={}
        
        # formulas taken from 9511344
        # running of charm and bottom mass
        self.c_fun=(lambda x, mu: (25.0/6.0*x)**(12.0/25.0)*(1.0+1.014*x+1.389*x**2) if mu<self.m_b else\
       (23.0/6.0*x)**(12.0/23.0)*(1.0+1.175*x+1.501*x**2))
        self.m_c_run=(lambda mu: self.m_c*self.c_fun(self.alphas(mu)/np.pi,mu)/self.c_fun(self.alphas(self.m_c)/np.pi,self.m_c))
        self.m_b_run=(lambda mu: self.m_b*self.c_fun(self.alphas(mu)/np.pi,mu)/self.c_fun(self.alphas(self.m_b)/np.pi,self.m_b))
        
        # NlO QCD corrections to ccbar and bbar widths (from 9511344)
        self.DeltaQCD=(lambda m: 1.0+5.67*self.alphas(m)/np.pi+(35.94-1.36*self.Nf_L(m))*(self.alphas(m)/np.pi)**2)
        self.Deltat=(lambda m, mQ: (self.alphas(m)/np.pi)**2*(1.57-2.0/3.0*np.log(m**2/self.m_t**2)\
                                                 +1.0/9.0*(np.log(mQ**2/m**2))**2))
      
        self.partial_width["c cbar"]=(lambda m: 0.0 if m<=self.pert_QCD else \
                         self.Nc*self.G_F*m/(4*np.sqrt(2.0)*np.pi)*self.m_c_run(m)**2\
                         *(self.DeltaQCD(m)+self.Deltat(m,self.m_c_run(m)))*(1-np.min([1.0,4*self.m_D**2/m**2]))**1.5)
        
        self.partial_width["b bbar"]=(lambda m: 0.0 if m<=self.pert_QCD else \
                         self.Nc*self.G_F*m/(4*np.sqrt(2.0)*np.pi)*self.m_b_run(m)**2\
                         *(self.DeltaQCD(m)+self.Deltat(m,self.m_b_run(m)))*(1-np.min([1.0,4*self.m_B**2/m**2]))**1.5)
        
        # strange is only LO, but always subdominant anyways
        self.partial_width["s sbar"]=(lambda m: 0.0 if m<=self.pert_QCD else self.Nc*self.G_F*m/(4*np.sqrt(2.0)*np.pi)*self.m_s**2*(1-np.min([1.0,4*self.m_K**2/m**2]))**1.5)
        
        # NLO QCD corrections to gluon width  
        self.Gamma_glue_LO=(lambda m: self.G_F*self.alphas(m)**2*m**3/(16.0*np.sqrt(2.0)*np.pi**3)*\
                              np.abs(self.loop_g((4*self.m_t**2)/m**2) + self.loop_g((4*self.m_b**2)/m**2)\
                                     +self.loop_g((4*self.m_c**2)/m**2))**2)
        self.glue_NLO=(lambda m: self.Nf_H(m)*(95.0/4.0-7.0/6.0*self.Nf_L(m))*self.alphas(m)/np.pi)

        self.partial_width["gluon gluon"]=(lambda m: 0.0 if m<=self.pert_QCD else self.Gamma_glue_LO(m)*(1.0+self.glue_NLO(m)))
        
        self.partial_width["e+ e-"] = (lambda m: self.G_F*m/(4*np.sqrt(2.0)*np.pi)*self.m_e**2*(1-np.min([1.0,4*self.m_e**2/m**2]))**1.5)
        
        self.partial_width["mu+ mu-"]=(lambda m: self.G_F*m/(4*np.sqrt(2.0)*np.pi)*self.m_mu**2*(1-np.min([1.0,4*self.m_mu**2/m**2]))**1.5)
        
        self.partial_width["tau+ tau-"]=(lambda m: self.G_F*m/(4*np.sqrt(2.0)*np.pi)*self.m_tau**2*(1-np.min([1.0,4*self.m_tau**2/m**2]))**1.5)
        
        # photon width
        A_f=(lambda t: 2.0*t*(1.0+(1.0-t)*self.loop_f(t)))
        A_W=(lambda t: -1.0*(2.0+3.0*t+3.0*t*(2.0-t)*self.loop_f(t)))

        self.partial_width["gamma gamma"]=(lambda m: self.G_F*self.alpha_EM**2*m**3/\
                              (128.0*np.sqrt(2.0)*np.pi**3)*\
                              np.abs(self.Nc*(1.0/3.0)**2*A_f(4.0*self.m_s**2/m**2)+\
                                     self.Nc*(1.0/3.0)**2*A_f(4.0*self.m_b**2/m**2)+\
                                     self.Nc*(2.0/3.0)**2*A_f(4.0*self.m_c**2/m**2)+\
                                     self.Nc*(2.0/3.0)**2*A_f(4.0*self.m_t**2/m**2)+\
                                     (1.0)**2*A_f(4.0*self.m_e**2/m**2)+\
                                     (1.0)**2*A_f(4.0*self.m_mu**2/m**2)+\
                                     (1.0)**2*A_f(4.0*self.m_tau**2/m**2)+\
                                     A_W(4.0*self.m_W**2/m**2))**2)

        
        
        datatemp_pi=np.loadtxt(data_dir+"higgs_portal/higgsportal_2pi.txt",skiprows=2).T
        datatemp_K=np.loadtxt(data_dir+"higgs_portal/higgsportal_2K.txt",skiprows=2).T
        interpol_pi=interp1d(datatemp_pi[0],datatemp_pi[1],fill_value=0.0,bounds_error=False,copy=True)
        interpol_K=interp1d(datatemp_K[0],datatemp_K[1],fill_value=0.0,bounds_error=False,copy=True)
        
        # 66% branching ratio to pi+ pi- vs pi0 pi0
        self.partial_width["pi0 pi0"]=(lambda m: 0.0 if m>self.pert_QCD else 0.333*interpol_pi(m))
        self.partial_width["pi+ pi-"]=(lambda m: 0.0 if m>self.pert_QCD else 0.667*interpol_pi(m))
        
        # 50% branching ratio to pi+ pi- vs pi0 pi0
        self.partial_width["K+ K-"]=(lambda m: 0.0 if m>self.pert_QCD else 0.5*interpol_K(m))
        self.partial_width["K0 K0bar"]=(lambda m: 0.0 if m>self.pert_QCD else 0.5*interpol_K(m))
        
        # (30) of 1809.01876, parametrizes unknown hadronic channels
        self.partial_width["other"]=(lambda m: 0.0 if m>self.pert_QCD else 3.34e-8* m**3*np.max([0.0,1-16.0*self.m_pi**2/m**2])**0.5)
        
        # total width
        self.total_width=(lambda m: sum([self.partial_width[channel](m) for channel in self.partial_width.keys()]))
        
        # branching ratios. Substract the unaccounted width, so that the branching ratios add up to 1.
        # this slightly inflates 2pi and 2K branching ratios
        self.branching_ratios={}
        for channel in self.partial_width.keys():
          self.branching_ratios[channel]=(lambda m, channel=channel: self.partial_width[channel](m)/(self.total_width(m)-self.partial_width["other"](m)))
        
        # lifetime
        self.ctau=(lambda m,stheta:1.0/self.total_width(m)/stheta**2*GeVtimescm)
                   
        # maximum mixing without tuning
        self.stheta_max=(lambda m,xi_Lambda=1.0: xi_Lambda**3*m**3/(self.vev*self.m_h**2))
        
        # minimum ctau without tuning
        self.ctau_min=(lambda m,xi_Lambda=1.0: self.ctau(m,self.stheta_max(m)))
        
        # Citation String
        #-----------------
        self.citation_string="""
        @article{Winkler:2018qyg,
        author = "Winkler, Martin Wolfgang",
        title = "{Decay and detection of a light scalar boson mixing with the Higgs boson}",
        eprint = "1809.01876",
        archivePrefix = "arXiv",
        primaryClass = "hep-ph",
        reportNumber = "NORDITA-2018-087",
        doi = "10.1103/PhysRevD.99.015018",
        journal = "Phys. Rev. D",
        volume = "99",
        number = "1",
        pages = "015018",
        year = "2019"
        }
        
        @article{Djouadi:1995gt,
        author = "Djouadi, A. and Spira, M. and Zerwas, P. M.",
        title = "{QCD corrections to hadronic Higgs decays}",
        eprint = "hep-ph/9511344",
        archivePrefix = "arXiv",
        reportNumber = "DESY-95-210, KA-TP-8-95",
        doi = "10.1007/s002880050120",
        journal = "Z. Phys. C",
        volume = "70",
        pages = "427--434",
        year = "1996"
        }
        """
        
        
    #######################
    # Universal functions
    #######################
    
    # calculates the weight for the decay to occur within the fiducial acceptance defined by Lxy_min and Lxy_max, while saturating 
    # the minimum lifetime bound.
    def weight(self,m,boost,eta,m_a=-1.0,xi_Lambda=1.0,Lxy_min=0.0,Lxy_max=100.0):
      # uses the default mediator mass for each model. Flag is not defined for higgs and darkphoton models
      if m_a<0.0 or self.portal=="darkphoton" or self.portal=="higgs": 
        return np.exp(-Lxy_min/np.abs(np.cosh(eta))/(boost*self.ctau_min(m,xi_Lambda=xi_Lambda)))\
              -np.exp(-Lxy_max/np.abs(np.cosh(eta))/(boost*self.ctau_min(m,xi_Lambda=xi_Lambda)))
      else: # specify the mediator mass
        return np.exp(-Lxy_min/np.abs(np.cosh(eta))/(boost*self.ctau_min(m,xi_Lambda=xi_Lambda,m_a=m_a)))\
              -np.exp(-Lxy_max/np.abs(np.cosh(eta))/(boost*self.ctau_min(m,xi_Lambda=xi_Lambda,m_a=m_a)))
    
    # function to write pythia card. Allow user to specify ctau, but print warning if lifetime is not realistic.
    # the parameter "m" is the mass of the dark meson which decays to the SM: \tilde omega for the vector portal, \tilde \eta otherwise.
    # xi_Lambda is the ratio of the confinement scale over the mass of the unstable meson.
    # xi_omega is m_omega/m_eta
    # xi_Ap is the m_A'/m_eta, only applicable for dark photon portal.
    def pythia_card(self,filename,m,xi_omega=1.0,xi_Lambda=1.0,xi_Ap=0.4,ctau=1.0,Nc=3,Nf=1,mH=125,mayDecay=True,production='ggf',userName="Simon"):
      if ctau < self.ctau_min(m,xi_Lambda) and mayDecay:
        print("Warning, requested lifetime may not be theoretically consistent.")
      
      if self.portal=="vector":
        m_eta=m/xi_omega
        m_omega=m
        Lambda=xi_Lambda*m
      else:
        m_eta=m
        m_omega=xi_omega*m
        Lambda=xi_Lambda*m
      if self.portal=="darkphoton":
        m_Ap=xi_Ap*m_eta
        
      today = date.today()
      
      file= open(filename,"w")
      file.write(
        f"""
        ##############################################################################################
        # Pythia 8 card for dark shower models, as defined in arXiv 2103.01238
        # Decay portal: {self.portal}
        # Input mass scale: {mH} GeV
        # Scalar meson mass: {m_eta} GeV
        # Vector meson mass: {m_omega} GeV
        # ctau: {ctau*10.0} mm
        # Number colors: {Nc}
        # Number flavors: {Nf}
        # Confinement scale: {Lambda} GeV
        # Generated by {userName} on {today.strftime("%B %d, %Y")}
        # lines starting with ! or # are commented out
        ##############################################################################################
        """)
      if(ctau < self.ctau_min(m,xi_Lambda) or m_eta<0.8 or (1.0<m_eta<2.0 and self.portal=="higgs")):
        file.write("""##############################################################################################
        # WARNINGS
        """)
        if(ctau < self.ctau_min(m,xi_Lambda)):
          file.write("""# * Requested lifetime may not be theoretically well motivated.
        """)
        if(m<0.8):
          file.write("""# * Dark pion mass < 2 x minimum dark quark mass allowed by pythia 8. 
        """)
        if(1.0<=m<=2.0 and self.portal=="higgs"):
          file.write("""# * 4pi, 2eta and 2rho branching ratios are unknown and were not included. 
        #   As a result, the branching ratios to Kaons and pions are likely slightly inflated. See 1809.01876 
        """)
        file.write("""##############################################################################################
        """)
      file.write(f"""  
        ! Settings used in the main program.
        Main:numberOfEvents = 10000                     ! number of events to generate
        Beams:eCM = 13000.                              ! CM energy of collision

        ! Settings related to output in init(), next() and stat().
        Init:showChangedSettings = on                   ! list changed settings
        Init:showChangedParticleData = off              ! list changed particle data
        Next:numberCount = 1000                         ! print message every n events
        Next:numberShowInfo = 0                         ! print event information n times
        Next:numberShowProcess = 0                      ! print process record n times
        Next:numberShowEvent = 0                        ! print event record n times

        ! For debugging purposes only
        ! PartonLevel:ISR = off
        ! PartonLevel:FSR = off
        ! PartonLevel:MPI = off
        ! HadronLevel:all = off

        ! Production settings
        """)
      if production=="wh" or production=="Wh" or production=="WH": 
        file.write("""
        HiggsSM:ffbar2HW = on # force leptonic W decay
        24:0:onMode=0 # bRatio=0.3213690 products=-1 2
        24:1:onMode=0 # bRatio=0.0164940 products=-1 4
        24:2:onMode=0 # bRatio=0.0165020 products=-3 2
        24:3:onMode=0 # bRatio=0.3206150 products=-3 4
        24:4:onMode=0 # bRatio=0.0000100 products=-5 2
        24:5:onMode=0 # bRatio=0.0005910 products=-5 4
        24:6:onMode=1 # bRatio=0.1081660 products=-11 12
        24:7:onMode=1 # bRatio=0.1081660 products=-13 14
        24:8:onMode=0 # bRatio=0.1080870 products=-15 16
        """)
      elif production=="zh" or production=="Zh" or production=="ZH": 
        file.write("""
        HiggsSM:ffbar2HZ = on # force leptonic Z decay
        23:0:onMode=0 # bRatio=0.1539950 products=1 -1 
        23:1:onMode=0 # bRatio=0.1194200 products=2 -2 
        23:2:onMode=0 # bRatio=0.1539840 products=3 -3 
        23:3:onMode=0 # bRatio=0.1192590 products=4 -4 
        23:4:onMode=0 # bRatio=0.1522720 products=5 -5 
        23:5:onMode=1 # bRatio=0.0335760 products=11 -11 
        23:6:onMode=0 # bRatio=0.0668060 products=12 -12 
        23:7:onMode=1 # bRatio=0.0335760 products=13 -13 
        23:8:onMode=0 # bRatio=0.0668060 products=14 -14 
        23:9:onMode=0 # bRatio=0.0335000 products=15 -15 
        23:10:onMode=0 # bRatio=0.0668060 products=16 -16 
        """)        
      elif production=="ggf" or production=="ggF" or production=="GGF" :  
        file.write("""HiggsSM:gg2H = on""")
      elif production=="vbf" or production=="VBF" :  
        file.write("""HiggsSM:ff2Hff = on""")
      elif production=="tth" or production=="ttH" or production=="TTH" :  
        file.write("""
          HiggsSM:gg2Httbar = on
          HiggsSM:qqbar2Httbar = on""") 
      else : # all 
        file.write("""HiggsSM:all = on""")
      
      file.write(f"""  
        ! decay the Higgs to two dark quark, turn off all SM branching ratios
        25:m0 ={mH}
        25:addChannel = 1 1.0 102 4900101 -4900101
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
        HiddenValley:Ngauge = {Nc}                          ! number of colors
        HiddenValley:nFlav = {Nf}                           ! number of flavors
        HiddenValley:fragment = on
        HiddenValley:FSR = on
        HiddenValley:alphaOrder = 1                      ! use running coupling
        HiddenValley:Lambda = {Lambda}                        ! dark confinement scale
        HiddenValley:pTminFSR = {1.1*Lambda}                      ! pT cut off on dark shower (IR regulator)
        HiddenValley:spinFv=0                            ! spin of bifundamentals, which are not used, but set for consistency
        HiddenValley:probVector={np.round(self.spin1_fraction(xi_Lambda,xi_omega),decimals=3)}
      
        ! dark sector meson mass spectrum
        4900101:m0 = {np.max([0.4,m_eta*0.4])}                                ! Dark Quark Mass, pythia 8 crashes if this is < 0.4 GeV
        4900111:m0 = {m_eta}                                ! Setting pi'0  Mass
        4900113:m0 = {m_omega}                                ! Setting omega'0 Mass
        """)
      if(Nf>1):
          file.write(f"""4900211:m0 = {m_eta}                                ! Setting pi'+  Mass
        4900213:m0 = {m_omega}                                ! Setting omega'+ Mass   
        """)
      if(self.portal == "darkphoton"):
          file.write(f"""999999:all = GeneralResonance void 1 0 0 {m_Ap} 0.001 0. 0. 0.        ! dark photon A'
        """) # mA'/mphi currently hardwired in 
      
      file.write("""
        ! dark meson decay chains
        """)
      
      
      # omega_0 decays, check kinematics. CHECK BRANCHING RATIOS!!!
      if(xi_omega>2.0):
        if(Nf>1):
          file.write(
          """4900113:addChannel = 1 0.5 91 4900111 4900111                       ! omega'0 -> pi'0 pi'0
        4900113:addChannel = 1 0.5 91 4900211 -4900211                      ! omega'0 -> pi'+ pi'-
        """
        )
        else:
          file.write(
          """4900113:addChannel = 1 1.0 91 4900111 4900111                       ! omega'0 -> pi'0 pi'0
        """
        ) 
      else:
        if(mayDecay and self.portal == "vector"):
          for key in self.decay_PDG_codes:
            if np.round(self.branching_ratios[key](m_omega),decimals=3)>0.0:
              file.write(f"""4900113:addChannel = 1 {np.round(self.branching_ratios[key](m_omega),decimals=3)} 91 {self.decay_PDG_codes[key]}                 ! pi0' -> {key}
        """)
          file.write(f"""4900113:tau0 = {ctau*10}                                         ! proper lifetime, in mm
        """)
        else:  
          file.write("""4900113:onMode = 0                                                ! omega'0, stable
        """)
        
      # omega+ and pi+, they only exist of Nf>1 
      if(Nf>1):
        # pi+ is always taken to be stable
        file.write("""4900211:onMode = 0                                                  ! pi'+, stable
        """)
        # if open, decay the dark rho to dark pions
        if(xi_omega>2.0): 
           
          file.write("""4900213:addChannel = 1 1.0 91 4900211 4900111                       ! omega'+ -> pi'+ pi'0
        """)
        else:
          file.write("""4900213:onMode = 0                                                ! omega'+, stable
        """)
      
      # pi0' decay
      if(self.portal == "darkphoton"):            
        file.write("""4900111:addChannel = 1 1.0 91 999999 999999                         ! pi0' -> A' A'
        """)
        if(mayDecay):
          for key in self.decay_PDG_codes:
            if np.round(self.branching_ratios[key](m_Ap),decimals=3)>0.0:
              file.write(f"""999999:addChannel = 1 {np.round(self.branching_ratios[key](m_Ap),decimals=3)} 91 {self.decay_PDG_codes[key]}              ! A' -> {key}
        """)
          file.write(f"""999999:tau0 = {ctau*10}                                          ! proper lifetime, in mm
        """)           
        else:
          file.write("""999999:onMode = 0
        """)            
      else:
        if(mayDecay and not self.portal == "vector"):
          for key in self.decay_PDG_codes:
            if np.round(self.branching_ratios[key](m_eta),decimals=3)>0.0:
              file.write(f"""4900111:addChannel = 1 {np.round(self.branching_ratios[key](m_eta),decimals=3)} 91 {self.decay_PDG_codes[key]}                 ! pi0' -> {key}
        """)
          file.write(f"""4900111:tau0 = {ctau*10}                                         ! proper lifetime, in mm
        """)       
        else:
          file.write("""4900111:onMode = 0                                                ! pi'0, stable
        """)  
      
      file.close()
    
  
  
  
    def citation(self):
      print("Branching ratios for this portal are based on:")
      print(self.citation_string)
      if self.portal=="vector":
        print("EWPT bounds were taken from:")
        print(self.citation_string_EWPT)
      
    
    
