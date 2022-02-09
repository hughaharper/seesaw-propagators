import numpy as np

SPMYR = 3.15576e13
GBAR = 9.81
PI = np.pi
PI2 = PI**2
RT = 8.3144621

class Litho(object):
    
    def __init__(self, **kwargs):
        '''
        Set Litho Defaults
        '''
        # param. for pl8 cooling
        self.ts = 0.0 # surface temp
        self.tm = 1350.0 # mantle temp
        self.diff = 8.0e-7 # thermal diffusivity
        self.cp = 1172.0 # heat capac.
        self.dp = 1.25e5 # plate thickness
        self.usp = 0.04 # FULL spr. rate (m / yr)
        
        # param. for thermal subsidence, sed. loading
        self.dref = 2600.0 # ridge crest depth
        self.alph = 3.2e-5 # vol. expansion coeff
        self.rw = 1025 # water density
        self.rs = 2300 # sed density
        self.rc = 2900 # crustal density
        self.rm = 3300 # mantle density
        # self.dr = # ridge axial depth
        self.dw = 4100 # mean sflr depth
        self.dc = 6000 # crustal thickness
        self.ds = 500 # sed thickness
        
        # param. for brittle YSE
        self.byerlpt = 0.786 # tension, byerlee coef for pressure, stress > 529.9 MPa
        self.byergpt = 0.679 # tension, byerlee coef for pressure, stress < 529.9 MPa
        self.byergst = 5.67e7 # tension, byerlee val for cohesion, stress < 529.9 MPa
        self.byergpc = 3.68 # compression, byerlee coef for pressure, stress < 113.2 MPa
        self.byerlpc = 2.12 # compression, byerlee coef for pressure, stress > 113.2 MPa
        self.byerlsc = 1.766e8 # compression, byerlee val for cohesion, stress > 113.2 MPa
        self.byerlpsh = 0.85 # shear, byerlee coef for pressure term, stress < 200 MPa
        self.byergcsh = 5.0e7 # shear, byerlee val for cohesion, stress > 200 MPa
        self.byergpsh = 0.6 # shear, byerlee coef for pressure, stress > 200 MPa
        self.phyd = 0.0 # pore pressure level
        
        # elastic params
        self.young = 6.5e10 # young's modulus
        self.pois = 0.25 # poisson's ratio
        self.telas = 600 # temp at base of elastic layer
        
        # param. for ductile flow law
        self.crust_flow = 0 # crustal flow law. 0 = dry olivine, 1 = wet olivine, 2 = diabase
        self.eps1 = 1.e-14 # strain rate
        self.str_exp = 3.5 # stress exponent
        self.str_pow = 2.4e-16 # stress amplitude factor for power law
        self.str_dor = 8.5e9 # stress const for Dorn law
        self.sren_dor = 5.7e11 # strain rate factor for Dorn law
        self.qp = 5.40e5 # activation energy for power law
        self.qd = 5.49e5 # activation energy for Dorn law
        
        # depth-related params
        self.z = np.linspace(0,2.5e4,1000,endpoint=False) # z coordinates
        
        # set any kwargs
        for key, value in kwargs.items():
            # should verify that the key is an attribute
            self.__setattr__(key,value)
        return
    
    def print_vals(self):
        '''
        Print litho parameter values
        '''
        print("Litho params: \n")
        return
    
    def get_depth_sflr(self, age):
        '''
        Compute seafloor depth given a plate age
        See T&S (3rd ed.) eq 4.211
        '''
        tage = age*SPMYR
        pref = self.rm * self.alph * self.dp \
            * (self.tm - self.ts)/(self.rm - self.rw)
        ii = np.arange(0,50)
        jj = (2*ii + 1)*(2*ii + 1)
        argex = (self.diff*jj*PI2*tage)/(self.dp**2)
        term = np.exp(-1.0*argex) / jj
        tsum = np.sum(term)
        depth = self.dref + pref*(0.5 - 4.0*tsum/PI2)
        return depth
    
    def get_obp(self, dsflr):
        '''
        Compute overburden pressure profile given
        depth of the seafloor
        '''
        # water column overburden
        pw = dsflr*self.rw*GBAR
        press = (self.z*self.rc*GBAR) - \
            (self.phyd*self.rw*GBAR)*(self.z + dsflr)
        return pw + press
    
    def get_temperature(self,age,model=1):
        '''
        Compute vertical temperature profile given
        age and a cooling model.
        0 = half space,
        1 = Plate, See T&S (3rd ed.) eq 4.130
        2 = Sleep (1975)
        '''
        tage = age*SPMYR
        if model == 1:
            jj = np.arange(1,51)
            argex = self.diff*jj*jj*PI2*tage/(self.dp**2)
            argsin = self.z[:,np.newaxis]*jj*PI/self.dp # broadcasts z and jj
            term = np.exp(-1.0*argex)*np.sin(argsin)/jj
            tsum = np.sum(term,axis=1)
            temp = self.ts + (self.tm - self.ts)*((self.z/self.dp) + (2.0*tsum/PI))
        elif model == 2:
            tu = self.usp/365/24/60/60 # m/yr to m/s
            x = tu*tage
            z_seg = 33.0e3
            latent_h = 1.028e9
            l_adiab = 1.0e-3
            d_adiab = 0.3e-3
            melt_grad = 3.0e-3
            conduct = 2.5104
            rhoc = conduct/self.diff
            gamma = 1 - (d_adiab*self.dp)/self.tm
            T_c = self.tm*(gamma - ((gamma*self.dc)/self.dp) + (self.dc/self.dp))
            T_seg = self.tm*(gamma - ((gamma*z_seg)/self.dp) + (z_seg/self.dp))
            R_p = (2*self.diff*PI)/(tu*self.dp)
            jj = np.arange(1,502)
            a_m = (tu/(2*self.diff))*(1 - np.sqrt(1 + (R_p*R_p*jj*jj)))
            B_m1 = np.cos((jj*PI*z_seg)/self.dp)*((1 - (z_seg/self.dp))*
                                                  self.tm*gamma - T_seg +
                                                  self.tm*(z_seg/self.dp))
            B_m2 = (1/(jj*PI))*np.sin((jj*PI*z_seg)/self.dp) * \
                (self.tm*gamma + melt_grad*self.dp - self.tm)
            B_m3 = np.cos((jj*PI*self.dc)/self.dp)*(T_seg - (z_seg - self.dc)*
                                                    melt_grad - (latent_h/rhoc) -
                                                    T_c)
            B_m4 = (1/(jj*PI))*np.sin((jj*PI*self.dc)/self.dp)* \
                (l_adiab*self.dp - melt_grad*self.dp) + (latent_h/rhoc) + \
                T_c - l_adiab*self.dc
            B_m = ((2*tu*rhoc)/(jj*PI))*(B_m1 + B_m2 + B_m3 + B_m4)
            A_m = 2/(1 + np.sqrt(1 + (R_p*R_p*jj*jj)))
            T_homo = A_m*B_m*np.sin((self.z[:,np.newaxis]*jj*PI)/self.dp)*np.exp(a_m*x)
            temp = (1/(tu*rhoc))*np.sum(T_homo,axis=1) + ((self.z*self.tm)/self.dp)
        return temp
    
    def get_ductile(self,temp):
        '''
        Get vertical ductile strength profile
        '''
        tempk = temp + 273.15
        # Dorn Law
        # ductstr = self.str_dor*(1-np.sqrt((tempk*RT/self.qd)*np.log(self.sren_dor/self.eps1)))
        # Power Law
        ductstr = np.empty_like(temp)
        if self.crust_flow != 0:
            if self.crust_flow == 1:
                self.str_exp = 3.0
                self.str_pow = 1.9e-15
                self.qp = 4.2e5            
            elif self.crust_flow == 2:
                self.str_exp = 4.7
                self.str_pow = 5.0e-28
                self.qp = 4.82e5
                
            idx = (self.z <= self.dc)
            ductstr[idx] = ((self.eps1/self.str_pow)*
                            np.exp(self.qp/(tempk[idx]*RT)))**(1.0/self.str_exp)
            self.str_exp = 3.5
            self.str_pow = 2.4e-16
            self.qp = 5.4e5
            idx = (self.z > self.dc)
            ductstr[idx] = ((self.eps1/self.str_pow)*
                            np.exp(self.qp/(tempk[idx]*RT)))**(1.0/self.str_exp)
        else:
            ductstr = ((self.eps1/self.str_pow)*np.exp(self.qp/(tempk*RT)))**(1.0/self.str_exp)
        return ductstr
    
    def get_byerlee(self,pressure,regime='c'):
        '''
        Get vertical brittle strength profile
        for a specified stress regime
        See Mueller & Phillips, 1995 appendix
        '''
        byerstr = np.empty_like(pressure)
        if 'c' in regime:
            # pressure vals <= 1.132e8
            idx = (pressure <= 1.132e8)
            byerstr[idx] = -1.0*(pressure[idx]*self.byergpc)
            # pressure vals > 1.132e8
            idx = (pressure > 1.132e8)
            byerstr[idx] = -1.0*(self.byerlsc + pressure[idx]*self.byerlpc)
            return byerstr
        if 't' in regime:
            # pressure vals <= 5.299e8
            idx = (pressure <= 5.299e8)
            byerstr[idx] = pressure[idx]*self.byerlpt
            # pressure vals > 5.299e8
            idx = (pressure > 5.299e8)
            byerstr[idx] = self.byergst + pressure[idx]*self.byergpt
            return byerstr
        if 's' in regime:
            # pressure vals <= 2.0e8
            idx = (pressure <= 2.0e8)
            byerstr[idx] = pressure[idx]*self.byerlpsh
            # pressure vals > 2.0e8
            idx = (pressure > 2.0e8)
            byerstr[idx] = self.byergcsh + pressure[idx]*self.byergpsh
            return byerstr
    
    def get_yse(self, age, thermal=1, **kwargs):
        '''
        Compute a yield strength envelope
        given a plate age
        ''' 
        dsf = self.get_depth_sflr(age)
        obp = self.get_obp(dsf)
        temp = self.get_temperature(age,model=thermal)
        dustr = self.get_ductile(temp)
        bystrC = self.get_byerlee(obp,'c') # compression
        bystrT = self.get_byerlee(obp,'t') # tension
        bystrS = self.get_byerlee(obp,'s') # shear
        yseC = np.fmax(bystrC,-1*dustr)
        yseT = np.fmin(bystrT,dustr)
        yseS = np.fmin(bystrS,dustr)
        return yseC, yseT, yseS