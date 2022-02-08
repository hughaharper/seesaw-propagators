class litho(object):
    
    def __init__(self,**kwargs):
        '''
        Set Litho Defaults
        '''
        # param. for pl8 cooling
        self.ts = 0.0 # surface temp
        self.tm = 1350.0 # mantle temp
        self.diff = 8.0e-7 # thermal diffusivity
        self.cp = 1172.0 # heat capac.
        self.dp = 1.25e5 # plate thickness
        self.usp = 4 # FULL spr. rate (cm / yr)
        
        # param. for thermal subsidence, sed. loading
        self.gbar = 9.81 # grav. accel.
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
        self.byergcsh = 0.5 # shear, byerlee val for cohesion, stress > 200 MPa
        self.byergpsh = 0.6 # shear, byerlee coef for pressure, stress > 200 MPa
        self.phyd = 1.0 # pore pressure level
        
        # elastic params
        self.young = 6.5e10 # young's modulus
        self.pois = 0.25 # poisson's ratio
        self.telas = 600 # temp at base of elastic layer
        
        # param. forductile flow law
        self.eps1 = 1.e-14 # strain rate
        self.str_exp = 3.0 # stress exponent
        self.str_pow = 7.0e-14 # stress amplitude factor for power law
        self.str_dor = 8.5e9 # stress const for Dorn law
        self.sren_dor = 5.7e11 # strain rate factor for Dorn law
        self.qp = 5.20e5 # activation energy for power law
        self.qd = 5.49e5 # activation energy for Dorn law
        
        # depth-related params
        self.zmt = -1.0 # mechanical thickness based on duct. str
        self.zn = -1.0 # depth of nodal plane
        self.zy = 1000 # depth to top of elastic layer
        return
    
    def print_vals(self):
        '''
        Print litho parameter values
        '''
        print("Litho defaults: \n")
        return
    
    def get_yse(self,age):
        '''
        Compute a yield strength envelope
        given a plate age
        '''
        return