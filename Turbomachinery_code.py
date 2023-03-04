# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 17:25:00 2023

@author: Nick de Jong
"""
#Libraries
import numpy as np
import math
from scipy.optimize import fsolve
from scipy.optimize import least_squares
from scipy.optimize import minimize
import rocketcea as rc
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer
import Nozzle_turbine as NT


#Values used during tests, corresponding to Global or software input
O_F_ = 5.0
Pa_ = 1.0e5
Tf_cool_ = 500.0
dptcool_ = 1.0e5
m_ = 2.0

#Placeholder for propellant class
class Propellant:
    #Oxidizer
    Ox_name = "LOX" #Oxidizer name for rocketCEA
    Ox_composition = "O 2" #Composition of oxidizer for rocketcea
    o_dens = 1141.0 #Oxidizer density
    ocp = 14307.0 #oxidizer cp
    h_ox = -12.979 #oxidizer enthalpy
    o_lamb = 1.0e-3
    omiu=1.0e-6
   
    #Fuel
    Fuel_name = "LH2" #Fuel name for rocketCEA
    Fuel_composition = "H 2" #Composition of fuel for rocketcea
    f_dens_l = 71.0 #liquid fuel density
    f_dens_g = 1.0 #gaseous fuel density
    f_gamma = 1.4 #fuel gamma
    fcp = 14307.0 #fuel cp
    h_fuel = -9.012 # fuel enthalpy
    R_f = 4.1573 #fuel gas constant
    f_lamb = 1.0e-3
    fmiu=1.0e-6
    
    Frozen_state=0
    
    #Propellant
    gama = 1.4
    tq = 0.9 #characteristic chemical time of propellant
    MR = 3 #mixture ratio

    def __init__(self,type): #Placeholder
        if(type==1):
            ox_dens=1141.0

prop_ = Propellant(0)

#Placeholder for default class
class Default:
    #Tolerances
    pres_tol = 0.01
    toll_c_star = 0.01
    toll_F_obj = 0.01
    Max_iterations_mass_flow = 10000
    toll_P_adapted = 0.01
    Safety_factor=1.3

    #Seeds
    Pres = 1e6
    inj_vel = 15

    #Injectors
    Cd = 0.7

    #Nozzle
    Nozzle_type = 0
    MR = 0
    De_max = 2.5
    De_turbine_noz_max = 2.5
    Theta_con = 60
    Theta_conical = 15
    Theta_bell = 55
    TH_exit_bell = 3
    R_u_ratio=1

    #Turbomachinery
    cycle_type = "EX"
    Eff_t = 0.6 #Turbine efficiency
    Eff_p = 0.6 #Pump efficiency
    Eff_m = 0.95 #Mechanical efficiency between turbine and pumps
    p_to = 1.0e5 #oxidizer tank storage pressure
    ptf = 1.0e5 #Fuel tank oxidizer pressure
    Wmotor = 1.0e6 #Power of the electric motor

    #Combustion chamber
    SF = 1.0

    #Cooling
    Dr = 0.01
    A=0.0003
    T_fuel_tanks = 20
    T_ox_tanks = 60

def0 = Default()

#Main function, which calls the function for the selected cycle type
def TurboM(Default : Default, prop : Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
    # Tf_cool and dptcool are given by cooling, corresponding respectively to temperature after cooling and pressure drop in cooling
    # m corresponds to the mass flow in the nozzle and it is given by the nozzle functions, as well as O_F, corresponding to mixture ratio
    # p_a is the ambient pressure, considered a global input
    match Default.cycle_type:
        case "GG": #Gas Generator
            turbo = GG(Default, prop, O_F, p_a, Tf_cool, dptcool, m)
        case "EX": #Expander cycle
            turbo = EX(Default, prop, O_F, Tf_cool, dptcool, m)
        case "SC": #Staged combustion cycle
            turbo = SC(Default, prop, O_F, p_a, Tf_cool, dptcool, m)
        case "CB": #Coolant bleed cycle
            turbo = CB(Default, prop, O_F, p_a, Tf_cool, dptcool, m)
        case "TO": #Combustion tap off cycle
            turbo = TO(Default, prop, O_F, p_a, Tf_cool, dptcool, m) #Will not be currently implmented
        case "NO": #No turbomachinery
            dptvalve = 0.0
            dptlines = 0.0
            return [Default.p_to-dptvalve-dptlines, Default.ptf-dptvalve-dptlines]
        case "EL": #Electric motor to dirve pumps
            turbo = EL(Default, prop, O_F, Tf_cool, dptcool, m)
        case _:
            print("Cycle Not recognized")
            return
    
    turbo.results();
    #return [turbo.dptop,turbo.dptfp,turbo.pt1,turbo.pt2,turbo.ptinj,turbo.Wt,turbo.Wop,turbo.Wfp]
    return turbo.ptinj


#Function that computes the expander cycle
class EX:
    #Global input
    ptanko : float
    ptankf : float
    prop : Propellant
    O_F : float
    eff_p : float
    eff_t : float
    eff_m : float

    #Software input
    Tf_cool : float
    dptcool : float
    m : float

    #Output
    dptop : float
    dptfp : float
    pt1 : float
    pt2 : float
    ptinj : float
    Wt : float
    Wop : float
    Wfp : float

    #Auxiliary
    dptvalve = 0.0
    dptlines = 0.0

    #Initialize values
    def __init__(self, DF : Default, prop : Propellant, O_F : float, Tf_cool : float, dptcool : float, m : float):
        print("Expander cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_p = DF.Eff_p
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m

    #Obtain results, calling optimization procedure and then computing variables of interest
    def results(self):
        res = minimize(self.opt, [1.0e6], method = 'Nelder-Mead', bounds=[[0.0, 10.0e12]])
        self.dptfp = res["x"][0]
        self.dptop, self.pt1, self.pt2, self.ptinj = fsolve(self.equations,[1.0e6,1.0e7,1.0e6,1.0e6],self.dptfp)
        self.Wop = self.O_F/(self.O_F+1.0) * self.m * self.dptop / (self.eff_p*self.prop.o_dens)
        self.Wfp = 1.0/(self.O_F+1.0) * self.m * self.dptfp / (self.eff_p*self.prop.f_dens_l)
        self.Wt = 1.0/(self.O_F+1.0) * self.m * self.eff_t * self.prop.fcp * self.Tf_cool * (1.0-(self.pt2/self.pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
        print(res)
        print([self.dptop, self.dptfp, self.pt1, self.pt2, self.ptinj])
        print([self.Wop,self.Wfp,self.Wt])

    #Optimize for maximum chamber pressure
    def opt(self,dptfp):
        root = least_squares(self.equations,[1.0e6,1.0e7,1.0e6,1.0e6], args = dptfp, bounds = ((0.0,0.0,0.0,0.0),(10.0e10,10.0e10,10.0e10,10.0e10)))
        return -1.0*root["x"][3]

    #System of equations to be solved
    def equations(self,vars,dptfp):
        dptop, p1t, p2t, pinj = vars
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - pinj,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - p1t,
            p2t - self.dptlines - pinj,
            self.O_F*dptop/(self.eff_p*self.prop.o_dens) + dptfp/(self.eff_p*self.prop.f_dens_l) - self.eff_m*self.eff_t*self.prop.fcp*self.Tf_cool*(1.0-(p2t/p1t)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
            ]

    
#Function that computes staged combustion
class SC:
    #Global input
    ptanko : float
    ptankf : float
    prop : Propellant
    O_F : float
    eff_p : float
    eff_t : float
    eff_m : float

    #Software input
    Tf_cool : float
    dptcool : float
    m : float

    #Output
    dptop : float
    dptfp : float
    pt1 : float
    pt2 : float
    ptinj : float
    Wt : float
    Wop : float
    Wfp : float

    #Auxiliary
    dptvalve = 0.0
    dptlines = 0.0

    def __init__(self, DF : Default, prop : Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Staged Combustion cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_p = DF.Eff_p
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m


#Function that computes coolant bleed
class CB:
    #Global input
    ptanko : float
    ptankf : float
    df : Default
    prop : Propellant
    O_F : float
    eff_p : float
    eff_t : float
    eff_m : float
    pa : float

    #Software input
    Tf_cool : float
    dptcool : float
    m : float

    #Output
    dptop : float
    dptfp : float
    pt1 : float
    pt2 : float
    ptinj : float
    Wt : float
    Wop : float
    Wfp : float
    l : float #fraction of fuel for bleed

    #Auxiliary
    dptvalve = 0.0
    dptlines = 0.0
    m_O = 1.0
    m_F = 1.0

    def __init__(self, DF : Default, prop : Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Coolant bleed cycle selected")
        self.df = DF
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_p = DF.Eff_p
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.m_O = O_F/(O_F+1.0) * m
        self.m_F = 1.0/(O_F+1.0) * m


    def results(self):
        res = minimize(self.opt, [1.0e5,0.01], method = 'Nelder-Mead', bounds=[[self.pa*1.2,10.0e12],[1.0e-5,1.0]])
        self.pt2 = res["x"][0]
        self.l = res["x"][1]
        self.dptop, self.pt1, self.dptfp = fsolve(self.equations,[1.0e6,1.0e7,1.0e7],(self.pt2,self.l))
        self.ptinj = self.pt1
        self.Wop = self.m_O * self.dptop / (self.eff_p*self.prop.o_dens)
        self.Wfp = (1.0/(1.0-self.l)) * self.m_F * self.dptfp / (self.eff_p*self.prop.f_dens_l)
        self.Wt = (self.l/(1.0-self.l)) * self.m_F * self.eff_t * self.prop.fcp * self.Tf_cool * (1.0-(self.pt2/self.pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
        print(res)
        print([self.dptop, self.dptfp, self.pt1, self.pt2, self.l])
        print([self.Wop,self.Wfp,self.Wt])

    def opt(self,vars):
        root = least_squares(self.equations,[1.0e6,1.0e7,1.0e7], args = vars, bounds = ((10.0,self.pa,10.0),(10.0e10,10.0e10,10.0e10)))
        Isp_m = self.get_Isp_m(root["x"][1])
        Isp_a = self.get_Isp_a(root["x"][1], vars[0], vars[1])
        return -(Isp_m + Isp_a*(vars[1]/(1.0-vars[1])) * 1.0/(self.O_F+1.0))/(1.0+(vars[1]/(1.0-vars[1])) * 1.0/(self.O_F+1.0))
    
    def get_Isp_m(self,pinj): #temporaty, needs modification
        return NT.Turbine_nozzle(self.m,pinj,self.prop,self.pa,self.df,self.prop.h_fuel,self.prop.h_ox,self.prop.f_dens_g,self.prop.o_dens)
    
    def get_Isp_a(self,pt1,pt2,l):
        T2t = self.Tf_cool*(pt2/pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma)
        return (l/(1.0-l))*self.m_F*math.sqrt(2.0*self.prop.R_f*T2t*self.prop.f_gamma/(self.prop.f_gamma-1.0) * (1.0 - (self.pa/pt2)**((self.prop.f_gamma-1.0)/self.prop.f_gamma)))

    def equations(self,vars,p2t,l):
        dptop, p1t, dptfp = vars
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - p1t,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - p1t,
            self.O_F*dptop/(self.eff_p*self.prop.o_dens) + (1.0/(1.0-l))*dptfp/(self.eff_p*self.prop.f_dens_l) - (l/(1.0-l))*self.eff_m*self.eff_t*self.prop.fcp*self.Tf_cool*(1.0-(p2t/p1t)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
            ]


#Function that computes Gas Generator
class GG:
    #Global input
    ptanko : float
    ptankf : float
    prop : Propellant
    O_F : float
    eff_p : float
    eff_t : float
    eff_m : float

    #Software input
    Tf_cool : float
    dptcool : float
    m : float

    #Output
    dptop : float
    dptfp : float
    pt1 : float
    pt2 : float
    ptinj : float
    Wt : float
    Wop : float
    Wfp : float

    #Auxiliary
    dptvalve = 0.0
    dptlines = 0.0

    def __init__(self, DF : Default, prop : Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Gas generator cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_p = DF.Eff_p
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m


#Function that computes Combustion Tap-Off cycle
class TO:
    #Global input
    ptanko : float
    ptankf : float
    prop : Propellant
    O_F : float
    eff_p : float
    eff_t : float
    eff_m : float

    #Software input
    Tf_cool : float
    dptcool : float
    m : float

    #Output
    dptop : float
    dptfp : float
    pt1 : float
    pt2 : float
    ptinj : float
    Wt : float
    Wop : float
    Wfp : float

    #Auxiliary
    dptvalve = 0.0
    dptlines = 0.0

    def __init__(self, DF : Default, prop : Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Tap-off cycle not supported")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_p = DF.Eff_p
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m


#Function that computes Electric driven pumps cycle
class EL:
     #Global input
    ptanko : float
    ptankf : float
    prop : Propellant
    O_F : float
    eff_p : float
    eff_m : float
    wmotor : float

    #Software input
    Tf_cool : float
    dptcool : float
    m : float

    #Output
    dptop : float
    dptfp : float
    ptinj : float
    Wop : float
    Wfp : float

    #Auxiliary
    dptvalve = 0.0
    dptlines = 0.0

    def __init__(self, DF : Default, prop : Propellant, O_F : float, Tf_cool : float, dptcool : float, m : float):
        print("Coolant bleed cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_p = DF.Eff_p
        self.eff_m = DF.Eff_m
        self.wmotor = DF.Wmotor
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m

    def results(self):
        self.dptop, self.dptfp, self.ptinj = fsolve(self.equations,[1.0e6,1.0e6,1.0e6])
        self.Wop = self.O_F/(self.O_F+1.0) * self.m * self.dptop / (self.eff_p*self.prop.o_dens)
        self.Wfp = 1.0/(self.O_F+1.0) * self.m * self.dptfp / (self.eff_p*self.prop.f_dens_l)
        print([self.dptop, self.dptfp, self.ptinj])
        print([self.Wop,self.Wfp,self.wmotor])

    def equations(self,vars):
        dptop, dptfp, pinj = vars
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - pinj,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - pinj,
            self.O_F*dptop/(self.eff_p*self.prop.o_dens) + dptfp/(self.eff_p*self.prop.f_dens_l) - self.eff_m*self.wmotor
            ]


#Main Function
if __name__ == '__main__':
    print('Loading...')
    print(TurboM(def0, prop_, O_F_, Pa_, Tf_cool_, dptcool_, m_))
    print('\nProcess Terminated')

    #Isp = NT.Turbine_nozzle(m_,3.0e5,prop_,Pa_,def0,prop_.h_fuel,prop_.h_ox,prop_.f_dens_g,prop_.o_dens)
    #print(Isp)