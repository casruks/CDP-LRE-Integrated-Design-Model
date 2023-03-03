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
Tf_cool_ = 600.0
dptcool_ = 1.0e5
m_ = 2.0

#Placeholder for propellant class
class Propellant:
    o_dens = 1141.0 #oxidizer density
    f_dens = 71.0 #fuel density
    f_gamma = 1.4 #fuel gamma
    fcp = 14307.0 #fuel cp
    R_f = 4.1573 #fuel gas constant

    def __init__(self,type): #Placeholder
        if(type==1):
            ox_dens=1141.0

prop_ = Propellant(0)

#Placeholder for default class
class Default:
    #Turbomachinery
    cycle_type = "EX"
    Eff_t = 0.6 #Turbine efficiency
    Eff_p = 0.6 #Pump efficiency
    Eff_m = 0.95 #Mechanical efficiency between turbine and pumps
    p_to = 1.0e5 #oxidizer tank storage pressure
    ptf = 1.0e5 #Fuel tank oxidizer pressure

def0 = Default()

#Main function, which calls the function for the selected cycle type
def TurboM(Default, prop, O_F, p_a, Tf_cool, dptcool, m):
    match Default.cycle_type:
        case "GG": #Gas Generator
            turbo = GG(Default.p_to, Default.ptf, prop, O_F, Default.Eff_p, Default.Eff_t, Default.Eff_m, p_a, Tf_cool, dptcool, m)
        case "EX": #Expander cycle
            turbo = EX(Default.p_to, Default.ptf, prop, O_F, Default.Eff_p, Default.Eff_t, Default.Eff_m, Tf_cool, dptcool, m)
        case "SC": #Staged combustion cycle
            turbo = SC(Default.p_to, Default.ptf, prop, O_F, Default.Eff_p, Default.Eff_t, Default.Eff_m, p_a, Tf_cool, dptcool, m)
        case "CB": #Coolant bleed cycle
            turbo = CB(Default.p_to, Default.ptf, prop, O_F, Default.Eff_p, Default.Eff_t, Default.Eff_m, p_a, Tf_cool, dptcool, m)
        case "TO": #Combustion tap off cycle
            turbo = TO(Default.p_to, Default.ptf, prop, O_F, Default.Eff_p, Default.Eff_t, Default.Eff_m, p_a, Tf_cool, dptcool, m) #Will not be currently implmented
        case "NO": #No turbomachinery
            dptvalve = 0.0
            dptlines = 0.0
            return [Default.p_to-dptvalve-dptlines, Default.ptf-dptvalve-dptlines]
        case "EL": #Electric motor to dirve pumps
            turbo = EL(Default.p_to, Default.ptf, prop, O_F, Default.Eff_pump, Default.Eff_turb, Default.eff_m, Tf_cool, dptcool, m)
        case _:
            print("Cycle Not recognized")
            return
    
    turbo.results();
    return [turbo.ptinj]
            

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

    def __init__(self, p_to : float, ptf : float, prop : Propellant, O_F : float, eff_pump : float, eff_turb : float, eff_m : float, Tf_cool : float, dptcool : float, m : float):
        print("Expander cycle selected")
        self.ptanko = p_to
        self.ptankf = ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_p = eff_pump
        self.eff_t = eff_turb
        self.eff_m = eff_m
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m


    def results(self):
        res = minimize(self.opt, [1.0e6], method = 'Nelder-Mead', bounds=[[0.0, 10.0e12]])
        self.dptfp = res["x"][0]
        self.dptop, self.pt1, self.pt2, self.ptinj = fsolve(self.equations,[1.0e6,1.0e7,1.0e6,1.0e6],self.dptfp)
        self.Wop = self.O_F/(self.O_F+1.0) * self.m * self.dptop / (self.eff_p*self.prop.o_dens)
        self.Wfp = 1.0/(self.O_F+1.0) * self.m * self.dptfp / (self.eff_p*self.prop.f_dens)
        self.Wt = 1.0/(self.O_F+1.0) * self.m * self.eff_t * self.prop.fcp * self.Tf_cool * (1.0-(self.pt2/self.pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
        print(res)
        print([self.dptop, self.dptfp, self.pt1, self.pt2, self.ptinj])
        print([self.Wop,self.Wfp,self.Wt])

    def opt(self,dptfp):
        root = least_squares(self.equations,[1.0e6,1.0e7,1.0e6,1.0e6], args = dptfp, bounds = ((0.0,0.0,0.0,0.0),(10.0e10,10.0e10,10.0e10,10.0e10)))
        return -1.0*root["x"][3]

    def equations(self,vars,dptfp):
        dptop, p1t, p2t, pinj = vars
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - pinj,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - p1t,
            p2t - self.dptlines - pinj,
            self.O_F*dptop/(self.eff_p*self.prop.o_dens) + dptfp/(self.eff_p*self.prop.f_dens) - self.eff_m*self.eff_t*self.prop.fcp*self.Tf_cool*(1.0-(p2t/p1t)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
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

    def __init__(self, p_to : float, ptf : float, prop : Propellant, O_F : float, eff_pump : float, eff_turb : float, eff_m : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Staged Combustion cycle selected")
        self.ptanko = p_to
        self.ptankf = ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_p = eff_pump
        self.eff_t = eff_turb
        self.eff_m = eff_m
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m


#Function that computes coolant bleed
class CB:
    #Global input
    ptanko : float
    ptankf : float
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

    def __init__(self, p_to : float, ptf : float, prop : Propellant, O_F : float, eff_pump : float, eff_turb : float, eff_m : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Coolant bleed cycle selected")
        self.ptanko = p_to
        self.ptankf = ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_p = eff_pump
        self.eff_t = eff_turb
        self.eff_m = eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.m_O = O_F/(O_F+1.0) * m
        self.m_F = 1.0/(O_F+1.0) * m


    def results(self):
        res = minimize(self.opt, [1.0e5,0.01], method = 'Nelder-Mead', bounds=[[self.pa,10.0e12],[1.0e-5,1.0]])
        self.pt2 = res["x"][0]
        self.l = res["x"][1]
        self.dptop, self.pt1, self.dptfp = fsolve(self.equations,[1.0e6,1.0e7,1.0e7],(self.pt2,self.l))
        self.Wop = self.m_O * self.dptop / (self.eff_p*self.prop.o_dens)
        self.Wfp = (1.0/(1.0-self.l)) * self.m_F * self.dptfp / (self.eff_p*self.prop.f_dens)
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
        return 600*(pinj/1.0e7)**(1.0/5.0)
    
    def get_Isp_a(self,pt1,pt2,l):
        T2t = self.Tf_cool*(pt2/pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma)
        return (l/(1.0-l))*self.m_F*math.sqrt(2.0*self.prop.R_f*T2t*self.prop.f_gamma/(self.prop.f_gamma-1.0) * (1.0 - (self.pa/pt2)**((self.prop.f_gamma-1.0)/self.prop.f_gamma)))

    def equations(self,vars,p2t,l):
        dptop, p1t, dptfp = vars
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - p1t,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - p1t,
            self.O_F*dptop/(self.eff_p*self.prop.o_dens) + (1.0/(1.0-l))*dptfp/(self.eff_p*self.prop.f_dens) - (l/(1.0-l))*self.eff_m*self.eff_t*self.prop.fcp*self.Tf_cool*(1.0-(p2t/p1t)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
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

    def __init__(self, p_to : float, ptf : float, prop : Propellant, O_F : float, eff_pump : float, eff_turb : float, eff_m : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Gas generator cycle selected")
        self.ptanko = p_to
        self.ptankf = ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_p = eff_pump
        self.eff_t = eff_turb
        self.eff_m = eff_m
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

    def __init__(self, p_to : float, ptf : float, prop : Propellant, O_F : float, eff_pump : float, eff_turb : float, eff_m : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("not supported")
        self.ptanko = p_to
        self.ptankf = ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_p = eff_pump
        self.eff_t = eff_turb
        self.eff_m = eff_m
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

    def __init__(self, p_to : float, ptf : float, prop : Propellant, O_F : float, eff_pump : float, eff_turb : float, eff_m : float, Tf_cool : float, dptcool : float, m : float):
        print("not supported")


#Main Function
if __name__ == '__main__':
    print('Loading...')
    print(TurboM(def0, prop_, O_F_, Pa_, Tf_cool_, dptcool_, m_))
    print('\nProcess Terminated')