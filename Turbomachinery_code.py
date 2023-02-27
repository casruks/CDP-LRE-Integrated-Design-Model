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

#Values used during tests, corresponding to Global or software input
P_entrance = 1.0e5
O_F = 5.0
eff_p = 0.6
eff_t = 0.6
eff_m = 0.95
Pa = 1.0e5
Tf_cool = 600.0
dptcool = 1.0e5
m = 2.0
CYTYPE = "EX"

#Placeholder for propellant class
class Propellant:
    o_dens = 1141.0
    f_dens = 71.0
    f_gamma = 1.4
    fcp = 14307.0
    R_f = 4.1573

    def __init__(self,type):
        if(type==1):
            ox_dens=1141.0

prop = Propellant(0)

#Main function, which calls the function for the selected cycle type
def TurboM(p_to, ptf, prop, O_F, eff_pump, eff_turb, eff_m, p_a, Tf_cool, dptcool, m, cycle_type):
    match cycle_type:
        case "GG": #Gas Generator
            turbo = GG(p_to, ptf, prop, O_F, eff_pump, eff_turb, eff_m, Tf_cool, dptcool, m)
        case "EX": #Expander cycle
            turbo = EX(p_to, ptf, prop, O_F, eff_pump, eff_turb, eff_m, Tf_cool, dptcool, m)
        case "SC": #Staged combustion cycle
            turbo = SC(p_to, ptf, prop, O_F, eff_pump, eff_turb, eff_m, Tf_cool, dptcool, m)
        case "CB": #Coolant bleed cycle
            turbo = CB(p_to, ptf, prop, O_F, eff_pump, eff_turb, eff_m, p_a, Tf_cool, dptcool, m)
        case "TO": #Combustion tap off cycle
            turbo = TO(p_to, ptf, prop, O_F, eff_pump, eff_turb, eff_m, Tf_cool, dptcool, m) #Will not be currently implmented
        case "NO": #No turbomachinery
            dptvalve = 0.0
            dptlines = 0.0
            return [p_to-dptvalve-dptlines, ptf-dptvalve-dptlines]
        case _:
            print("Cycle Not recognized")
            return
    
    turbo.results();
    return [turbo.dptop,turbo.dptfp,turbo.pt1,turbo.pt2,turbo.ptinj,turbo.Wt,turbo.Wop,turbo.Wfp]
            

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
    dptop = 0.0
    dptfp = 0.0
    pt1 = 0.0
    pt2 = 0.0
    ptinj = 0.0
    Wt = 0.0
    Wop = 0.0
    Wfp = 0.0

    #Auxiliary
    dptvalve = 0.0
    dptlines = 0.0

    def __init__(self, p_to, ptf, prop, O_F, eff_pump, eff_turb, eff_m, Tf_cool, dptcool, m):
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
        dptfp = res["x"][0]
        dptop, pt1, pt2, ptinj = fsolve(self.equations,[1.0e6,1.0e7,1.0e6,1.0e6],dptfp)
        Wop = O_F/(O_F+1.0) * m * dptop / (eff_p*prop.o_dens)
        Wfp = 1.0/(O_F+1.0) * m * dptfp / (eff_p*prop.f_dens)
        Wt = 1.0/(O_F+1.0) * m * eff_t * prop.fcp * Tf_cool * (1.0-(pt2/pt1)**((prop.f_gamma-1.0)/prop.f_gamma))
        print(res)
        print([dptop, dptfp, pt1, pt2, ptinj])
        print([Wop,Wfp,Wt])

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
    dptop: float
    dptfp: float
    pt1: float
    pt2: float
    ptinj: float
    Wt: float
    Wop: float
    Wfp: float

    #Auxiliary
    dptvalve = 0.0
    dptlines = 0.0

    def __init__(self, p_to, ptf, prop, O_F, eff_pump, eff_turb, eff_m, Tf_cool, dptcool, m):
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
    dptop = 0.0
    dptfp = 0.0
    pt1 = 0.0
    pt2 = 0.0
    ptinj = 0.0
    Wt = 0.0
    Wop = 0.0
    Wfp = 0.0
    l = 0.0 #fraction of fuel for bleed

    #Auxiliary
    dptvalve = 0.0
    dptlines = 0.0
    m_O = 1.0
    m_F = 1.0

    def __init__(self, p_to, ptf, prop, O_F, eff_pump, eff_turb, eff_m, p_a, Tf_cool, dptcool, m):
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
        pt2 = res["x"][0]
        l = res["x"][1]
        dptop, pt1, dptfp = fsolve(self.equations,[1.0e6,1.0e7,1.0e7],(pt2,l))
        Wop = self.m_O * dptop / (eff_p*prop.o_dens)
        Wfp = (1.0/(1.0-l)) * self.m_F * dptfp / (eff_p*prop.f_dens)
        Wt = (l/(1.0-l)) * self.m_F * eff_t * prop.fcp * Tf_cool * (1.0-(pt2/pt1)**((prop.f_gamma-1.0)/prop.f_gamma))
        print(res)
        print([dptop, dptfp, pt1, pt2, l])
        print([Wop,Wfp,Wt])

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
    dptop: float
    dptfp: float
    pt1: float
    pt2: float
    ptinj: float
    Wt: float
    Wop: float
    Wfp: float

    #Auxiliary
    dptvalve = 0.0
    dptlines = 0.0

    def __init__(self, p_to, ptf, prop, O_F, eff_pump, eff_turb, eff_m, Tf_cool, dptcool, m):
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
    dptop: float
    dptfp: float
    pt1: float
    pt2: float
    ptinj: float
    Wt: float
    Wop: float
    Wfp: float

    #Auxiliary
    dptvalve = 0.0
    dptlines = 0.0

    def __init__(self, p_to, ptf, prop, O_F, eff_pump, eff_turb, eff_m, Tf_cool, dptcool, m):
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


#Main Function
if __name__ == '__main__':
    print('Loading...')
    TurboM(P_entrance, P_entrance, prop, O_F, eff_p, eff_t, eff_m, Pa, Tf_cool, dptcool, m, CYTYPE)
    print('\nProcess Terminated')