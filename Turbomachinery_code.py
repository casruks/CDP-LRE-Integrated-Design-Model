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
from rocketcea.cea_obj_w_units import CEA_Obj
import Nozzle_turbine as NT
import Aux_classes as aux

# Errors
    # 0: Cycle not recognized
    # 1: Optimization did not converge
    # 2: system could not be solved
    # 3: pump potency does not equal turbine potency
    # 4: pump or turbine potency is negative

# Warnings
    # Currently no warnings emitted


#aux function, which calls the function for the selected cycle type
def TurboM(Default : aux.Default, prop : aux.Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
    # Default : [Default class object] contains information on the default values employed in the application
    # prop : [Propellant class object] contains the propellant data
    # p_a : [Pa] atmospheric pressure
    # Tf_cool : [K] [escalar] fuel temperature after cooling channel
    # dptcool : [Pa] [escalar] pressure losses in cooling channel
    # m : [kg/s] mass flow in the nozzle
    match Default.cycle_type:
        case 0: #Expander cycle
            turbo = EX(Default, prop, O_F, Tf_cool, dptcool, m)
        case 1: #Coolant bleed cycle
            turbo = CB(Default, prop, O_F, p_a, Tf_cool, dptcool, m)
        case 2: #Gas Generator
            turbo = GG(Default, prop, O_F, p_a, Tf_cool, dptcool, m)
        case 3: #Staged combustion cycle (currently unsupported)
            turbo = SC(Default, prop, O_F, p_a, Tf_cool, dptcool, m)
            return [0, 0, 0, 0, 1] #Returns 1 in last position, indicating aux to break
        case 4: #Electric motor to dirve pumps
            turbo = EL(Default, prop, O_F, Tf_cool, dptcool, m)
        case 5: #No turbomachinery
            dptvalve = Default.v_loss
            dptlines = Default.line_loss
            dptmixing = 0.0
            return [((Default.p_to-dptvalve-dptlines)*O_F+Default.ptf-dptvalve-dptlines)/(1.0+O_F) - dptmixing, 0.0, 0.0, 0.0, 0.0, m, 0]
        case 6: #Combustion tap off cycle (currently unsuported)
            turbo = TO(Default, prop, O_F, p_a, Tf_cool, dptcool, m) #Will not be currently implmented
            return [0, 0, 0, 0, 0, 0, 1] #Returns 1 in last position, indicating aux to break
        case _:
            print("Cycle Not recognized")
            return [0, 0, 0, 0, 0, 0, 1] #Returns 1 in last position, indicating aux to break
    
    turbo.results();
    return [turbo.ptinj, turbo.Wop, turbo.Wfp, turbo.Wt, turbo.l, turbo.mt, turbo.br] #Returns 0 in last position, indicating aux to continue


#Function that computes the expander cycle
class EX:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_t : float #[-] efficiecny of turbine
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    pt1 : float #[Pa] Total pressure at turbine inlet
    pt2 : float #[Pa] Total pressure at turbine outlet
    ptinj : float #[Pa] Total pressure at injector inlet
    Wt : float #[W] Turbine power
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power
    l = 1.0 #[-] fraction of fuel for bleed
    mt = 1.0 #Total mass flow

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines

    #Flags
    br = 0 #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, Tf_cool : float, dptcool : float, m : float):
        print("Expander cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.dptvalve = DF.v_loss
        self.dptlines = DF.line_loss

    #Obtain results, calling optimization procedure and then computing variables of interest
    def results(self):
        #Minimization
        res = minimize(self.opt, [1.0e6], method = 'Nelder-Mead', bounds=[[0.0, 10.0e12]])
        self.dptfp = res["x"][0]

        #Computation of results
        self.dptop, self.pt1, self.pt2, self.ptinj = fsolve(self.equations,[1.0e6,1.0e7,1.0e6,1.0e6],self.dptfp)
        self.Wop = self.O_F/(self.O_F+1.0) * self.m * self.dptop / (self.eff_po*self.prop.o_dens)
        self.Wfp = 1.0/(self.O_F+1.0) * self.m * self.dptfp / (self.eff_pf*self.prop.f_dens_l)
        self.Wt = 1.0/(self.O_F+1.0) * self.m * self.eff_t * self.prop.fcp * self.Tf_cool * (1.0-(self.pt2/self.pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
        self.mt = self.m
        
        #Check to see if results are coherent
        if(abs(self.Wop+self.Wfp-self.Wt*self.eff_m) > 0.01 ):
            self.br = self.br | 1<<3
        if(self.Wop < 0.0 or self.Wfp < 0.0):
            self.br = self.br | 1<<4
        if(not res["success"]):
            self.br = self.br | 1<<1

        #Debugging information
        print([self.dptop, self.dptfp, self.pt1, self.pt2, self.ptinj])
        print([self.Wop,self.Wfp,self.Wt])

    #Optimize for maximum chamber pressure
    def opt(self,dptfp):
        root = least_squares(self.equations,[1.0e6,1.0e7,1.0e6,1.0e6], args = dptfp, bounds = ((0.0,0.0,0.0,0.0),(10.0e10,10.0e10,10.0e10,10.0e10)))

        if(not root["success"] or abs(sum(root["fun"])) > 0.01):
            self.br = self.br | 1<<2
            return abs(sum(root["fun"]))*9999
        else:
            self.br = self.br & ~(1<<2)

        return -1.0*root["x"][3]

    #System of equations to be solved
    def equations(self,vars,dptfp):
        dptop, p1t, p2t, pinj = vars
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - pinj,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - p1t,
            p2t - self.dptlines - pinj,
            self.O_F*dptop/(self.eff_po*self.prop.o_dens) + dptfp/(self.eff_pf*self.prop.f_dens_l) - self.eff_m*self.eff_t*self.prop.fcp*self.Tf_cool*(1.0-(p2t/p1t)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
            ]

    
#Function that computes staged combustion
class SC:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_t : float #[-] efficiecny of turbine
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling
    pa : float #[Pa] Ambient pressure

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    pt1 : float #[Pa] Total pressure at turbine inlet
    pt2 : float #[Pa] Total pressure at turbine outlet
    ptinj : float #[Pa] Total pressure at injector inlet
    Wt : float #[W] Turbine power
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power
    l : float #[-] fraction of oxidizer for pre-combustor
    T1t : float #[K] Temperature after pre-combustor
    mt = 1.0 #Total mass flow

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines
    dptmix = 0.0 #[Pa] mixing total pressure losses
    dptcomb = 0.0 #[Pa] pressure losses in combustion

    #Flags
    br = 0 #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Staged Combustion cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.dptvalve = DF.v_loss
        self.dptlines = DF.line_loss

    #Obtain results, calling optimization procedure and then computing variables of interest
    def results(self):
        #Minimization
        res = minimize(self.opt, [1.0e5,0.01], method = 'Nelder-Mead', bounds=[[self.pa*1.2,10.0e10],[1.0e-5,0.9]])
        self.pt2 = res["x"][0]
        self.l = res["x"][1]

        #Computation of results
        self.dptop, self.pt1, self.dptfp, self.ptinj = fsolve(self.equations,[1.0e6,1.0e7,1.0e7,1.0e7],(self.pt2,self.l))
        self.Wop = self.O_F/(self.O_F+1.0) * self.m * self.dptop / (self.eff_po*self.prop.o_dens)
        self.Wfp = 1.0/(self.O_F+1.0) * self.m * self.dptfp / (self.eff_pf*self.prop.f_dens_l)
        self.Wt = (1.0+self.l*self.O_F)/(self.O_F+1.0) * self.m * self.eff_t * self.prop.fcp * self.T1t * (1.0-(self.pt2/self.pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
        
        #Check to see if results are coherent
        if(abs(self.Wop+self.Wfp-self.Wt*self.eff_m) > 0.01 ):
            self.br = self.br | 1<<3
        if(self.Wop < 0.0 or self.Wfp < 0.0):
            self.br = self.br | 1<<4
        if(not res["success"]):
            self.br = self.br | 1<<1

        #Debugging information
        print([self.dptop, self.dptfp, self.pt1, self.pt2, self.ptinj])
        print([self.Wop,self.Wfp,self.Wt])

    #Optimize for maximum chamber pressure
    def opt(self,vars):
        root = least_squares(self.equations,[1.0e6,1.0e7,1.0e7,1.0e7], args = vars, bounds = ((10.0,self.pa,10.0,self.pa*1.5),(10.0e10,10.0e10,10.0e10,10.0e10)))

        if(not root["success"] or abs(sum(root["fun"])) > 0.01):
            self.br = self.br | 1<<2
            return sum(root["fun"])*9999
        else:
            self.br = self.br & ~(1<<2)

        Isp_m = self.get_Isp_m(root["x"][1]) #Need to pass on the modified composition, not yet implemented....
        return -Isp_m
    
    def get_Isp_m(self,pinj): #temporaty, needs modification
        return NT.Turbine_nozzle(self.m,pinj,self.prop,self.pa,self.df,self.prop.h_fuel,self.prop.h_ox,self.prop.f_dens_g,self.prop.o_dens,self.O_F)

    #System of equations to be solved
    def equations(self,vars,p2t,l):
        dptop, p1t, dptfp, pinj = vars
        pcomb = self.ptankf - self.dptvalve + dptfp - self.dptcool
        self.T1t = 1.0 #Function to get combustion temperature in pre-combustor
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - pinj,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - self.dptmix - self.dptcomb - p1t,
            p2t - self.dptlines - pinj,
            self.O_F*dptop/(self.eff_po*self.prop.o_dens) + dptfp/(self.eff_pf*self.prop.f_dens_l) - (1.0+l*self.O_F)*self.eff_m*self.eff_t*self.prop.fcp*self.T1t*(1.0-(p2t/p1t)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
            ]


#Function that computes coolant bleed
class CB:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_t : float #[-] efficiecny of turbine
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling
    pa : float #[Pa] Ambient pressure

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    pt1 : float #[Pa] Total pressure at turbine inlet
    pt2 : float #[Pa] Total pressure at turbine outlet
    ptinj : float #[Pa] Total pressure at injector inlet
    Wt : float #[W] Turbine power
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power
    l : float #[-] fraction of fuel for bleed
    mt = 1.0 #Total mass flow

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines
    m_O = 1.0 #[kg/s] Oxidizer mass flow
    m_F = 1.0 #[kg/s] Fuel mass flow

    #Flags
    br = 0 #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Coolant bleed cycle selected")
        self.df = DF
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.m_O = O_F/(O_F+1.0) * m
        self.m_F = 1.0/(O_F+1.0) * m
        self.dptvalve = DF.v_loss
        self.dptlines = DF.line_loss

    #Obtain results, calling optimization procedure and then computing variables of interest
    def results(self):
        #Minimization
        res = minimize(self.opt, [5.0e5,5.0e6], method = 'Nelder-Mead', bounds=[[self.pa*1.1,1.0e12],[self.pa*1.4,5.0e12]],)
        self.pt2 = res["x"][0]
        self.pt1 = res["x"][1]

        #Computation of results
        self.dptop, self.l, self.dptfp = fsolve(self.equations,[1.0e6,0.05,1.0e7],(self.pt2,self.pt1))
        self.ptinj = self.pt1
        self.Wop = self.m_O * self.dptop / (self.eff_po*self.prop.o_dens)
        self.Wfp = (1.0/(1.0-self.l)) * self.m_F * self.dptfp / (self.eff_pf*self.prop.f_dens_l)
        self.Wt = (self.l/(1.0-self.l)) * self.m_F * self.eff_t * self.prop.fcp * self.Tf_cool * (1.0-(self.pt2/self.pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
        self.mt = (self.l/(1.0-self.l)) * self.m_F + self.m
        
        #Check to see if results are coherent
        if(abs(self.Wop+self.Wfp-self.Wt*self.eff_m) > 0.01 ):
            self.br = self.br | 1<<3
        if(self.Wop < 0.0 or self.Wfp < 0.0):
            self.br = self.br | 1<<4
        if(not res["success"]):
            self.br = self.br | 1<<1

        #Debugging information
        print([self.dptop, self.dptfp, self.pt1, self.pt2, self.l])
        print([self.Wop,self.Wfp,self.Wt])

    #Optimize for maximum chamber pressure
    def opt(self,vars):
        root = least_squares(self.equations,[1.0e6,0.05,1.0e7], args = vars, bounds = ((1000.0,1.0e-5,1000.0),(10.0e10,0.8,10.0e10)), max_nfev=9999)

        if(not root["success"] or abs(sum(root["fun"])) > 0.01):
            self.br = self.br | 1<<2
            return abs(sum(root["fun"]))*9999
        else:
            self.br = self.br & ~(1<<2)

        Isp_m = self.get_Isp_m(vars[1])
        Isp_a = self.get_Isp_a(vars[1], vars[0])
        return -(Isp_m + Isp_a*(root["x"][1]/(1.0-root["x"][1]))/(self.O_F+1.0))/(1.0+(root["x"][1]/(1.0-root["x"][1]))/(self.O_F+1.0))
    
    #Get Isp of aux nozzle
    def get_Isp_m(self,pinj): #temporaty, needs modification
        return NT.Turbine_nozzle(self.m,pinj,self.prop,self.pa,self.df,self.prop.h_fuel,self.prop.h_ox,self.prop.f_dens_g,self.prop.o_dens,self.O_F)
    
    #get Isp of open cycle auxiliary nozzle
    def get_Isp_a(self,pt1,pt2):
        T2t = self.Tf_cool*(pt2/pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma)
        return math.sqrt(2.0*self.prop.R_f*T2t*self.prop.f_gamma/(self.prop.f_gamma-1.0) * (1.0 - (self.pa/pt2)**((self.prop.f_gamma-1.0)/self.prop.f_gamma)))

    #System of equations to be solved
    def equations(self,vars,p2t,p1t):
        dptop, l, dptfp = vars
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - p1t,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - p1t,
            self.O_F*dptop/(self.eff_po*self.prop.o_dens) + (1.0/(1.0-l))*dptfp/(self.eff_pf*self.prop.f_dens_l) - (l/(1.0-l))*self.eff_m*self.eff_t*self.prop.fcp*self.Tf_cool*(1.0-(p2t/p1t)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
            ]


#Function that computes Gas Generator
class GG:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_t : float #[-] efficiecny of turbine
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling
    pa : float #[Pa] Ambient pressure

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    pt1 : float #[Pa] Total pressure at turbine inlet
    pt2 : float #[Pa] Total pressure at turbine outlet
    ptinj : float #[Pa] Total pressure at injector inlet
    Wt : float #[W] Turbine power
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power
    l : float #[-] fraction of fuel for bleed
    mt = 1.0 #Total mass flow

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines
    dptmix = 0.0 #[Pa] mixing total pressure losses
    dptcomb = 0.0 #[Pa] pressure losses in combustion
    m_O = 1.0 #[kg/s] Oxidizer mass flow
    m_F = 1.0 #[kg/s] Fuel mass flow

    #Flags
    br = 0 #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Gas generator cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.df = DF
        self.dptvalve = DF.v_loss
        self.dptlines = DF.line_loss
        self.l = DF.l_def
        self.ispObj = CEA_Obj( oxName=self.prop.Ox_name, fuelName=self.prop.Fuel_name,cstar_units='m/s',pressure_units='bar',temperature_units='K',isp_units='sec',density_units='kg/m^3',specific_heat_units='J/kg-K',viscosity_units='poise',thermal_cond_units='W/cm-degC')

    #Obtain results, calling optimization procedure and then computing variables of interest
    def results(self):
        #Here the main nozzle function should be called, for a given chamber pressure to get the mass flow required, then the total mass flow for that pressure is computed and it is minimized
        #Minimization
        res = minimize(self.opt, [5.0e5], method = 'Nelder-Mead', bounds=[[self.pa*1.1,1.0e8]])
        self.pt2 = res["x"][0]

        #Computation of results
        self.ispObj = CEA_Obj( oxName=self.prop.Ox_name, fuelName=self.prop.Fuel_name,cstar_units='m/s',pressure_units='bar',temperature_units='K',isp_units='sec',density_units='kg/m^3',specific_heat_units='J/kg-K',viscosity_units='poise',thermal_cond_units='W/cm-degC')
        root = least_squares(self.equations,[1.0e7,1.0e7,1.0e7,1.0e7], args = [self.pt2], bounds = ((10.0,self.pa*1.1,10.0,self.pa*1.5),(10.0e10,10.0e10,10.0e10,10.0e10)))
        self.dptop = root["x"][0]; self.pt1 = root["x"][1]; self.dptfp = root["x"][2]; self.ptinj = root["x"][3];
        if(not root["success"] or abs(sum(root["fun"])) > 0.01):
            self.br = self.br | 1<<2
        else:
            self.br = self.br & ~(1<<2)
        
        self.mt = (1.0/(1.0-self.l)) * self.m
        self.Wop = self.mt*self.O_F/(self.O_F+1.0) * self.dptop / (self.eff_po*self.prop.o_dens)
        self.Wfp = self.mt/(self.O_F+1.0) * self.dptfp / (self.eff_pf*self.prop.f_dens_l)
        self.T1t = self.ispObj.get_Tcomb(Pc=(self.pt1 - self.dptmix - self.dptcomb)/1.0e5,MR=self.O_F)
        self.Wt = self.l * self.mt * self.eff_t * (self.prop.fcp+self.prop.ocp)*0.5 * self.T1t * (1.0-(self.pt2/self.pt1)**((1.2-1.0)/1.2))
        
        #Check to see if results are coherent
        if(abs(self.Wop+self.Wfp-self.Wt*self.eff_m) > 0.01 ):
            self.br = self.br | 1<<3
        if(self.Wop < 0.0 or self.Wfp < 0.0):
            self.br = self.br | 1<<4
        if(not res["success"]):
            self.br = self.br | 1<<1

        #Debugging information
        print([self.dptop, self.dptfp, self.pt1, self.pt2, self.ptinj])
        print([self.Wop,self.Wfp,self.Wt])

    #Optimize for maximum chamber pressure
    def opt(self,vars):
        self.ispObj = CEA_Obj( oxName=self.prop.Ox_name, fuelName=self.prop.Fuel_name,cstar_units='m/s',pressure_units='bar',temperature_units='K',isp_units='sec',density_units='kg/m^3',specific_heat_units='J/kg-K',viscosity_units='poise',thermal_cond_units='W/cm-degC')
        root = least_squares(self.equations,[1.0e7,1.0e7,1.0e7,1.0e7], args = vars, bounds = ((1000.0,self.pa*1.5,1000.0,self.pa*1.5),(10.0e10,10.0e10,10.0e10,10.0e10)))
        if(not root["success"] or abs(sum(root["fun"])) > 0.01):
            self.br = self.br | 1<<2
            return np.Inf
        else:
            self.br = self.br & ~(1<<2)
        
        Isp_m = self.get_Isp_m(root["x"][3])
        Isp_a = self.get_Isp_a(root["x"][1], vars[0], self.l)
        return -(Isp_m + Isp_a*(self.l/(1.0-self.l)))/(1.0+(self.l/(1.0-self.l)))
    
    def get_Isp_m(self,pinj): #temporaty, needs modification
        return NT.Turbine_nozzle(self.m,pinj,self.prop,self.pa,self.df,self.prop.h_fuel,self.prop.h_ox,self.prop.f_dens_g,self.prop.o_dens,self.O_F)
    
    #get Isp of open cycle auxiliary nozzle
    def get_Isp_a(self,pt1,pt2,l):
        T2t = self.Tf_cool*(pt2/pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma)
        return NT.Turbine_nozzle(self.m*l/(1.0-l),pt2,self.prop,self.pa,self.df,self.prop.h_fuel,self.prop.h_ox,self.prop.f_dens_g,self.prop.o_dens,self.O_F)

    #System of equations to be solved
    def equations(self,vars,p2t):
        dptop, p1t, dptfp, pinj = vars
        self.T1t = self.ispObj.get_Tcomb(Pc=pinj/1.0e5,MR=self.O_F) # Function that returns the combustion chamber temperature
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - pinj,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - pinj,
            pinj + self.dptmix + self.dptcomb - p1t,
            self.O_F*dptop/(self.eff_po*self.prop.o_dens) + dptfp/(self.eff_pf*self.prop.f_dens_l) - self.l*(self.O_F+1.0)*self.eff_m*self.eff_t*(self.prop.fcp+self.prop.ocp)*0.5*self.T1t*(1.0-(p2t/p1t)**((1.2-1.0)/1.2))
            ]


#Function that computes Combustion Tap-Off cycle
class TO:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_t : float #[-] efficiecny of turbine
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling
    pa : float #[Pa] Ambient pressure

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    pt1 : float #[Pa] Total pressure at turbine inlet
    pt2 : float #[Pa] Total pressure at turbine outlet
    ptinj : float #[Pa] Total pressure at injector inlet
    Wt : float #[W] Turbine power
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power
    l : float #[-] fraction of fuel for bleed

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines
    m_O = 1.0 #[kg/s] Oxidizer mass flow
    m_F = 1.0 #[kg/s] Fuel mass flow

    #Flags
    br = 0 #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Tap-off cycle not supported")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.dptvalve = DF.v_loss
        self.dptlines = DF.line_loss


#Function that computes Electric driven pumps cycle
class EL:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling
    wmotor : float #[W] Power of electrical motor
    Wt : float #[W] Equal to wmotor, used for conveniance

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    ptinj : float #[Pa] Total pressure at injector inlet
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power
    l = 1.0 #[-] fraction of fuel for bleed
    mt = 1.0 #Total mass flow

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines

    #Flags
    br = 0 #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, Tf_cool : float, dptcool : float, m : float):
        print("Electrical motor cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_m = DF.Eff_m
        self.wmotor = DF.Wmotor
        self.Wt = DF.Wmotor
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.dptvalve = DF.v_loss
        self.dptlines = DF.line_loss

    #Obtain results, no optimization procedure needed for this cycle
    def results(self):
        #Computation of results
        root = least_squares(self.equations,[1.0e6,1.0e6,1.0e6])
        self.dptop, self.dptfp, self.ptinj = root["x"][0], root["x"][1], root["x"][2]
        if(not root["success"] or abs(sum(root["fun"])) > 0.01):
            self.br = self.br | 1<<2
            return 
        else:
            self.br = self.br & ~(1<<2)
        
        self.Wop = self.O_F/(self.O_F+1.0) * self.m * self.dptop / (self.eff_po*self.prop.o_dens)
        self.Wfp = 1.0/(self.O_F+1.0) * self.m * self.dptfp / (self.eff_pf*self.prop.f_dens_l)
        self.mt = self.m

        #Check to see if results are coherent
        if(abs(self.Wop+self.Wfp-self.Wt*self.eff_m) > 0.01 ):
            self.br = self.br | 1<<3
        if(self.Wop < 0.0 or self.Wfp < 0.0):
            self.br = self.br | 1<<4

        #Debugging information
        print([self.dptop, self.dptfp, self.ptinj])
        print([self.Wop,self.Wfp,self.wmotor])

    #Equations to be solved
    def equations(self,vars):
        dptop, dptfp, pinj = vars
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - pinj,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - self.dptlines - pinj,
            (self.O_F*dptop/(self.eff_po*self.prop.o_dens) + dptfp/(self.eff_pf*self.prop.f_dens_l))*self.m/(1.0+self.O_F) - self.eff_m*self.wmotor
            ]


#main Function
if __name__ == '__main__':
    print('Loading...')
    default = aux.Default(0)
    default.cycle_type = 2 # 0:EX (expander) - 1:CB (coolant bleed) - 2:GG (gas generator) - 3:SC (staged combustion) - 4:EL (electrical) - 5:PF (pressure fed)

    #Values used during tests, corresponding to Global or software input
    O_F_ = 6.04 #[-]
    Pa_ = 1.0e3 #[Pa] ambient pressure
    Tf_cool_ = 322.0 #[K] temperature after cooling
    dptcool_ = 711 #[Pa] pressure drop in cooling channel
    m_ = 456.0 #[kg/s] mass flow in the nozzle
    default = aux.Default(0)
    prop = aux.Propellant(0)
    print(TurboM(default, prop, O_F_, Pa_, Tf_cool_, dptcool_, m_))
    
    #counter = 0
    #for m_ in np.linspace(0.1,400,100):
    #    for Tf_cool_ in np.linspace(100,900,100):
    #        ptinj_, Wop_, Wfp_, Wt_, l_, mt_, br_ = TurboM(default, prop, O_F_, Pa_, Tf_cool_, dptcool_, m_)
    #        if(br_):
    #            print("error for m=" + str(m_) + ", Tf=" + str(Tf_cool_))
    #            counter = counter + 1

    #print(counter)
    print('\nProcess Terminated')