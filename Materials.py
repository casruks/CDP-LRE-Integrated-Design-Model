import math as mth
import Aux_classes as aux

class Materials:
    def __init__(self, material, density, yieldstress_l, Emod, OpTemp_u, k, cost, heat_cap, ultstress, mu, thr_exp):
        self.material = material #Description of the materials
        self.density = density #Density of the material [kg/m^3]
        self.yieldstress_l = yieldstress_l #Yieldstress [Pa]
        self.Emod = Emod #Young's Modulus [Pa]
        self.OpTemp_u = OpTemp_u #Maximum Operational Temperature in [K]
        self.k = k #Thermal Conductivity of the Material in [W/m*K]
        self.cost = cost #Cost per kilogram [Eur/kg] 
        self.heat_cap = heat_cap #Heat capacity [J/kg*K]
        self.ulstress = ultstress #Ultimate Tensile Strength [Pa]
        self.mu = mu #Poisson Ratio [dimensionless]  
        self.thr_exp = thr_exp  #Coefficient of Thermal Expansion [mm/m*K]
       


Custom                  =       Materials('Custom Material',                    21000.0,  2300.0e6,   471.0e9, 2200.0,   48.0,   938.0,  145.0,  2500.0e6,   0.265,  7.25e-6)
Rhenium                 =       Materials('Rhenium-Iridium alloy',              21000.0,  2300.0e6,   471.0e9, 2200.0,   48.0,   938.0,  145.0,  2500.0e6,   0.265,  7.25e-6)
Al_7075_T6              =       Materials('Aluminium 7075 T6',                  2770.0,   530.0e6,    72.0e9,  373.15,   130.0,  3.97,   929.0,  580.0e6,    0.335,  24.1e-6)
Ti6Al4V                 =       Materials('Tungsten & Aluminium alloy',         4428.0,   898.0e6,    115.0e9, 693.15,   6.7,    23.8,   570.0,  980.0e6,    0.370,  9.1e-6)
Haynes_188              =       Materials('Nickel Alloy',                       9010.0,   419.0e6,    244.0e9, 1423.0,   10.8,   30.6,   436.0,  953.0e6,    0.323,  12.6e-6)
Inc_X_750               =       Materials('Inconel-750, Nickel-Chromium Alloy', 8290.0,   763.0e6,    222.0e9, 1090.0,   12.4,   20.4,   480.0,  1140.0e6,   0.312,  14.0e-6)
Inc_600                 =       Materials('Inconel-600, Nickel-Chromium Alloy', 8470.0,   290.0e6,    206.0e9, 982.0,    15.9,   20.2,   480.0,  586.0e6,    0.302,  13.1e-6)
Inc_718                 =       Materials('Inconel-718, Nickel-Chromium Alloy', 8190.0,   1036.0e6,   211.0e9, 825.0,    11.4,   16.6,   458.0,  914.0e6,    0.302,  13.4e-6)
Inc_A_286               =       Materials('A-286_Nickel-Chromium Alloy',        7920.0,   7920.0e6,   201.0e9, 700.0,    23.9,   5.43,   430.0,  1160.0e6,   0.330,  13.7e-6)
Columbium_c103          =       Materials('Niobium (Colombium) - cold rolled',  8950.0,   730.0e6,    92.0e9,  1255.0,   44.0,   225.0,  350.0,  800.0e6,    0.390,   8.4e-6)
Copper_structural       =       Materials('Copper',                             8940.0,   40.0e6,     128.0e9, 573.0,    398.0,  5.60,   387.0,  160.0e6,    0.350,  16.9e-6)
D6AC_Steel              =       Materials('D6AC Steel Alloy',                   7870.0,   1450.0e6,   210.0e9, 1180.0,   52.0,   1.12,   480.0,  2360.0e6,   0.295,  13.1e-6)
default                 =       Materials('default_coating',                    0.0,      0.0,        0.0,     0.0,      1,      0.0,    0.0,    0.0,        0.0,    0.0)

#Coating Materials:
Coat_Custom             =       Materials('Carbon-Carbon Matrix coating',       1950.0,   0,          0,          2400.0,   37.4,   6.4,   2093.0,   0,   0,   18.0e-6)
Copper                  =       Materials('Copper coating',                     8940.0,   40.0e6,     128.0e9,    573.0,    390.0,  6.4,   388.0,    0,   0,   18.0e-6)
Narloy_Z                =       Materials('Copper Alloy Coating',               9130.0,   0,          0,          740.0,    290.0,  6.4,   388.0,    0,   0,   22.0e-6)
GRCop_84                =       Materials('GRPCop_84 Copper alloy coating',     8756.0,   0,          0,          1000.0,   351.0,  6.4,   388.0,    0,   0,   20.0e-6)
Silica                  =       Materials('Scilica Coating',                    1700.0,   0,          0,          1200.0,   0.55,   6.4,   1256.0,   0,   0,   20.0e-6)
Carbon                  =       Materials('Carbon-Carbon Matrix coating',       1950.0,   0,          0,          2400.0,   37.4,   6.4,   2093.0,   0,   0,   20.0e-6)


##Computing Mass nozzle:
class ReferenceEngine:
    def __init__(self, pc, thrust, arear, rt, mprop, rhoprop, FS, Material_NCG, Material_P, Material_V, mfrac_tube, mfrac_manifold, mfrac_jacket, mfrac_chamber, mfrac_pump, mfrac_valve, RefMass):
        self.pc = pc #Combustion chamber pressure of the reference engine
        self.thrust = thrust #Thrust of the reference engine
        self.arear = arear #Area ratio 
        self.rt = rt #Throat area
        self.mprop = mprop #Mass flow rate
        self.rhoprop = rhoprop #Density of the propellant
        self.FS = FS #Material safety factor
        self.Material_NCG = Material_NCG #Nozzle material
        self.Material_P = Material_P #Turbopump material
        self.Material_V = Material_V #Valve material
        self.mfrac_tube = mfrac_tube #Mass fraction of regenerative nozzle tubes
        self.mfrac_manifold = mfrac_manifold #Mass fraction of regenerative nozzle manifold
        self.mfrac_jacket = mfrac_jacket #Mass fraction of regenerative nozzle jacket
        self.mfrac_chamber = mfrac_chamber #Mass fraction of combustion chamber + gas generator(if applicable)
        self.mfrac_pump = mfrac_pump #Mass fraction of turbopump
        self.mfrac_valve = mfrac_valve #Mass fraction of valve
        self.RefMass = RefMass #Total mass of the reference engine

RL10    = ReferenceEngine(3.20e6, 7.34e4, 61.1,  0.076, 16.85, 351.91, 1.1, Inc_718, D6AC_Steel, D6AC_Steel, 0.0903, 0.2575, 0.0484, 0.2281, 0.3096, 0.0661, 138) #Expander Cycle Reference
LE5     = ReferenceEngine(3.65e6, 1.03e5, 140.0, 0.068, 23.33, 343.83, 1.1, Inc_718, Al_7075_T6, Al_7075_T6, 0.1419, 0.0578, 0.1021, 0.2804, 0.3042, 0.1136, 255)#Gas Generator Cycle Reference #fix rhoprop
SSME    = ReferenceEngine(2.04e7, 2.28e6, 77.5,  0.138, 512.6, 361.89, 1.2, Inc_718, D6AC_Steel, D6AC_Steel, 0.0907, 0.1984, 0.0737, 0.1955, 0.2763, 0.1654, 3177)#Stage Combustion Cycle Reference #fixrhoprop

#Mass estimation function Nozzle Tubes:
def Mass(Pc, material_N, material_V, arear, rt, mprop, FS, cycle, x, R, t, O_prop,F_prop,O_F):
    
    total_surf = 0
    Mass_warnings = 0
    Mass_error = 0

    reference = RL10

    if cycle == 3:
        reference = SSME
    elif cycle == 2:
        reference == LE5
    else:
        reference == RL10

    #Warnings:
    for i in range(len(t)):
        if t[i] < 0.001:
            Mass_warnings=Mass_warnings|(1<<0)
        else:
            Mass_warnings = Mass_warnings & (~(1<<0))
    #Density of the Propellant:
    rho_prop = ((O_prop*F_prop)*(1+OF))/(F_prop*OF+O_prop)
    
    #Nozzle Mass:
    for i in range(len(x)-1):
        if t[i] < 0.001:
            total_surf += mth.dist([x[i+1],R[i+1]],[x[i],R[i]])*t[i]*R[i]
        else:
            total_surf += mth.dist([x[i+1],R[i+1]],[x[i],R[i]])*0.001*R[i]
    
    NozzleMass = total_surf*2*mth.pi*material_N.density

    if NozzleMass < 0:
        Mass_error=Mass_error|(1<<0)
    
    #Cooling Mass:
    TubeMass = reference.mfrac_tube*(((Pc/reference.pc)**1)*((material_N.density/reference.Material_NCG.density)**1)*(((material_N.yieldstress_l/FS)/(reference.Material_NCG.yieldstress_l/reference.FS))**(-1))*((arear/reference.arear)**2)*((rt/reference.rt)**2)) #Dimensionless Nozzle Tube Mass
    ManifoldMass = reference.mfrac_manifold*(((Pc/reference.pc)**1)*((material_N.density/reference.Material_NCG.density)**1)*((mprop/reference.mprop)**1)*((material_N.yieldstress_l/reference.Material_NCG.yieldstress_l)**(-1))*((rho_prop/reference.rhoprop)**(-1))*((rt/reference.rt)**2)) #Dimensionless Nozzle Manifold Mass

    if (TubeMass + ManifoldMass) < 0:
        Mass_error=Mass_error|(1<<1)
        return 0, Mass_error,Mass_warnings

    #TurboPump Mass
    #PumpMass = reference.mfrac_pump*(((Pc/reference.pc)**0.15)*((material_P.density/reference.material_P.density)**1)*(((material_P.yieldstress_l/FS)/(reference.Material_p.yieldstress_l/reference.FS)*(-1)))*((mprop/reference.mprop)**0.9)*((rho_prop/reference.rhoprop)**(-0.45))*((Ns/Ns0)**(-0.6)))
    PumpMass = reference.mfrac_pump*(((Pc/reference.pc)**0.53)*((mprop/reference.mprop)**0.53)*(rho_prop/reference.rhoprop)**(-0.53)) #Dimensionless mass of turbopump (historic data method)
    
    if PumpMass < 0:
        Mass_error=Mass_error|(1<<2)
        return 0, Mass_error,Mass_warnings


    #Valves
    ValveMass = reference.mfrac_valve*(((Pc/reference.pc)**1.0)*((material_V.density/reference.Material_V.density)**1)*(((material_V.yieldstress_l/FS)/(reference.Material_V.yieldstress_l/reference.FS))**(-1))*((mprop/reference.mprop)**1)*((rho_prop/reference.rhoprop)**(-1)))#Dimensionless mass of Valves
    
    if ValveMass < 0:
        Mass_error=Mass_error|(1<<3)
        return 0, Mass_error,Mass_warnings
    
    #Total:
    Total_Mass = (ManifoldMass + TubeMass + PumpMass + ValveMass)*reference.RefMass + NozzleMass
    return Total_Mass, Mass_error, Mass_warnings


def RhoProp(O_prop, F_prop, OF):
    rho_prop = 0
    rho_prop = ((O_prop*F_prop)*(1+OF))/(F_prop*OF+O_prop)
    return rho_prop


##Cost function 
def Cost(m_engine,f1,f3,R,n,cycle):
    TotalCost_MY = 0
    cost_error = 0
    cost_warning = 0

    if R < 0.8:
        cost_warning = cost_warning|(1<<0)

    if cycle == 5:
        f2 = 0.25+1.1096*R**51.968
    elif cycle == 4:
        f2 = 0.27+1.0750*R**46.994
    else:
        f2 = 0.33+1.0769*R**51.360

    f4  = -0.0553*mth.log(n) + 1.0011
    C_D = 1.1*f1*f2*f3*162*m_engine**0.58
    F_E = 5*f4*m_engine**0.46
    TotalCost_MY = C_D + F_E

    if TotalCost_MY < 0:
        cost_error = cost_error|(1<<0)
        return 0, cost_warning, cost_error
    return TotalCost_MY, cost_warning, cost_error

def RhoProp(O_prop, F_prop, OF):
        rho_prop = 0
        rho_prop = ((O_prop*F_prop)*(1+OF))/(F_prop*OF+O_prop)
        return rho_prop

##Reuseability: Prediction of low cycle fatigue life of the thrust chamber
def Reuseability(material, sigma_T, Twg, Twc, H, l, w, p):
    
    ReuseabilityError = 0
    ReuseabilityWarning = 0

    DT = Twg - Twc #Twg = Temperature on the gas side of the chamber, Twc = Temperature on the coolant channel wall temperature

    if Twg < Twc:
        ReuseabilityError == ReuseabilityError|(1<<0)
        return 0,ReuseabilityError,ReuseabilityWarning

    w = l
    
    #Inelastic Strain:
    if sigma_T > material.yieldstress_l:
        ep_pl1 =  (sigma_T - material.yieldstress_l)/material.Emod
    else:
        ep_pl1 == 0
        ReuseabilityWarning = ReuseabilityWarning|(1<<0)
    
    ep_pl2 = (material.Emod*(material.thr_exp*DT)**2)/(12*material.yieldstress_l*(1-material.mu)**2)
    ep_1 = (ep_pl1+ep_pl2)

    # Deflection:
    def1 = 2*((H/ep_1) - mth.sqrt(((H/ep_1)**2) - ((l/4)**2)))
    def2 = (ep_1*p*l**2)/(4*H*material.yieldstress_l)
    def_tot = def1+def2

    #Ligament Deformation:
    # t_N = (N*w*def_tot)/(l+w)
    # t_N_min = (2*H*(l+w) - N*w*def_tot)/(l+w)
    # t_N_max = (2*H*(l+w)**2 + N*w*l*def_tot)/(l+w)**2


    # #Critical Thickness:
    q = material.ultstress/material.yieldstress_l
    tcr = 2*H*(1-mth.e**(-q))

    # #Instability life:
    Nf = (tcr*(l+w))/(def_tot*w) 
    return Nf, ReuseabilityError, ReuseabilityWarning
