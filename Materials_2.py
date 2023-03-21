import math as mth
import Aux_classes as aux

class Materials:
    def __init__(self, material, density, yieldstress_l, Emod, OpTemp_u, k, cost, heat_cap, ultstress, mu, thr_exp):
        self.material = material
        self.density = density
        self.yieldstress_l = yieldstress_l
        self.Emod = Emod
        self.OpTemp_u = OpTemp_u
        self.k = k
        self.cost = cost
        self.heat_cap = heat_cap
        self.ulstress = ultstress
        self.mu = mu
        self.thr_exp = thr_exp 
       
##Properties of Material Class:
#material: Description of the materials
#density: in [kg/m^3]
#YieldStress_l: Yieldstress of the material in [Pa]
#Emod: Young's Modulus of the material in [Pa]
#OpTemp_u: Maximum Operational Temperature in [K]
#k: Thermal Conductivity of the Material in [W/m*K]
#cost: Cost per kilogram [Eur/kg] 
#heat_cap: Heat capacity [J/kg*K]
# ulstress: Ultimate Tensile Strength [Pa]
# mu: Poisson Ratio [dimensionless]    
#thr_exp: Coefficient of Thermal Expansion [strain/K]

Custom                  =       Materials('Custom Material',                    21000.0,  2300.0e6,   471.0e9, 2200.0,   48.0,   938.0,  145.0,  2500.0e6,   0.265,  7.25)
Rhenium                 =       Materials('Rhenium-Iridium alloy',              21000.0,  2300.0e6,   471.0e9, 2200.0,   48.0,   938.0,  145.0,  2500.0e6,   0.265,  7.25)
Al_7075_T6              =       Materials('Aluminium 7075 T6',                  2770.0,   530.0e6,    72.0e9,  373.15,   130.0,  3.97,   929.0,  580.0e6,    0.335,  24.1)
Ti6Al4V                 =       Materials('Tungsten & Aluminium alloy',         4428.0,   898.0e6,    115.0e9, 693.15,   6.7,    23.8,   570.0,  980.0e6,    0.370,  9.1)
Haynes_188              =       Materials('Nickel Alloy',                       9010.0,   419.0e6,    244.0e9, 1423.0,   10.8,   30.6,   436.0,  953.0e6,    0.323,  12.6)
Inc_X_750               =       Materials('Inconel-750, Nickel-Chromium Alloy', 8290.0,   763.0e6,    222.0e9, 1090.0,   12.4,   20.4,   480.0,  1140.0e6,   0.312,  14.0)
Inc_600                 =       Materials('Inconel-600, Nickel-Chromium Alloy', 8470.0,   290.0e6,    206.0e9, 982.0,    15.9,   20.2,   480.0,  586.0e6,    0.302,  13.1)
Inc_718                 =       Materials('Inconel-718, Nickel-Chromium Alloy', 8190.0,   1036.0e6,   211.0e9, 825.0,    11.4,   16.6,   458.0,  914.0e6,    0.302,  13.4)
Inc_A_286               =       Materials('A-286_Nickel-Chromium Alloy',        7920.0,   7920.0e6,   201.0e9, 700.0,    23.9,   5.43,   430.0,  1160.0e6,   0.330,  13.7)
Columbium_c103          =       Materials('Niobium (Colombium) - cold rolled',  8950.0,   730.0e6,    92.0e9,  1255.0,   44.0,   225.0,  350.0,  800.0e6,    0.390,   8.4)
Copper_structural       =       Materials('Copper',                             8940.0,   40.0e6,     128.0e9, 573.0,    398.0,  5.60,   387.0,  160.0e6,    0.350,  16.9) 
D6AC_Steel              =       Materials('D6AC Steel Alloy',                   7870.0,   1450.0e6,   210.0e9, 1180.0,    52.0,  1.12,   480.0,  2360.0e6,   0.295,  13.1)
default                 =       Materials('default_coating',                       0.0,       0.0,        0.0,   0.0,        1,   0.0,     0.0,     0.0,       0.0,   0.0)

#Coating Materials:
Coat_Custom             =       Materials('Carbon-Carbon Matrix coating',       1950.0,   0,          0,          2400.0,   37.4,   6.4,   2093.0,   0,   0,   18.0)
Copper                  =       Materials('Copper coating',                     8940.0,   0,          0,          573.0,    390.0,  6.4,   388.0,    0,   0,   18.0)
Narloy_Z                =       Materials('Copper Alloy Coating',               9130.0,   0,          0,          740.0,    290.0,  6.4,   388.0,    0,   0,   22.0)
GRCop_84                =       Materials('GRPCop_84 Copper alloy coating',     8756.0,   0,          0,          1000.0,   351.0,  6.4,   388.0,    0,   0,   20.0)
Silica                  =       Materials('Scilica Coating',                    1700.0,   0,          0,          1200.0,   0.55,   6.4,   1256.0,   0,   0,   20.0)
Carbon                  =       Materials('Carbon-Carbon Matrix coating',       1950.0,   0,          0,          2400.0,   37.4,   6.4,   2093.0,   0,   0,   20.0)


##Computing Mass nozzle: 
#Nozzle Surface Fucntion:
def Nozzle_mass(x,R,t,material):
    total_surf = 0
    # nozmass_error = 0
    # nozmass_warnings = 0
    
    # Warning 0: "The thickness is lower than 1mm, a default thickness of 1mm has been used"
    # Warning 1: "The material density is lower than 0"
    # Warning 2: "The nozzle mass is lower than 0
    
    #Example: error implementation
    # if Pc>300:
    #     if (warnings & (1<<0))==False:
    #         warnings=warnings|(1<<0)
    # elif warnings & (1<<0):
    #     warnings=warnings & (~(1<<0))


    # if t < 0.001:
    #     if (nozmass_warnings & (1<<0)) == False:
    #         nozmass_warnings=nozmass_warnings|(1<<0)
    # else:
    #     nozmass_warnings = nozmass_warnings & (~(1<<0))
    
    # if material.density < 0:
    #     if (nozmass_warnings & (1<<0)) == False:
    #         nozmass_warnings=nozmass_warnings|(1<<1)
    # else:
    #     nozmass_warnings = nozmass_warnings & (~(1<<1))
    #     return 0,nozmass_error, nozmass_warnings
     

    for i in range(len(x)-1):
        if t[i] < 0.001:
            total_surf += mth.dist([x[i+1],R[i+1]],[x[i],R[i]])*t[i]*R[i]
        else:
            total_surf += mth.dist([x[i+1],R[i+1]],[x[i],R[i]])*0.001*R[i]
    NozzleMass = total_surf*2*mth.pi*material.density

    # if NozzleMass < 0:
    #     nozmass_error=nozmass_error|(1<<0)
    #     return 0,nozmass_error,nozmass_warnings
    
    return NozzleMass 

class ReferenceEngine:
    def __init__(self, pc, thrust, arear, rt, mprop, rhoprop, FS, Material_NCG, Material_P, Material_V, mfrac_tube, mfrac_manifold, mfrac_jacket, mfrac_chamber, mfrac_pump, mfrac_valve, RefMass):
        self.pc = pc
        self.thrust = thrust
        self.arear = arear
        self.rt = rt
        self.mprop = mprop
        self.rhoprop = rhoprop
        self.FS = FS
        self.Material_NCG = Material_NCG
        self.Material_P = Material_P
        self.Material_V = Material_V
        self.mfrac_tube = mfrac_tube
        self.mfrac_manifold = mfrac_manifold
        self.mfrac_jacket = mfrac_jacket
        self.mfrac_chamber = mfrac_chamber
        self.mfrac_pump = mfrac_pump
        self.mfrac_valve = mfrac_valve
        self.RefMass = RefMass

RL10    = ReferenceEngine(3.20e6, 7.34e4, 61.1,  0.076, 16.85, 351.91, 1.1, Inc_718, D6AC_Steel, D6AC_Steel, 0.0903, 0.2575, 0.0484, 0.2281, 0.3096, 0.0661, 138) #Expander Cycle Reference
LE5     = ReferenceEngine(3.65e6, 1.03e5, 140.0, 0.068, 23.33, 343.83, 1.1, Inc_718, Al_7075_T6, Al_7075_T6, 0.1419, 0.0578, 0.1021, 0.2804, 0.3042, 0.1136, 255)#Gas Generator Cycle Reference #fix rhoprop
SSME    = ReferenceEngine(2.04e7, 2.28e6, 77.5,  0.138, 512.6, 361.89, 1.2, Inc_718, D6AC_Steel, D6AC_Steel, 0.0907, 0.1984, 0.0737, 0.1955, 0.2763, 0.1654, 3177)#Stage Combustion Cycle Reference #fixrhoprop




#Mass estimation function Nozzle Tubes:
def Mass(Pc, material_N, material_V, arear, rt, mprop, FS, rhoprop, cycle):
    
    Mass_warnings = 0
    # Mass_error = 0

    if Pc < 0:
        Mass_warnings=Mass_warnings|(1<<0)
    
    if material_N.density < 0:
        Mass_warnings=Mass_warnings|(1<<1)    
    
    if material_V.density < 0:
        Mass_warnings=Mass_warnings|(1<<2)
    
    if arear < 0:
        Mass_warnings=Mass_warnings|(1<<3)
    
    if rt < 0:
        Mass_warnings=Mass_warnings|(1<<4)
    
    if mprop < 0:
        Mass_warnings=Mass_warnings|(1<<5)
    
    if FS < 0:
        Mass_warnings=Mass_warnings|(1<<6)
    
    if rhoprop < 0:
        Mass_warnings=Mass_warnings|(1<<7)
  

    reference = RL10

    if cycle == 'SC':
        reference = SSME
    elif cycle == 'EX':
        reference == RL10
    else:
        reference == LE5


    #Cooling Mass:
    TubeMass = reference.mfrac_tube*(((Pc/reference.pc)**1)*((material_N.density/reference.Material_NCG.density)**1)*(((material_N.yieldstress_l/FS)/(reference.Material_NCG.yieldstress_l/reference.FS))**(-1))*((arear/reference.arear)**2)*((rt/reference.rt)**2)) #Dimensionless Nozzle Tube Mass
    ManifoldMass = reference.mfrac_manifold*(((Pc/reference.pc)**1)*((material_N.density/reference.Material_NCG.density)**1)*((mprop/reference.mprop)**1)*((material_N.yieldstress_l/reference.Material_NCG.yieldstress_l)**(-1))*((rhoprop/reference.rhoprop)**(-1))*((rt/reference.rt)**2)) #Dimensionless Nozzle Manifold Mass

    # if TubeMass + ManifoldMass < 0:
        # Mass_warnings=Mass_warnings|(1<<10)
        # return 0, Mass_error,Mass_warnings

    #TurboPump Mass
    #PumpMass = reference.mfrac_pump*(((Pc/reference.pc)**0.15)*((material_P.density/reference.material_P.density)**1)*(((material_P.yieldstress_l/FS)/(reference.Material_p.yieldstress_l/reference.FS)*(-1)))*((mprop/reference.mprop)**0.9)*((rhoprop/reference.rhoprop)**(-0.45))*((Ns/Ns0)**(-0.6)))
    PumpMass = reference.mfrac_pump*(((Pc/reference.pc)**0.53)*((mprop/reference.mprop)**0.53)*(rhoprop/reference.rhoprop)**(-0.53)) #Dimensionless mass of turbopump (historic data method)
    
    # if PumpMass < 0:
        # Mass_warnings=Mass_warnings|(1<<11)
        # return 0, Mass_error,Mass_warnings


    #Valves
    ValveMass = reference.mfrac_valve*(((Pc/reference.pc)**1.0)*((material_V.density/reference.Material_V.density)**1)*(((material_V.yieldstress_l/FS)/(reference.Material_V.yieldstress_l/reference.FS))**(-1))*((mprop/reference.mprop)**1)*((rhoprop/reference.rhoprop)**(-1)))#Dimensionless mass of Valves
    
    # if ValveMass < 0:
        # Mass_warnings=Mass_warnings|(1<<11)
        # return 0, Mass_error,Mass_warnings
    
    #Total:
    Total_Mass = (ManifoldMass + TubeMass + PumpMass + ValveMass)*reference.RefMass 
    return Total_Mass, Mass_warnings


def RhoProp(O_prop, F_prop, OF):
    rho_prop = 0
    rho_prop = ((O_prop*F_prop)*(1+OF))/(F_prop*OF+O_prop)
    return rho_prop


##Cost function 
def Cost(m_engine,f1,R,f3,n,a_m):
    TotalCost_MY = 0
    f2  = 0.27+1.075*R**46.994
    f4  = -0.0553*mth.log(n) + 1.0011
    C_D = 1.1*f1*f2*f3*162*m_engine**0.58
    F_E = a_m*f4*m_engine**0.46
    TotalCost_MY = C_D + F_E
    
    #a_m is correction factor 4 in the GUI
    
    return TotalCost_MY


def RhoProp(O_prop, F_prop, OF):
        rho_prop = 0
        rho_prop = ((O_prop*F_prop)*(1+OF))/(F_prop*OF+O_prop)
        return rho_prop

##Reuseability:
#Thinning of the coolant pipe wall after 1 cycle

def Reuseability(material, Twg_max, Twc_max, Twg_min, Twc_min, DT, H, Chanel_width, P_chamber, P_channel, mu, N):
    Ti_max = (Twg_max + Twc_max)/2
    T0_max = 0.35*Ti_max #The 0.35 comes from the assumption that the cooling channels are rectangular in shape. This needs to be changed for any other cross section. 
    
    Ti_min = (Twg_min + Twc_min)/2
    T0_min = 0.35*Ti_min 
    
    DT_ep_pl2 = DT
    l = Chanel_width
    w = l
    p = P_channel - P_chamber
    Reuseability = 'Not Reuseable'

    #Inelastic Strain:
    ep_pl1 = (material.k*(Ti_max - T0_max) - (Ti_min - T0_min)) - (2*material.yieldstress_l/material.Emod)
    ep_pl2 = ((material.k*DT_ep_pl2)**2)/(12*material.yieldstress_l*(1-mu)**2)
    ep_1 = 2*(ep_pl1+ep_pl2)

    #Deflection:
    def1 = 2*((H/ep_1) - mth.sqrt(((H/ep_1)**2) - ((l/4)**2)))
    def2 = (ep_1*p*l**2)/(4*H*material.yieldstress_l)
    def_tot = def1 + def2

    #Ligament Deformation:
    t_N = (N*w*def_tot)/(l+w)
    # t_N_min = (2*H*(l+w) - N*w*def_tot)/(l+w)
    #t_N_max = (2*H*(l+w)**2 + N*w*l*def_tot)/(l+w)**2

    # #Fatigue and Creep Rupture Damage: Funtion does not seem to work when creep is 
    q = 0.2*((160e6 - material.yieldstress_l)/material.yieldstress_l)**(0.6)
    # ep_c1avg = material.k*(T1 - T0)
    # ep_c1 = ep_c1avg*(((q-1)/q)*((t_N_min/t_N_max) - 1))*((t_N_min/t_N_max)**((q-1)/q) - 1)**(-1)
    # ep_c2 = ep_c1avg
    # ep_c_tot = (2/mth.sqrt(3))*mth.sqrt((ep_c1**2 + ep_c1*ep_c2 + ep_c2**2))
    
    #Critical Thickness:
    tcr = 2*H*(1-mth.e**(-q))

    if t_N > tcr:
        Reuseability = 'Not Reuseable after selected amount of cycles'
    else:
        Reuseability = 'Reuseable System'

    return Reuseability

x = 1,2,3
y = 3,4,5
t = 0.001,0.001,0.001

print(Mass(20e6, Rhenium, Rhenium, 60.1, 0.3, 50, 1.1, 300, 'SC'))