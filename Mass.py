import math as mth

class Materials:
    def __init__(self, material, density, yieldstress_l, Emod, OpTemp_u, k, cost):
        self.material = material
        self.density = density
        self.yieldstress_l = yieldstress_l
        self.Emod = Emod
        self.OpTemp_u = OpTemp_u
        self.k = k
        self.cost = cost
        #self.reusability = reusability

Rhenium                 =       Materials('Rhenium',                            21000.0,  2300.0e6,   471.0e9, 2200.0,   48.0,   938.0)#, 1190.0e6)
Aluminium_7075_T6       =       Materials('Aluminium 7075 T6',                  2810.0,   505.0e6,    72.0e9,  373.15,   130.0,  3.97)#, 129.0e6)
Ti6Al4V                 =       Materials('Tungsten & Aluminium alloy',         4428.0,   110.0e6,    114.0e9, 693.15,   6.7,    23.8)#, 409.0e6)
Haynes_188              =       Materials('Nickel Alloy',                       8980.0,   464.0e6,    244.0e9, 1423.0,   10.8,   30.6)#, 0)
Inconel_X_750           =       Materials('Inconel-750, Nickel-Chromium Alloy', 8280.0,   850.0e6,    222.0e9, 1090.0,   12.0,   20.4)#, 672.0e6)
Inconel_600             =       Materials('Inconel-600, Nickel-Chromium Alloy', 8470.0,   290.0e6,    206.0e9, 982.0,    15.9,   20.2)#, 640.0e6)
Inconel_718             =       Materials('Inconel-718, Nickel-Chromium Alloy', 8190.0,   1036.0e6,   211.0e9, 825.0,    11.4,   16.6)#, 774.0e6)
Inconel_A_286           =       Materials('A-286_Nickel-Chromium Alloy',        7920.0,   7920.0e6,   201.0e9, 700.0,    23.9,   5.43)#, 0)
Columbium_c103          =       Materials('Niobium (Colombium) - cold rolled',  8600.0,   550.0e6,    130.0e9, 1255.0,   0.54,   225.0)#, 109e6)
Copper_structural       =       Materials('Copper',                             8940.0,   40.0e6,     128.0e9, 573.0,    398.0,  5.60)#, 104e6)   

#Coating Materials:
Copper                  =       Materials('Copper coating',                     8940.0,   0,          0,          573.0,    390.0,  6.4)#, 0)
Narloy_Z                =       Materials('Copper Alloy Coating',               9130.0,   0,          0,          740.0,    290.0,  6.4)#, 0)
GRCop_84                =       Materials('GRPCop_84 Copper alloy coating',     8756.0,   0,          0,          1000.0,   351.0,  6.4)#, 0)
Silica                  =       Materials('Scilica Coating',                    1700.0,   0,          0,          1200.0,   0.55,   6.4)#, 0)
Carbon                  =       Materials('Carbon-Carbon Matrix coating',       1950.0,   0,          0,          2400.0,   37.4,   6.4)#, 0)


##Computing Mass nozzle: 
#Nozzle Surface Fucntion:
def Nozzle_surface(x,R):
    total_dist = 0
    for i in range(len(x)-1):
        total_dist += mth.dist([x[i+1],R[i+1]],[x[i],R[i]])
    vol = total_dist*2*mth.pi
    return vol

class ReferenceEngine:
    def __init__(self, pc, thrust, arear, rt, mprop, rhoprop, FS, TotalMass):
        self.pc = pc
        self.thrust = thrust
        self.arear = arear
        self.rt = rt
        self.mprop = mprop
        self.rhoprop = rhoprop
        self.FS = FS
        self.TotalMass = TotalMass

RL10    = ReferenceEngine(3.20e6, 7.34e4, 61.1,  0.076, 16.85, 351.91, 1.1, 138) #Expander Cycle Reference
LE5     = ReferenceEngine(3.65e6, 1.03e5, 140.0, 0.068, 23.33, 5.5, 1.1, 255)#Gas Generator Cycle Reference #fix rhoprop
SSME    = ReferenceEngine(2.04e7, 2.28e6, 77.5,  0.138, 512.6, 6, 1.2, 3177)#Stage Combustion Cycle Reference #fixrhoprop


##Cycle:
#Pressure Ratio             = (Pc/reference.pc)
#Density Ratio Inconel 718  = (material.density/Inconel_718.density)
#Density Ratio Al_7075      = (material.density/Aluminium_7075_T6.density)
#Stress Ratio Inconel 718   = ((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/RL10.FS)
#Stress Ratio Al_7075       = ((material.yieldstress_l/FS)/(Aluminium_7075_T6.yieldstress_l/RL10.FS)
#AreaRatio Ratio            = (arear/RL10.arear)
#ThroatArea Ratio           = (rt/RL10.rt)
#Propelant mass Ratio       = (mprop/RL10.mprop)
#Propelant Density Ratio    = (rhoprop/x)



#Mass estimation function Nozzle Tubes:
# def Mass_PressureFed(Pc, material, arear, rt, mprop, FS, rhoprop,Ns):
#     x=5
#     y=5
#     #Nozzle Mass
#     TubeMass = 0.0629*(((Pc/RL10.pc)**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/RL10.FS))**(-1))*((arear/RL10.arear)**1)*((rt/RL10.rt)**2)) #Dimensionless Nozzle Tube Mass
#     ManifoldMass = 0.2301*(((Pc/RL10.pc)**1)*((material.density/Inconel_718.density)**1)((mprop/RL10.mprop)**1)((material.yieldstress_l/Inconel_718.yieldstress_l)**(-1))*((rhoprop/x)**(-1))*((rt/RL10.rt)**1)) #Dimensionless Nozzle Manifold Mass
#     JacketMass = 0.0210*(((Pc/RL10.pc)**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/RL10.FS))**(-1))*((arear/RL10.arear)**1.5)*((rt/RL10.rt)**3)) #Dimensionless Nozzle Jacket Mass
    
#     #Chamber Mass
#     ChamberMass = 0.2007*(((Pc/RL10.pc)**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/RL10.FS)**(-1)))*((rt/RL10.rt)**1))
    
    
#     #Valves
#     ValveMass = (0.0202+0.185)*(((Pc/RL10.pc)**1)*((material.density/Aluminium_7075_T6.density)**1)*(((material.yieldstress_l/FS)/(Aluminium_7075_T6.yieldstress_l/RL10.FS)**(-1)))*((mprop/RL10.mprop)**1)*((rhoprop/x)**(-1)))
  
#     return TubeMass, ManifoldMass, JacketMass, ChamberMass, ValveMass


#Mass estimation function Nozzle Tubes:
def Mass_Expander_Regenerative(Pc, material, arear, rt, mprop, FS, rhoprop): #include Ns into input if design method is chosen
    x = reference.rhoprop
    #y= 

    reference = RL10
    mfrac_tube = 0.0629+0.0274
    mfrac_manifold = 0.2301+0.0274
    mfrac_jacket = 0.0210+0.0274
    mfrac_chamber = 0.2007+0.0274
    mfrac_pump = 0.1411+0.1411+0.0274
    mfrac_Valve = 0.0202+0.0185+0.0274

    #Nozzle Mass
    TubeMass = mfrac_tube*(((Pc/reference.pc)**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/reference.FS))**(-1))*((arear/reference.arear)**1)*((rt/reference.rt)**2))#*0.0629 #Dimensionless Nozzle Tube Mass
    ManifoldMass = mfrac_manifold*(((Pc/reference.pc)**1)*((material.density/Inconel_718.density)**1)((mprop/reference.mprop)**1)((material.yieldstress_l/Inconel_718.yieldstress_l)**(-1))*((rhoprop/x)**(-1))*((rt/reference.rt)**1))#*0.2301 #Dimensionless Nozzle Manifold Mass
    JacketMass = mfrac_jacket*(((Pc/reference.pc)**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/reference.FS))**(-1))*((arear/reference.arear)**1.5)*((rt/reference.rt)**3))#*0.0210 #Dimensionless Nozzle Jacket Mass
    
    #Chamber Mass
    ChamberMass = mfrac_chamber*((Pc/reference.pc)**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/reference.FS)**(-1)))*((rt/reference.rt)**1)#*(0.1411+0.1411) 
    
    #TurboPump Mass
    #PumpMass = mfrac_pump*(((Pc/reference.pc)**0.15)*((material.density/Aluminium_7075_T6.density)**1)*(((material.yieldstress_l/FS)/(Aluminium_7075_T6.yieldstress_l/reference.FS)*(-1)))*((mprop/reference.mprop)**0.9)*((rhoprop/x)**(-0.45))*((Ns/y)**(-0.6)))#*(0.1411+0.1411)
    PumpMass = mfrac_pump*(((Pc/reference.pc)**0.53)*((mprop/reference.mprop)**0.53)*(rhoprop/reference.rhoprop)**(-0.53))

    #Valves
    ValveMass = mfrac_Valve*(((Pc/reference.pc)**1)*((material.density/Aluminium_7075_T6.density)**1)*(((material.yieldstress_l/FS)/(Aluminium_7075_T6.yieldstress_l/reference.FS)**(-1)))*((mprop/reference.mprop)**1)*((rhoprop/x)**(-1)))#*(0.0202+0.185)
    

    #Total:
    Total_Mass_Expander = (TubeMass + ManifoldMass + JacketMass + ChamberMass + PumpMass + ValveMass)*reference.TotalMass 

    return Total_Mass_Expander

#Mass estimation function Nozzle Tubes:
def Mass_GasGenerator_Regenerative(Pc, material, arear, rt, mprop, FS, rhoprop): #include Ns into input if design method is chosen
    x = reference.rhoprop
    #y= 

    reference = SSME
    mfrac_tube = 0.0640+0.0267
    mfrac_manifold = 0.1717+0.0267
    mfrac_jacket = 0.0470+0.0267
    mfrac_chamber = 0.1137+0.0551+0.0267 #2nd massfraction of gas generator
    mfrac_pump = 0.1089+0.1407+0.0267
    mfrac_Valve = 0.0698+0.069+0.0266

    #Nozzle Mass
    TubeMass = mfrac_tube*(((Pc/reference.pc)**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/reference.FS))**(-1))*((arear/reference.arear)**1)*((rt/reference.rt)**2))#*0.0629 #Dimensionless Nozzle Tube Mass
    ManifoldMass = mfrac_manifold*(((Pc/reference.pc)**1)*((material.density/Inconel_718.density)**1)((mprop/reference.mprop)**1)((material.yieldstress_l/Inconel_718.yieldstress_l)**(-1))*((rhoprop/x)**(-1))*((rt/reference.rt)**1))#*0.2301 #Dimensionless Nozzle Manifold Mass
    JacketMass = mfrac_jacket*(((Pc/reference.pc)**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/reference.FS))**(-1))*((arear/reference.arear)**1.5)*((rt/reference.rt)**3))#*0.0210 #Dimensionless Nozzle Jacket Mass
    
    #Chamber Mass
    ChamberMass = mfrac_chamber*((Pc/reference.pc)**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/reference.FS)**(-1)))*((rt/reference.rt)**1)#*(0.1411+0.1411) 
    
    #TurboPump Mass
    #PumpMass = mfrac_pump*(((Pc/reference.pc)**0.15)*((material.density/Aluminium_7075_T6.density)**1)*(((material.yieldstress_l/FS)/(Aluminium_7075_T6.yieldstress_l/reference.FS)*(-1)))*((mprop/reference.mprop)**0.9)*((rhoprop/x)**(-0.45))*((Ns/y)**(-0.6)))#*(0.1411+0.1411)
    PumpMass = mfrac_pump*(((Pc/reference.pc)**0.53)*((mprop/reference.mprop)**0.53)*(rhoprop/reference.rhoprop)**(-0.53))

    #Valves
    ValveMass = mfrac_Valve*(((Pc/reference.pc)**1)*((material.density/Aluminium_7075_T6.density)**1)*(((material.yieldstress_l/FS)/(Aluminium_7075_T6.yieldstress_l/reference.FS)**(-1)))*((mprop/reference.mprop)**1)*((rhoprop/x)**(-1)))#*(0.0202+0.185)
    

    #Total:
    Total_Mass_GasGenerator = (TubeMass + ManifoldMass + JacketMass + ChamberMass + PumpMass + ValveMass)*reference.TotalMass 

    return Total_Mass_GasGenerator

##Cost function 
def Cost(m_engine, R, n):
    TotalCost = 0
    f1  = 1.0
    f2  = 0.27+1.075*R**46.994
    f3  = 1
    f4  = -0.0553*mth.log(n) + 1.0011
    a_m = 4.0

    C_D = 1.1*f1*f2*f3*162*m_engine**0.58

    F_E = a_m*f4*m_engine**0.46

    TotalCost = C_D + F_E
    return TotalCost


#print(mfrac_chamber+mfrac_jacket+mfrac_manifold+mfrac_tube+mfrac_pump+mfrac_Valve)