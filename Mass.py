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
Ti6Al4V                 =       Materials('Tungsten & Aluminium alloy',         4428.0,   110.0e6,    114.0e9, 693.15,   6.7,  23.8)#, 409.0e6)
Haynes_188              =       Materials('Nickel Alloy',                       8980.0,   464.0e6,    244.0e9, 1423.0,   10.8, 30.6)#, 0)
Inconel_X_750           =       Materials('Inconel-750, Nickel-Chromium Alloy', 8280.0,   850.0e6,    222.0e9, 1090.0,   12.0,   20.4)#, 672.0e6)
Inconel_600             =       Materials('Inconel-600, Nickel-Chromium Alloy', 8470.0,   290.0e6,    206.0e9, 982.0,    15.9, 20.2)#, 640.0e6)
Inconel_718             =       Materials('Inconel-718, Nickel-Chromium Alloy', 8190.0,   1036.0e6,   211.0e9, 825.0,    11.4, 16.6)#, 774.0e6)
Inconel_A_286           =       Materials('A-286_Nickel-Chromium Alloy',        7920.0,   7920.0e6,   201.0e9, 700.0,    23.9, 5.43)#, 0)
Columbium_c103          =       Materials('Niobium (Colombium) - cold rolled',  8600.0,   550.0e6,    130.0e9, 1255.0,   0.54, 225.0)#, 109e6)
Copper_structural       =       Materials('Copper',                             8940.0,   40.0e6,     128.0e9, 573.0,    398.0, 5.60)#, 104e6)   

#Coating Materials:
Copper                  =       Materials('Copper coating',                     8940.0,   0,          0,          573.0,    390.0,  6.4)#, 0)
Narloy_Z                =       Materials('Copper Alloy Coating',               9130.0,   0,          0,          740.0,    290.0,  6.4)#, 0)
GRCop_84                =       Materials('GRPCop_84 Copper alloy coating',     8756.0,   0,          0,          1000.0,   351.0,  6.4)#, 0)
Silica                  =       Materials('Scilica Coating',                    1700.0,   0,          0,          1200.0,   0.55, 6.4)#, 0)
Carbon                  =       Materials('Carbon-Carbon Matrix coating',       1950.0,   0,          0,          2400.0,   37.4, 6.4)#, 0)


##Computing Mass nozzle: 
#Nozzle Surface Fucntion:
def Nozzle_surface(x,R):
    total_dist = 0
    for i in range(len(x)-1):
        total_dist += mth.dist([x[i+1],R[i+1]],[x[i],R[i]])
    vol = total_dist*2*mth.pi
    return vol

##Expander Cycle:
#Pressure Ratio             = (Pc/(3.20e3))
#Density Ratio Inconel 718  = (material.density/Inconel_718.density)
#Density Ratio Al_7075      = (material.density/Aluminium_7075_T6.density)
#Stress Ratio Inconel 718   = ((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/1.1)
#Stress Ratio Al_7075       = ((material.yieldstress_l/FS)/(Aluminium_7075_T6.yieldstress_l/1.1)
#AreaRatio Ratio            = (arear/61.1)
#ThroatArea Ratio           = (rt/0.076)
#Propelant mass Ratio       = (mprop/16.85)
#Propelant Density Ratio    = (rhoprop/x)



#Mass estimation function Nozzle Tubes:
def Mass_Expander(Pc, material, arear, rt, mprop, FS, rhoprop,Ns):
    x=5
    y=5
    #Nozzle Mass
    TubeMass = ((Pc/(3.20e3))**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/1.1))**(-1))*((arear/61.1)**1)*((rt/0.076)**2) #Dimensionless Nozzle Tube Mass
    ManifoldMass = ((Pc/(3.20e3))**1)*((material.density/Inconel_718.density)**1)((mprop/16.85)**1)((material.yieldstress_l/Inconel_718.yieldstress_l)**(-1))*((rhoprop/x)**(-1))*((rt/0.076)**1) #Dimensionless Nozzle Manifold Mass
    JacketMass = ((Pc/3.20e3)**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/1.1))**(-1))*((arear/61.1)**1.5)*((rt/0.076)**3) #Dimensionless Nozzle Jacket Mass
    
    #Chamber Mass
    ChamberMass = ((Pc/(3.20e3))**1)*((material.density/Inconel_718.density)**1)*(((material.yieldstress_l/FS)/(Inconel_718.yieldstress_l/1.1)**(-1)))*((rt/0.076)**1) 
    
    #TurboPump Mass
    PumpMass = ((Pc/(3.20e3))**0.15)*((material.density/Aluminium_7075_T6.density)**1)*(((material.yieldstress_l/FS)/(Aluminium_7075_T6.yieldstress_l/1.1)*(-1)))*((mprop/16.85)**0.9)*((rhoprop/x)**(-0.45))*((Ns/y)**(-0.6))
    
    #Valves
    ValveMass = ((Pc/(3.20e3))**1)*((material.density/Aluminium_7075_T6.density)**1)*(((material.yieldstress_l/FS)/(Aluminium_7075_T6.yieldstress_l/1.1)**(-1)))*((mprop/16.85)**1)*((rhoprop/x)**(-1))
  
    return TubeMass, ManifoldMass, JacketMass, ChamberMass, PumpMass, ValveMass



