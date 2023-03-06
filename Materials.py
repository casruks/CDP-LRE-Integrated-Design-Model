import math as mth

class Materials:
    def __init__(self, material, density, yieldstress, Emod, OpTemp_u, Conductivity, cost):
        self.material = material
        self.density = density
        self.yieldstress = yieldstress
        self.Emod = Emod
        self.OpTemp_u = OpTemp_u
        self.Conductivity = Conductivity
        self.cost = cost

Rhenium                 =       Materials('Rhenium',                            21000,  2300*10**6, 471*10**9, 1270,   48,   938)
Aluminium_7075_T6       =       Materials('Aluminium 7075 T6',                  2810,   505*10**6,  72*10**9,  373.15, 130,  3.97)
Ti6Al4V                 =       Materials('Tungsten & Aluminium alloy',         4428,   110*10**6,  114*10**9, 693.15, 6.7,  23.8)
Haynes_188              =       Materials('Nickel Alloy',                       8980,   464*10**6,  244*10**9, 1423,   10.8, 30.6)
Inconel_X_750           =       Materials('Inconel-750, Nickel-Chromium Alloy', 8280,   850*10**6,  222*10**9, 1090,   12,   20.4)
Inconel_600             =       Materials('Inconel-600, Nickel-Chromium Alloy', 8470,   290*10**6,  206*10**9, 982,    15.9, 20.2)
Inconel_718             =       Materials('Inconel-718, Nickel-Chromium Alloy', 8190,   1036*10**6, 211*10**9, 825,    11.4, 16.6)
Inconel_A_286           =       Materials('A-286, Nickel-Chromium Alloy',       7920,   7920*10**6, 201*10**9, 700,    23.9, 5.43)
Columbium_c103          =       Materials('Niobium (Colombium) - cold rolled',  8600,   550*10**6,  130*10**9, 1255,   0.54, 225)

# #Coating Materials:
# Copper
# Narloy_Z
# GRCop_84
# Silica
# Alumina




        
##Selecting Materials:
def Material_Select(pressure,safety,temp):
    Select = [Rhenium,Aluminium_7075_T6,Ti6Al4V,Haynes_188,Inconel_X_750,Inconel_600,Inconel_718,Inconel_A_286,Columbium_c103]
    o = []
    for i in Select:
        if i.yieldstress/safety > pressure and i.OpTemp_u > temp:
            o.append((i.material))
        else:
            print('No Suitable Material')
    return o

##Computing Mass nozzle: 
#Form Coordinates (integration of nozzle & material)
def coordinate(a,b):
    if len(a) != len(b):
        raise Exception("#x-coordinates not equal to #y-coordiates for nozzle") 
    coordinate = []
    for i in range(len(a)):
        coordinate.append((a[i],b[i]))
    return coordinate

#Mass estimation function: 
def Mass(x,t,material):
    total_dist = 0
    for i in range(len(x)-1):
        total_dist += mth.dist(x[i+1],x[i])
    vol = total_dist*t*2*mth.pi
    Mass = vol*material.density
    return Mass




print(Material_Select(3000*10**6,2,3000))


