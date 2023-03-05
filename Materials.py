import math as mth

class Materials:
    def __init__(self, material, density, yieldstress_l, yieldstress_u, Emod, OpTemp_l, OpTemp_u, cost):
        self.material = material
        self.density = density
        self.yieldstress_l = yieldstress_l
        self.yieldstress_u = yieldstress_u
        self.Emod = Emod
        self.OpTemp_l = OpTemp_l
        self.OpTemp_u = OpTemp_u
        self.cost = cost

Rhenium = Materials('Rhenium', 21000, 1060*10**6, 2300*10**6,471*10**9, 0, 997, 938)
Rhenium2 = Materials('Rhenium2', 21000, 1070*10**6, 2300*10**6,471*10**9, 0, 997, 938)


        
##Selecting Materials:
def Material_Select(pressure):
    Select = [Rhenium, Rhenium2]
    o = []
    for i in Select:
        if i.yieldstress_l > pressure:
            o.append((i.material))
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




print(Material_Select(25))

