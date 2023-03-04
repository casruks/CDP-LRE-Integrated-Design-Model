import math as mth


class Materials:
    def __init__(self, material, density, yieldstress_l, yieldstress_u, Emod, OpTemp_l, OpTemp_u, cost,k):
        self.material = material
        self.density = density
        self.yieldstress_l = yieldstress_l
        self.yieldstress_u = yieldstress_u
        self.Emod = Emod
        self.OpTemp_l = OpTemp_l
        self.OpTemp_u = OpTemp_u
        self.cost = cost
        self.k=k

Rhenium = Materials('Rhenium', 21000, 1060*10**6, 2300*10**6,471*10**9, 0, 997, 938,39600)
Rhenium2 = Materials('Rhenium2', 21000, 1060*10**6, 2300*10**6,471*10**9, 0, 997, 938,39600)

def Material_Select(pressure):
    a = [Rhenium, Rhenium2]
    for i in a:
        if i.yieldstress_l > pressure:
            a = i.material
    return a


#Computing Mass: 
x = [(1,2),(3,4),(5,6)]
 
def Mass(x,t,material):
    total_dist = 0
    for i in range(len(x)-1):
        total_dist += mth.dist(x[i+1],x[i])
    
    vol = total_dist*t*2*mth.pi
    Mass = vol*material.density
    return Mass


Mass(x,0.1,Rhenium)
Material_Select(25)





    
