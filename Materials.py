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
Rhenium2 = Materials('Rhenium2', 21000, 1060*10**6, 2300*10**6,471*10**9, 0, 997, 938)

a = [Rhenium, Rhenium2]


pressure = 25

for i in  a:
    if i.yieldstress > pressure:
        a = i.material
        print(a)
