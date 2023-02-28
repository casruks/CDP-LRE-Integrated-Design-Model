class Subs:
    def __init__(self):
        self.heatingvalues = '42 is the answer to the whole universe'
        self.Compound = 0
        self.mass = 0



def Igniters (m,Hc,H0):

    Power = m *(Hc-H0)

    heatingvalues = [19.8*10**6,13.5*10**6,9.2*10**6,6.5*10**6,2.9*10**6,2.5*10**6]
    heatingcompounds = ['hydrogenoxygen','methaneoxygen','magnesiumteflonviton','boronpotassiumnitratewax','blackpowder','hydrogenperoxide']
    time = default.ignburntime  # Igniter burn time
    Compound = [Subs() for i in range(n)]

    for i in range(n):
        Compound[i].heatingvalues = heatingvalues[i]
        Compound[i].compounds = heatingcompounds[i]
        Compound[i].mass = Power/Compound[i].heatingvalues*time/14


    return Compound

#n = 6
#nH2 = 1/7
#n0 = 6/7
#ho = 0 * nH2 + -6*10**3 *n0 * 32*10**(-3)
#hc = 110.1*10**3 * n0/(32*10**(-3)) + 99.97*10**3* nH2 /( 2*10**(-3))
#Compound = Igniters(467,hc,ho)
#for i in range(n):
#   print(Compound[i].compounds,Compound[i].mass)