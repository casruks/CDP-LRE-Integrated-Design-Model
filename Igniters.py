#Igniter function uses the following inputs:
#m - mass flow
#Hc - enthalpy of propellants at the end of the combustion chamber
#H0 - enthalpy of propellants at the beggining of the combustion chamber
#The output of the igniters function consists of a vector of the Subs class (this is defined right after this comment). For each element of the vector,
#a different combination of substances is used to compute the mass needed for the igniter

class Subs:
    def __init__(self):
        self.heatingvalues = '42 is the answer to the whole universe'
        self.Compound = 0
        self.mass = 0
def Enthalpy (Propellant,Tc,of)
    Hc_ox = (Propellant.o_nist_enthalpy_coef[0] * Tc + Propellant.o_nist_enthalpy_coef[1] * Tc ** 2 / 2
             + Propellant.o_nist_enthalpy_coef[2] * Tc ** 3 / 3 + Propellant.o_nist_enthalpy_coef[3] * Tc ** 4 / 4
             - Propellant.o_nist_enthalpy_coef[4] / Tc + Propellant.o_nist_enthalpy_coef[5]
             - Tc + Propellant.o_nist_enthalpy_coef[7])
    Hc_fuel = (Propellant.f_nist_enthalpy_coef[0] * Tc + Propellant.f_nist_enthalpy_coef[1] * Tc ** 2 / 2
               + Propellant.f_nist_enthalpy_coef[2] * Tc ** 3 / 3 + Propellant.f_nist_enthalpy_coef[3] * Tc ** 4 / 4
               - Propellant.f_nist_enthalpy_coef[4] / Tc + Propellant.f_nist_enthalpy_coef[5]
               - Tc + Propellant.f_nist_enthalpy_coef[7])
    #nfuel+nox = 1 Nox/Nfuel=of nfuel = 1/(of+1) + nox = of/(of+1)
    H = Hc_ox*1/(of+1)+ Hc_fuel*of/(of+1)
    return H

def Igniters (m,Propellant,default,Tc,of):

    Hc = Enthalpy(Propellant,Tc,of)
    H0 = Enthalpy(Propellant,298,of)

    Power = m *(Hc-H0)

    heatingvalues = [19.8*10**6,13.5*10**6,9.2*10**6,6.5*10**6,2.9*10**6,2.5*10**6]
    heatingcompounds = ['hydrogenoxygen','methaneoxygen','magnesiumteflonviton','boronpotassiumnitratewax','blackpowder','hydrogenperoxide']
    time = default.ignburntime  # Igniter burn time
    Compound = [Subs() for i in range(n)]

    for i in range(n):
        Compound[i].heatingvalues = heatingvalues[i]
        Compound[i].compounds = heatingcompounds[i]
        Compound[i].mass = Power/Compound[i].heatingvalues*time/14


    return Compound[0]


#n = 6
#nH2 = 1/7
#n0 = 6/7
#ho = 0 * nH2 + -6*10**3 *n0 * 32*10**(-3)
#hc = 110.1*10**3 * n0/(32*10**(-3)) + 99.97*10**3* nH2 /( 2*10**(-3))
#Compound = Igniters(467,hc,ho)
#for i in range(n):
#   print(Compound[i].compounds,Compound[i].mass)