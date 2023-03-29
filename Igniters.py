#Igniter function uses the following inputs:
#m - mass flow
#Propellant - propellant being used and its properties
#default
#tc - chamber temperature
#of - oxidizer to fuel ratio
#type - type of igniter specified by the user
#The output of the igniters function consists of the igniter fuel mass, igniter oxidizer mass, igniter propellant total mass and warning.

#class Subs:
    #def __init__(self):
        #self.heatingvalues = '42 is the answer to the whole universe'
        #self.Compound = 0
        #self.mass = 0


def Enthalpy (Propellant,Tc,of):
    i = 0
    Hc_ox = (Propellant.o_nist_enthalpy_coef[0+i] * Tc + Propellant.o_nist_enthalpy_coef[1+i] * Tc ** 2 / 2
             + Propellant.o_nist_enthalpy_coef[2+i] * Tc ** 3 / 3 + Propellant.o_nist_enthalpy_coef[3+i] * Tc ** 4 / 4
             - Propellant.o_nist_enthalpy_coef[4+i] / Tc + Propellant.o_nist_enthalpy_coef[5+i]
             - Tc + Propellant.o_nist_enthalpy_coef[7+i])
    #print(Hc_ox)
    Hc_ox = Hc_ox / Propellant.ox_M

    Hc_fuel = (Propellant.f_nist_enthalpy_coef[0+i] * Tc + Propellant.f_nist_enthalpy_coef[1+i] * Tc ** 2 / 2
               + Propellant.f_nist_enthalpy_coef[2+i] * Tc ** 3 / 3 + Propellant.f_nist_enthalpy_coef[3+i] * Tc ** 4 / 4
               - Propellant.f_nist_enthalpy_coef[4+i] / Tc + Propellant.f_nist_enthalpy_coef[5+i]
               - Tc + Propellant.f_nist_enthalpy_coef[7+i])

    Hc_fuel=Hc_fuel/ Propellant.f_M
  #  print(Hc_fuel)
    if (Propellant.Fuel_name == "Ethanol"):
        if (Tc <= 2):
            Hc_fuel = 252.26
            Hc_fuel = Hc_fuel / Propellant.f_M

        elif (Tc <= 2.25):
            Hc_fuel = 252.26 + 180 * (Tc*1000-2000)/1000
            Hc_fuel = Hc_fuel / Propellant.f_M

        elif (Tc <= 2.5):
            Hc_fuel = 297.26 + 183 * (Tc*1000-2250)/1000
            Hc_fuel = Hc_fuel / Propellant.f_M

        elif (Tc <= 2.75):
            Hc_fuel = 343.01 + 186 * (Tc*1000-2500)/1000
            Hc_fuel = Hc_fuel / Propellant.f_M

        elif (Tc <= 3):
            Hc_fuel = 389.51 + 188 * (Tc*1000-2750)/1000
            Hc_fuel = Hc_fuel / Propellant.f_M

        else:
            Hc_fuel = 436.51 + 188 * (Tc * 1000 - 3000) / 1000
            Hc_fuel = Hc_fuel / Propellant.f_M

    #print("Hc_fuel: "+str(Hc_fuel))

    H = Hc_fuel*(1/(of+1))+Hc_ox*(1-1/(of+1))
  #  print(H)

    return H


def Igniters (m,Propellant,default,Tc,of,type):
    Tc = Tc / 1000
    wr = 0
    match type:
        case '00':
            Hc = Enthalpy(Propellant,Tc,of)*10**3
            Power = m * (Hc)
            heatingvalue = Propellant.heatingvalue*(1/(default.ign_o_f+1))
            mass_igniter = Power/heatingvalue * default.ignburntime /(default.fudgefactor)
            mfuel = mass_igniter/(1+default.ign_o_f)
            mox = default.ign_o_f *mass_igniter/ (1+default.ign_o_f)
            if mass_igniter > 100 and Propellant.Fuel_name == "LH2":
                wr = wr | (1<<1)
            elif mass_igniter > 100 and Propellant.Fuel_name != "LH2":
                wr = wr | (1<<0)
            return mass_igniter,mfuel,mox,wr
        case '10':
            Hc = Enthalpy(Propellant, Tc, of) * 10 ** 3
            Power = m * (Hc)
            heatingvalue = 50*10**6 * (1 / (default.ign_o_f + 1))
            mass_igniter = Power / heatingvalue * default.ignburntime / (default.fudgefactor)
            mfuel = mass_igniter / (1 + default.ign_o_f)
            mox = default.ign_o_f * mass_igniter / (1 + default.ign_o_f)
            if mass_igniter > 100:
               wr = wr | (1<<0)
            return mass_igniter, mfuel, mox,wr
        case '11':
            Hc = Enthalpy(Propellant, Tc, of) * 10 ** 3
            Power = m * (Hc)
            heatingvalue = 119.96*10**6 * (1 / (default.ign_o_f + 1))
            mass_igniter = Power / heatingvalue * default.ignburntime / (default.fudgefactor)
            mfuel = mass_igniter / (1 + default.ign_o_f)
            mox = default.ign_o_f * mass_igniter / (1 + default.ign_o_f)
            if mass_igniter > 100:
               wr = wr | (1<<1)
            return mass_igniter, mfuel, mox,wr
        case '20':
            Hc = Enthalpy(Propellant, Tc, of) * 10 ** 3
            Power = m * (Hc)
            heatingvalue = 2.5*10**6 #* (1 / (default.ign_o_f + 1))
            mass_igniter = Power / heatingvalue * default.ignburntime / (default.fudgefactor)
            mfuel = mass_igniter / (1 + default.ign_o_f)
            mox = default.ign_o_f * mass_igniter / (1 + default.ign_o_f)
            if mass_igniter > 100:
               wr = wr | (1<<0)
            return mass_igniter, mfuel, mox,wr
        case '30':
            Hc = Enthalpy(Propellant, Tc, of) * 10 ** 3
            Power = m * (Hc)
            heatingvalue = 22.6*10**6 * (1 / (default.ign_o_f + 1))
            mass_igniter = Power / heatingvalue * default.ignburntime / (default.fudgefactor)
            mfuel = mass_igniter / (1 + default.ign_o_f)
            mox = default.ign_o_f * mass_igniter / (1 + default.ign_o_f)
            if mass_igniter > 100:
               wr = wr | (1<<0)
            return mass_igniter, mfuel, mox,wr
        case '31':
            Hc = Enthalpy(Propellant, Tc, of) * 10 ** 3
            Power = m * (Hc)
            heatingvalue = 36.4*10**6 * (1 / (default.ign_o_f + 1))
            mass_igniter = Power / heatingvalue * default.ignburntime / (default.fudgefactor)
            mfuel = mass_igniter / (1 + default.ign_o_f)
            mox = default.ign_o_f * mass_igniter / (1 + default.ign_o_f)
            if mass_igniter > 100:
               wr = wr | (1<<0)
            return mass_igniter, mfuel, mox,wr
        case '32':
            Hc = Enthalpy(Propellant, Tc, of) * 10 ** 3
            Power = m * (Hc)
            heatingvalue = 28.326 * 10 ** 6 * (1 / (default.ign_o_f + 1))
            mass_igniter = Power / heatingvalue * default.ignburntime / (default.fudgefactor)
            mfuel = mass_igniter / (1 + default.ign_o_f)
            mox = default.ign_o_f * mass_igniter / (1 + default.ign_o_f)
            if mass_igniter > 100:
                wr = wr | (1 << 0)
            return mass_igniter, mfuel, mox, wr
        case '33':
            Hc = Enthalpy(Propellant, Tc, of) * 10 ** 3
            Power = m * (Hc)
            heatingvalue = 32.94 * 10 ** 6 * (1 / (default.ign_o_f + 1))
            mass_igniter = Power / heatingvalue * default.ignburntime / (default.fudgefactor)
            mfuel = mass_igniter / (1 + default.ign_o_f)
            mox = default.ign_o_f * mass_igniter / (1 + default.ign_o_f)
            if mass_igniter > 100:
                wr = wr | (1 << 0)
            return mass_igniter, mfuel, mox, wr

    #heatingvalues = [119.96*10**6*(1/(default.ign_o_f+1)),13.5*10**6,9.2*10**6,6.5*10**6,2.9*10**6,2.5*10**6]
    #heatingcompounds = ['hydrogenoxygen_LHV','methaneoxygen','magnesiumteflonviton','boronpotassiumnitratewax','blackpowder','hydrogenperoxide']
    #time = default.ignburntime  # Igniter burn time
    #Compound = [Subs() for i in range(n)]

    #for i in range(n):
     #   Compound[i].heatingvalues = heatingvalues[i]
      #  Compound[i].compounds = heatingcompounds[i]
       # Compound[i].mass = Power/Compound[i].heatingvalues/default.fudgefactor*time


    #return Compound[0].mass

#class Default:
 #   ignburntime = 3.5
  #  fudgefactor = 20
   # ign_o_f = 0.7

#class Propellant:
    #f_nist_enthalpy_coef = [121.4010,4.816880,-0.763012,0.043232,-40.78650,-11.38110,305.3440,95.35340]
    #f_nist_enthalpy_coef = [
            #42.31,
            #-4.393,
            #1.28243,
            #-0.106876,
            #-21.5339,
            #-37.5151,
            #163.08,
            #0,
            #34.066,
            #-11.263,
           # 12.4328,
         #   -2.773,
          #  -0.15856,
        #    -10.081,
       #     173.71,
      #      0,
     #   ]
    #o_nist_enthalpy_coef = [20.91,10.72,-2.02,0.1464,9.2457,5.338,237.62,0,
    #                        31.33,-20.235,57.87,-36.51,-0.007374,-8.9035,246.79,0]
    #ox_M = 32*10**(-3)
   # f_M = 46.7*10**(-3)
  #  heatingvalue = 119.96 * 10 ** 6
  #  enthalpy_298_ox = 0*10**6
 #   enthalpy_298_f = 5.01766 * 10**6
#    Fuel_name = "Ethanol"


#prop = Propellant
#default = Default
#Tc = 2250
#of = 6
#ofign = 0.7

# mox +mfuel = 467 mox/mfuel = 6 mox = mfuel *6, mfuel = 1/7
# mox+mfuel = m mox/mfuel = of mox = mfuel *of, mfuel = m/(1+of),

#mox + mfuel = m
#mox/mfuel = o/f_ign
#mox = o/f_ign * mfuel
#mfuel ( 1 + o/f_ign) = m
#mfuel = m/(1+o/f_ign)
#mox = o/f_ign*m / (1+o/f_ign)

#type = '00'
#ignitermass,mox,mfuel,wr_ign = Igniters(467,prop,default,Tc,of,type)
#y = mox+mfuel
#print(ignitermass,mox,mfuel,wr_ign)



