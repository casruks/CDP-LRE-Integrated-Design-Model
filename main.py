import Turbomachinery_code as Turbo
import CombustionChamber as Comb
import Injector as Inj
import Igniters as Ign

#Default values
class Default:
    #Tolerances
    pres_tol = 0.01

    #Seeds
    Pres = 1e6
    inj_vel = 15

    #Injectors
    Cd = 1.0

    #Nozzle


    #Turbomachinery
    Eff_t = 0.5
    Eff_p = 0.5
    Eff_m = 0.95

    #Combustion chamber
    SF = 1.0


    #Cooling


    #Material
    material = "This"

    #Init
    def __init__(self,type):
        if(type==1):
            self.Cd = 1.0

default = Default(0)

#Propellant class
class Propellant:
    #Oxidizer
    o_name = "LOX"
    o_dens = 1141.0

    #Fuel
    f_name = "LH"
    f_dens_l = 71.0
    f_dens_g = 1.0
    f_gamma = 1.4
    fcp = 14307.0
    R_f = 4.1573

    #Propellant
    lamb = 0.0
    tq = 0.0

    def __init__(self,type):
        if(type==1):
            ox_dens=1141.0

prop = Propellant(0)

#Main Function
if __name__ == '__main__':
    Thrust = 15000 #= input("Introduce thrust")
    Thrust_time = 30 #= input("Introduce thrust time")

    p_new = 0.0
    p_old = default.Pres
    inj_vel = default.inj_vel
    while abs(p_new-p_old)/p_old > default.pres_tol:
        p_new = p_old
        #Compute nozzle
        At, Ae = nozzle()

        #Compute chamber
        h_comb = Comb.CombustionChamber(p_new, At, prop, default.material, default.SF, inj_vel, D0)

        #Compute regenerative


        #Compute Turbo
        Turbo.TurboM()

        #Cmpute Injector



    #Compute Ignitor

    #Compute Masses

    #Compute reliability

    #Compute costs

    print("Starting...")