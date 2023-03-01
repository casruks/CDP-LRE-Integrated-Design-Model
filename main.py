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

    #Igniters
    ignburntime = 4

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
    ocp = 14307.0
    o_lamb = 0.0
   
    #Fuel
    f_name = "LH"
    f_dens_l = 71.0
    f_dens_g = 1.0
    f_gamma = 1.4
    fcp = 14307.0
    R_f = 4.1573
    f_lamb = 0.0
    
    #Propellant
    gama = 1.4
    tq = 0.0

    def __init__(self,type):
        match type:
            case 0:
                f_name = "LH"
                o_name = "LOX"

        
            case 1:
                f_name = "CH4"

prop = Propellant(0)
bool = 0 #this variable is used to show the combustor function we are in the first loop
#Main Function
if __name__ == '__main__':
    Thrust = 15000 #= input("Introduce thrust")
    Thrust_time = 30 #= input("Introduce thrust time")

    p_new = 0.0
    p_old = default.Pres
    inj_vel = default.inj_vel
    while abs(p_new-p_old)/p_old > default.pres_tol:
        p_new = p_old
        #Compute nozzle (1)
        At, Ae = nozzle()

        #COmpute injector (1)

        #Compute chamber - needs Chamber temperature + oxider to fuel ratio from previous functions (Tc and of)
        h_comb, Achamber, ThicknessChamber = Comb.CombustionChamber(p_new, At, prop, default.material, default.SF, inj_vel, D0,Tc,of,bool)

        #COmpute nozzle (2)

        #Compute regenerative


        #Compute Turbo
        Turbo.TurboM()

        #Cmpute Injector (2)


    bool = 1 #Shows the combustor it is out of the loop in order to compute mass!
    #Compute Ignitor - m is the mass flow, Hc is enthalpy of propelants at chamber exit, H0 is enthalpy of propelants at chamber entry
    #For further information on igniter output, see comments on first line of the igniters functions
    igniter_results = Igniters(m,Hc,H0)
    #Compute Masses

    #Compute reliability

    #Compute costs

    print("Starting...")