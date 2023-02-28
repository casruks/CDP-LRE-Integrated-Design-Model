import Turbomachinery_code as Turbo
import Injector as Inj

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

#Main Function
if __name__ == '__main__':
    print("Starting...")