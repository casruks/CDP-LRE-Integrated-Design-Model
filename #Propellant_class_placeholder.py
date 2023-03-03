#Propellant class
class Propellant:
    #Oxidizer
    Ox_name = "LOX" #Oxidizer name for rocketCEA
    o_dens = 1141.0 #Oxidizer density
    ocp = 14307.0 #oxidizer cp
    o_lamb = 0.0
   
    #Fuel
    Fuel_name = "LH" #Fuel name for rocketCEA
    f_dens_l = 71.0 #liquid fuel density
    f_dens_g = 1.0 #gaseous fuel density
    f_gamma = 1.4 #fuel gamma
    fcp = 14307.0 #fuel cp
    R_f = 4.1573 #fuel gas constant
    f_lamb = 0.0
    
    Frozen_state=0
    
    #Propellant
    gama = 1.4
    tq = 0.0 #characteristic chemical time of propellant

    def __init__(self,type):
        match type:
            case 0:
                f_name = "LH"
                o_name = "LOX"

        
            case 1:
                f_name = "CH4"

prop = Propellant(0)