#inputs: pc, At, Propellant, Material, safety factor, velocity,d0,Tc,of,bool
#outputs: Thickness, Area of the Chamber, Conductive Heat transfer factor, Mass
import math

#Pc - Chamber Pressure
#At - Throat Area
#Propellant - Class with Propellant qualities
#Material - Class with Material qualities
#Safety factor - Factor to compute thickness/etc
#Velocity - Velocity of propellants after injector
#d0 - diameter of drops after injection
#Tc - Chamber Temperature
#of - oxidizer/fuel ratio
#bool - variable to say if this is on the pressure iteration loop


def CombustionChamber (Pc,At,Propellant,Material,Safety,velocity,d0,Tc,of,bool):

    Vi = 4/3*math.pi*(d0/2)**3
    Vf = Vi * 0.1 #random value for now
    d = (3/4*Vf)**(1/3)*2
    time_f = -(d-d0)/Propellant.lamb_f #dquadrado
    time_o = -(d-d0)/Propellant.lamb_o


    nfuel = 1/(of+1)
    nox = 1-nfuel

    gama = nox * Propellant.o_gama + nfuel * Propellant.f_gama
    GAMA = math.sqrt(gama)*(2/(gama+1))**((gama+1)/(2*gama-2))

    if time_f>time_o:
        time = time_f
    else:
        time = time_o

    LengthChamber = velocity/time
    lstar = Propellant.tq * GAMA * math.sqrt(propellant.f_gamma*Tc)



    Vchamber = lstar * At
    Achamber = Vchamber / LengthChamber
    Rchamber = (Achamber/math.pi)**(1/2)



    Thickness = Pc * Rchamber * Safety / Material.yieldstress

   if bool == 1:
        Mass = Achamber * Thickness * Material.density

    a = 0.023
    #weighted averages for the several Parameters
    ro = nox*Propellant.o_density + nfuel*Propellant.f_density
    Pr = 4*gama/(9*gama-5)
    cp = nox*Propellant.ocp + nfuel*Propellant.fcp
    niu = nox*Propellant.omiu + nox*Propellant.fmiu
    k = (Pr/(niu*cp))

    heattransfer = a * ro**0.8 * velocity**0.8 * (1/(Rchamber*2))**0.2 * (k*Pr**0.33/niu**0.8)
    if bool == 0:
        return (heattransfer,Achamber,Thickness)
    else:
        return (heattransfer,Achamber,Thickness,Mass)
