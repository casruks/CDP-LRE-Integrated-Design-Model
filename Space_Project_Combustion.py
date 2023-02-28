#inputs: pc, At, type of propellant, material, safety factor
#outputs: Length - check, thickness - check, Vchamber - check, h (heatflow) - check, Achamber - check
import math

#d = d0 - lambda *t
#define d0 as something that comes from injectors
#vinj = sqrt(1/ksi)*sqrt(2*deltap/ro) -> need delta P and ksi from inject type
# considering 90% oxidizer vaporized in order to compute t
# t + velocity - > chamber length
# Lstar = Vc / At - > Vc
# chamber length + Vc - > Achamber
#h = a * ro^0.8 * v^0.8 * (1/Dc)^0.2 * (k*Pr^0.33/niu^0.8)
#Pr = 4*gama/(9*gama-5)
#k = Pr / (niu*cp)
# combustion chamber - shapes - cylindrical so far
# t = pc  * Rchamber * safety factor / sigma material
# Propellant dude needs to compute gama, cp, niu for the average o/f ratio
import numpy as np

def Lambda(Propellant):
    if Propellant == 'blabla':
        return 2
    elif Propellant == 'blabla1':
        return 3
    else:
        return 4

def Lstar(Propellant):
    if Propellant == 'blabla':
        return 2
    elif Propellant == 'blabla1':
        return 3
    else:
        return 4

def PropertiesMaterial(Material):
    if Material == 'blabla':
        return 2
    elif Material == 'blabla1':
        return 3
    else:
        return 4

def CombustionChamber (Pc,At,Propellant,Material,Safety,deltap,ksi,d0):
    lamb = Lambda(Propellant)
    Vi = 4/3*math.pi*(d0/2)**3
    Vf = Vi * 0.1 #random value for now
    d = (3/4*Vf)**(1/3)*2
    time = -(d-d0)/lamb #dquadrado
    ro = 1000
    velocity = math.sqrt(1/ksi)*math.sqrt(2*deltap/ro)

    LengthChamber = velocity/time
    lstar = Lstar(Propellant)

    Vchamber = lstar * At
    Achamber = Vchamber / LengthChamber
    Rchamber = (Achamber/math.pi)**(1/2)
    sigma = 1000000
    sigma = PropertiesMaterial(Material)
    Thickness = Pc * Rchamber * Safety / sigma
    k = 1000
    Mass = Achamber * Thickness * k #k is a value to be determined

    a = 0.023
    ro = 5
    Pr = 5
    cp = 5
    niu = 5
    k = (Pr/(niu*cp))
    # ro,niu,Pr comes from propellant as well, so I'm using placeholders
    h = a * ro**0.8 * velocity**0.8 * (1/(Rchamber*2))**0.2 * (k*Pr**0.33/niu**0.8)
    return (Mass,Thickness,h,LengthChamber,Vchamber,Achamber)

[Mass,Thickness,h_chamber,LengthChamber,Vchamber,Achamber]=CombustionChamber(1*10**6,0.02,'blabla','blabla',1.5,2*10^5,1,0.001)
print(Mass,Thickness,h_chamber,LengthChamber,Vchamber,Achamber)