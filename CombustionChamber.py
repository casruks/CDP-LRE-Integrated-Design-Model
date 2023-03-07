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
class Propellant:
    o_lamb = 1*10**(-6)
    f_lamb = 1.41*10**(-5)
    gama = 1.4
    tq = 0.75*10**(-3)
    o_dens = 1
    f_dens_g = 1
    ocp=1
    fcp=1
    omiu=1
    fmiu=1
class Material:
    density = 1
    yieldstress_l =1



def CombustionChamber (Pc,At,Propellant,Material,Safety,velocity,d0,Tc,of,bool,rho_c,cp_c,mu_c,k_c,Pr_c):

    Vi = 4/3*math.pi*(d0/2)**3
    Vf = Vi * 0.3 #random value for now
    d = (3/4*Vf/math.pi)**(1/3)*2
    time_f = -(d**2-d0**2)/Propellant.f_lamb #dquadrado
    time_o = -(d**2-d0**2)/Propellant.o_lamb
    kloads = 1

    nfuel = 1/(of+1)
    nox = 1-nfuel

    gama = Propellant.gama
    GAMA = math.sqrt(gama)*(2/(gama+1))**((gama+1)/(2*gama-2))

    if time_f>time_o:
        time = time_f
    else:
        time = time_o

    LengthChamber = velocity*time
    #lstar = Propellant.tq * GAMA * math.sqrt(gama*Tc)

    lstar =  Propellant.lstar



    Vchamber = lstar * At
    Achamber = Vchamber / LengthChamber
    print(Achamber)
    Rchamber = (Achamber/math.pi)**(1/2)
    dchamber = Rchamber * 2



    Thickness = Pc * Rchamber * Safety / Material.yieldstress_l

    if bool == 1:
        Mass = kloads *(1/(LengthChamber/dchamber)+2)*Material.density*Safety/Material.yieldstress_l*Vchamber*Pc

    a = 0.023
    ro = rho_c
    Pr = Pr_c
    cp = cp_c
    niu = mu_c
    k = k_c
    #weighted averages for the several Parameters
    #ro = nox*Propellant.o_dens + nfuel*Propellant.f_dens_g
    #Pr = 4*gama/(9*gama-5)
    #cp = nox*Propellant.ocp + nfuel*Propellant.fcp
    #niu = nox*Propellant.omiu + nox*Propellant.fmiu
    #k = (Pr/(niu*cp))

    heattransfer = a * ro**0.8 * velocity**0.8 * (1/(Rchamber*2))**0.2 * (k*Pr**0.33/niu**0.8)

    Re=velocity*ro*dchamber/niu
    if bool == 0:
        return (heattransfer,dchamber,Thickness,LengthChamber,Re)
    else:
        return (heattransfer,dchamber,Thickness,LengthChamber,Re,Mass)


prop = Propellant
mat = Material

ht,dc,t,lc,re = CombustionChamber(20300000,0.053,prop,mat,2,30,1.5*10**-4,3400,6,0)
#print(dc,t,lc,)
#At = 0.053
#Ac = 2.96*At
#y=math.sqrt(Ac/math.pi)*2
#print(y)