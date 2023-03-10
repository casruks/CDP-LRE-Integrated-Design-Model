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
<<<<<<< Updated upstream
#class Propellant:
 #   o_lamb = 1*10**(-6)
  ## gama = 1.4
   # tq = 0.75*10**(-3)
   # o_dens = 1
   # f_dens_g = 1
   # ocp=1
   # fcp=1
   # omiu=1
   # fmiu=1
#class Material:
 #   density = 1
  #  yieldstress_l =1
    
=======
class Propellant:
    o_lamb = 1*10**(-6)
    f_lamb = 1*10**-3
    lstar = 0.95
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

class Default:
    SF = 1.0
    D_0 = 1.35 * 10 ** -4
    kloads = 1
    inj_velocity = 15

>>>>>>> Stashed changes


def CombustionChamber (Pc,At,Propellant,Material,default,Tc,of,bool,rho_c,cp_c,mu_c,k_c,Pr_c):

    Safety = default.SF
    velocity = default.inj_velocity
    d0 = default.D_0
    #Input Sanity Check

    #if(d0 < 1.35*10**-5):
        #quit("Diameter of droplet coming from injector is too small")
    #if(d0 > 0.5*10**-3):
        #quit("Diameter of droplet coming from injector is too large")

    #if(velocity < 15):
        #quit("Velocity of droplet coming from injector is too small")
    #if(velocity > 60):
        #quit("Velocity of droplet coming from injector is too large")


    #setting boundaries for chamber diameter based on expected convergence rate
    Ac_low = At * 1.5
    Ac_high = At * 3.5
    d_low = math.sqrt(Ac_low / math.pi) * 2
    d_high = math.sqrt(Ac_high / math.pi) * 2

    #first iteration for length and diameter

    Vi = 4.0/3.0*math.pi*(d0/2.0)**3.0
    Vf = Vi * 0.3 #random value for now
    d = (3.0/4.0*Vf/math.pi)**(1.0/3.0)*2.0
    time_f = -(d**2.0-d0**2.0)/Propellant.f_lamb #dquadrado
    time_o = -(d**2.0-d0**2.0)/Propellant.o_lamb

    if time_f>time_o:
        time = time_f
    else:
        time = time_o

    LengthChamber = velocity * time

    lstar = Propellant.lstar
    Vchamber = lstar * At
    Achamber = Vchamber / LengthChamber
    Rchamber = (Achamber / math.pi) ** (1.0 / 2.0)
    dchamber = Rchamber * 2.0
    print("dchamber initial: " + str(dchamber))

    #sanity check for the outputs
    if (dchamber>d_high):
        i = 1
        print("AHAHAHAHAH")
        while(dchamber>d_high):

            #factor measures the ratio between initial and final droplet volume
            factor = 0.3 - 0.01*i
            Vf = Vi * factor  # random value for now
            d = (3.0 / 4.0 * Vf / math.pi) ** (1.0 / 3.0) * 2.0
            if (factor < 0):
                quit("Cannot calculate chamber dimensions in the accepted boundaries")
            time_f = -(d ** 2.0 - d0 ** 2.0) / Propellant.f_lamb  # dquadrado
            time_o = -(d ** 2.0 - d0 ** 2.0) / Propellant.o_lamb
            if time_f > time_o:
                time = time_f
            else:
                time = time_o

            #computation of new chamber length and diameter
            LengthChamber = velocity * time
            Achamber = Vchamber / LengthChamber
            print("dchamber: " + str(dchamber))
            Rchamber = (Achamber / math.pi) ** (1.0 / 2.0)
            dchamber = Rchamber * 2.0
            i=i+1
    factor = 0.3

    if (dchamber < d_low):
        i = 1
        while (dchamber < d_low and factor<0.4):
            print("EHEHEHEHEHEHE")
            # factor measures the ratio between initial and final droplet volume
            factor = 0.3 + 0.01 * i
            Vf = Vi * factor  # random value for now
            d = (3.0 / 4.0 * Vf / math.pi) ** (1.0 / 3.0) * 2.0
            if (factor > 0.4):
                quit("Cannot calculate chamber dimensions in the accepted boundaries")
            time_f = -(d ** 2.0 - d0 ** 2.0) / Propellant.f_lamb  # dquadrado
            time_o = -(d ** 2.0 - d0 ** 2.0) / Propellant.o_lamb
            if time_f > time_o:
                time = time_f
            else:
                time = time_o

            # computation of new chamber length and diameter
            LengthChamber = velocity * time
            Achamber = Vchamber / LengthChamber
            print("dchamber loop 2: " + str(dchamber))
            Rchamber = (Achamber / math.pi) ** (1.0 / 2.0)
            dchamber = Rchamber * 2.0
            i = i + 1

    #gama = Propellant.gama
    #GAMA = math.sqrt(gama)*(2.0/(gama+1.0))**((gama+1.0)/(2.0*gama-2.0))



    #lstar = Propellant.tq * GAMA * math.sqrt(gama*Tc)




    Thickness = Pc * Rchamber * Safety / Material.yieldstress_l

    if bool == 1:
        kloads = default.kloads
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
default = Default

ht,dc,t,lc,re = CombustionChamber(20300000, 0.053, prop, mat, default, 3400, 6, 0, 1, 1, 1, 1, 1)
print(dc,lc)
Ac_low = 0.053 * 1.5
Ac_high = 0.053 * 3.5
d_low = math.sqrt(Ac_low/math.pi)*2
d_high = math.sqrt(Ac_high/math.pi)*2
print(d_low,d_high)
#At = 0.053
#Ac = 2.96*At
#y=math.sqrt(Ac/math.pi)*2
#print(y)