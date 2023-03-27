#inputs: pc, At, Propellant, Material, safety factor, velocity,d0,Tc,of,bool
#outputs: Thickness, Area of the Chamber, Conductive Heat transfer factor, Mass
import math
import warnings


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


    

#class Propellant:
 #   o_lamb = 1*10**(-6)
  #  f_lamb = 1.5*10**-6
   # lstar = [1,4]
   # o_dens = 1
    #f_dens_g = 1
    #ocp=1
    #fcp=1
    #omiu=1
    #fmiu=1
#class Material:
 #   density = 1
  #  yieldstress_l =1

#class Default:
 #   SF = 1.0
  #  kloads = 1
   # inj_velocity = 20
    #ConvergenceRatio_l = 1.5
    #ConvergenceRatio_h = 3.5
    #factor = 0.3
    #a = 0.023




def CombustionChamber (conv,Pc,At,Propellant,Material,default,velocity_f,velocity_ox,d_f,d_o,bool,rho_c,cp_c,mu_c,k_c,Pr_c,A_est,of):



    wr = 0
    er = 0
    A_minimum = A_est
    factor = default.factor #this is the factor that correlates initial droplet volume to final droplet volume.
    #final droplet Volume = initial droplet volume * factor
    Control_cycle_one = 0 #variable used to specify whether diminishing final volume of droplet is enough to reach convergence ratio
    Control_cycle_two = 0 #same thing but increasing final volume of droplet
    test = 0
    Safety = default.SF
    #velocity = default.inj_velocity
    #Input Sanity Check

    ##errors:
    if(default.ConvergenceRatio_l>default.ConvergenceRatio_h):
        er = er | (1 << 6)
        return(0,0,0,0,0,wr,er)

    if (d_f < 0) or (d_o < 0):
        er = er | (1 << 0)
        return(0,0,0,0,0,wr,er)
    if (d_o > 1 * 10 ** -3) or (d_f > 1*10 **-3):
        er = er | (1 << 1)
        return (0, 0, 0, 0, 0,wr,er)
        # quit("Diameter of droplet coming from injector is too large")

    if (velocity_ox < 0):
        er = er | (1 << 2)
        return (0, 0, 0, 0, 0, wr,er)
        # quit("Velocity of oxidizer droplet coming from injector is too small")

    if (velocity_f < 0):
        er = er | (1 << 3)
        return (0, 0, 0, 0, 0, wr,er)
        # quit("Velocity of fuel droplet coming from injector is too small")

    ## warnings :
    if(d_o < 10*10**-6) or (d_f < 10*10**-6) :
        wr = wr |(1<<0)
        #quit("Diameter of droplet coming from injector is too small")
    if(d_o > 250*10**-6) or (d_f > 250*10**-6):
        wr = wr | (1 << 1)
        #quit("Diameter of droplet coming from injector is too large")

    if(velocity_ox < 10):
        wr = wr |(1<<2)
        #quit("Velocity of oxidizer droplet coming from injector is too small")
    if(velocity_ox > 50):
        wr = wr | (1 << 3)
        #quit("Velocity of oxidizer droplet coming from injector is too large")

    if (velocity_f < 10):
        wr = wr | (1 << 4)
        #quit("Velocity of fuel droplet coming from injector is too small")
    if (velocity_f > 300):
        wr = wr | (1 << 5)
        #quit("Velocity of fuel droplet coming from injector is too large")

    #setting boundaries for chamber diameter based on expected convergence rate
    Ac_low_1 = A_minimum
    Ac_low_2=At *default.ConvergenceRatio_l

    if Ac_low_1 > Ac_low_2:
        Ac_low = Ac_low_1
    else:
        Ac_low = Ac_low_2

    Ac_high = At * default.ConvergenceRatio_h

    if Ac_high < A_est:
        er = er |(1 << 4)
        return (0, 0, 0, 0, 0, wr,er)

    d_low = math.sqrt(Ac_low / math.pi) * 2
    d_high = math.sqrt(Ac_high / math.pi) * 2

    #print("D_low:" +str(d_low))
    #print("D_high."+ str(d_high))

    #first iteration for length and diameter

    Vi_o = 4.0/3.0*math.pi*(d_o/2.0)**3.0
    Vf_o = Vi_o * factor #random value for now
    d_1 = (3.0/4.0*Vf_o/math.pi)**(1.0/3.0)*2.0
    time_o = -(d_1**2.0-d_o**2.0)/Propellant.o_lamb

    Vi_f = 4.0 / 3.0 * math.pi * (d_f / 2.0) ** 3.0
    Vf_f = Vi_f * factor  # random value for now
    d_2 = (3.0 / 4.0 * Vf_f / math.pi) ** (1.0 / 3.0) * 2.0
    time_f = -(d_2 ** 2.0 - d_f ** 2.0) / Propellant.f_lamb  # dquadrado

    LengthChamber_f = time_f * velocity_f
    LengthChamber_ox= velocity_ox * time_o

    if LengthChamber_ox>LengthChamber_f:
        LengthChamber = LengthChamber_ox
    else:
        LengthChamber = LengthChamber_f


    #print("LengthChamber:" + str(LengthChamber))

    lstar = (Propellant.lstar[0]+Propellant.lstar[1])/2
    Vchamber = lstar * At
    Achamber = Vchamber / LengthChamber
    Rchamber = (Achamber / math.pi) ** (1.0 / 2.0)
    dchamber = Rchamber * 2.0
    #print('dchamber_in ' + str(dchamber))
    #print("dchamber initial: " + str(dchamber))

    #sanity check for the outputs
    if (dchamber>d_high):
      #  print("AHAHAHAHAH")
        while(dchamber>d_high):

            #factor measures the ratio between initial and final droplet volume
            factor = factor - 0.01
            if (factor < 0):
                Control_cycle_one = 1
                break
            Vi_o = 4.0 / 3.0 * math.pi * (d_o / 2.0) ** 3.0
            Vf_o = Vi_o * factor  # random value for now
            d_1 = (3.0 / 4.0 * Vf_o / math.pi) ** (1.0 / 3.0) * 2.0
            time_o = -(d_1 ** 2.0 - d_o ** 2.0) / Propellant.o_lamb

            Vi_f = 4.0 / 3.0 * math.pi * (d_f / 2.0) ** 3.0
            Vf_f = Vi_f * factor  # random value for now
            d_2 = (3.0 / 4.0 * Vf_f / math.pi) ** (1.0 / 3.0) * 2.0
            time_f = -(d_2 ** 2.0 - d_f ** 2.0) / Propellant.f_lamb  # dquadrado

            LengthChamber_f = time_f * velocity_f
            LengthChamber_ox = velocity_ox * time_o

            if LengthChamber_ox > LengthChamber_f:
                LengthChamber = LengthChamber_ox
            else:
                LengthChamber = LengthChamber_f

            Achamber = Vchamber / LengthChamber
            #print("dchamber: " + str(dchamber))
            Rchamber = (Achamber / math.pi) ** (1.0 / 2.0)
            dchamber = Rchamber * 2.0
            #print("Chamber Length" + str(LengthChamber))
    if(Control_cycle_one == 1):
        #print("MUAHAHAHAHAHAHHAHAHAHAHAH")
        while (dchamber > d_high):
            lstar = lstar - 0.01
            #print(lstar)
            if (lstar < Propellant.lstar[0]):
                test = 1
                break

            Vchamber = lstar * At
            Achamber = Vchamber / LengthChamber
            Rchamber = (Achamber / math.pi) ** (1.0 / 2.0)
            dchamber = Rchamber * 2.0
            #print(dchamber)

    if(test == 1):
        lstar = Propellant.lstar[0]
        Vchamber = lstar * At
        Achamber = Ac_high
        Rchamber = (Achamber / math.pi) ** (1.0 / 2.0)
        dchamber = Rchamber * 2.0
        LengthChamber = Vchamber/ Achamber
        wr = wr | (1 << 6)
        #warnings.warn("Chamber Length was increased more than what it takes to fully atomize propellant, so that chamber diameter fits convergence ratio parameters")

    factor = default.factor

    if (dchamber < d_low):
        while (dchamber < d_low):
            #print("EHEHEHEHEHEHE")
            # factor measures the ratio between initial and final droplet volume
            factor = factor + 0.01
            if (factor > 0.4):
                Control_cycle_two = 1
                break
            Vi_o = 4.0 / 3.0 * math.pi * (d_o / 2.0) ** 3.0
            Vf_o = Vi_o * factor  # random value for now
            d_1 = (3.0 / 4.0 * Vf_o / math.pi) ** (1.0 / 3.0) * 2.0
            time_o = -(d_1 ** 2.0 - d_o ** 2.0) / Propellant.o_lamb

            Vi_f = 4.0 / 3.0 * math.pi * (d_f / 2.0) ** 3.0
            Vf_f = Vi_f * factor  # random value for now
            d_2 = (3.0 / 4.0 * Vf_f / math.pi) ** (1.0 / 3.0) * 2.0
            time_f = -(d_2 ** 2.0 - d_f ** 2.0) / Propellant.f_lamb  # dquadrado

            LengthChamber_f = time_f * velocity_f
            LengthChamber_ox = velocity_ox * time_o

            if LengthChamber_ox > LengthChamber_f:
                LengthChamber = LengthChamber_ox
            else:
                LengthChamber = LengthChamber_f

            Achamber = Vchamber / LengthChamber
            #print("dchamber loop 2: " + str(dchamber))
            Rchamber = (Achamber / math.pi) ** (1.0 / 2.0)
            dchamber = Rchamber * 2.0

    if (Control_cycle_two == 1):
        #print("MUAHIHIHIHIHIHIH")
        while (dchamber < d_low):
            lstar = lstar + 0.01

            Vchamber = lstar * At
            Achamber = Vchamber / LengthChamber
            Rchamber = (Achamber / math.pi) ** (1.0 / 2.0)
            dchamber = Rchamber * 2.0
            #print("dchamber:" + str(dchamber))

    if (lstar > Propellant.lstar[1]):
        wr = wr | (1 << 7)
        #warnings.warn(
            #"Cannot calculate chamber dimensions with current Lstar range, final diameter of the chamber is too small. As such program increased"
            #"Lstar to " + str(lstar))

    Thickness = Pc * Rchamber * Safety / Material.yieldstress_l

    if bool == 1:
        kloads = default.kloads
        Mass = kloads *(1/(LengthChamber/dchamber)+2)*Material.density*Safety/Material.yieldstress_l*Vchamber*Pc
        if Mass < 0:
            er = er | (1<<7)
            return(0,0,0,0,0,wr,er)

    a = default.a
    ro = rho_c
    Pr = Pr_c
    cp = cp_c
    niu = mu_c
    k = k_c

    velocity = velocity_ox*(1-1/(1+of)) + velocity_f*(1/(1+of))
    heattransfer = a * ro**0.8 * velocity**0.8 * (1/(Rchamber*2))**0.2 * (k*Pr**0.33/niu**0.8)
    if heattransfer < 0:
        er = er | (1<<8)
        return(0,0,0,0,0,wr,er)
    Re=velocity*ro*dchamber/niu

    if LengthChamber > 1.2:
        wr = wr | (1 << 8)

    if bool == 0:
        return (heattransfer,dchamber,Thickness,LengthChamber,Re,wr,er)
    else:
        return (heattransfer,dchamber,Thickness,LengthChamber,Re,Mass,wr,er)


#prop = Propellant
#mat = Material
#default = Default
#A_t = (0.0106451400000000,	0.0139612624000000,	0.0283225240000000,	0.0349934784000000,	0.0133870700000000,	0.0110515908000000,	0.0784998430000000,	0.0378063760000000,	0.117780409600000,	0.0421999156000000,	0.132902960000000,	0.132902960000000	,0.132257800000000,	0.117780409600000,	0.109419136000000,	0.620256824000000,	2.83870400000000e-06,	5.54837600000000e-06,	7.21288880000000e-05	,6.50321280000000e-05	,4.98063520000000e-05,	0.000148386800000000,	9.67740000000000e-05,	0.000381934720000000,	0.000255418844000000,	0.000255418844000000,	0.000547740840000000,	0.00102580440000000,0.000929030400000000,0.000723224360000000,0.00134838440000000)
#conv = (2.90000000000000	,3.27000000000000	,2.54000000000000,	2.52000000000000	,4,	5.33500000000000,	2	,4.10000000000000	,2.04000000000000,	2.51000000000000,	1.67000000000000	,1.67000000000000,	1.62000000000000	,2.04000000000000	,1.58000000000000	,1.30700000000000,	0,	0	,0	,3.80000000000000,	0,	0,	0,	4.20000000000000,	3.10000000000000,	3.10000000000000,	0,	0,	4.90000000000000,	0,	3)
#for i in range(len(A_t)):
#ht, dc, t, lc, re, wr, er = CombustionChamber(1.5, 20300000, 0.042, prop, mat, default, 200, 30, 100* 10 ** -6,68 * 10 ** -6, 0, 1, 1, 1, 1, 1, 0.00000000000000053,6)
#print(dc,lc,wr,er)
#print(dc,lc,er,wr)
#Ac_low = 0.053 * 1.5
#Ac_high = 0.053 * 3.5
#d_low = math.sqrt(Ac_low/math.pi)*2
#d_high = math.sqrt(Ac_high/math.pi)*2
#print(d_low,d_high)
#At = 0.053
#Ac = 2.96*At
#y=math.sqrt(Ac/math.pi)*2
#print(y