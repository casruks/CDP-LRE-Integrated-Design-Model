import Turbomachinery_code as Turbo
import CombustionChamber as Comb
import Injector as Inj
import Igniters as Ign
import Reliability as Rel
import Nozzle_loop_1 as Nz_1
import Nozzle_loop_2 as Nz_2
import Nozzle_turbine as Nz_t
import Cooling
import Materials as Ms
import numpy as np
import Aux_classes as aux
import math
import matplotlib.pyplot as plt

Thrust_ = 1860000 #= input("Introduce thrust")
Thrust_time_ = 150 #= input("Introduce thrust time")
Pamb_ = 1e3 #= input("Introudce ambient pressure (Pa)")
prop = aux.Propellant(0)
default = aux.Default(0)
dat = [aux.Data(Thrust_, Thrust_time_, Pamb_)]

#Main Function
def Main(d : aux.Data):
    p_old = 0.0
    p_new = default.Pres
    bool = 0 #this variable is used to show the combustor function we are in the first loop
    regCool=Cooling.RegenerativeCool(); #inicialise cooling
    errors=[]
    warnings=[]
    #First loop, pressure convergence
    while abs(p_new-p_old)/p_new > default.pres_tol:
        p_old = p_new

        #Compute nozzle (1)
        #Inputs: 
        ## Chamber pressure in bars
        ## Thrust in Newton
        ## Ambient pressure in bars
        ## Propellants class
        ## Default class
        d[-1].m_nozz,d[-1].Tc,d[-1].O_F,d[-1].At,d[-1].eps,d[-1].Isp,rho_c,cp_c,mu_c,k_c,Pr_c,errors_nz1,warnings_nz1 = Nz_1.Nozzle_loop_1(p_new/100000.0,d[-1].Thrust,d[-1].Pa/100000.0,prop,default,default.Nozzle_type)
        # Outputs:
        ## mass flow rate in kg/s
        ## Chamber temperatures in K
        ## Oxidizer to fuel ratio
        ## Throat area in m^2
        ## Expansion ratio
        ## Specific impulse in s
        ## Density of gas in comustion chamber kg/m^3
        ## Specific heat at constant pressure in chamber in J/kgK
        ## Viscosity (mu) in combustion chamber in Pa*s
        ## Conduction constant in combustion chamber in W/(m*degC), it's in degC but as long as it is used with temperature differences it shouldn't make a difference
        ## Prandtl number in combustion chamber
        ## Error value: 1== there has been a fatal error
        #if(error==1): return False
        
        errors.append((errors_nz1))
        warnings.append((warnings_nz1))

        for i in range(len(errors)):
            if errors[i]!=0:
                return False, errors, warnings;
            
        #Compute injector (1)
        ## Added A_est, representing estimated total orifice area (A_ox+A_f) for sanity check with combustion chamber dimensions..
        d[-1].v_iox, d[-1].v_if, d[-1].D_f, d[-1].D_ox, d[-1].dp, d[-1].eta_s, d[-1].m_ox, d[-1].m_f, d[-1].n_ox, d[-1].n_f, A_est, er, wr = Inj.injector1(default, prop, p_new, d[-1].m_nozz, d[-1].O_F)

        errors.append((er))
        warnings.append((wr))

        for i in range(len(errors)):
            if errors[i]!=0:
                return False, errors, warnings;
        #Compute Chamber
        #inputs
        ## Pressure at chamber in bars
        ## Throat area in m^2
        ## Material
        ## default
        ## fuel injection velocity
        ## oxidizer injection velocity
        ## Chamber temperature in K
        ## oxidizer to fuel ratio
        ## bool variable that specifies if inside loop or not
        ## density, cp, miu,k, prandlt - this all comes from nozzle just leave it like that
        ## IMPORTANT - this function is currently using hardcoded droplet diameter, bc droplet diameter coming from injector does not make sense.
        d[-1].h_comb, d[-1].Dc, d[-1].ThicknessChamber, d[-1].Chamber_L, d[-1].Re_c,wr_comb,er_comb= Comb.CombustionChamber(p_new, d[-1].At, prop, Ms.Rhenium, default, d[-1].v_if, d[-1].v_iox, d[-1].D_f, d[-1].D_ox, bool,rho_c,cp_c,mu_c/10,k_c,Pr_c,A_est)
        
        #errors.append((er_comb))
        #warnings.append((wr_comb))

        #for i in range(len(errors)):
        #    if errors[i]!=0:
        #        return False, errors, warnings;
        #outputs
        ## conductive heat transfer coefficient
        ## chamber diameter in m
        ## chamber thickness in m
        ## chamber length in m
        ## reynolds number

        #COmpute nozzle (2)
        #Inputs
        ## Chamber pressure in bars
        ## Combustion chamber temperature in K
        ## Propellants class
        ## Material properties through material class
        ## Nozzle type (from default class)
        ## O_F ratio
        ## Expansion ratio
        ## Throat area in m^2
        ## Mass flow rate in kg/s
        ## Combustion chamber diameter in m
        ## Default class
        t_noz,x_noz,y_noz,Tw_ad_noz,h_c_noz,d[-1].D_t,d[-1].D_e,d[-1].L_nozzle_con,d[-1].L_nozzle_div,d[-1].L_tot,x_noz_cool,y_noz_cool,errors_nz2,warnings_nz2=Nz_2.Nozzle_loop(p_new/100000.0, d[-1].Tc, prop, Ms.Rhenium, default.Nozzle_type, d[-1].O_F, d[-1].eps, d[-1].At, d[-1].m_nozz, d[-1].Dc, default)
        #Outputs
        ## Array of thickness at the discretized points in the nozzle (corresponding to x_noz) in m
        ## Array of discretized positions in the nozzle where all variables and properties are calcualted (positions in m)
        ## Array of radius in the nozzle at the various discretization points (corresponding to x_noz) in m
        ## Array of adiabatic wall temperature at the various discretization points (corresopnding to x_noz) in K
        ## Array of coefficients of convective heat transfer at the various discretization points (corresponding to x_noz) in W/(m^2K)
        ## Diameter of the throat in m
        ## Exit diameter in m
        ## Length of nozzle convergent in m
        ## Length of nozzle divergent in m
        ## Total length of the nozzle in m
        ## Array of x coordinates used for the cooling properties
        ## Array of y coordinates used for the cooling properties

        errors.append((errors_nz2))
        warnings.append((warnings_nz2))

        for i in range(len(errors)):
            if errors[i]!=0:
                return False, errors, warnings;
        #Compute regenerative
       
        #Tf_cool,dptcool=regCool.Run_for_Toperating1D(Tw_ad_noz, h_c_noz, t_noz,prop,Ms.Rhenium,default.A,default.T_fuel_tanks,d[-1].m_nozz/(1.0+d[-1].O_F)/default.n,x_noz_cool[-1],y_noz_cool)
        #Tf_cool,dptcool_c=regCool.Run_for_Toperating0D(d[-1].Tc, d[-1].h_comb, d[-1].ThicknessChamber,prop,Ms.Rhenium,d[-1].Chamber_L*default.Dr,Tf_cool,d[-1].m_nozz/(1.0+d[-1].O_F)/default.n,d[-1].Chamber_L)
        #dptcool=dptcool+dptcool_c
        Coolobj = Cooling.CoolingClass()
        Coolobj_c = Cooling.CoolingClass()
        nozzle_mass = Ms.Nozzle_mass(x_noz,y_noz,t_noz,Ms.D6AC_Steel)
        chamber_mass= Comb.CombustionChamber(p_new, d[-1].At, prop, Ms.D6AC_Steel, default, d[-1].v_if, d[-1].v_iox, d[-1].Tc, d[-1].O_F, 1,rho_c,cp_c,mu_c/10,k_c,Pr_c,A_est)
        #Inputs 
        #Inicial temperature in K
        # heat capacitance of the material in J/K
        #operation time in seconds
        #mass of the nozzle in kg
        #emissivity
        #adiabatic wall temperature of the nozzle np.array
        #convection coefficient in W/(m2 K) np.array
        #Thickness of the nozzle in m np.array
        #coating thickness in m 
        #propellant
        #material of the nozzle
        #material of the coating
        #Hydralic diamter of the cooling system piping in m 
        #Segment area of contact in m2
        #inicial coolant temperature in K
        #coolant mass flow in kg/s
        #length of the nozzle in m
        #Radius of the nozzle in m
        # Option to overwrite area with input area bool
        # which regenerative function to call (default 0, do not change!)
        #
        c=420 # HEAT CAPACITY, TEMPORARY, UNTIL IMPLEMENTED IN MATERIALS
        (Tf_cool,Tw_wall_nozzle_calculated, dptcool,_,type_variable_nozzle,err_nozzle_cooling,warn_nozzle_cooling) = Coolobj.Run_cooling(
            default.T0,
            c,
            default.operationtime,
            nozzle_mass,
            default.eps,
            Tw_ad_noz,
            h_c_noz,
            t_noz,
            np.array([default.default_coating_thickness for i in range(len(t_noz))]),
            prop,
            Ms.Rhenium,
            default.default_coating,
            default.Dr,
            default.A,
            default.T_fuel_tanks,
            d[-1].m_nozz / (1.0 + d[-1].O_F) / default.n,
            x_noz_cool[-1],
            y_noz_cool,
            default.overwriteA,
            default.regenerative_case,
        )
        check=0

        #errors.append((err_nozzle_cooling))
        #warnings.append((warn_nozzle_cooling))

        #for i in range(len(errors)):
        #    if errors[i]!=0:
        #        return False, errors, warnings;
        #Outputs
        #end temperature of the coolant in K
        #wall temperature array
        #pressure loss
        #_
        #errors
        #warnings

        #Inputs
        #Inicial temperature in K
        # heat capacitance of the material in J/K
        #operation time in seconds
        #mass of the chamber in kg
        #emissivity
        #adiabatic wall temperature of the chamber np.array
        #convection coefficient in W/(m2 K) np.array
        #Thickness of the chamber in m np.array
        #default thickness of the coating to
        #propellant
        #material of the chamber
        #deafult material of the coating
        #Hydralic diamter of the cooling system piping in m 
        #Segment area of contact in m2
        #inicial coolant temperature in K
        #coolant mass flow in kg/s
        #length of the chamber in m
        #Radius of the chamber in m~
        # Option to overwrite area with input area bool
        # which regenerative function to call (default 0, do not change!)
        (
            Tf_cool,
            Tw_wall_chamber_calculated,
            dptcool_c,
            _,
            type_variable_chamber,
            err_chamber_cooling,
            warn_chamber_cooling,
        ) = Coolobj_c.Run_cooling(
            default.T0,
            c,
            default.operationtime,
            chamber_mass,
            default.eps,
            np.array([d[-1].Tc]),
            np.array([d[-1].h_comb]),
            np.array([d[-1].ThicknessChamber]),
            np.array([default.default_coating_thickness]),
            prop,
            Ms.Rhenium,
            default.default_coating,
            default.Dr,
            np.array([d[-1].Chamber_L * d[-1].Dc * math.pi]),
            Tf_cool,
            d[-1].m_nozz / (1.0 + d[-1].O_F) / default.n,
            d[-1].Chamber_L,
            np.array([d[-1].Dc / 2]),
            True,
            default.regenerative_case,
        )
        #Outputs
        #end temperature of the coolant in K
        #wall temperature array
        #pressure loss
        #_
        #errors
        #warnings

        #errors.append((err_chamber_cooling))
        #warnings.append((warn_chamber_cooling))

        #for i in range(len(errors)):
        #    if errors[i]!=0:
        #        return False, errors, warnings;
        #dptcool = dptcool + dptcool_c

        #Compute Turbo
        dp_cool=np.max(dptcool)
        #dp_cool=0.004
        #Tf_cool=400
        d[-1].ptinj, d[-1].W_Opump, d[-1].W_Fpump, d[-1].W_turb, d[-1].fuel_frac, error_t = Turbo.TurboM(default, prop, d[-1].O_F, d[-1].Pa, Tf_cool, dp_cool, d[-1].m_nozz)
        errors.append((error_t))
        

        for i in range(len(errors)):
            if errors[i]!=0:
                return False, errors, warnings;
        
        #Compute Injector (2)
        p_new, dp_ox, dp_f, er_in, wr_in = Inj.injector2(default, prop, d[-1].v_iox, d[-1].v_if, d[-1].ptinj, d[-1].eta_s)

        errors.append((er_in))
        warnings.append((wr_in))

        for i in range(len(errors)):
            if errors[i]!=0:
                return False, errors, warnings;
        print("P_new: " + str(p_new))
              
    d[-1].Pc = p_new
    bool = 1 #Shows the combustor it is out of the loop in order to compute mass!
    
    # Plot nozzle contour
    #x_plot1=np.linspace(y_noz[0],y_noz[-1],len(y_noz))
    #y_plot=np.linspace(y_noz[0],y_noz[-1],len(y_noz))
    #x_plot2=np.linspace(y_noz[0],y_noz[-1],len(y_noz))
    
    #for i in range(len(x_noz)):
    #    x_plot1[i]=y_noz[i]
    #    y_plot[i]=x_noz[-1]-x_noz[i]
    #    x_plot2[i]=-x_plot1[i]
    
    #x_h1=np.linspace(x_plot2[0],x_plot1[0],100)
    #y_h1=np.linspace(x_plot2[0],x_plot1[0],100)
    #for i in range(len(x_h1)):
    #    y_h1[i]=y_plot[0]
    
    #x_h2=np.linspace(x_plot2[-1],x_plot1[-1],100)
    #y_h2=np.linspace(x_plot2[-1],x_plot1[-1],100)
    #for i in range(len(x_h2)):
    #    y_h2[i]=y_plot[-1]
    #plt.plot(x_plot1,y_plot,color='k')
    #plt.plot(x_plot2,y_plot,color='k')
    #plt.plot(x_h1,y_h1,color='k')
    #plt.plot(x_h2,y_h2,color='k')
    #plt.xlim(0, x_plot[-1]*1.1)
    #plt.ylim(0,y_plot[-1]*1.1)
    #ax = plt.gca()
    #ax.set_aspect('equal', adjustable='box')
    #plt.show()


    #Computes igniters mass
    #inputs
    ## massflow
    ## propellant
    ## default
    ## chamber temperature
    ## oxidizer to fuel ratio of main chamber
    Ign_propellant_mass, Ign_fuel_mass, Ign_ox_mass,wr_ign = Ign.Igniters(d[-1].m_nozz,prop,default,d[-1].Tc,d[-1].O_F,default.type)
    #outputs
    ## mass used for igniter
    ## By this order: igniter propellant total mass, igniter fuel mass, igniter oxidizer mass, warnings

    #Compute reliability 
    Reliability = Rel.Reliability(default, d[-1].time, d[-1].Thrust, d[-1].Thrust, default.val)
    
    #Compute Mass:
    NozzleMass = Ms.Nozzle_mass(x_noz,y_noz,t_noz,Ms.D6AC_Steel)
    print('Nozzle mass =', NozzleMass)
    ChamberMass = Comb.CombustionChamber(p_new, d[-1].At, prop, Ms.D6AC_Steel, default, d[-1].v_if, d[-1].v_iox, d[-1].Tc, d[-1].O_F, 1,rho_c,cp_c,mu_c/10,k_c,Pr_c,A_est)[5]
    IgnitorMass, mfuel, mox, wr = Ign.Igniters(d[-1].m_nozz,prop,default,d[-1].Tc,d[-1].O_F,default.type)
    #! Ns- pump specific speed not yet implemented
    #Mass = ChamberMass + IgnitorMass #+ Ms.Mass(p_new,Ms.Rhenium,Ms.Rhenium,Ms.Rhenium,d[-1].Eps,d[-1].A_t,0,aux.Default.Safety_factor,0,Turbo.Ns)
    rho_prop = Ms.RhoProp(prop.f_dens_l,prop.o_dens,aux.Data.O_F)
    #Mass = ChamberMass + IgnitorMass +NozzleMass + Ms.Mass(p_new,Ms.Rhenium,Ms.Rhenium,d[-1].Eps,d[-1].A_t, d[-1].m_nozz ,aux.Default.Safety_factor,rho_prop,'SC','try')
    #print(NozzleMass, ChamberMass, IgnitorMass) #kg
     # #Computing costs:
    # n_engine = 1
    # #Cost = Ms.Cost(Mass, Reliability, n_engine)

    print("Calculations finished")
    return 0,0,0,0,0,0,0,7,0,0,0,0,0,0

if __name__ == '__main__':
    Main(dat)