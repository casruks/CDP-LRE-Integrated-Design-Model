import Turbomachinery_code as Turbo
import CombustionChamber as Comb
import Injector as Inj
import Igniters as Ign
import Reliability as Rel
import Nozzle_loop_1 as Nz_1
import Nozzle_loop_2 as Nz_2
import Nozzle_turbine as Nz_t
import Cooling
import Materials as Mt
import numpy as np
import Aux_classes as aux

Thrust_ = 1350000 #= input("Introduce thrust")
Thrust_time_ = 180 #= input("Introduce thrust time")
Pamb_ = 1000 #= input("Introudce ambient pressure (Pa)")
prop = aux.Propellant(0)
default = aux.Default(0)
dat = aux.Data(Thrust_, Thrust_time_, Pamb_)

#Main Function
def Main(d : aux.Data):
    p_old = 0.0
    p_new = default.Pres
    bool = 0 #this variable is used to show the combustor function we are in the first loop
    regCool=Cooling.RegenerativeCool(); #inicialise cooling

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
        d.m_nozz,d.Tc,d.O_F,d.At,d.eps,d.Isp,rho_c,cp_c,mu_c,k_c,Pr_c = Nz_1.Nozzle_loop_1(p_new/100000.0,d.Thrust,d.Pa/100000.0,prop,default,default.Nozzle_type)
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
        

        #Compute injector (1)
        d.v_iox, d.v_if, d.D_f, d.D_ox, d.dp, d.eta_s, d.m_ox, d.m_f, d.n_ox, d.n_f, d.P_D = Inj.injector1(default, prop, p_new, d.m_nozz, d.O_F)


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
        d.h_comb, d.Dc, d.ThicknessChamber, d.Chamber_L, d.Re_c= Comb.CombustionChamber(p_new, d.At, prop, Mt.Rhenium, default, d.v_if, d.v_iox, d.Tc, d.O_F, bool,rho_c,cp_c,mu_c/10,k_c,Pr_c)
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
        t_noz,x_noz,y_noz,Tw_ad_noz,h_c_noz,d.D_t,d.D_e,d.L_nozzle_con,d.L_nozzle_div,d.L_tot,x_noz_cool,y_noz_cool=Nz_2.Nozzle_loop(p_new/100000.0, d.Tc, prop, Mt.Rhenium, default.Nozzle_type, d.O_F, d.eps, d.At, d.m_nozz, d.Dc, default)
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

        
        #Compute regenerative
        Re=10^5
        Tf_cool,dptcool=regCool.Run_for_Toperating1D(Tw_ad_noz, h_c_noz, t_noz,prop,Mt.Rhenium,default.A,default.T_fuel_tanks,d.m_nozz/(1.0+d.O_F)/default.n,x_noz[-1],y_noz)
        Tf_cool,dptcool_c=regCool.Run_for_Toperating0D(d.Tc, h_comb, ThicknessChamber,prop,Mt.Rhenium,Chamber_L*default.Dr,Tf_cool,d.m_nozz/(1.0+d.O_F)/default.n,Chamber_L)
        dptcool=dptcool+dptcool_c

        #Compute Turbo
        dp_cool=np.max(dptcool)
        d.ptinj, d.W_Opump, d.W_Fpump, d.W_turb, error = Turbo.TurboM(default, prop, d.O_F, d.Pa, Tf_cool, dp_cool, d.m_nozz)
        if(error): return False
        
        #Cmpute Injector (2)
        p_new, dp_ox, dp_f = Inj.injector2(default, prop, d.v_iox, d.v_if, d.D_f, d.D_ox, d.ptinj, d.eta_s)
        print("P_new: " + str(p_new))
              
    d.Pc = p_new
    bool = 1 #Shows the combustor it is out of the loop in order to compute mass!

    #Computes igniters mass
    #inputs
    ## massflow
    ## propellant
    ## default
    ## chamber temperature
    ## oxidizer to fuel ratio of main chamber
    igniter_compound = Ign.Igniters(d.m_nozz,prop,default,d.Tc,d.O_F)
    #outputs
    ## mass used for igniter


    #Compute reliability 
    #Reliability=Rel.Reliability(default, d.time, d.Thrust, d.Thrust, 0)


    #Compute masses
    chamber_material = Mt.Rhenium
    nozzle_material = Mt.Rhenium
    nozzlemass = Mt.Mass(x_noz,y_noz,t_noz,nozzle_material)
    h_comb, Dc, ThicknessChamber, Chamber_L,Re_c= Comb.CombustionChamber(p_new, d.At, prop, Mt.Rhenium, default,d.v_if,d.v_iox, d.Tc, d.O_F, bool,rho_c,cp_c,mu_c/10,k_c,Pr_c)
    #chambermass = Comb.Mass
    #Data.Total_mass=totalmass
    #totalmass = nozzlemass + chambermass
    

    #Computing costs:
    #Data.Total_cost=cost
    #cost = chamber_material.cost*chambermass + nozzle_material.cost*nozzlemass
    
    
    #Reusability:
    #Reuseability_chamber = Mt.Reusability(comb.Pc,chamber_material)
    
    print("Calculations finished")
    return True

if __name__ == '__main__':
    Main(dat)
