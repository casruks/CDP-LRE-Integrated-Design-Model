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

p_a = 1.0e5
Thrust = 15000
Thurst_time = 60




#Default values
class Default:
    #Tolerances
    pres_tol = 0.01
    toll_c_star = 0.01
    toll_F_obj = 0.01
    Max_iterations_mass_flow = 10000
    toll_P_adapted = 0.01
    Safety_factor=1.3

    #Seeds
    Pres = 1e6
    inj_vel = 15

    #Injectors
    Cd = 0.7

    #Nozzle
    Nozzle_type = 0
    MR = 0
    De_max = 2.5
    Theta_con = 60
    Theta_conical = 15
    Theta_bell = 55
    TH_exit_bell = 3
    R_u_ratio=1

    #Turbomachinery
    cycle_type = "EX"
    Eff_t = 0.6 #Turbine efficiency
    Eff_p = 0.6 #Pump efficiency
    Eff_m = 0.95 #Mechanical efficiency between turbine and pumps
    p_to = 1.0e5 #oxidizer tank storage pressure
    p_tf = 1.0e5 #Fuel tank oxidizer pressure

    #Combustion chamber
    SF = 1.0

    #Cooling
    Dr = 0.01
    A=0.0003

    T_fuel_tanks = 20
    T_ox_tanks = 60

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
    Ox_name = "LOX" #Oxidizer name for rocketCEA
    o_dens = 1141.0 #Oxidizer density
    ocp = 14307.0 #oxidizer cp
    o_lamb = 1.0e-3
    omiu=1.0e-6
   
    #Fuel
    Fuel_name = "LH2" #Fuel name for rocketCEA
    f_dens_l = 71.0 #liquid fuel density
    f_dens_g = 1.0 #gaseous fuel density
    f_gamma = 1.4 #fuel gamma
    fcp = 14307.0 #fuel cp
    R_f = 4.1573 #fuel gas constant
    f_lamb = 1.0e-3
    fmiu=1.0e-6
    
    Frozen_state=0
    
    #Propellant
    gama = 1.4
    tq = 0.9 #characteristic chemical time of propellant

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
    Pamb = 1.01325 #= input("Introudce ambient pressure (bar)")

    p_new = 0.0
    p_old = default.Pres
    inj_vel = default.inj_vel

    
    T_w_after_cooling = 0;#Temperature after cooling
    regCool=Cooling.RegenerativeCool();#inicialise cooling
    while abs(p_new-p_old)/p_old > default.pres_tol:
        p_new = p_old
        #Compute nozzle (1)
        m,Tc,O_F,At,eps,Isp = Nz_1.Nozzle_loop_1(p_new,Thrust,Pamb,prop,default)

        #Compute injector (1)
            # placeholders for propellant reference factor K_prop =1
        mu_prop = 2.69e-3    # [lbm/(ft-s)]
        sig_prop = 17        # [dynes/cm]    
        rho_prop = 47.7      # [lbm/ft3]
        Cd = 0.7
        v_iox, v_if, D_f, D_o = Inj.injector1(Cd, m, O_F, prop.o_dens, prop.f_dens_l, mu_prop, sig_prop, rho_prop)
        #Compute chamber - needs Chamber temperature + oxider to fuel ratio from previous functions (Tc and of)
        h_comb, Dc, ThicknessChamber = Comb.CombustionChamber(p_new, At, prop, Mt.Rhenium, default.SF, inj_vel, D_o, Tc, O_F, bool)

        #COmpute nozzle (2)
        t_noz,x_noz,y_noz,Tw_ad_noz,h_c_noz,P_noz,T_noz,Re_t=Nz_2.Nozzle_loop(p_new, Tc, prop, Mt.Rhenium, default.Nozzle_type, O_F, eps, At, m, Dc, default)
        
        #Compute regenerative
        Tf_cool, T_w_after_cooling,dptcool=regCool.Run(Tw_ad_noz[0], h_c_noz, t_noz[0],prop,Mt.Rhenium,default.Dr,default.A,default.T_fuel_tanks,Re_t,m/(1+O_F),default.L)

        #Compute Turbo
        ptinj = Turbo.TurboM(default, prop, O_F, p_a, Tf_cool, dptcool, m)

        #Cmpute Injector (2)
        p_new, dp_ox, dp_f = Inj.injector2(v_iox, v_if, D_f, D_o, ptinj, Cd, prop.o_dens, prop.f_dens_l)
        print(p_new)
    bool = 1 #Shows the combustor it is out of the loop in order to compute mass!
    #Compute Ignitor - m is the mass flow, Hc is enthalpy of propelants at chamber exit, H0 is enthalpy of propelants at chamber entry
    #For further information on igniter output, see comments on first line of the igniters functions
    # igniter_results = Igniters(m,Hc,H0,default)
    #Compute Masses
    print(p_new)
    #Compute reliability
    ## cycle = ['D_FR_SC', 'D_FF_SC', 'S_FR_SC', 'S_OR_SC', 'S_FR_GG', 'SP_EX']
    ## Prop = ['LOX_LH2', 'LOX_RP1']
    # Rel.Reliability(t, cycle, Fnom, Fop, N, prop, 0)

    #Compute costs

    print("Starting...")
