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

Thrust_ = 1500 #= input("Introduce thrust")
Thrust_time_ = 180 #= input("Introduce thrust time")
Pamb_ = 100000 #= input("Introudce ambient pressure (Pa)")



#Default values
class Default:
    #Tolerances
    pres_tol = 0.01
    toll_c_star = 0.01
    toll_F_obj = 0.01
    Max_iterations_mass_flow = 10000
    toll_P_adapted = 0.01
    Safety_factor=1.3
    noz_res=0.01

    #Seeds
    Pres = 1e6
    inj_vel = 15

    #Injectors
    Cd = 0.7

    #Nozzle
    Nozzle_type = 0
    MR = 6.03
    De_max = 2.5
    De_turbine_noz_max = 2.5
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
    ptf = 1.0e5 #Fuel tank oxidizer pressure
    Wmotor = 1.0e6 #Power of the electric motor

    #Combustion chamber
    SF = 1.0

    #Cooling
    Dr = 0.01
    A=0.0003
    T_fuel_tanks = 20
    T_ox_tanks = 60
    n=1
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
    Ox_composition = "O 2" #Composition of oxidizer for rocketcea
    o_dens = 1141.0 #Oxidizer density
    ocp = 14307.0 #oxidizer cp
    h_ox = -12.979 #oxidizer enthalpy
    o_lamb = 1.0e-3
    o_nist_enthalpy_coef = [1, 1, 1, 1, 1, 1, 1, 1]  # for shomate equation
    omiu=1.0e-6
   
    #Fuel
    Fuel_name = "LH2" #Fuel name for rocketCEA
    Fuel_composition = "H 2" #Composition of fuel for rocketcea
    f_dens_l = 71.0 #liquid fuel density
    f_dens_g = 1.0 #gaseous fuel density
    f_gamma = 1.4 #fuel gamma
    fcp = 14307.0 #fuel cp
    h_fuel = -9.012 # fuel enthalpy
    R_f = 4.1573 #fuel gas constant
    f_lamb = 1.0e-6
    fmiu=1.0e-6
    f_nist_enthalpy_coef = [1, 1, 1, 1, 1, 1, 1, 1]  # for shomate equation
    MR = 3 #mixture ratio
    
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
#Main Function
def Main(Thrust, Thrust_time, Pamb):
    p_old = 0.0
    p_new = default.Pres
    inj_vel = default.inj_vel
    bool = 0 #this variable is used to show the combustor function we are in the first loop
    
    regCool=Cooling.RegenerativeCool();#inicialise cooling
    while abs(p_new-p_old)/p_new > default.pres_tol:
        p_old = p_new
        #Compute nozzle (1)
        m,Tc,O_F,At,eps,Isp,rho_c,cp_c,mu_c,k_c,Pr_c = Nz_1.Nozzle_loop_1(p_new/100000,Thrust,Pamb/100000,prop,default)

        #Compute injector (1)
            # placeholders for propellant reference factor K_prop =1
        mu_prop = 2.69e-3    # [lbm/(ft-s)]
        sig_prop = 17        # [dynes/cm]    
        rho_prop = 47.7      # [lbm/ft3]
        Cd = 0.7
        v_iox, v_if, D_f, D_o = Inj.injector1(Cd, m, O_F, prop.o_dens, prop.f_dens_l, mu_prop, sig_prop, rho_prop)
        #Compute chamber - needs Chamber temperature + oxider to fuel ratio from previous functions (Tc and of)
        h_comb, Dc, ThicknessChamber, Chamber_L,Re_c= Comb.CombustionChamber(p_new, At, prop, Mt.Rhenium, default.SF, inj_vel, D_o, Tc, O_F, bool,rho_c,cp_c,mu_c,k_c,Pr_c)

        #COmpute nozzle (2)
        t_noz,x_noz,y_noz,Tw_ad_noz,h_c_noz,P_noz,T_noz,Re_t=Nz_2.Nozzle_loop(p_new/100000, Tc, prop, Mt.Rhenium, default.Nozzle_type, O_F, eps, At, m, Dc, default)
        
        #Compute regenerative

        #0D, nozzle (do not delete, in case 1D doesnt work)
        #Tf_cool, T_w_after_cooling,dptcool=regCool.Run(Tw_ad_noz[0], h_c_noz[0], t_noz[0],prop,Mt.Rhenium,default.Dr,default.A,default.T_fuel_tanks,Re_t,m/(1+O_F),x_noz[-1])

        #1D Nozzle + 0D chamber
        #placeholder Re
        Re=10^5
        Pr=1
        #400-600K
        #1200 nozzle
        #2500 combustion~

        Tf_cool,dptcool=regCool.Run_for_Toperating1D(Tw_ad_noz, h_c_noz, t_noz,prop,Mt.Rhenium,default.A,default.T_fuel_tanks,m/(1+O_F)/default.n,x_noz[-1],y_noz)
        Tf_cool,dptcool_c=regCool.Run_for_Toperating0D(Tc, h_comb, ThicknessChamber,prop,Mt.Rhenium,Chamber_L*default.Dr,Tf_cool,m/(1+O_F)/default.n,Chamber_L)
        dptcool=dptcool+dptcool_c
        #Tf_cool=450
        #dptcool=1000000
        #Compute Turbo
        ptinj = Turbo.TurboM(default, prop, O_F, Pamb, Tf_cool[-1], dptcool, m)

        #Cmpute Injector (2)
        p_new, dp_ox, dp_f = Inj.injector2(v_iox, v_if, D_f, D_o, ptinj, Cd, prop.o_dens, prop.f_dens_l)
        print(p_new)
        
        
    bool = 1 #Shows the combustor it is out of the loop in order to compute mass!
    #Compute Ignitor - m is the mass flow, Hc is enthalpy of propelants at chamber exit, H0 is enthalpy of propelants at chamber entry
    #For further information on igniter output, see comments on first line of the igniters functions
    # igniter_results = Igniters(m,prop,default,Tc,O_F)
    #Compute Masses
    print(p_new)
    print(Isp)
    print(eps)
    print(At)
    print(D_f)
    print(D_o)
    print("Colling D: ",regCool.D)
    #Compute reliability
    cycle = ['D_FR_SC', 'D_FF_SC', 'S_FR_SC', 'S_OR_SC', 'S_FR_GG', 'SP_EX']
    Prop = ['LOX_LH2', 'LOX_RP1']
    Reliability=Rel.Reliability(t, cycle, Fnom, Fop, N, prop, 0)

    #Compute masses
    chamber_material = Mt.Rhenium
    nozzle_material = Mt.Rhenium
    nozzlemass = Mt.Mass(x_noz,y_noz,t_noz,nozzle_material)
    chambermass = Comb.Mass
    totalmass = nozzlemass + chambermass
    
    #Computing costs:
    cost = Mt.chamber_material.cost*chambermass + Mt.nozzle_material.cost*nozzlemass
    
    print("Starting...")
    return p_new,Isp,m,150,Tc,Chamber_L

if __name__ == '__main__':
    Main(Thrust_, Thrust_time_, Pamb_)
