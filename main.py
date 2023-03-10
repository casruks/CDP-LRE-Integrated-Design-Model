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

Thrust_ = 1350000 #= input("Introduce thrust")
Thrust_time_ = 180 #= input("Introduce thrust time")
Pamb_ = 1000 #= input("Introudce ambient pressure (Pa)")


#Default values
class Default:
    
    #Seeds
    Pres = 1e6
    inj_vel = 15

    #Injectors
    Cd = 0.7
    d_ox = d_f = (0.0135 +0.281)*0.0254/2.0      # 3.74 mm http://libgen.rs/book/index.php?md5=3D236B9BDD4070690CA83952058D9A1F p.113
    InjTypes = ['like', 'unlike', 'pintle']
    InjType = InjTypes[2]
    mu_prop = 2.69e-3       # [lbm/(ft-s)], 1 lbm/(ft-s) = 1.4881639 Pa.s
    sig_prop = 17.0         # [dynes/cm], 1 dyn/cm = 1e-7 N/m     
    rho_prop = 47.7         # [lbm/ft3], 1 lbm/ft3 = 16.0185 kg/m3
    p_center = p_j = 1  #measured centerline pressure, measured mean jet pressure
    
    #Nozzle
    Nozzle_type = 0 # Type of nozzle, 0=conical, 1=bell
    MR = 0 # O_F ratio, 0=optimize for c*
    De_max = 2.5 # Maximum exit diameter of the nozzle
    De_turbine_noz_max = 2.5 # Maximum exit diameter for the turbine exhaust nozzle
    Theta_con = 60 # Angle of the convergent part of the nozzle in degrees
    Theta_conical = 15 # Angle of the divergent part for the conical nozzle, in degrees
    Theta_bell = 55 # Angle of the divergent part for the bell nozzle, in degrees
    TH_exit_bell = 3 # Exit angle for the bell nozzle, in degrees
    R_u_ratio=1 # Ratio between curvature radius and throat radius (for convergent throat section in bell, and for whole throat section in conical)
    R_u_bell=0.382 # Ratio between curvature radius and throat radius for divergent throat section in bell nozzle
    #Tolerances (For the nozzle)
    toll_c_star = 0.01 # Tollerance on the absolute difference between two iteration values of O_F ratio (Nozzle_1)
    toll_F_obj = 0.01  # Tollerance on the normalized difference between thrust calculated in iteration and target thrust (Nozzle_1)
    Max_iterations_mass_flow = 10000 # Maximum iteration for the third part of Nozzle_1 code (Nozzle_1)
    toll_P_adapted = 0.01 # Tollerance on the normalized difference between exit pressure and ambient pressure (Nozzle_1)
    noz_res=150 # Number of points in the discretization of the whole nozzle (Nozzle_2)

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
    lstar=0.8

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
    Safety_factor=1.3

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
    o_lamb = 1.0e-6
    o_nist_enthalpy_coef = [1, 1, 1, 1, 1, 1, 1, 1]  # for shomate equation
    omiu=1.0e-6
    lstar=0.9
   
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
    MR = 6.1 #mixture ratio
    
    Frozen_state=0 # Frozen state of the propellant 0=chemical equilibrium flow, 1=frozen flow (from throat onwards)
    
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


#Data class
class Data:
    #Global
    Thrust = 0.0
    time = 0.0
    Pa = 0.0
    O_F = 0.0
    Total_mass=0.0 # Total mass of the engine [kg]
    Total_cost=0.0 #Total cost of the engine [Eur]

    #Nozzle
    Pc = 0.0 
    Isp_noz = 0.0 #[s]
    m_nozz = 0.0 #[kg/s]
    L_div=0.0 # Length of the divergent part [m]
    L_con=0.0 # Length of the convergent part [m]
    A_t=0.0 # Throat area [m^2]
    Eps=0.0 # Expansion ratio [-]
    Dt=0.0 # Throat diameter [m]
    De=0.0 # Exit diameter [m]


    #Turbo
    W_Opump = 0.0
    W_Fpump = 0.0
    W_turb = 0.0
    fuel_frac = 0.0
    Ptinj = 0.0
    dptop = 0.0
    dptfp = 0.0

    #Combustion

    #Cooling

    #Injector

    #Ignitor

    #Material

    def __init__(self, Th, t, p):
        self.Thrust = Th
        self.time = t
        self.Pa = p

dat = Data(Thrust_, Thrust_time_, Pamb_)


#Main Function
def Main(d : Data):
    p_old = 0.0
    p_new = default.Pres
    inj_vel = default.inj_vel
    bool = 0 #this variable is used to show the combustor function we are in the first loop

    regCool=Cooling.RegenerativeCool();#inicialise cooling
    while abs(p_new-p_old)/p_new > default.pres_tol:
        p_old = p_new
        #Compute nozzle (1)

        #Inputs: 
        ## Chamber pressure in bars
        ## Thrust in Newton
        ## Ambient pressure in bars
        ## Propellants class
        ## Default class
        m,Tc,O_F,At,eps,Isp,rho_c,cp_c,mu_c,k_c,Pr_c = Nz_1.Nozzle_loop_1(p_new/100000.0,d.Thrust,d.Pa/100000.0,prop,default)

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
        
        Data.m_nozz=m 
        Data.O_F=O_F
        Data.A_t=At
        Data.Eps=eps
        Data.Isp_noz=Isp
        #Compute injector (1)
        v_iox, v_if, D_f, D_ox, dp, eta_s = injector1(default, prop, p_c, m, OF)
        
        #Compute chamber - needs Chamber temperature + oxider to fuel ratio from previous functions (Tc and of)
        h_comb, Dc, ThicknessChamber, Chamber_L,Re_c= Comb.CombustionChamber(p_new, At, prop, Mt.Rhenium, default.SF, inj_vel, D_ox, Tc, O_F, bool,rho_c,cp_c,mu_c/10.0,k_c,Pr_c)

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

        t_noz,x_noz,y_noz,Tw_ad_noz,h_c_noz=Nz_2.Nozzle_loop(p_new/100000.0, Tc, prop, Mt.Rhenium, default.Nozzle_type, O_F, eps, At, m, Dc, default)
        
        #Outputs
        ## Array of thickness at the discretized points in the nozzle (corresponding to x_noz) in m
        ## Array of discretized positions in the nozzle where all variables and properties are calcualted (positions in m)
        ## Array of radius in the nozzle at the various discretization points (corresponding to x_noz) in m
        ## Array of adiabatic wall temperature at the various discretization points (corresopnding to x_noz) in K
        ## Array of coefficients of convective heat transfer at the various discretization points (corresponding to x_noz) in W/(m^2K)
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

        Tf_cool,dptcool=regCool.Run_for_Toperating1D(Tw_ad_noz, h_c_noz, t_noz,prop,Mt.Rhenium,default.A,default.T_fuel_tanks,m/(1.0+O_F)/default.n,x_noz[-1],y_noz)
        Tf_cool,dptcool_c=regCool.Run_for_Toperating0D(Tc, h_comb, ThicknessChamber,prop,Mt.Rhenium,Chamber_L*default.Dr,Tf_cool,m/(1.0+O_F)/default.n,Chamber_L)
        dptcool=dptcool+dptcool_c
        #Tf_cool=450
        #dptcool=1000000
        #Compute Turbo
        dp_cool:float
        dp_cool=np.max(dptcool)
        ptinj = Turbo.TurboM(default, prop, O_F, d.Pa, Tf_cool, dp_cool, m)
        
        #Cmpute Injector (2)
        p_new, dp_ox, dp_f = Inj.injector2(default, prop, v_iox, v_if, D_f, D_ox, ptinj, eta_s)
        print("P_new: " + str(p_new))
        
        
    bool = 1 #Shows the combustor it is out of the loop in order to compute mass!
    #Compute Ignitor - m is the mass flow, Hc is enthalpy of propelants at chamber exit, H0 is enthalpy of propelants at chamber entry
    #For further information on igniter output, see comments on first line of the igniters functions
    # igniter_results = Igniters(m,prop,default,Tc,O_F)
    #Compute Masses
    print(p_new)
    print(Isp)
    print(eps)
    print(At)
    print(m)
    print(D_f)
    print(D_ox)
    print("Colling D: ",regCool.D)
    #Compute reliability
    cycle = ['D_FR_SC', 'D_FF_SC', 'S_FR_SC', 'S_OR_SC', 'S_FR_GG', 'SP_EX']
    Prop = ['LOX_LH2', 'LOX_RP1']
    t=d.time
    Fnom=d.Thrust
    #Reliability=Rel.Reliability(t, cycle, Fnom, Fop, N, prop, 0)

    #Compute masses
    chamber_material = Mt.Rhenium
    nozzle_material = Mt.Rhenium
    nozzlemass = Mt.Mass(x_noz,y_noz,t_noz,nozzle_material)
    h_comb, Dc, ThicknessChamber, Chamber_L,Re_c,chambermass= Comb.CombustionChamber(p_new, At, prop, Mt.Rhenium, default.SF, inj_vel, D_o, Tc, O_F, bool,rho_c,cp_c,mu_c/10,k_c,Pr_c)
    #chambermass = Comb.Mass
    Data.Total_mass=totalmass
    totalmass = nozzlemass + chambermass
    
    #Computing costs:
    Data.Total_cost=cost
    cost = chamber_material.cost*chambermass + nozzle_material.cost*nozzlemass
    
    #Mass Estimation Funtion for Hydro-lox engines:
    #M = 0.00051*Thrust_**0.92068
    
    #Reusability:
    #Reuseability_chamber = Mt.Reusability(comb.Pc,chamber_material)
    
    print("Starting...")
    return p_new,Isp,m,m*d.time,Tc,Chamber_L

if __name__ == '__main__':
    Main(dat)
