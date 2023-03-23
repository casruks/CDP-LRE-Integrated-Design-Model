import Materials as Mt

#Default values
class Default:
    #Seeds
    Pres = 1e6
    inj_vel = 20

    #Tolerances
    pres_tol = 0.01 #Tolerance for pressure convergence

    #Injectors
    Cd = 0.7
    d_ox = d_f = (0.0135 +0.281)*0.0254/2.0      # 3.74 mm http://libgen.rs/book/index.php?md5=3D236B9BDD4070690CA83952058D9A1F p.113
    mu_prop = 2.69e-3       # [lbm/(ft-s)], 1 lbm/(ft-s) = 1.4881639 Pa.s
    sig_prop = 17.0         # [dynes/cm], 1 dyn/cm = 1e-7 N/m     
    rho_prop = 47.7         # [lbm/ft3], 1 lbm/ft3 = 16.0185 kg/m3
    p_center = p_j = 1      #measured centerline pressure, measured mean jet pressure
     # Should be implemented for GUI user input
    InjTypes = ['like', 'unlike', 'pintle'] #Please keep this one in default as well. Used for Injector1() via Default.InjTypes 
    InjType = InjTypes[2] #This one is placeholder
    dp_user = 0.0 #bar, used when user specifies a pressure drop.
    dp_state = False #If False, user did not specify dp, if True user did specifiy dp.
    #Nozzle
    Nozzle_type = 0 # Type of nozzle, 0=conical, 1=bell

    MR = 6.04 # O_F ratio, 0=optimize for c* LIMIT: >=0
    De_max = 2.5 # Maximum exit diameter of the nozzle LIMIT: >0
    De_turbine_noz_max = 2.5 # Maximum exit diameter for the turbine exhaust nozzle LIMIT: >0
    Theta_con = 45 # Angle of the convergent part of the nozzle in degrees  LIMITS: 0<THETA<90
    Theta_conical = 15 # Angle of the divergent part for the conical nozzle, in degrees LIMITS: 0<THETA<90
    Theta_bell = 55 # Angle of the divergent part for the bell nozzle, in degrees LIMITS: 0<THETA<90
    TH_exit_bell = 3 # Exit angle for the bell nozzle, in degrees LIMITS: 0<THETA<90
    R_u_ratio=1 # Ratio between curvature radius and throat radius (for convergent throat section in bell, and for whole throat section in conical) LIMITS: RU>0
    R_u_bell=0.382 # Ratio between curvature radius and throat radius for divergent throat section in bell nozzle LIMITS: RU>0
    Eps_max=75 # Maximum expansion ratio. LIMITS: >1

    #Tolerances (For the nozzle)
    toll_c_star = 0.01 # Tollerance on the absolute difference between two iteration values of O_F ratio (Nozzle_1) LIMIT: >0
    toll_F_obj = 0.01  # Tollerance on the normalized difference between thrust calculated in iteration and target thrust (Nozzle_1) LIMIT: >0
    Max_iterations_mass_flow = 10000 # Maximum iteration for the third part of Nozzle_1 code (Nozzle_1) LIMIT:>0
    toll_P_adapted = 0.01 # Tollerance on the normalized difference between exit pressure and ambient pressure (Nozzle_1) LIMIT: >0
    noz_res=150 # Number of points in the discretization of the whole nozzle (Nozzle_2) LIMIT: >0
    n_cool=90 # Number of points used for cooling properties in the divergent part of the nozzle (Nozzle_2) LIMIT: >0
    SF_noz=1.3 # Safety factor for nozzle thickness

    #Turbomachinery
    cycle_type = 0 # 0:EX (expander) - 1:CB (coolant bleed) - 2:GG (gas generator) - 3:SC (staged combustion) - 4:EL (electrical) - 5:PF (pressure fed)
    Eff_t = 0.6 #[-] Turbine efficiency
    Eff_po = 0.6 #[-] oxidizer Pump efficiency
    Eff_pf = 0.6 #[-] fuel Pump efficiency
    Eff_m = 0.95 #[-] Mechanical efficiency between turbine and pumps
    p_to = 1.0e5 #[Pa] oxidizer tank storage pressure
    ptf = 1.0e5 #[Pa] Fuel tank oxidizer pressure
    Wmotor = 1.0e6 #[W] Power of the electric motor
    v_loss = 100 #[Pa] valve losses
    l_def = 0.1 #[-] Default percentage of total mass flow for GG

    #Combustion chamber
    SF = 1.0 # Safety factor for thickness estimation, put on advanced inputs
    a = 0.023 # Constant for Cornellise equation DOES NOT need to be put on advanced inputs (it is fixed)
    D_0 = 1.50*10**-4 #this is just a placeholder, bc the values from the injectors were not okay, do not put on inputs
    # DO NOT DELETE D_0, I've been reading literature, we cannot assume injector diameter == to droplet diameter as we have been doing so far. they have a very big difference
    kloads = 1 #Correction factor for mass, put on advanced inputs
    ConvergenceRatio_l = 1.5 #Minimum acceptable Convergence ratio
    ConvergenceRatio_h = 3.5 #Maximum acceptable Convergence ratio
    factor = 0.3  # this is the factor that correlates initial droplet volume to final droplet volume. final droplet Volume = initial droplet volume * factor

    # Cooling
    Dr = 0.01  # [m] hydralic diameter of the coolant channel
    A = 0.0003  # [m2] area of contact for each segment of the cooling
    T_fuel_tanks = 20  # [K] temperature of the fuel tanks, considered the inicial coolant temperature
    # T_ox_tanks = 60 #[K] temperature of the oxidiser tanks
    n = 1  # number of coolant chanels
    default_coating = Mt.Materials(
        "default_coating", 0.0, 0.0, 0.0, 0.0, 1, 0
    )  # default coolant
    default_coating_thickness = 0  # default coolant thickness
    T0 = 293.5  # [k] default inicial temperature
    eps = 0.85  # default emissivity
    overwriteA = False  # option to overwrite the surface area calculated by the program with the input variable A, given by the user or default class
    regenerative_case = 0  # option of which function to use in regenerative cooling; 0 corresponds to the explicit function Run1D()
    operationtime = (
        10000000000000000  # [s] default operation time (large to imply infinite time)
    )
    perimeter_percentage = 1  # percentage of the perimeter used for cooling chambers

    #Igniters

    ignburntime = 4 #Put on advanced inputs, it is the ignition burn time.
    ign_o_f=0.7 #fuel ratio of the igniter propellant
    fudgefactor = 20 #factor to correct for mass overestimation
    type = '00'

    #Material
    #material = Ms.Rhenium #Default Material selected OCCHIO ALLE STRONZATE
    Safety_factor = 1.3 #Default safety factor
    material_list = [Mt.Custom,Mt.Rhenium, Mt.Al_7075_T6, Mt.Ti6Al4V, Mt.Haynes_188, Mt.Inc_X_750, Mt.Inc_600, Mt.Inc_718, Mt.Inc_A_286,Mt.Columbium_c103,
                    Mt.Copper_structural, Mt.D6AC_Steel]
    coating_list = [Mt.Coat_Custom,Mt.Copper, Mt.Narloy_Z, Mt.GRCop_84, Mt.Silica, Mt.Carbon]
    noz_mat_select = material_list[0]
    chamber_mat_select = material_list[0]
    
    
    #Mass: All other inputs are called from other variables already here    
    valv_mat_select = Mt.Al_7075_T6 #Valve material --> Does not need to appear in GUI
    t = 1e-3 #Default thickness for nozzle when thickness is below 1mm

    #Cost:
    tech_ready = 1.25
    exp_factor = 1.3
    learn_factor = 0 

    #Reliabiliy
     # Should be implemented for GUI user input
     # Fop should also be included, but commented out for now as to not mess with code
    cycles = ["Expander Cycle", "Staged Combustion Cycle", "Gas Generator Cycle"]
    Prop = ['LOX_LH2', 'LOX_RP1'] #Only used for de-rating or up-rating, so does not impact general reliability.
    N = 1   # number of engines, input range: as long as it is >0.
    #Fop = d.Thrust # operating thrust, so if theres de-rating or up-rating. 0.4*d.Thrust < Fop < 1.1*d.Thrust
     # Defaults
    delta = 0.1017
    Fref =  2278e3
        #Other parameters for reliablity [Fernandez, 2022] can be found in :
        # - Reliability_Data/Cycle_Data.csv 
        # - Reliability_Data/Propellant_Uprating_Data.csv
    val = False     # set to true for validation of figure with [fernandez, 2022]
    #Init
    def __init__(self,type):
        if(type==1):
            self.Cd = 1.0


# Propellant class
class Propellant:
    # Oxidizer
    def Oxigen(self):
        self.Ox_name = "LOX"  # Oxidizer name for rocketCEA
        self.Ox_composition = "O 2"  # Composition of oxidizer for rocketcea
        self.o_dens = 1141.0  # Oxidizer density
        self.ocp = 14307.0  # oxidizer cp
        self.h_ox = -12.979  # oxidizer enthalpy
        self.o_lamb = 1.0e-6
        self.o_nist_enthalpy_coef = [
            20.91,
            10.72,
            -2.02,
            0.1464,
            9.2457,
            5.338,
            237.62,
            0,
            31.33,
            -20.235,
            57.87,
            -36.51,
            -0.007374,
            -8.9035,
            246.79,
            0,
        ]  # for shomate equation
        self.omiu = 1.0e-6
        self.molarmass = 15.999  # g/mol

    # Fuel
    def LH(self):
        self.Fuel_name = "LH2"  # Fuel name for rocketCEA
        self.Fuel_composition = "H 2"  # Composition of fuel for rocketcea
        self.f_dens_l = 71.0  # liquid fuel density
        self.f_dens_g = 1.0  # gaseous fuel density
        self.f_gamma = 1.4  # fuel gamma
        self.fcp = 14307.0  # fuel cp
        self.h_fuel = -9.012  # fuel enthalpy
        self.R_f = 4.1573  # fuel gas constant
        self.f_lamb = 9.6e-6
        self.fmiu = 1.0e-6
        self.f_nist_enthalpy_coef = [
            43.31,
            -4.293,
            1.27243,
            -0.096876,
            -20.5339,
            -38.5151,
            162.08,
            0,
            33.066,
            -11.363,
            11.4328,
            -2.773,
            -0.15856,
            -9.981,
            172.71,
            0,
        ]  # for shomate equation
        self.heatingvalue = 119.96 * 10**6  # for the fuel only!
        self.molarmass = 2.016  # g/mol
        self.enthalpy_298_f = 0 #J/kh
        # Propellant

    gama = 1.4
    tq = 0.9  # characteristic chemical time of propellant
    Frozen_state = 0  # Frozen state of the propellant 0=chemical equilibrium flow, 1=frozen flow (from throat onwards)
    lstar = [0.76, 1.02]

    def __init__(self, type):
        match type:
            case 0:
                f_name = "LH"
                o_name = "LOX"
                self.LH()
                self.Oxigen()

            case 1:
                f_name = "CH4"
                o_name = "LOX"
                self.methane()
                self.Oxigen()
            case 2:
                f_name = "RP1"
                o_name = "LOX"
                self.Rp1()
                self.Oxigen()
            case 3:
                f_name = "Ethanol"
                o_name = "LOX"
                self.Ethanol()
                self.Oxigen()
            case 4:
                f_name = "Ethanol"
                o_name = "LOX"
                self.UDMH()
                self.Oxigen()

    def Rp1(self):
        self.Fuel_name = "RP1"  # Fuel name for rocketCEA
        self.Fuel_composition = "C 1 H 1.95"  # Composition of fuel for rocketcea
        self.f_dens_l = 804.59  # liquid fuel density
        self.f_dens_g = 10**5 / (
            298 * 8.31455 / self.molarmass
        )  # gaseous fuel density with ideal gas for standard conditions
        self.f_gamma = 1.24  # fuel gamma
        self.fcp = 1.88  # fuel cp
        self.h_fuel = -1758.456  # fuel enthalpy
        self.R_f = 4.75  # fuel gas constant
        self.f_lamb = 115e-3
        self.fmiu = 2.166e-6
        self.f_nist_enthalpy_coef = []  # for shomate equation
        self.heatingvalue = -43100  # for the fuel only!
        self.molarmass = 175  # g/mol
        self.enthalpy_298_f = 0.137899*10**6#J/kh

    def Ethanol(self):
        self.Fuel_name = "Ethanol"  # Fuel name for rocketCEA
        self.Fuel_composition = "C 2 H 6 O 1"  # Composition of fuel for rocketcea
        self.f_dens_l = 789  # liquid fuel density
        self.f_dens_g = 1.59  # gaseous fuel density
        self.f_gamma = 1.13  # fuel gamma
        self.fcp = 111.7 / 46.0684 * 10**3  # fuel cp
        self.h_fuel = -234  # fuel enthalpy
        self.R_f = 8.314 / 46.0684 * 10**3  # fuel gas constant
        self.f_lamb = 166.4e-3
        self.fmiu = 1199.4e-6
        self.f_nist_enthalpy_coef = []  # for shomate equation
        self.heatingvalue = -29672  # for the fuel only!
        self.molarmass = 46.07
        self.enthalpy_298_f = 5.01766*10**6#J/kh

    def UDMH(self):
        self.Fuel_name = "UDMH"  # Fuel name for rocketCEA
        self.Fuel_composition = "C 2 H 8 N 2"  # Composition of fuel for rocketcea

        self.f_dens_l = 793  # IDK YET  # liquid fuel density

        self.f_dens_g = 10**5 / (
            298 * 8.31455 / self.molarmass
        )  # gaseous fuel density with ideal gas for standard conditions
        self.f_gamma = 1.152  # fuel gamma
        self.fcp = 164.05 / 60.098 * 10**3  # fuel cp
        self.h_fuel = 83.3 * 10**3 / 60.0983  # fuel enthalpy
        self.R_f = 8.314 / 60.0983 * 10**3  # fuel gas constant
        self.f_lamb = 168e-3
        self.fmiu = 1.1071e-05
        self.f_nist_enthalpy_coef = []  # for shomate equation
        self.heatingvalue = -32928  # for the fuel only!
        self.molarmass = 60.1
        self.enthalpy_298_f = 0.886280*10**6#J/kh

    def methane(self):
        self.Fuel_name = "Methane"  # Fuel name for rocketCEA
        self.Fuel_composition = "C 1 H 4"  # Composition of fuel for rocketcea

        self.f_dens_l = 451.13  # IDK YET  # liquid fuel density

        self.f_dens_g = 0.64861  # gaseous fuel density
        self.f_gamma = 2.2313 / 1.7082  # fuel gamma
        self.fcp = 2231  # fuel cp
        self.h_fuel = -74.87  # fuel enthalpy
        self.R_f = 8.314 / 16.0425 * 10**3  # fuel gas constant
        self.f_lamb = 0.033937
        self.fmiu = 480e-6
        self.f_nist_enthalpy_coef = [
            85.81217,
            11.26467,
            -2.114146,
            0.13819,
            -26.42221,
            -153.5327,
            224.4143,
            -74.8731,
        ]  # for shomate equation
        self.heatingvalue = -55511  # for the fuel only!
        self.molarmass = 16.04  # g/mol
        self.enthalpy_298_f = 4.56*10**6#J/kh

    def NTO(self):
        self.Ox_name = "NTO"  # Oxidizer name for rocketCEA
        self.Ox_composition = "N 2 O 4"  # Composition of oxidizer for rocketcea
        self.o_dens = 1442.46  # Oxidizer density
        self.ocp = 1.549 * 10**3  # oxidizer cp
        self.h_ox = (
            -19.56 / 92.0110 * 10**3
        )  # oxidizer enthalpy in liquid phase. Gas phase is +9.08/92.0110*10**3
        self.o_lamb = 131 * 10 ** (-3)
        self.o_nist_enthalpy_coef = [
            128.6220,
            2.524345,
            -0.520883,
            0.03663,
            -11.55704,
            -59.22619,
            417.0444,
            9.078988,
        ]  # for shomate equation
        self.omiu = 393 * 10 ** (-6)
        self.molarmass = 92.010  # g/mol
        self.enthalpy_298_f = 1.036*10**6 #J/kh


#Data class
class Data:
    #Global in
    Thrust = 0.0 #[N]
    time = 0.0 #[s]
    Pa = 0.0 #[Pa] atmospheric pressure
    O_F = 0.0 #[-] mixture ratio

    #Global out
    Isp = 0.0 #[s]
    cstar = 0.0 #[m/s]
    m = 0.0 #[kg/s] mass flow
    Mprop = 0.0 #[kg] Total propellant mass
    Mtot = 0.0 #[kg] Total engine mass
    cost = 0.0 #[EUR] Total engine cost
    rel = 0.0 #[%] Reliability of entire system

    #Nozzle
    Pc = 0.0 
    m_nozz = 0.0 #[kg/s]
    L_div=0.0 # Length of the divergent part [m]
    L_con=0.0 # Length of the convergent part [m]
    L_total=0.0 # Total length of the nozzle [m]
    A_t=0.0 # Throat area [m^2]
    Eps=0.0 # Expansion ratio [-]
    Dt=0.0 # Throat diameter [m]
    De=0.0 # Exit diameter [m]

    #Turbo
    W_Opump = 0.0 #[W] Oxidizer pump power
    W_Fpump = 0.0 #[W] Fuel pump power
    W_turb = 0.0 #[W] Turbine power
    fuel_frac = 0.0 #[-] Fraction of fuel discarded in open cycles
    ptinj = 0.0 #[Pa] Total pressure at injector inlet
    dptop = 0.0 #[Pa] Total pressure rise over oxidizer pump
    dptfp = 0.0 #[Pa] Total pressure rise over fuel pump
    turbo_m = 0.0 #[kg/s] Total mass flow

    #Combustion
    h_comb = 0.0 #Conductive heat transfer coefficient in chamber
    Dc = 0.00 #Diameter of Combustion Chamber
    ThicknessChamber = 0.0  #Thickness of CC
    Chamber_L = 0.0 # Length of CC
    chambermass = 0.0 #mass of the CC
    Tc = 0.0 # [K] Combustion temperature
    Mnoz = 0.0 # [kg] Total nozzle mass

    #Cooling
    type_variable_nozzle = 0
    type_variable_chamber = 0
    
    maximum_thermal_stress=0 #Maximum thermal stress in Pa
    safety_factor_cooling=0 #Safety factor in respect to the stress
    max_temperature_inner=0 #maximum temperature of the inner wall [k]
    max_temperature_outer=0 #maximum temperature of the outer wall [k]
    dptcool_cooling=0   #pressure drop due to the cooling system [Pa]


    #Injector
    v_iox = 0.0     # Injection velocity oxidizer [m/s]
    v_if = 0.0      # Injection velocity fuel [m/s]
    D_f = 0.0       # Droplet size fuel [m]
    D_ox = 0.0      # Droplet size oxidizer [m]
    dp = 0.0        # Pressure drop over injector [Pa]
    eta_s = 0.0     # Stability criteria factor [-]
    m_ox = 0.0      # Mass flow oxidizer per orifice [kg/s]
    m_f = 0.0       # Mass flow fuel per orifice [kg/s]
    n_ox = 0.0      # Number of oxidizer orifices 
    n_f = 0.0       # Number of fuel orifices       
    P_D = 0.0       # Momentum ratio [-]
    
    #Ignitor
    Igniter_compound = 0.0 #mass

    #Material

    
    #Reliability
    Reliability = [] # List with floats, to account for multiple SG cyle variants 
    cycle = ""       # Data cycle selected
    Prop = ""        # Propellant selected (only relevant for derating)
    N = 0.0          # Number of engines

    def __init__(self, Th, t, p):
        self.Thrust = Th
        self.time = t
        self.Pa = p

