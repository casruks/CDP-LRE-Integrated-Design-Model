import Cooling
import Materials as Mt

p_a = 1.0e5
Thrust = 15000
Thurst_time = 60


# Default values
class Default:
    # Tolerances
    pres_tol = 0.01
    toll_c_star = 0.01
    toll_F_obj = 0.01
    Max_iterations_mass_flow = 10000
    toll_P_adapted = 0.01
    Safety_factor = 1.3

    # Seeds
    Pres = 1e6
    inj_vel = 15

    # Injectors
    Cd = 0.7

    # Nozzle
    Nozzle_type = 0
    MR = 0
    De_max = 2.5
    De_turbine_noz_max = 2.5
    Theta_con = 60
    Theta_conical = 15
    Theta_bell = 55
    TH_exit_bell = 3
    R_u_ratio = 1

    # Turbomachinery
    cycle_type = "EX"
    Eff_t = 0.6  # Turbine efficiency
    Eff_p = 0.6  # Pump efficiency
    Eff_m = 0.95  # Mechanical efficiency between turbine and pumps
    p_to = 1.0e5  # oxidizer tank storage pressure
    ptf = 1.0e5  # Fuel tank oxidizer pressure
    Wmotor = 1.0e6  # Power of the electric motor

    # Combustion chamber
    SF = 1.0

    # Cooling
    Dr = 0.01
    A = 0.0003
    T_fuel_tanks = 20
    T_ox_tanks = 60

    # Igniters
    ignburntime = 4

    # Material
    material = "This"

    # Init
    def __init__(self, type):
        if type == 1:
            self.Cd = 1.0


default = Default(0)

# Propellant class
class Propellant:
    # Oxidizer
    Ox_name = "LOX"  # Oxidizer name for rocketCEA
    Ox_composition = "O 2"  # Composition of oxidizer for rocketcea
    o_dens = 1141.0  # Oxidizer density
    ocp = 14307.0  # oxidizer cp
    h_ox = -12.979  # oxidizer enthalpy
    o_lamb = 1.0e-3
    omiu = 1.0e-6

    # Fuel
    Fuel_name = "LH2"  # Fuel name for rocketCEA
    Fuel_composition = "H 2"  # Composition of fuel for rocketcea
    f_dens_l = 71.0  # liquid fuel density
    f_dens_g = 1.0  # gaseous fuel density
    f_gamma = 1.4  # fuel gamma
    fcp = 14307.0  # fuel cp
    h_fuel = -9.012  # fuel enthalpy
    R_f = 4.1573  # fuel gas constant
    f_lamb = 1.0e-3
    fmiu = 1.0e-6
    MR = 3  # mixture ratio

    Frozen_state = 0

    # Propellant
    gama = 1.4
    tq = 0.9  # characteristic chemical time of propellant

    def __init__(self, type):
        match type:
            case 0:
                f_name = "LH"
                o_name = "LOX"

            case 1:
                f_name = "CH4"


prop = Propellant(0)
bool = (
    0  # this variable is used to show the combustor function we are in the first loop
)
# Main Function
if __name__ == "__main__":
    Q = 0
    Tw_ad_noz = [3000, 3000, 3000]
    h_c_noz = [1300000, 1400000, 1500000]
    t_noz = [0.001, 0.001, 0.001]
    L = 1
    Re_t = 1000000
    m = 100
    O_F = 2
    A = default.A
    #default.A = [A for i in range(len(Tw_ad_noz))]
    # Temperature after cooling
    T_w_after_cooling = 0

    # inicialise cooling
    regCool = Cooling.RegenerativeCool()

    Tf_cool, T_w_after_cooling, dptcool = regCool.Run1D(
        Tw_ad_noz,
        h_c_noz,
        t_noz,
        prop,
        Mt.Rhenium,
        default.Dr,
        default.A,
        default.T_fuel_tanks,
        Re_t,
        m / (1 + O_F),
        L,
    )
    print("Tf_cool: ", Tf_cool)
    print("T_w_after_cooling: ", T_w_after_cooling)
    print("Q: ", regCool.Q)