import Cooling

# import Materials as Mt
import numpy as np
import math

# import Mass as Ms
import statistics

# import Materials_2 as Mt
import Aux_classes
import Materials as Ms
import itertools

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
    Dr = 0.005  # [m] hydralic diameter of the coolant channel
    A = 0.0003  # [m2] area of contact for each segment of the cooling
    T_fuel_tanks = 20  # [K] temperature of the fuel tanks, considered the inicial coolant temperature
    # T_ox_tanks = 60 #[K] temperature of the oxidiser tanks
    n = 1  # number of coolant chanels
    # default_coating = Mt.Materials(
    #   "default_coating", 0.0, 0.0, 0.0, 0.0, 1, 0
    #  )  # default coolant
    default_coating_thickness = 0  # default coolant thickness
    T0 = 293.5  # [k] default inicial temperature
    eps = 0.85  # default emissivity
    overwriteA = False  # option to overwrite the surface area calculated by the program with the input variable A, given by the user or default class
    regenerative_case = 0  # option of which function to use in regenerative cooling; 0 corresponds to the explicit function Run1D()
    operationtime = (
        10000000000000000  # [s] default operation time (large to imply infinite time)
    )

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
    n = 40
    Tw_ad_noz = np.array([200 for i in range(n)])
    h_c_noz = [1200 for i in range(n)]
    t_noz = np.array([0.005 for i in range(n)])
    # Tw_ad_noz = 9000
    # h_c_noz = 1200
    # t_noz = 0.001
    L = 400000.3
    # Re_t = 1000000
    m = 2500000
    O_F = 600000
    # A = default.A
    # default.A = [A for i in range(len(Tw_ad_noz))]
    # Temperature after cooling
    T_w_after_cooling = 0
    y = np.array([(200000.30 - 0.26) / (n - 1) * (n - 1 - i) + 0.26 for i in range(n)])

    # y = [(2.30 - 0.26) / 2 for i in range(n)]

    # inicialise cooling
    regCool = Cooling.RegenerativeCool()

    Ms.Rhenium.OpTemp_u = 7
    Ms.Rhenium.k = 4500000
    Aux_classes.Default.T_fuel_tanks = 20
    c = 340
    operationtime = 3000000000000000
    m_casing = 1000
    eps = 1.2
    overwriteA = False
    case_run = 0

    Coolobj = Cooling.CoolingClass()
    Coolobj_c = Cooling.CoolingClass()
    Tc = 3000
    h_comb = 1500
    ThicknessChamber = 0.003
    m_nozz = 25
    Chamber_L = 0.8
    Dc = 0.4
    O_F = 2000000000000000000000000000000000000
    Tf_cool = 1

    # chamber_mass = Ms.Chamber_mass(Dc, Chamber_L, ThicknessChamber, Ms.Rhenium)
    err_chamber_cooling = 0
    warn_err_chamber_cooling = 0

    y_noz_cool = np.array(
        [(2.30 - 0.26) / (n - 1) * (n - 1 - i) + 0.26 for i in range(n)]
    )
    x_noz_cool = np.array([3])
    # cool = Cooling.CoolingClass()
    # nozzle_mass = Ms.Nozzle_mass(
    #    x_noz_cool, y_noz_cool, statistics.mean(t_noz), Ms.Rhenium
    # )
    Ms.Rhenium.OpTemp_u = 2
    err_out = 0

    chamber_mass = 500
    nozzle_mass = 500

    # Ranges
    T0_range = [0.001, 100, 100000]  # range(1000)
    c_range = [1, 10, 100000]  # range(1, 1001)
    time_range = [0.01, 1, 1000000]  # range(1, 200)
    nozzle_mass_range = [1, 50, 500, 10000]  # range(1, 5000)
    eps_range = [0.1, 0.5, 1]  # [i / 10 for i in range(1, 10)]
    number_points = 110
    Tw_ad_range = [
        # np.array([j for i in range(number_points)]) for j in range(1, 5000, 100)
        np.array([1 for i in range(number_points)]),
        np.array([100 for i in range(number_points)]),
        np.array([900 for i in range(number_points)]),
        np.array([2000 for i in range(number_points)]),
        np.array([200000 for i in range(number_points)])
    ]
    h_c_range = [
        # np.array([j for i in range(number_points)]) for j in range(1, 50000, 100)
        np.array([0.01 for i in range(number_points)]),
        np.array([10 for i in range(number_points)]),
        np.array([100 for i in range(number_points)]),
        np.array([10000 for i in range(number_points)]) ,
        np.array([100000 for i in range(number_points)])
    ]
    t_range = [
        # np.array([j / 1000 for i in range(number_points)]) for j in range(1, 4000)
        np.array([0.0001 for i in range(number_points)]),
        np.array([0.01 for i in range(number_points)]),
        np.array([10 for i in range(number_points)])
    ]
    coating_thickness_range = [
        0.0001,
        0.001,
        0.1,
        10,
        1000,
    ]  # [j / 1000 for j in range(1, 4000)]
    prop_range = [Aux_classes.Propellant(i) for i in range(5)]
    material_range = [Ms.Rhenium, Ms.Al_7075_T6]
    coating_range = [Ms.Copper]
    Dr_range = [0.01, 1, 10, 100, 1000]  # [j / 1000 for j in range(1, 1000)]
    A = [0.01, 1, 10, 100, 1000]  # [j / 10 for j in range(1, 300, 5)]
    T_fuel_range = range(0, 300, 20)
    m_range = [1, 10, 100, 500, 1000, 5000]  # [j / 10 for j in range(1, 5000, 10)]
    O_F_range = [1, 4, 10, 100]  # range(1, 40)
    n_range = [1, 10, 100, 1000]  # range(1, 100)
    x_len_range = [0.001, 0.1, 5, 100]  # [j / 10 for j in range(1, 50)]
    y_range = [
        0.001,
        0.1,
        5,
        100,
    ]  # [np.array([j / 100 for i in range(number_points)]) for j in range(1, 700)]

    (Aux_classes.Default.overwriteA) = True
    Aux_classes.Default.regenerative_case = 0

    Range = itertools.product(
        T0_range,
        c_range,
        time_range,
        nozzle_mass_range,
        eps_range,
        Tw_ad_range,
        h_c_range,
        t_range,
        coating_thickness_range,
        prop_range,
        material_range,
        coating_range,
        Dr_range,
        A,
        T_fuel_range,
        m_range,
        O_F_range,
        n_range,
        x_len_range,
        y_range,
    )

    for (
        Aux_classes.Default.T0,
        c,
        Aux_classes.Default.operationtime,
        nozzle_mass,
        Aux_classes.Default.eps,
        Tw_ad_noz,
        h_c_noz,
        t_noz,
        Aux_classes.Default.default_coating_thickness,
        prop,
        Ms.Rhenium,
        Aux_classes.Default.default_coating,
        Aux_classes.Default.Dr,
        Aux_classes.Default.A,
        Aux_classes.Default.T_fuel_tanks,
        m_nozz,
        O_F,
        Aux_classes.Default.n,
        x_noz_cool[-1],
        y_noz_cool,
    ) in Range:
        (
            Tf_cool,
            Tw_wall_nozzle_calculated,
            dptcool,
            _,
            type_variable_nozzle,
            T_outer_wall_nozzle,
            err_nozzle_cooling,
            warn_nozzle_cooling,
        ) = Coolobj.Run_cooling(
            Aux_classes.Default.T0,
            c,
            Aux_classes.Default.operationtime,
            nozzle_mass,
            Aux_classes.Default.eps,
            Tw_ad_noz,
            h_c_noz,
            t_noz,
            np.array(
                [
                    Aux_classes.Default.default_coating_thickness
                    for i in range(len(t_noz))
                ]
            ),
            prop,
            Ms.Rhenium,
            Aux_classes.Default.default_coating,
            Aux_classes.Default.Dr,
            Aux_classes.Default.A,
            Aux_classes.Default.T_fuel_tanks,
            m_nozz / (1.0 + O_F) / Aux_classes.Default.n,
            x_noz_cool[-1],
            y_noz_cool,
            Aux_classes.Default.overwriteA,
            Aux_classes.Default.regenerative_case,
        )

    print("\nRun", case_run)
    print("Tf_cool: ", Tf_cool)
    print(
        "p loss",
        dptcool,
    )
    print("Tw_wall_calculated", Tw_wall_nozzle_calculated[-1])
    print("Outer wall temperature", T_outer_wall_nozzle[-1])
    print("errors:", err_nozzle_cooling)
    print("warning", warn_nozzle_cooling)
    print("type: ", type_variable_nozzle)
"""  
    Aux_classes.Default.n = 20
    perimeter_percentage = 0.4
    alpha = 2 * math.pi * perimeter_percentage / Aux_classes.Default.n
    # A_chamber=[Chamber_L * Dc * 2*math.pi]
    y_for_cooling_channel = np.amin(y_noz_cool)
    A_chamber = [L * y_for_cooling_channel * alpha]
    print("Running for chamber")
    (
        Tf_cool,
        Tw_wall_chamber_calculated,
        dptcool_c,
        _,
        type_variable_chamber,
        T_outer_wall_chamber,
        err_nozzle_cooling,
        warn_chamber_cooling,
    ) = Coolobj_c.Run_cooling(
        Aux_classes.Default.T0,
        c,
        Aux_classes.Default.operationtime,
        chamber_mass,
        Aux_classes.Default.eps,
        np.array([Tc]),
        np.array([h_comb]),
        np.array([ThicknessChamber]),
        np.array([Aux_classes.Default.default_coating_thickness]),
        prop,
        Ms.Rhenium,
        Aux_classes.Default.default_coating,
        Aux_classes.Default.Dr,
        np.array(A_chamber),
        Tf_cool,
        m_nozz / (1.0 + O_F) / Aux_classes.Default.n,
        Chamber_L,
        np.array([Dc / 2]),
        True,
        Aux_classes.Default.regenerative_case,
    )
    print("\nRun", case_run)
    print("Tf_cool: ", Tf_cool)
    print(
        "p loss",
        dptcool_c,
    )
    print("Tw_wall_calculated", Tw_wall_chamber_calculated[-1])
    print("Outer wall temperature", T_outer_wall_chamber[-1])
    print("errors:", err_chamber_cooling)
    print("warning", warn_chamber_cooling)

    max_temperature_inner = np.maximum(
        np.max(T_outer_wall_chamber), np.max(T_outer_wall_nozzle)
    )
    max_temperature_outer = np.maximum(
        np.max(Tw_wall_chamber_calculated), np.max(Tw_wall_nozzle_calculated)
    )
    maximum_thermal_stress = (
        6.12
        * 10**-6
        * Ms.Rhenium.Emod
        * np.maximum(
            np.max(
                np.array(Tw_wall_chamber_calculated) - np.array(T_outer_wall_chamber)
            ),
            np.max(np.array(Tw_wall_nozzle_calculated) - np.array(T_outer_wall_nozzle)),
        )
    )
    safety_factor_cooling = Ms.Rhenium.yieldstress_l / maximum_thermal_stress

    print(
        "maximum deltaT",
        np.maximum(
            np.max(
                np.array(Tw_wall_chamber_calculated) - np.array(T_outer_wall_chamber)
            ),
            np.max(np.array(Tw_wall_nozzle_calculated) - np.array(T_outer_wall_nozzle)),
        ),
    )

chamber_D = 0.8
ThicknessChamber = 0.1
Aux_classes.Default.default_coating_thickness = 1 * 10**-3
Ms.Rhenium.thr_exp = 7.25
Aux_classes.Default.default_coating.thr_exp = 2
Aux_classes.Default.default_coating.Emod = 2 * 10**6
maximum_thermal_stress = 0

(
    maximum_thermal_stress,
    safety_factor_cooling,
    max_temperature_inner,
    max_temperature_outer,
    err_out,
) = Cooling.outputs(
    T_outer_wall_chamber,
    T_outer_wall_nozzle,
    Tw_wall_chamber_calculated,
    Tw_wall_nozzle_calculated,
    Ms.Rhenium,
    Ms.Rhenium,
    Aux_classes.Default.default_coating,
    np.array([ThicknessChamber for i in range(len(t_noz))]),
    t_noz,
    np.array(
        [Aux_classes.Default.default_coating_thickness for i in range(len(t_noz))]
    ),
    y,
    chamber_D / 2,
    err_out,
)
print("max_temperature_inner", max_temperature_inner)
print("max_temperature_outer", max_temperature_outer)
print("maximum_thermal_stress", maximum_thermal_stress)
print("safety_factor_cooling", safety_factor_cooling)
print("errors: ", err_chamber_cooling, " ", err_nozzle_cooling, " ", err_out)

"""
"""     
    (
        T_co_calculated,
        Tw_wall_calculated,
        ploss,
        m_flow_fuel,
        err,
        warn,
    ) = cool.Run_cooling(
        275,
        c,
        operationtime,
        m_casing,
        eps,
        Tw_ad_noz,
        h_c_noz,
        t_noz,
        default.default_coating_thickness,
        prop,
        Ms.Rhenium,
        default.default_coating,
        default.Dr,
        default.A,
        default.T_fuel_tanks,
        m,
        L,
        y,
        overwriteA,
        case_run,
    )

     (
        T_co_calculated,
        Tw_wall_calculated,
        ploss,
        m_flow_fuel,
    ) = regCool.Main_regenerative_run_function(
        Tw_ad_noz,
        h_c_noz,
        t_noz,
        default.default_coating_thickness,
        prop,
        Ms.Rhenium,
        default.default_coating,
        default.Dr,
        default.A,
        default.T_fuel_tanks,
        m,
        L,
        y,
        overwriteA,
        case_run,
    ) """

""" Tf_cool, dptcool = regCool.Run_for_Toperating1D(
        Tw_ad_noz,
        h_c_noz,
        t_noz,
        default.default_coating_thickness,
        prop,
        Ms.Rhenium,
        default.default_coating,
        default.A,
        default.T_fuel_tanks,
        m,
        L,
        y,
    )
    print("Toperating1D")
    print("D: ", regCool.D)
    print("Tf_cool: ", Tf_cool)
    # print("T_w_after_cooling: ", T_w_after_cooling)
    print("Q: ", regCool.Q)
    print("p loss", dptcool[len(dptcool) - 1])

    # Tf_cool, Twall, dptcool = regCool.Run1D(
    # Tw_ad_noz,
    # h_c_noz,
    # t_noz,
    #  prop,
    # Ms.Rhenium,
    # 0.02,
    # default.A,
    # m,
    #  m * 0.02 / prop.fmiu * 1 / (math.pi * (0.02 / 2) ** 2),
    #  m,
    #  L,
# )

# print("D: ", regCool.D)
# print("Tf_cool: ", Tf_cool)
# print("T_w_after_cooling: ", T_w_after_cooling)
# print("Q: ", regCool.Q)
# print("Twall: ", Twall)
print("\nToperating0D")
default.A = sum([L * math.pi * 2 * y[i] for i in range(len(y))])
Tf_cool, dptcool = regCool.Run_for_Toperating0D(
    Tw_ad_noz[0],
    h_c_noz[0],
    t_noz[0],
    prop,
    Ms.Rhenium,
    default.A,
    default.T_fuel_tanks,
    m,
    L,
)

print("D: ", regCool.D)
print("Tf_cool: ", Tf_cool)
# print("T_w_after_cooling: ", T_w_after_cooling)
print("Q: ", regCool.Q)


Dr = 0.01

default.A = [L / len(y) * math.pi * 2 * y[i] for i in range(len(y))]
# print(default.A)
Re = 4 * m / (prop.fmiu) * 4 / (math.pi * Dr)
Tf_cool, Tw_wall_calculated, dptcool = regCool.Run1D(
    Tw_ad_noz,
    h_c_noz,
    t_noz,
    default.default_coating_thickness,
    prop,
    Ms.Rhenium,
    default.default_coating,
    Dr,
    default.A,
    default.T_fuel_tanks,
    m,
    L,
    y
)
print("\nRun1D")
print("Tf_cool: ", Tf_cool)
print("p loss", dptcool)
print("Tw_wall_calculated", Tw_wall_calculated[-1])
# print("y", y)


m1, Tf_cool, Tw_wall_calculated, ploss = regCool.Run1D_iterative_for_m(
    Tw_ad_noz,
    h_c_noz,
    t_noz,
    default.default_coating_thickness,
    prop,
    Ms.Rhenium,
    default.default_coating,
    Dr,
    default.T_fuel_tanks,
    L,
    y
)

print("\nRun for mass flow")
print("m", m1)
print("Tf_cool: ", Tf_cool[-1])
print("p loss", dptcool)
print("Tw_wall_calculated", Tw_wall_calculated[-1]) """
