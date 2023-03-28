import Cooling


# import Materials as Mt
import numpy as np
import math

# import Mass as Ms
import statistics

# import Materials_2 as Mt
import Aux_classes
import Materials as Ms

# import itertools

default = Aux_classes.Default(0)


# Main Function
if __name__ == "__main__":
    # inicialise cooling
    regCool = Cooling.RegenerativeCool()

    Coolobj = Cooling.CoolingClass()
    Coolobj_c = Cooling.CoolingClass()

    Q = 0
    n = 2000
    prop = Aux_classes.Propellant(0)
    Tw_ad_noz = np.array([3000 for i in range(n)])
    h_c_noz = [6000 for i in range(n)]
    t_noz = np.array([0.003 for i in range(n)])
    L = 4.3
    T_w_after_cooling = 0
    # Ms.Rhenium.OpTemp_u = 700
    Aux_classes.Default.T_fuel_tanks = 20
    # c = 340
    operationtime = 30000
    m_casing = 1000
    eps = 0.85
    overwriteA = False
    case_run = 0
    Tc = 3000  # Tc = 3588
    h_comb = 34000
    ThicknessChamber = 0.0508
    m_nozz = 515
    Chamber_L = 2
    Dc = 0.4
    O_F = 6.24
    err_chamber_cooling = 0
    warn_err_chamber_cooling = 0
    y_noz_cool = np.array(
        [(2.30 - 0.26) / (n - 1) * (n - 1 - i) + 0.26 for i in range(n)]
    )
    x_noz_cool = np.array([L])
    err_out = 0

    chamber_mass = 1
    nozzle_mass = 1

    alpha = 2 * math.pi * 1 / 1
    y_for_cooling_channel = np.amin(y_noz_cool)
    # A_nozzle = x_noz_cool[-1] * y_for_cooling_channel * alpha
    A_nozzle = sum(
        [
            (alpha * y_noz_cool[i]) * x_noz_cool[-1] / len(y_noz_cool)
            for i in range(len(y_noz_cool))
        ]
    )
    # print("A_nozzle", A_nozzle)
    # Ms.Inc_A_286.k = 4000000000000
    Aux_classes.Default.default_coating_thickness = 0.12e-3
    Coolobj.Q = 0
    Coolobj_c.Q = 0
    Aux_classes.Default.Dr = 0.01
    type_variable_nozzle = -1
    main_material = Ms.Haynes_188
    coating = Ms.GRCop_84
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
        main_material.heat_cap,
        operationtime,
        nozzle_mass,
        Aux_classes.Default.eps,
        Tw_ad_noz,
        h_c_noz,
        t_noz,
        np.array(
            [Aux_classes.Default.default_coating_thickness for i in range(len(t_noz))]
        ),
        prop,
        main_material,
        coating,
        0.0005,  # Aux_classes.Default.Dr,
        A_nozzle,
        Aux_classes.Default.T_fuel_tanks,
        m_nozz / (1.0 + O_F) / Aux_classes.Default.n,
        x_noz_cool[-1],
        y_noz_cool,
        True,
        default.regenerative_case,
    )

    print("\nRun Nozzle")
    print("Tf_cool: ", Tf_cool)
    print(
        "p loss",
        dptcool,
    )
    print("Tw_wall_calculated", max(Tw_wall_nozzle_calculated))
    print("Outer wall temperature", max(T_outer_wall_nozzle))
    print("errors:", err_nozzle_cooling)
    print("warning", warn_nozzle_cooling)
    print("type: ", type_variable_nozzle)
    print("Q: ", Coolobj.Q)
    print()

    Aux_classes.Default.n = 1
    perimeter_percentage = 1
    alpha = 2 * math.pi * perimeter_percentage / Aux_classes.Default.n
    # A_chamber=[Chamber_L * Dc * 2*math.pi]
    y_for_cooling_channel = np.amin(y_noz_cool)
    print("Running for chamber")
    A_chamber = Chamber_L * Dc / 2 * alpha
    ThicknessChamber = 0.02
    (
        Tf_cool,
        Tw_wall_chamber_calculated,
        dptcool_c,
        _,
        type_variable_chamber,
        T_outer_wall_chamber,
        err_chamber_cooling,
        warn_chamber_cooling,
    ) = Coolobj_c.Run_cooling(
        default.T0,
        Ms.Columbium_c103.heat_cap,
        operationtime,
        chamber_mass,
        eps,
        np.array([Tc]),
        np.array([h_comb]),
        np.array([ThicknessChamber]),
        np.array([0]),
        prop,
        Ms.Columbium_c103,
        Ms.Rhenium,
        0.0005,  # Aux_classes.Default.Dr,
        np.array([A_chamber]),
        Tf_cool,
        m_nozz / (1.0 + O_F) / Aux_classes.Default.n,
        Chamber_L,
        np.array([Dc / 2]),
        True,
        Aux_classes.Default.regenerative_case,
    )
    Coolobj_c.Q = Coolobj_c.Q * Aux_classes.Default.n
    print("\nRun Chamber")
    print("Tf_cool: ", Tf_cool)
    print(
        "p loss",
        dptcool_c,
    )
    print("Tw_wall_calculated", Tw_wall_chamber_calculated[-1])
    print("Outer wall temperature", T_outer_wall_chamber)
    print("errors:", err_chamber_cooling)
    print("warning", warn_chamber_cooling)
    print("type: ", type_variable_chamber)
    print()

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
    y_noz_cool,
    chamber_D / 2,
    err_out,
)
print("max_temperature_inner", max_temperature_inner)
print("max_temperature_outer", max_temperature_outer)
print("maximum_thermal_stress", maximum_thermal_stress)
print("safety_factor_cooling", safety_factor_cooling)
print("errors: ", err_nozzle_cooling, " ", err_chamber_cooling, " ", err_out)
