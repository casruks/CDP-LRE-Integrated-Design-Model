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

import matplotlib.pyplot as plt

default = Aux_classes.Default(0)


# Main Function
if __name__ == "__main__":
    # inicialise cooling
    regCool = Cooling.RegenerativeCool()

    Coolobj = Cooling.CoolingClass()
    Coolobj_c = Cooling.CoolingClass()

    Q = 0
    n = 130

    prop = Aux_classes.Propellant(0)
    # prop.fcp=40000
    Tw_ad_noz = np.array([3000 for i in range(n)])
    h_c_noz = [20000 for i in range(n)]
    t_noz = np.array([0.7112e-3 for i in range(n)])
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

    y_noz_cool = y_noz_cool / 2

    y_noz_cool = np.flip(y_noz_cool)

    x_noz_cool = np.array([L])
    err_out = 0

    chamber_mass = 1
    nozzle_mass = 1
    Aux_classes.Default.n = 1
    alpha = 2 * math.pi * 1 / Aux_classes.Default.n
    y_for_cooling_channel = np.amin(y_noz_cool)
    # A_nozzle = x_noz_cool[-1] * y_for_cooling_channel * alpha
    A_nozzle = [
        (alpha * y_noz_cool[i]) * x_noz_cool[-1] / len(y_noz_cool)
        for i in range(len(y_noz_cool))
    ]
    print(A_nozzle)
    print()
    # A_nozzle = [0.07 for i in range(len(y_noz_cool))]

    A_nozzle = [
        (alpha * y_noz_cool[i])
        * (x_noz_cool[-1] / len(y_noz_cool))
        / (
            math.cos(
                math.atan(
                    abs(y_noz_cool[i + 1] - y_noz_cool[i])
                    / (x_noz_cool[-1] / len(y_noz_cool))
                )
            )
        )
        for i in range(len(y_noz_cool) - 1)
    ]

    A_nozzle.append(
        (
            (alpha * y_noz_cool[-1])
            * (x_noz_cool[-1] / len(y_noz_cool))
            / (
                math.cos(
                    math.atan(
                        abs(y_noz_cool[-1] - y_noz_cool[-2])
                        / (x_noz_cool[-1] / len(y_noz_cool))
                    )
                )
            )
        )
    )
    print(A_nozzle)
    # A_nozzle = [x * 0.2 for x in A_nozzle]

    # print("A_nozzle", A_nozzle)
    # Ms.Inc_A_286.k = 4000000000000
    Aux_classes.Default.default_coating_thickness = 0
    Coolobj.Q = 0
    Coolobj_c.Q = 0
    Aux_classes.Default.Dr = 0.01
    type_variable_nozzle = -1
    main_material = Ms.Narloy_Z  # Ms.Inc_A_286
    coating = Ms.GRCop_84

    # main_material.k = 40000

    h_noz = np.array(
        [
            8733.65604295,
            8854.05554601,
            8978.12019254,
            9106.01881387,
            9237.93062674,
            9374.04602696,
            9514.56745446,
            9659.7103371,
            9809.70412108,
            9964.79339688,
            10125.23913032,
            10291.32000966,
            10463.33392034,
            10641.59956045,
            10826.45821091,
            11018.27567582,
            11217.44440931,
            11424.38584676,
            11639.5529588,
            11863.43304739,
            12096.5508033,
            12096.5508033,
            12426.21609851,
            12736.64207645,
            13028.41445151,
            13301.73020355,
            13556.52605556,
            13792.56436969,
            14009.49367689,
            14206.89329497,
            14384.30747069,
            14541.27228717,
            14677.33732899,
            14792.08336154,
            14885.1368356,
            14956.18175226,
            15004.96924782,
            15031.32514328,
            15035.15562496,
            15016.45116614,
            39612.26348024,
            38148.51344142,
            36922.6743845,
            35825.08053871,
            34811.85052942,
            33862.04075575,
            32963.02586085,
            32106.58429016,
            31287.0736599,
            30500.35293133,
            29743.28563261,
            29013.42281717,
            28308.77468082,
            27627.71284859,
            26968.85693341,
            26331.02038483,
            25713.15333411,
            25114.34732547,
            24533.75795865,
            24035.42846263,
            23000.14999869,
            22007.21407856,
            21068.03339276,
            20183.76987132,
            19352.45540089,
            18571.21438575,
            17836.90979833,
            17146.06150518,
            16496.34660973,
            15884.55713129,
            15308.03496408,
            14764.31935813,
            14251.13230056,
            13766.36986918,
            13308.0946454,
            12874.5188912,
            12463.99269175,
            12074.99323532,
            11706.1147484,
            11356.06044831,
            11023.62733357,
            10707.70288258,
            10407.25944707,
            10121.33836426,
            9849.05583614,
            9589.58663011,
            9342.19911779,
            9106.10883994,
            8880.6939408,
            8665.33295727,
            8459.44299049,
            8262.49839192,
            8073.99007365,
            7893.44799504,
            7720.43235437,
            7554.53234292,
            7395.3626827,
            7242.56292059,
            7095.79452619,
            6954.74045326,
            6819.10261405,
            6688.70138692,
            6563.07233165,
            6442.07052109,
            6325.46338335,
            6213.03233564,
            6104.57173627,
            5999.88805739,
            5898.79841039,
            5801.13089113,
            5706.72301873,
            5615.42193744,
            5527.08279352,
            5441.56932664,
            5358.752394,
            5278.51034996,
            5200.72791758,
            5125.29612902,
            5052.11176703,
            4981.07729989,
            4912.07894244,
            4845.07194302,
            4779.95185229,
            4716.64000511,
            4655.06186386,
            4595.14656527,
            4536.82681409,
            4480.0389078,
            4424.72211987,
            4370.81888135,
        ]
    )

    T_noz = np.array(
        [
            3421.8003346,
            3418.73181735,
            3415.71796259,
            3412.7638408,
            3409.87508551,
            3407.05796637,
            3404.31947309,
            3401.66741215,
            3399.11051842,
            3396.65858426,
            3394.32260935,
            3392.11497484,
            3390.04964645,
            3388.14241194,
            3386.41115942,
            3384.87620476,
            3383.56067767,
            3382.49097853,
            3381.69732073,
            3381.21437662,
            3381.08204976,
            3381.08204976,
            3381.23759976,
            3381.51434817,
            3381.86827386,
            3382.2473836,
            3382.59268539,
            3382.8392831,
            3382.91765872,
            3382.75516593,
            3382.27772784,
            3381.41170822,
            3380.08590507,
            3378.23359913,
            3375.79457791,
            3372.71704826,
            3368.95934917,
            3364.49138086,
            3359.29567648,
            3353.36805919,
            3516.67773468,
            3568.16551476,
            3610.71641827,
            3648.55844544,
            3683.42242176,
            3716.14990486,
            3747.25990354,
            3777.09809961,
            3805.90159668,
            3833.85076021,
            3861.08079477,
            3887.69769282,
            3913.78720985,
            3939.41607115,
            3964.63940212,
            3989.50174177,
            4014.04165724,
            4038.28606799,
            4062.26220442,
            4094.87668835,
            4059.40676617,
            4035.38605011,
            4017.38217014,
            4003.20302027,
            3991.77530542,
            3982.44341324,
            3974.77248416,
            3968.49134867,
            3963.30538675,
            3959.07472029,
            3955.66186588,
            3952.95412801,
            3950.85860486,
            3949.29705435,
            3948.20213918,
            3947.51550989,
            3947.18673034,
            3947.17106291,
            3947.42917887,
            3947.92544509,
            3948.6284228,
            3949.51022446,
            3950.54434734,
            3951.70894357,
            3952.98234583,
            3954.34629989,
            3955.75545044,
            3957.25227157,
            3958.79288119,
            3960.36532351,
            3961.95953585,
            3963.56379314,
            3965.16991592,
            3966.77017136,
            3968.35768774,
            3969.92616637,
            3971.47018322,
            3972.98492747,
            3974.46664366,
            3975.91159738,
            3977.31700489,
            3978.63719359,
            3979.95610931,
            3981.22961918,
            3982.45655376,
            3983.63615059,
            3984.76794512,
            3985.85137634,
            3986.88671806,
            3987.87421983,
            3988.81413089,
            3989.70706981,
            3990.55397012,
            3991.35552036,
            3992.11290894,
            3992.82704265,
            3993.49919971,
            3994.13062401,
            3994.72269405,
            3995.27665759,
            3995.8292072,
            3996.31166097,
            3996.76031925,
            3997.17650442,
            3997.56157279,
            3997.91703294,
            3998.24436588,
            3998.54481982,
            3998.81993442,
            3999.07089961,
        ]
    )

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
    A_chamber = [Chamber_L * Dc / 2 * alpha]
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
        0.00000005,  # Aux_classes.Default.Dr,
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

# plt.plot(y_noz_cool)
# plt.show()

print(len(y_noz_cool))
