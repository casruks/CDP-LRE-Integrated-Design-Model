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
    t_noz = np.array([0.7e-3 for i in range(n)])
    L = 3
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

    # print("sum(A_nozzle)", sum(A_nozzle))
    # Ms.Inc_A_286.k = 4000000000000
    Aux_classes.Default.default_coating_thickness = 0
    Coolobj.Q = 0
    Coolobj_c.Q = 0
    Aux_classes.Default.Dr = 0.01
    type_variable_nozzle = -1
    main_material = Ms.Narloy_Z
    coating = Ms.Narloy_Z

    # main_material.k = 40000

    h_noz = np.array(
        [
            9112.84407274,
            9238.42418163,
            9367.82658003,
            9501.22728499,
            9638.81313392,
            9780.78261031,
            9927.34674379,
            10078.73009169,
            10235.17181059,
            10396.92682666,
            10564.26711516,
            10737.48309986,
            10916.88518457,
            11102.80543003,
            11295.59939066,
            11495.64812657,
            11703.36040784,
            11919.17512866,
            12143.56395018,
            12377.03419109,
            12620.13198489,
            12620.13198489,
            12963.89129632,
            13287.55824096,
            13591.74387741,
            13876.65402494,
            14142.22337017,
            14388.20513545,
            14614.23426702,
            14819.87400457,
            15004.65150051,
            15168.08586849,
            15309.71073893,
            15429.09263247,
            15525.84599948,
            15599.64548685,
            15650.23581174,
            15677.43950372,
            15681.16269527,
            15661.39908264,
            41214.004097,
            39549.41015167,
            38169.73277153,
            36940.76882078,
            35810.82761475,
            34755.64777315,
            33760.62330712,
            32816.32756204,
            31916.2356384,
            31055.53864328,
            30230.58959743,
            29438.49269805,
            28676.87226804,
            27943.76801196,
            27237.48779124,
            26556.57073847,
            25899.7110301,
            25265.74522215,
            24653.60674279,
            24126.9491134,
            23013.89597886,
            21951.29712285,
            20952.11232904,
            20016.54092898,
            19142.44501846,
            18325.56411209,
            17561.85588693,
            16847.37331786,
            16178.38333967,
            15551.41275564,
            14963.2605164,
            14410.97942057,
            13891.87572013,
            13403.47749415,
            12943.52513897,
            12509.94809038,
            12100.8521732,
            11714.49898275,
            11349.2929961,
            11003.77087106,
            10676.58528931,
            10366.49745408,
            10072.40236199,
            9793.16853927,
            9527.87054268,
            9275.59973427,
            9035.54069793,
            8806.9172401,
            8589.01683162,
            8381.17930848,
            8182.7924563,
            7993.28631793,
            7812.13266133,
            7638.83899403,
            7472.94780296,
            7314.03230948,
            7161.69457305,
            7015.66337166,
            6875.39160924,
            6740.65752094,
            6611.15873917,
            6486.6127003,
            6366.75453503,
            6251.3374059,
            6140.12884279,
            6032.9118137,
            5929.48214632,
            5829.64877077,
            5733.2319874,
            5640.06294148,
            5549.98325619,
            5462.84372625,
            5378.50399113,
            5296.83190941,
            5217.7032848,
            5141.00073192,
            5066.59121112,
            4994.41541934,
            4924.35252554,
            4856.30960309,
            4790.19878294,
            4725.93718295,
            4663.44625938,
            4602.65184834,
            4543.48361034,
            4485.87494845,
            4429.76280049,
            4375.08730069,
            4321.79181881,
            4269.82242796,
        ]
    )

    T_noz = np.array(
        [
            3436.99184066,
            3433.87148968,
            3430.80661803,
            3427.80239459,
            3424.86456394,
            3421.99952134,
            3419.21439874,
            3416.51716395,
            3413.91673513,
            3411.42311325,
            3409.04753578,
            3406.80265551,
            3404.70274898,
            3402.76396037,
            3401.00458734,
            3399.44541739,
            3398.1101246,
            3397.02573913,
            3396.22320482,
            3395.73804339,
            3395.61114874,
            3395.61114874,
            3395.7864735,
            3396.09506971,
            3396.49274268,
            3396.92705331,
            3397.33826652,
            3397.66042833,
            3397.82263993,
            3397.75055382,
            3397.36808657,
            3396.59931798,
            3395.37052584,
            3393.61228833,
            3391.26157316,
            3388.26372458,
            3384.57425692,
            3380.16036742,
            3375.00209088,
            3369.09303497,
            3532.79463673,
            3590.06392726,
            3637.1449321,
            3679.03572081,
            3717.70586409,
            3754.09601117,
            3788.7842385,
            3822.15034874,
            3854.4554115,
            3885.89652478,
            3916.61952789,
            3946.73863496,
            3976.34516115,
            4005.50873635,
            4034.28636317,
            4062.72231704,
            4090.85381875,
            4118.7062316,
            4146.3033588,
            4182.8496551,
            4142.96394773,
            4116.47740285,
            4096.82874483,
            4081.49918063,
            4069.18961529,
            4059.18512843,
            4050.99314435,
            4044.26650788,
            4038.74805426,
            4034.2404746,
            4030.5867072,
            4027.65961352,
            4025.35368752,
            4023.58049145,
            4022.26476755,
            4021.34187518,
            4020.75517281,
            4020.45589591,
            4020.40120621,
            4020.55284631,
            4020.87759102,
            4021.3461447,
            4021.90333258,
            4022.58523142,
            4023.34023845,
            4024.15200391,
            4025.00258648,
            4025.87890674,
            4026.76891884,
            4027.66220616,
            4028.54906097,
            4029.42217866,
            4030.27433948,
            4031.10009614,
            4031.89439317,
            4032.65344353,
            4033.37417627,
            4034.00994103,
            4034.64614269,
            4035.2380331,
            4035.78490577,
            4036.28624406,
            4036.74223603,
            4037.1526215,
            4037.51834833,
            4037.83985567,
            4038.11839918,
            4038.35493213,
            4038.55086597,
            4038.70770386,
            4038.82677664,
            4038.9097835,
            4038.95835848,
            4038.97429104,
            4038.95914733,
            4038.91472308,
            4038.87900668,
            4038.78194106,
            4038.6608028,
            4038.51728184,
            4038.35310571,
            4038.16973638,
            4037.96893865,
            4037.75209828,
            4037.52079476,
            4037.27646625,
            4037.02052871,
            4036.75439226,
            4036.47922538,
            4036.1964055,
        ]
    )

    m_co = m_nozz / (1.0 + O_F) / Aux_classes.Default.n
    # m_co = 180000000000000000000000
    # print("sum(A_nozzle)", sum(A_nozzle))
    # A_nozzle = [x * 10000 for x in A_nozzle]
    A_nozzle = Cooling.Nozzle_area_calculation(alpha, y_noz_cool, x_noz_cool)
    print("sum(A_nozzle)", sum(A_nozzle))

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
        T_noz,
        h_noz,
        t_noz,
        np.array(
            [Aux_classes.Default.default_coating_thickness for i in range(len(t_noz))]
        ),
        prop,
        main_material,
        coating,
        0.00005,  # Aux_classes.Default.Dr,
        A_nozzle,
        Aux_classes.Default.T_fuel_tanks,
        m_co,
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
        m_co,
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
x_plot = np.array([x_noz_cool[-1] / len(T_noz) * x for x in range(len(T_noz))])
# print(x_plot)
# x_plot = np.flip(x_plot)
plt.plot(x_plot, T_noz)


plt.plot(x_plot, Tw_wall_nozzle_calculated)
plt.show()
