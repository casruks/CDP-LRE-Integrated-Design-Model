import math
from scipy.integrate import quad
import scipy.optimize
import scipy.constants
import numpy as np

# import Materials as Mt

# import Materials_2 as Mt
from statistics import mean
import Aux_classes


# To do: calculate the mass flow required in order to end temperature to be equal to operating temperature!


# Holds all the functions to calculate the Heat tranfer in the rocket engine.

# Goal:
#   Calculate the Temperature at the wall, to check if it melts
#   Find the end coolant temperature
#   Find the coolant pressure loss

# Organization:
#   As classes. Each class holds the specific type of cooling (heatsink, radiation,etc...)
#   Each class holds functions to calculate T at the wall:
# Heatsink: Tcalculation()
# RadiationCool: Tcalculation()
# RegenerativeCool: Main_regenerative_run_function()
#   Other functions are defined in each class, which can be used to calculate the heat transmited
#   Each class holds the T at the wall and the heat transmited as a class variable

# All units are in SI


# ------------------Classes-------------------

# General cooling class
# Can be used to organise/hold all types of cooling for a specific application
class CoolingClass:
    def _init_(self):
        self.heatsink = Heatsink()
        self.radiationcool = RadiationCool()
        self.regencool = RegenerativeCool()

    def Run_cooling(
        self,
        T0: float,
        c: float,
        operationtime: float,
        m_casing: float,
        eps: float,
        Tr: np.array,
        hg: float,
        thickness: float,
        coating_thickness: float,
        prop: Aux_classes.Propellant,
        Mater,
        coating,
        Dr: float,
        A: float,
        Ti_co: float,
        m: float,
        L: float,
        y: np.array,
        overwriteA: bool,
        case_run: int,
        thickness_for_regenerative: float,
    ):
        self.heatsink = Heatsink()
        self.radiationcool = RadiationCool()
        self.regencool = RegenerativeCool()
        # Inicialise variables, namely output ones
        err = 0
        warn = 0
        T_co_calculated = Ti_co
        Tw_wall_calculated = np.array([1])
        ploss = 0
        m_flow_fuel = 0
        type_variable = -5
        T_outer_wall = np.array([1])
        # Check if input variables are positive
        if (
            check_positive_args(
                T0,
                c,
                operationtime,
                m_casing,
                eps,
                Tr,
                hg,
                thickness,
                coating_thickness,
                Dr,
                A,
                Ti_co,
                m,
                L,
                y,
            )
            == False
        ):
            err = err | (1 << 0)
            return (
                T_co_calculated,
                Tw_wall_calculated,
                ploss,
                m_flow_fuel,
                0,
                T_outer_wall,
                err,
                warn,
            )

        # Calculate total contact area for cooling system
        if overwriteA == False:
            A_total = sum(
                [self.regencool.FindA(y, L, len(y), i) for i in range(len(y))]
            )
        else:
            A_total = A
            # print("A: ",A)

        A_rad = sum([self.regencool.FindA(y, L, len(y), i) for i in range(len(y))])

        # Determine which operating temperature to use
        if coating.OpTemp_u == 0 or np.any(coating_thickness == 0):
            TestTemp = Mater.OpTemp_u
        else:
            TestTemp = coating.OpTemp_u

        # Calculate the heatsink solution
        err = self.heatsink.Tcalculation(
            T0, np.amax(Tr), np.amax(hg), A_rad, c, operationtime, m_casing, err
        )
        type_variable = 0
        Tw_wall_calculated = [self.heatsink.T_calculated]
        T_outer_wall = Tw_wall_calculated
        self.Q = self.heatsink.Q

        if err != 0:
            return (
                T_co_calculated,
                Tw_wall_calculated,
                ploss,
                m_flow_fuel,
                type_variable,
                T_outer_wall,
                err,
                warn,
            )

        if self.heatsink.T_calculated > TestTemp:

            type_variable = 1
            # Calculate the radiation cooling solution
            if np.any(x <= 0 for x in coating_thickness):
                k_arg = Mater.k
                t_Arg = mean(thickness)
            else:
                k_arg = Mater.k * coating.k
                t_Arg = mean(thickness) * coating.k + mean(coating_thickness) * Mater.k
            err = self.radiationcool.Tcalculation(
                2000,
                100,
                np.amax(Tr),
                eps,
                k_arg,
                t_Arg,
                np.amax(hg),
                err,
            )
            self.Q = self.radiationcool.q * A_rad
            T_co_calculated = Ti_co
            Tw_wall_calculated = [self.radiationcool.T_calculated]
            T_outer_wall = [self.radiationcool.T_outer_wall]
            # print("T_outer_wall radiation:", T_outer_wall)
            # print("Tw_wall_calculated radiation:", Tw_wall_calculated)
            ploss = 0
            m_flow_fuel = 0
            if err != 0:
                return (
                    T_co_calculated,
                    Tw_wall_calculated,
                    ploss,
                    m_flow_fuel,
                    type_variable,
                    T_outer_wall,
                    err,
                    warn,
                )

            if self.radiationcool.T_calculated > TestTemp:
                type_variable = 2
                # Calculate the Regenerative Cooling Solution
                (
                    T_co_calculated,
                    Tw_wall_calculated,
                    ploss,
                    m_flow_fuel,
                    err,
                ) = self.regencool.Main_regenerative_run_function(
                    Tr,
                    hg,
                    thickness_for_regenerative,
                    coating_thickness,
                    prop,
                    Mater,
                    coating,
                    Dr,
                    A_total,
                    Ti_co,
                    m,
                    L,
                    y,
                    overwriteA,
                    case_run,
                    err,
                )
                T_outer_wall = self.regencool.T_out_wall
                # print("Tw_wall_calculated: ",Tw_wall_calculated)
                # Update heat extracted by cooling
                self.Q = self.regencool.Q
        # print("type_variable", type_variable)
        # Check if output variables are within reason
        if check_positive_args(T_co_calculated) == False:
            err = err | (1 << 7)

        if check_positive_args(Tw_wall_calculated) == False or np.any(
            [x > TestTemp for x in Tw_wall_calculated]
            or np.any([x > TestTemp for x in T_outer_wall])
        ):
            # print("Tw_wall_calculated: ",Tw_wall_calculated)
            # err = err | (1 << 8)
            warn = warn | (1 << 0)
        if check_positive_args(ploss) == False or ploss > 10**10:
            err = err | (1 << 9)
        if check_positive_args(m_flow_fuel) == False:
            err = err | (1 << 10)

        if check_positive_args(self.Q) == False:
            err = err | (1 << 11)
        if check_positive_args(T_outer_wall) == False:
            err = err | (1 << 12)

        # if(m_flow_fuel > 6000 or Tw_wall_calculated[-1] > 2000 or T_co_calculated > 1000 or ploss > 10**5):
        # warn=warn|(1<<1)
        # if T_co_calculated > Tw_wall_calculated[-1] and type_variable == 2:
        #    err = err | (1 << 6)
        return (
            T_co_calculated,
            Tw_wall_calculated,
            ploss,
            m_flow_fuel,
            type_variable,
            T_outer_wall,
            err,
            warn,
        )


# Heat sink model
# Calculates the Temperature at the wall after a time interval of operation time
# Q is the total energy transfered, and is in J
class Heatsink:
    def _init_(self, Q=-1, m=-1):
        self.Q = Q
        # heat
        self.m = m
        # nozzle mass required
        self.T_calculated = -1

    # Final T after a certain operation time dt, assuming a nozzle mass
    def Tcalculation(self, T0, Tf, h, A, c, dt, m, err):
        # if(check_positive_args(T0, Tf, h, A, c, dt, m)==False):
        # print(A)
        self.T_calculated = (T0 - Tf) * float(math.e ** (-h * A / (m * c) * dt)) + Tf
        self.Q = (self.T_calculated - T0) * c / dt  # Watts
        if check_positive_args(self.T_calculated) == False:
            err = err | (1 << 1)
        return err


# Ratiation cool
# Calculates the equilibrium temperature, taking into account radiated heat
# Q is a power, and is in W
class RadiationCool:
    def _init_(self, Q=-1, t=-1):
        self.Q = Q
        # heat
        self.t = t
        # thickness
        self.T_calculated = -1

    # Defines system of equations required in order to find the end temperature. Auxiliary function
    def Tcalculation_system(self, x, Tr, eps, k, t, h):
        Tout, Ti = x
        return [
            h * (Tr - Ti) - (Ti - Tout) * k / t,
            eps * scipy.constants.sigma * Tout**4 - h * (Tr - Ti),
        ]

    # Calculates the end equilibrium temperature
    # Requires temperature at the wall guesses for the iterative method
    def Tcalculation(self, Toutguess, Tinguess, Tr, eps, k, t, h, err):
        # if check_positive_args(Tr, eps, k, t) == False:

        x0 = [Toutguess, Tinguess]
        sol = scipy.optimize.least_squares(
            self.Tcalculation_system, x0, args=(Tr, eps, k, t, h)
        )
        if not sol["success"] or abs(sum(sol["fun"])) > 0.01:
            self.T_calculated = 10**10
            self.T_outer_wall = 10**10
            # err = err | (1 << 13)
        else:
            self.T_calculated = sol["x"][1]
            self.T_outer_wall = sol["x"][0]
        # print("h * (Tr - Ti) - (Ti - Tout) * k / t: ",h * (Tr - sol[1]) - (sol[1] - sol[0]) * k / t)

        # print("self.T_calculated",sol[1])
        # print("self.T_outer_wall",sol[0])
        self.q = (self.T_calculated - self.T_outer_wall) * k / t
        if check_positive_args(self.T_calculated) == False:
            err = err | (1 << 2)
        return err


# Regenerative Cooling
#
# Main Goal:Get equilibrium wall temperature, end coolant temperature and pressure drop
#
# Options:
# Run(): gives the end temperature of the coolant, the wall temperature and the pressure loss, for a 0D case
# Run1D(): gives the end temperature of the coolant, the wall temperature and the pressure loss, for a 1D case
# Run_for_Toperating0D(): sizes the hydralic diamter to get the wall cooled to operating temperature. Returns the coolant temperature and the pressure loss
# Run_for_Toperating1D(): same as #Run_for_Toperating0D(), but in 1D
# Run1D_iterative_for_m(): Find the massflow required to get the wall cooled to operating temperature. Returns the massflow, the coolant temperature, the wall temperature, and the pressure loss
#
# Note on units
# Q is a power, and is in W
class RegenerativeCool:
    def _init_(self):
        self.Q = 0
        # heat
        self.t = 0
        # thickness
        self.T_calculated = -1

    # ----------------------------------------------------------------------------------------------------------------

    # Calculates and returns equilibrium wall temperature, and end coolant temperature

    def Tcalculation1D(
        self, Taw: float, Ti_co: float, A: float, hg: float, ArrayCounter: int
    ):
        q = (Taw - Ti_co) / (
            1 / (hg) + self.t[ArrayCounter] / self.Mater.k + 1 / self.hco
        )

        # print("self.t[ArrayCounter] / self.Mater.k: ",self.t[ArrayCounter] / self.Mater.k)
        # print("self.t[ArrayCounter] / self.Mater.k: ",self.t[ArrayCounter] / self.Mater.k)
        # print("1 / hg: ",1 / hg)
        # print("1 / self.hco: ", 1 / self.hco)
        # print()
        # print("I AM HERE")
        # q = 10**8
        # print("A: ",sum(A))
        self.Q = self.Q + q * A[ArrayCounter]

        # self.Q += q * A
        # print(A[ArrayCounter])
        # self.Prop.fcp=283.73e3
        Tinext_co = Ti_co + q * A[ArrayCounter] / (self.Prop.fcp * self.m_flow_fuel)
        # Tinext_co = Ti_co + q * A/ (self.Prop.fcp * self.m_flow_fuel)
        # print("self.hco", self.hco)
        # T_wall = self.t[ArrayCounter] / self.Mater.k * q + Ti_co + q / self.hco
        T_wall = Taw - q / hg
        # print("T_wall : ", T_wall)
        # print("hg: ", hg)
        self.T_outer_Wall_loop_val = T_wall - self.t[ArrayCounter] / self.Mater.k * q
        # print("self.t[ArrayCounter]",self.t[ArrayCounter])
        # print("T_outer by wall: ", T_wall - self.t[ArrayCounter] / self.Mater.k * q)
        # print("T_out by coolant: ", Ti_co - q / self.hco)
        # Tinext_co: end coolant temperature
        # T_wall: wall temperature
        check_positive_args(Tinext_co, T_wall)
        return Tinext_co, T_wall

    # Calculates pressure loss
    # auxiliary function of Run()
    def auxiliary_calculation_f_for_pressureloss(self, local_f, Dr, local_Re):
        roughness = 3.5  # Default from a paper
        return (
            local_f
            - (
                1
                / (
                    -2
                    * math.log10(
                        roughness / Dr / 3.7 + 2.51 / (local_Re * math.sqrt(local_f))
                    )
                )
            )
            ** 2
        )

    def pressureloss(self, m_flow_fuel: float, Dr: float, L: float):
        local_Re = Dr * m_flow_fuel / self.Prop.fmiu
        local_f = scipy.optimize.fsolve(
            self.auxiliary_calculation_f_for_pressureloss, 1, args=(Dr, local_Re)
        )
        delta_p = local_f * m_flow_fuel**2 / (2 * self.Prop.f_dens_g) * L / Dr
        return delta_p

    def Inicialise(
        self,
        t: float,
        Prop: Aux_classes.Propellant,
        Mater,
        Dr: float,
        Re: float,
        m_flow_fuel: float,
    ):
        self.Q = 0
        # self.Pr = 4 * Prop.f_gamma / (9 * Prop.f_gamma - 5)
        self.Pr = 1  # PLACEHOLDER
        self.f = (1.82 * math.log10(Re) - 1.64) ** (-2)
        # print("Re: ", Re)
        self.Nu = (
            self.f
            / 8
            * (Re - 1000)
            * self.Pr
            / (1 + 12.7 * math.sqrt(self.f / 8) * (self.Pr ** (2 / 3) - 1))
        )
        self.hco = self.Nu * Mater.k / Dr
        # print("self.hco: ", self.hco)

        self.Prop = Prop
        self.m_flow_fuel = m_flow_fuel

        self.t = t
        self.Mater = Mater

    # One of the Main functions for regenerative cooling
    # Takes material properties, flow properties, and coolant pipe properties
    # Returns end wall temperature, end coolant temperature, pressure drop
    # Saves/Updates Q

    # Run function, but for 1D case
    def Run1D(
        self,
        Tr: np.array,
        hg: float,
        t: float,
        t_coating: float,
        Prop: Aux_classes.Propellant,
        Mater,
        Mater_coating,
        Dr,
        A,
        Ti_co,
        m_flow_fuel,
        L,
        err,
    ):
        # if check_positive_args(Tr, hg, t, Dr, A, Ti_co, m_flow_fuel, L) == False:

        # self.t = t * Mater_coating.k + t_coating * Mater.k
        self.Mater = Mater
        if np.any(t_coating <= 0):
            self.t = t
            self.Mater.k = Mater.k
        else:
            self.Mater.k = Mater.k * Mater_coating.k
            self.t = t * Mater_coating.k + t_coating * Mater.k

        Re = m_flow_fuel / Prop.fmiu * 4 / (math.pi * Dr)
        # Inicialise variables
        self.Inicialise(self.t, Prop, self.Mater, Dr, Re, m_flow_fuel)
        self.D = Dr

        # T_co_calcualted, T_wall_calcualted = self.Tcalculation(Tr, Ti_co, A, hg)

        # Calculate pressure loss
        ploss = self.pressureloss(m_flow_fuel, Dr, L)

        # Inicialise variables with placeholders
        T_co_calcualted = [0 for i in range(len(Tr) + 1)]
        T_co_calcualted[0] = Ti_co
        T_wall_calcualted = [0 for i in range(len(Tr))]
        self.T_out_wall = [0 for i in range(len(Tr))]
        self.Q = 0
        # Calculate the wall temperature and the coolant temperature for each point along the wall
        for i in range(len(Tr)):
            T_co_calcualted[i + 1], T_wall_calcualted[i] = self.Tcalculation1D(
                Tr[i], T_co_calcualted[i], A, hg[i], i
            )
            # print(T_co_calcualted[i])
            self.T_out_wall[i] = self.T_outer_Wall_loop_val
            # self.Q=self.Q+(Tr[i]-T_wall_calcualted[i])*hg[i]

        # print(T_co_calcualted)
        if check_positive_args(T_co_calcualted, T_wall_calcualted, ploss) == False:
            err = err | (1 << 3)
        return T_co_calcualted[-1], T_wall_calcualted, ploss, err

    # ----------------------------------------------------------------------------------------------------------------

    def FindA(self, y, L, l, i):
        return (2 * math.pi * y[i]) * L / l

    def OverrideA(self, A):
        self.A = A

    def Main_regenerative_run_function(
        self,
        Tr: np.array,
        hg: float,
        t: float,
        t_coating: float,
        Prop: Aux_classes.Propellant,
        Mater,
        Mater_coating,
        Dr: float,
        A: float,
        Ti_co: float,
        m_flow_fuel: float,
        L: float,
        y: float,
        overwriteA: bool,
        case_run: int,
        err,
    ):
        if overwriteA == False:
            A_arg = [self.FindA(y, L, len(y), i) for i in range(len(y))]
        else:
            # A_arg = [A / len(y) for i in range(len(y))]
            A_arg = A

        match case_run:
            case 0:
                T_co_calculated, Tw_wall_calculated, ploss, err = self.Run1D(
                    Tr,
                    hg,
                    t,
                    t_coating,
                    Prop,
                    Mater,
                    Mater_coating,
                    Dr,
                    A_arg,
                    Ti_co,
                    m_flow_fuel,
                    L,
                    err,
                )
            case _:
                raise ValueError("Non existent regenerative cooling case called")

        return T_co_calculated, Tw_wall_calculated, ploss, m_flow_fuel, err


# Sanitisaion function


def check_positive_args(*args):
    for arg in args:
        if isinstance(arg, (int, float, np.int32, np.generic)):
            if arg < 0:
                # print("All numerical arguments must be positive")
                return False
        elif isinstance(arg, (list, tuple, np.ndarray)):
            if any(x < 0 for x in arg):
                # print("All elements of numerical arguments must be positive")
                return False
        else:
            raise ValueError("Unsupported argument type: {type(arg)._name_}")
    return True


def thermal_stress(
    main_material,
    coating_material,
    main_material_thickness,
    coating_thickness,
    y,
    delta_T,
):
    if np.all(coating_thickness != 0):

        Ri = y + coating_thickness / 2
        Ro = y + coating_thickness + main_material_thickness / 2

        delta_Ri_indp = 2 * math.pi * delta_T * coating_material.thr_exp
        p_contact = delta_Ri_indp / (
            Ri**2 / (coating_material.Emod * coating_thickness)
            + Ro**2 / (main_material.Emod * main_material_thickness)
        )

        stress_main = np.amax(p_contact * Ro / main_material_thickness)
        stress_coating = np.amax(p_contact * Ri / coating_thickness)

        return stress_main, stress_coating
    else:
        v = main_material.mu
        stress_main = (
            delta_T * main_material.thr_exp * main_material.Emod / (2 * (1 - v))
        )
        # print("stress main: ", stress_main)

        stress_coating = 0
        return stress_main, stress_coating


def outputs(
    T_outer_wall_chamber,
    T_outer_wall_nozzle,
    Tw_wall_chamber_calculated,
    Tw_wall_nozzle_calculated,
    nozzle_material,
    chamber_material,
    coating_material,
    chamber_thickness,
    nozzle_thickness,
    coating_thickness,
    y,
    chamber_Radius,
    err,
):
    # pLACE_HOLDER_THERMAL_EXPANSION_COEFFICIENT = 6.12 * 10**-6
    T_outer_wall_chamber = np.array(T_outer_wall_chamber)
    T_outer_wall_nozzle = np.array(T_outer_wall_nozzle)
    Tw_wall_chamber_calculated = np.array(Tw_wall_chamber_calculated)
    Tw_wall_nozzle_calculated = np.array(Tw_wall_nozzle_calculated)

    max_temperature_outer = np.maximum(
        np.max(T_outer_wall_chamber), np.max(T_outer_wall_nozzle)
    )
    max_temperature_inner = np.maximum(
        np.max(Tw_wall_chamber_calculated), np.max(Tw_wall_nozzle_calculated)
    )

    Tw_wall_chamber_calculated = np.array(Tw_wall_chamber_calculated)
    Tw_wall_nozzle_calculated = np.array(Tw_wall_nozzle_calculated)
    delta_T = np.maximum(
        np.max((Tw_wall_chamber_calculated + T_outer_wall_chamber) / 2)
        - Aux_classes.Default.T0,
        np.max((Tw_wall_nozzle_calculated + T_outer_wall_nozzle) / 2)
        - Aux_classes.Default.T0,
    )
    # print("DeltaT: ", delta_T)
    if check_positive_args(delta_T) == False:
        err = err | (1 << 12)

    if delta_T == np.max(
        ((Tw_wall_chamber_calculated + T_outer_wall_chamber) / 2)
        - Aux_classes.Default.T0
    ):
        main_material = chamber_material
        main_material_thickness = chamber_thickness
        Radius = chamber_Radius
    else:
        main_material = nozzle_material
        main_material_thickness = nozzle_thickness
        Radius = y[
            np.argmax(
                ((Tw_wall_chamber_calculated + T_outer_wall_chamber) / 2)
                - Aux_classes.Default.T0
            )
        ]
    stress_main, stress_coating = thermal_stress(
        main_material,
        coating_material,
        main_material_thickness,
        coating_thickness,
        Radius,
        delta_T,
    )

    safety_factor_cooling_main = main_material.yieldstress_l / stress_main
    if stress_coating != 0:
        safety_factor_cooling_coating = coating_material.yieldstress_l / stress_coating
    else:
        safety_factor_cooling_coating = (
            0  # so it doesn't appear as the max value in the return
        )

    return (
        max(stress_main, stress_coating),
        max(safety_factor_cooling_main, safety_factor_cooling_coating),
        max_temperature_inner,
        max_temperature_outer,
        err,
    )


def Nozzle_area_calculation(alpha, y_noz_cool, x_noz_cool):
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
    return A_nozzle
