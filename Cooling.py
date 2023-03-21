import math
from scipy.integrate import quad
import scipy.optimize
import scipy.constants
import numpy as np
import Materials as Mt
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
    def __init__(self):
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
        Mater: Mt.Materials,
        coating: Mt.Materials,
        Dr: float,
        A: float,
        Ti_co: float,
        m: float,
        L: float,
        y: np.array,
        overwriteA: bool,
        case_run: int,
    ):

        # Inicialise variables, namely output ones
        err = 0
        warn = 0
        T_co_calculated = Ti_co
        Tw_wall_calculated = [-1]
        ploss = 0
        m_flow_fuel = 0
        type_variable = 0
        T_outer_wall=0
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
                overwriteA,
                case_run,
            )
            == False
        ):
            err = err | (1 << 0)
            return T_co_calculated, Tw_wall_calculated, ploss, m_flow_fuel, T_outer_wall,err, warn

        # Calculate total contact area for cooling system
        if overwriteA == False:
            A_total = sum(
                [self.regencool.FindA(y, L, len(y), i) for i in range(len(y))]
            )
        else:
            A_total = A

        # Determine which operating temperature to use
        if coating.OpTemp_u == 0:
            TestTemp = Mater.OpTemp_u
        else:
            TestTemp = coating.OpTemp_u

        # Calculate the heatsink solution
        err = self.heatsink.Tcalculation(
            T0, mean(Tr), mean(hg), A_total, c, operationtime, m_casing, err
        )

        Tw_wall_calculated = [self.heatsink.T_calculated]
        T_outer_wall=[Tw_wall_calculated]
        self.Q = self.heatsink.Q

        if err != 0:
            return T_co_calculated, Tw_wall_calculated, ploss, m_flow_fuel, T_outer_wall,err, warn

        if self.heatsink.T_calculated > TestTemp:

            type_variable = 1
            # Calculate the radiation cooling solution
            err = self.radiationcool.Tcalculation(
                700,
                500,
                mean(Tr),
                eps,
                Mater.k * coating.k,
                mean(thickness) * coating.k + mean(coating_thickness) * Mater.k,
                mean(hg),
                err,
            )
            self.Q = self.radiationcool.Q
            T_co_calculated = Ti_co
            Tw_wall_calculated = [self.radiationcool.T_calculated]
            T_outer_wall=[self.radiationcool.T_outer_wall]
            ploss = 0
            m_flow_fuel = 0
            if err != 0:
                return (
                    T_co_calculated,
                    Tw_wall_calculated,
                    ploss,
                    m_flow_fuel,
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
                    thickness,
                    coating_thickness,
                    prop,
                    Mater,
                    coating,
                    Dr,
                    A,
                    Ti_co,
                    m,
                    L,
                    y,
                    overwriteA,
                    case_run,
                    err,
                )
                T_outer_wall=self.regencool.T_out_wall
        # Update heat extracted by cooling
        self.Q = self.regencool.Q

        # Check if output variables are within reason
        if check_positive_args(T_co_calculated) == False or T_co_calculated > 1000:
            err = err | (1 << 7)

        if check_positive_args(Tw_wall_calculated) == False or any(
            x > TestTemp for x in Tw_wall_calculated
        ):
            err = err | (1 << 8)
        if check_positive_args(ploss) == False or ploss > 10**6:
            err = err | (1 << 9)
        if check_positive_args(m_flow_fuel) == False or m_flow_fuel > 30:
            err = err | (1 << 10)

        if check_positive_args(self.Q) == False:
            err = err | (1 << 11)

        if check_positive_args(T_outer_wall) == False:
            err = err | (1 << 12)

        # if(m_flow_fuel > 6000 or Tw_wall_calculated[-1] > 2000 or T_co_calculated > 1000 or ploss > 10**5):
        # warn=warn|(1<<1)
        if T_co_calculated > Tw_wall_calculated[-1]:
            err = err | (1 << 6)
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
    def __init__(self, Q=-1, m=-1):
        self.Q = Q
        # heat
        self.m = m
        # nozzle mass required
        self.T_calculated = -1

    # Final T after a certain operation time dt, assuming a nozzle mass
    def Tcalculation(self, T0, Tf, h, A, c, dt, m, err):
        # if(check_positive_args(T0, Tf, h, A, c, dt, m)==False):

        self.T_calculated = (T0 - Tf) * math.e ** (-h * A / (m * c) * dt) + Tf
        if check_positive_args(self.T_calculated) == False:
            err = err | (1 << 1)
        return err


# Ratiation cool
# Calculates the equilibrium temperature, taking into account radiated heat
# Q is a power, and is in W
class RadiationCool:
    def __init__(self, Q=-1, t=-1):
        self.Q = Q
        # heat
        self.t = t
        # thickness
        self.T_calculated = -1


    # Defines system of equations required in order to find the end temperature. Auxiliary function
    def Tcalculation_system(self, x, Tr, eps, k, t, h):
        Ti, Tout = x
        return [
            h * (Tr - Ti) - (Ti - Tout) * k / t,
            eps * scipy.constants.sigma * Tout**4 - h * (Tr - Ti),
        ]

    # Calculates the end equilibrium temperature
    # Requires temperature at the wall guesses for the iterative method
    def Tcalculation(self, Toutguess, Tinguess, Tr, eps, k, t, h, err):
        # if check_positive_args(Tr, eps, k, t) == False:

        x0 = [Toutguess, Tinguess]
        sol = scipy.optimize.fsolve(
            self.Tcalculation_system, x0, args=(Tr, eps, k, t, h)
        )
        self.T_calculated = sol[1]
        self.T_outer_wall= sol[0]

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
    def __init__(self):
        self.Q = 0
        # heat
        self.t = 0
        # thickness
        self.T_calculated = -1

    # ----------------------------------------------------------------------------------------------------------------

    # Calculates and returns equilibrium wall temperature, and end coolant temperature
    

    def Tcalculation1D(
        self, Tr: float, Ti_co: float, A: float, hg: float, ArrayCounter: int
    ):
        q = (Tr - Ti_co) / (1 / hg + self.t[ArrayCounter] / self.Mater.k + 1 / self.hco)
        self.Q += q * A[ArrayCounter]
        # self.Q += q * A
        # print(A[ArrayCounter])
        Tinext_co = Ti_co + q * A[ArrayCounter] / (self.Prop.fcp * self.m_flow_fuel)
        # Tinext_co = Ti_co + q * A/ (self.Prop.fcp * self.m_flow_fuel)
        # print("self.hco", self.hco)
        # T_wall = self.t[ArrayCounter] / self.Mater.k * q + Ti_co + q / self.hco
        T_wall = Tr - q / hg
        self.T_outer_Wall_loop_val=T_wall-self.t[ArrayCounter]/self.Mater.k*q
        # Tinext_co: end coolant temperature
        # T_wall: wall temperature
        check_positive_args(Tinext_co, T_wall)
        return Tinext_co, T_wall

    # Calculates pressure loss
    # auxiliary function of Run()
    def auxiliary_calculation_f_for_pressureloss(self,local_f,Dr,local_Re):
        roughness=3.5 #Default from a paper
        return local_f-(1/(-2*math.log10(roughness/Dr/3.7+2.51/(local_Re*local_f))))**2
    def pressureloss(self, m_flow_fuel: float, Dr: float, L: float):
        local_Re=Dr*m_flow_fuel/self.Prop.fmiu
        local_f=scipy.optimize.fsolve(self.auxiliary_calculation_f_for_pressureloss,1,args=(Dr,local_Re))
        delta_p = local_f * m_flow_fuel**2 / (2 * self.Prop.f_dens_l) * L / Dr
        return delta_p*10**6

    def Inicialise(
        self,
        t: float,
        Prop: Aux_classes.Propellant,
        Mater: Mt.Materials,
        Dr: float,
        Re: float,
        m_flow_fuel: float,
    ):
        self.Q = 0
        # self.Pr = 4 * Prop.f_gamma / (9 * Prop.f_gamma - 5)
        self.Pr = 0.69  # PLACEHOLDER
        self.f = (1.82 * math.log10(Re) - 1.64) ** (-2)
        self.Nu = (
            self.f
            / 8
            * (Re - 1000)
            * self.Pr
            / (1 + 12.7 * math.sqrt(self.f / 8) * (self.Pr**2 / 3 - 1))
        )
        self.hco = self.Nu * Mater.k / Dr

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
        Mater: Mt.Materials,
        Mater_coating: Mt.Materials,
        Dr,
        A,
        Ti_co,
        m_flow_fuel,
        L,
        err,
    ):
        # if check_positive_args(Tr, hg, t, Dr, A, Ti_co, m_flow_fuel, L) == False:

        self.t = t * Mater_coating.k + t_coating * Mater.k
        self.Mater = Mater
        self.Mater.k = Mater.k * Mater_coating.k
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
        # Calculate the wall temperature and the coolant temperature for each point along the wall
        for i in range(len(Tr)):
            T_co_calcualted[i + 1], T_wall_calcualted[i] = self.Tcalculation1D(
                Tr[i], T_co_calcualted[i], A, hg[i], i
            )
            self.T_out_wall[i]=self.T_outer_Wall_loop_val
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
        Mater: Mt.Materials,
        Mater_coating: Mt.Materials,
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
            A = [self.FindA(y, L, len(y), i) for i in range(len(y))]
        else:
            A = [A for i in range(len(y))]

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
                    A,
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
                print("All numerical arguments must be positive")
                return False
        elif isinstance(arg, (list, tuple, np.ndarray)):
            if any(x < 0 for x in arg):
                print("All elements of numerical arguments must be positive")
                return False
        else:
            raise ValueError("Unsupported argument type: {type(arg).__name__}")
    return True
