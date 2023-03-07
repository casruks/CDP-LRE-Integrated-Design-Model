import math
from scipy.integrate import quad
import scipy.optimize
import scipy.constants


#To do: calculate the mass flow required in order to end temperature to be equal to operating temperature!



# Holds all the functions to calculate the Heat tranfer in the rocket engine.

# Goal:
#   Calculate the Temperature at the wall, to check if it melts

# Organization:
#   As classes. Each class holds the specific type of cooling (heatsink, radiation,etc...)
#   Each class holds functions to calculate T at the wall:
# Heatsink: Tcalculation()
# RadiationCool: Tcalculation()
# RegenerativeCool: Run()
#   Other functions are defined in each class, which can be used to calculate the heat transmited
#   Each class holds the T at the wall and the heat transmited as a class variable

# All units are in SI


# ------------------Classes-------------------

# General cooling class
# Can be used to organise/hold all types of cooling for a specific application
class Cooling:
    def __init__(self):
        self.heatsink = Heatsink()
        self.radiationcool = RadiationCool()
        self.regencool = RegenerativeCool()


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
    def Tcalculation(self, T0, Tf, h, A, c, dt, m):
        self.T_calculated = (T0 - Tf) * math.e ** (-h * A / (m * c) * dt) + Tf

    # Heat absorved assuming a certain nozzle mass
    def Qcalculation(self, T0, Tf, h, A, c, dt, m):
        Q_int_arg = (
            lambda time: h
            * (Tf - ((T0 - Tf) * math.e ** (-h * A / (m * c) * time) + Tf))
            * A
        )
        self.Q = quad(Q_int_arg, 0, dt)

    # mass such that the nozzle doesn't melt, assiming Tmelt as the wall temperature
    def mcalculation(self, T0, Tf, h, A, Tmelt, c, dt):
        if Tmelt == Tf:
            raise ValueError(
                "T_melt==Tf, error in mass calculation, for this equality only works after infinite time has passed"
            )
        self.m = h * A * dt / (-math.log((Tmelt - Tf) / (T0 - Tf)) * c)


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

    # Calculate the necessary thickness, assuming that the end temperature inside is Tmelt
    def thickcalculation(self, Tmelt, Tr, eps, k, h):
        self.t = (
            (h * (Tr - Tmelt) / (eps * scipy.constants.sigma) ** (1 / 4) - Tmelt)
            / (Tr - Tmelt)
            * k
            / h
        )
        self.Q = h(Tr - Tmelt)

    # Defines system of equations required in order to find the end temperature. Auxiliary function
    def Tcalculation_system(self, x, Tr, eps, k, t, h):
        Ti, Tout = x
        return [
            h * (Tr - Ti) - (Tout - Ti) * k / t,
            eps * scipy.constants.sigma * Tout**4 - h * (Tr - Ti),
        ]

    # Calculates the end equilibrium temperature
    # Requires temperature at the wall guesses for the iterative method
    def Tcalculation(self, Toutguess, Tinguess, Tr, eps, k, t):
        x0 = [Toutguess, Tinguess]
        sol = scipy.optimize.fsolve(self.Tcalculation_system, x0, args=(Tr, eps, k, t))
        self.T_calculated = sol[0]
        return self.T_calculated


# Regenerative Cooling
# Currently 0D
# Use function Run() to get equilibrium wall temperature, end coolant temperature and pressure drop
# in the coolant
# Run() Also inicialises all the components required
# Q is a power, and is in W
class RegenerativeCool:
    def __init__(self):
        self.Q = 0
        # heat
        self.t = 0
        # thickness
        self.T_calculated = -1

    # Calculates and returns equilibrium wall temperature, and end coolant temperature
    # auxiliary function of Run()
    def Tcalculation(self, Tr, Ti_co, A, hg):

        q = (Tr - Ti_co) / (1 / hg + self.t / self.Mater.k + 1 / self.hco)
        self.Q += q * A
        Tinext_co = Ti_co + q * A / (self.Prop.fcp * self.m_flow_fuel)
        T_wall = self.t / self.Mater.k * q + Ti_co + q / self.hco

        # Tinext_co: end coolant temperature
        # T_wall: wall temperature
        return Tinext_co, T_wall

    def Tcalculation1D(self, Tr, Ti_co, A, hg, ArrayCounter):

        q = (Tr - Ti_co) / (1 / hg + self.t[ArrayCounter] / self.Mater.k + 1 / self.hco)
        self.Q += q * A
        Tinext_co = Ti_co + q * A / (self.Prop.fcp * self.m_flow_fuel)
        T_wall = self.t[ArrayCounter] / self.Mater.k * q + Ti_co + q / self.hco

        # Tinext_co: end coolant temperature
        # T_wall: wall temperature
        return Tinext_co, T_wall

    # Calculates pressure loss
    # auxiliary function of Run()
    def pressureloss(self, m_flow_fuel, Dr, L):
        delta_p = self.f * m_flow_fuel**2 / (2 * self.Prop.f_dens_l) * L / Dr
        return delta_p

    def Inicialise(self, t, Prop, Mater, Dr, Re, m_flow_fuel):
        self.Q = 0
        #self.Pr = 4 * Prop.f_gamma / (9 * Prop.f_gamma - 5)
        self.Pr=0.69 #PLACEHOLDER
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

    # Main function for regenerative cooling
    # Takes material properties, flow properties, and coolant pipe properties
    # Returns end wall temperature, end coolant temperature, pressure drop
    # Saves/Updates Q
    def Run(self, Tr, hg, t, Prop, Mater, Dr, A, Ti_co, Re, m_flow_fuel, L):
        self.Inicialise(t, Prop, Mater, Dr, Re, m_flow_fuel)

        T_co_calcualted, T_wall_calcualted = self.Tcalculation(Tr, Ti_co, A, hg)
        ploss = self.pressureloss(m_flow_fuel, Dr, L)

        # T_co_list=[0 for i in range(len(Tr))]
        # T_co_list[0]=Ti_co
        # for i in Tr:
        #   T_co_calcualted[i], T_wall_calcualted[i] = self.Tcalculation(Tr[i],T_co_list[i],A[i])
        return T_co_calcualted, T_wall_calcualted, ploss

    # Run function, but for 1D case
    def Run1D(self, Tr, hg, t, Prop, Mater, Dr, A, Ti_co, Re, m_flow_fuel, L):
        self.Inicialise(t, Prop, Mater, Dr, Re, m_flow_fuel)

        # T_co_calcualted, T_wall_calcualted = self.Tcalculation(Tr, Ti_co, A, hg)
        ploss = self.pressureloss(m_flow_fuel, Dr, L)

        T_co_calcualted = [0 for i in range(len(Tr) + 1)]
        T_co_calcualted[0] = Ti_co
        T_wall_calcualted = [0 for i in range(len(Tr))]
        for i in range(len(Tr)):
            T_co_calcualted[i + 1], T_wall_calcualted[i] = self.Tcalculation1D(
                Tr[i], T_co_calcualted[i], A, hg[i], i
            )
        #print(T_co_calcualted)
        return T_co_calcualted, T_wall_calcualted, ploss

    def Run_for_Toperating0D(self, Tr, hg, t, Prop, Mater, Dr, A, Ti_co, Re, m_flow_fuel, L):
        self.Inicialise(t, Prop, Mater, Dr, Re, m_flow_fuel)

        q = (Tr - Ti_co) / (1 / hg + self.t / self.Mater.k + 1 / self.hco)
        self.Q += q * A
        Tinext_co = Ti_co + q * A / (self.Prop.fcp * self.m_flow_fuel)
        T_wall = self.t[ArrayCounter] / self.Mater.k * q + Ti_co + q / self.hco

