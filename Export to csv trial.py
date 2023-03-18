import csv
from dataclasses import dataclass

@dataclass
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

    #Combustion
    h_comb = 0.0 #Conductive heat transfer coefficient in chamber
    Dc = 0.00 #Diameter of Combustion Chamber
    ThicknessChamber = 0.0  #Thickness of CC
    Chamber_L = 0.0 # Length of CC
    chambermass = 0.0 #mass of the CC
    Tc = 0.0 # [K] Combustion temperature
    Mnoz = 0.0 # [kg] Total nozzle mass

    #Cooling


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


data = Data()

def export_to_csv(filename, data):
    members = [attr for attr in dir(data) if not callable(getattr(data, attr)) and not attr.startswith("__")]
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(members)
        writer.writerow([getattr(data, attr) for attr in members])

# Export data to CSV
export_to_csv('Data.csv', data)







