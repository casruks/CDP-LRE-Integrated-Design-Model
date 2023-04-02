import csv
import numpy as np
from matplotlib import pyplot as plt

def Reliability(default, t, Fnom, Fop, val):
    '''
    Determines the reliability of the LRE, based on LH2/LOX. \n
    default: loads default values \n
    t: Thrust time [s] \n
    Fnom: Nominal designed thrust level [N] \n
    val: True or False, for validation.
    '''

    cycles = default.cycles
    cycle = default.cycle_type # 0:EX (expander) - 1:CB (coolant bleed) - 2:GG (gas generator) - 3:SC (staged combustion) - 4:EL (electrical) - 5:PF (pressure fed)
    #cycle = default.cycle
    prop = default.Prop[0]
    N = default.N
    # Data for Cycle Impact (effect of engine cycle on reliability)
    CyclesData = {}
    with open('Reliability_Data/Cycle_Data.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for i, row in enumerate(reader):
            if i != 0:
                # Cycle Code = (lambda, Cf, explanation of cycle type)
                # [Cf is fraction of catastrophic failures, as explained in Fernández et al. (2022)]
                CyclesData[row[0]] = (float(row[1]), float(row[2]), row[4])
    
    # Data for Thrust Impact (effect of total thrust on reliability)
    delta = default.delta
    Fref = default.Fref #SSME
    
    # Data for De-rating/Up-rating Impact (effect of operating at non-nominal thrust on reliability)
    RatingData = {}
    with open('Reliability_Data/Propellant_Uprating_Data.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for i, row in enumerate(reader):
            if i != 0:
                # Propellant Code = (p, q, explanation of propellant)
                # [p and q are parameters as used in Fernández et al. (2022)]
                RatingData[row[0]] = (float(row[1]), float(row[2]), row[4])
    
    # Calculates a first order estimate for LRE reliability based on basic engine parameters
    def reliability_general(t, cycle, Fnom, Fop, N, prop): 
        """t = total burn time in seconds \n
        cycle = engine cycle type code \n
        Fnom = total nominal thrust \n
        Fop = total operating thrust \n
        N = number of engines \n
        prop = propellant code \n"""

        wr = 0
        ## Warning codes
        #01 - Engine cycle not included in model, setting to SC cycle reliability.
        if cycle == 1 or cycle >= 4:
            wr = wr|(1<<0)
            print('Warning,', cycle, 'is not modelled for reliability. returning staged combustion reliability.')
        else:
            wr = wr&(~(1<<0))

        if cycle == 0: #"Expander Cycle"
            cycle = 'SP_EX'
        elif cycle == 2: #"Gas Generator Cycle"
            cycle = 'S_FR_GG'
        elif cycle == 3: #"Staged Combustion Cycle"
            cycle = 'D_FR_SC'
        else: 
            cycle = 'D_FR_SC'
        
        # Cycle Impact
        lambda_ref = CyclesData[cycle][0]
        # Nominal Thrust Impact
        Fnom_single = Fnom/N
        lambda_1 = lambda_ref*(Fnom_single/Fref)**delta
        # De-rating/Up-rating Impact
        Fop_single = Fop/N
        alpha = Fop_single/Fnom_single
        if not prop in RatingData.keys():
            raise Exception("Submitted Propellant Type is not valid. Valid Propellant Codes are: " + str(RatingData.keys()))
        p, q, _ = RatingData[prop]
        lambda_2 = lambda_1*((1 - p) + p*np.exp(-q*(1 - alpha)))
        R = np.exp(-N*lambda_2*t)
        return R, wr  

    # VALIDATION --> Yields the same figure as Fig.3 in Fernández et al. (2022)
    def validate(Fref):
        N = [4, 5, 6, 7, 8, 9]
        t = np.arange(0, 500, 0.1)
    
        plt.figure()
        ax = plt.axes()
        for Ni in N:
            R = reliability_general(t, "D_FR_SC", Fref, Fref, Ni, "LOX_LH2")
            R = R[0]
            ax.plot(t, R, label = "N = " + str(Ni))
    
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Reliability')
        ax.legend()
        plt.show()
        
    if val == True:
        validate(Fref)
    return reliability_general(t, cycle, Fnom, Fop, N, prop)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    