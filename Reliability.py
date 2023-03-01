import csv
import numpy as np
from matplotlib import pyplot as plt

def Reliability(t, cycle, Fnom, Fop, N, prop, val):
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
    delta = 0.1017
    Fref = 2278*1000
    
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
        
        if not cycle in CyclesData.keys():
            raise Exception("Submitted Cycle Type is not valid. Valid Cycle Type Codes are: " + str(CyclesData.keys()))
        
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
        
        # Number of Engines and Total Burn Time Impact
        return np.exp(-N*lambda_2*t)
    
    if val == 1:
    # VALIDATION --> Yields the same figure as Fig.3 in Fernández et al. (2022)
        def validate(Fref):
            N = [4, 5, 6, 7, 8, 9]
            t = np.arange(0, 500, 0.1)
        
            plt.figure()
            ax = plt.axes()
            for Ni in N:
                R = reliability_general(t, "D_FR_SC", Fref, Fref, Ni, "LOX_LH2")
                ax.plot(t, R, label = "N = " + str(Ni))
        
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Reliability')
            ax.legend()
            plt.show()
        return validate(Fref)
    return reliability_general(t, cycle, Fnom, Fop, N, prop)

