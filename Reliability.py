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
        return R    

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
        
#########################################################################################################
from datetime import date
import csv 
import Aux_classes as aux
default = aux.Default(0)
def Validation1(default):
    # cycles = default.cycles
    # lst = ['R']

    # with open('Verification_rel_'+str(date.today())+'.csv', mode='w', newline='') as file:
    #     writer = csv.writer(file)
    #     #writer.writerow(inj1lst)
    #     writer.writerows(map(lambda x: [x], lst))

    # lst_t = np.arange(1, 1e3, 150)
    # lst_F = np.arange(1e3, 1e10, 1e7)
    # val = False
    # for cycle in cycles:
    #     for t in lst_t:
    #         for F in lst_F:
    #             lst.append(Reliability(default, t, F, F, val))

    # with open('Verification_rel_'+str(date.today())+'.csv', mode='w', newline='') as file:
    #     writer = csv.writer(file)
    #     #writer.writerow(inj1lst)
    #     writer.writerows(map(lambda x: [x], lst))
    # print('done')

    #### singular SSME val - 'D_FR_SC'
    t_SSME = 510
    F_SSME = 2278e3
    
    t_RS68 = 411
    F_RS68 = 3300e3

    val = False
    
    ## SSME
    default.cycle_type = 3
    R_DBFRSC = Reliability(default, t_SSME, F_SSME, F_SSME, val)
    print('R_SSME =', R_DBFRSC)
    R_SSME = 0.9983     # [B. K. Wood (Boeing & Rocketdyne), 2002]
    R_SSME = 0.997
    print('error =', (abs(R_DBFRSC-R_SSME)/R_SSME)*100,'%')

    ## RS68 - !PREDICTED!
    default.cycle_type = 2
    R_RS68 = Reliability(default, t_RS68, F_RS68, F_RS68, val)
    print('R_RS68 =', R_RS68)
    R_ref = 0.9987     # [B. K. Wood (Boeing & Rocketdyne), 2002]
    print('error =', (abs(R_RS68-R_ref)/R_ref)*100,'%')
#Validation1(default)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     