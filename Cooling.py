import math
from scipy.integrate import quad
import scipy.optimize
import scipy.constants 
class Cooling:
    def __init__(self):
        self.heatsink = Heatsink();
        self.radiationcool = RadiationCool();


class Heatsink():
    def __init__(self, Q=-1, ,m=-1):
        self.Q=Q; #heat
        self.m=m; #nozzle mass required
        self.T_calculated=-1;

    #Final T after a certain operation time dt, assuming a nozzle mass
    def Tcalculation(self, T0, Tf,h,A,m=self.m,c,dt):
        self.T_calculated = (T0-Tf)*math.e**(-*h*A/(m*c)*dt)+Tf

    #Heat absorved assuming a certain nozzle mass
    def Qcalculation(self,T0, Tf,h,A,m=self.m,c,dt)
       Q_int_arg =  lambda time: h*(Tf-((T0-Tf)*math.e**(-*h*A/(m*c)*time)+Tf))*A
       self.Q=quad(Q_int_arg,0,dt)

    #mass such that the nozzle doesn't melt
    def mcalculation(self,T0, Tf,h,A,Tmelt,c,dt)):
        if T_melt==Tf:
            raise ValueError('T_melt==Tf, error in mass calculation, for this equality only works after infinite time has passed')
        self.m=h*A*dt/(-math.log((Tmelt-Tf)/(T0-Tf))*c)

class RadiationCool:
    def __init__(self, Q=-1, t=-1):
        self.Q=Q; #heat
        self.t=t; #thickness
        self.T_calculated=-1;

    #Calculate the necessary thickness, assuming that the end temperature inside is Tmelt
    def thickcalculation(self,Tmelt,Tr,eps,k,h):
        self.t=(h*(Tr-Tmelt)/(eps*scipy.constants.sigma)**(1/4)-Tmelt)/(Tr-Tmelt)*k/h
        self.Q=h(Tr-Tmelt)
    
    def Tcalculation_system(x,Tr, eps,k,t=self.t):
        Ti, Tout=x
        return [ h*(Tr-Ti)-(Tout-Ti)*k/t,
            eps*scipy.constants.sigma*Tout**4- h*(Tr-Ti)
        ]


    def Tcalculation(self,Toutguess,Tinguess,Tr, eps,k,t=self.t):
        x0=[Toutguess,Tinguess]
        sol=scipy.optimize.fsolve(Tcalculation_system,x0,args=(Tr, eps,k,t))
        self.T_calculated=sol[0]

class RegenerativeCool:
    def __init__(self, Q=-1,t=-1):
        self.Q=Q; #heat
        self.t=t; #thickness
        self.T_calculated=-1;
    
    #minimum thickness so that it does not melt
    def thicknesscalculation():




