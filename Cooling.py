import math
from scipy.integrate import quad
import scipy.optimize
import scipy.constants 





class Cooling:
    def __init__(self):
        self.heatsink = Heatsink();
        self.radiationcool = RadiationCool();


class Heatsink():
    def __init__(self, Q=-1, m=-1):
        self.Q=Q; #heat
        self.m=m; #nozzle mass required
        self.T_calculated=-1;

    #Final T after a certain operation time dt, assuming a nozzle mass
    def Tcalculation(self, T0, Tf,h,A,c,dt,m=self.m):
        self.T_calculated = (T0-Tf)*math.e**(-h*A/(m*c)*dt)+Tf

    #Heat absorved assuming a certain nozzle mass
    def Qcalculation(self,T0, Tf,h,A,c,dt,m=self.m):
       Q_int_arg =  lambda time: h*(Tf-((T0-Tf)*math.e**(-h*A/(m*c)*time)+Tf))*A
       self.Q=quad(Q_int_arg,0,dt)

    #mass such that the nozzle doesn't melt
    def mcalculation(self,T0, Tf,h,A,Tmelt,c,dt):
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
    
    def Tcalculation_system(self,x,Tr, eps,k,t=self.t):
        Ti, Tout=x
        return [ h*(Tr-Ti)-(Tout-Ti)*k/t,
            eps*scipy.constants.sigma*Tout**4- h*(Tr-Ti)
        ]


    def Tcalculation(self,Toutguess,Tinguess,Tr, eps,k,t=self.t):
        x0=[Toutguess,Tinguess]
        sol=scipy.optimize.fsolve(self.Tcalculation_system,x0,args=(Tr, eps,k,t))
        self.T_calculated=sol[0]

class RegenerativeCool:
    def __init__(self):
        self.Q=0; #heat
        self.t=0; #thickness
        self.T_calculated=-1;

    def Tcalculation(self,Tr,Ti_co,A):

            q=(Tr-Ti_co)/(1/hg+t/self.Mater.k+1/self.hco)
            self.Q += q*A
            Tinext_co=Ti_co+q*A/(self.Prop.fcp*self.m_flow_fuel)
            T_wall=self.t/self.Mater.k*q+Ti_co+q/self.hco

            return Tinext_co, T_wall

    def pressureloss(m_flow_fuel,Dr,L):
        delta_p=self.f*m_flow_fuel**2/(2*self.Prop.density)*L/Dr

    def Run(self,Tr, hg, t, Prop ,Mater  ,Dr,A,Ti_co,Re,m_flow_fuel,L):
        self.Q=0
        self.Pr=4*Prop.f_gamma/(9*Prop.f_gamma-5)
        self.f=(1.82*math.log10(Re)-1.64)**(-2)
        self.Nu=self.f/8*(Re-1000)*self.Pr/(1+12.7*math.sqrt(self.f/8)*(self.Pr^2/3-1))
        self.hco=self.Nu*self.Mater.k/Dr
        
        self.Prop = Prop
        self.m_flow_fuel = m_flow_fuel

        self.t = t
        self.Mater = Mater()
        
        T_co_calcualted, T_wall_calcualted = self.Tcalculation(Tr,Ti_co,A)
        ploss=self.pressureloss(m_flow_fuel,Dr,L)

        #T_co_list=[0 for i in range(len(Tr))]
        #T_co_list[0]=Ti_co
        #for i in Tr:
         #   T_co_calcualted[i], T_wall_calcualted[i] = self.Tcalculation(Tr[i],T_co_list[i],A[i])
        return T_co_calcualted, T_wall_calcualted, ploss

    


        



