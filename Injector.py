import numpy as np
import matplotlib.pyplot as plt

def injector1(default, propellant, p_c, m, OF):
    
    '''
    Computes injection velocities, droplet size abd pressure drop over injector. \n
    default = default class containing variables such as Cd, d_o, InjType..; \n
    propellant = propellant class containing variables such as rho_ox, rho_f..; \n
    p_c = initial chamber pressure; \n
    m = mass flow toward injector; \n
    OF = mass ratio. \n
    '''
    
    #Default variables
    d_ox = default.d_ox
    d_f = default.d_f
    C_d = default.Cd
    InjType = default.InjType
    mu_prop = default.mu_prop    # [lbm/(ft-s)], 1 lbm/(ft-s) = 1.4881639 Pa.s
    sig_prop = default.sig_prop  # [dynes/cm], 1 dyn/cm = 1e-7 N/m 
    rho_prop = default.rho_prop  # [lbm/ft3], 1 lbm/ft3 = 16.0185 kg/m3 
    p_center = default.p_center
    p_j = default.p_j
    
    #Propellant densities
    rho_ox = propellant.o_dens
    rho_f = propellant.f_dens_l    
    
    def Massflow(m, OF):
        m_ox = (OF/(OF+1.0))*m
        m_f = m - m_ox
        return m_ox, m_f
     
    #! Maybe remove this and just replace it for eta_s alone, allow user to 
    #  change it accordingly.
    ## For throttled 0.2, for Unthrottled engines 0.3
    if InjType == 'like':
        eta_s = 0.1 #10-15% [Humble et al.]
    elif InjType == 'unlike':
        eta_s = 0.2 #20-25% [Humble et al.]
    elif InjType == 'pintle':
        eta_s = 0.1 # [Humble et al.], [DARE, Sparrow]
        
    dp = eta_s * p_c
    v_iox = C_d * (2*dp/rho_ox)**0.5
    v_if = C_d * (2*dp/rho_f)**0.5
    m_ox, m_f = Massflow(m, OF)    
    n_ox = m_ox / (rho_ox * v_iox * (np.pi/4) * d_ox**2)
    n_f = m_f / (rho_f * v_if * (np.pi/4) * d_f**2) 
    print('n_f =', n_f, 'n_ox =', n_ox)
    
    mu_wax = 2.69e-3        # [lbm/(ft-s)], 1 lbm/(ft-s) = 1.4881639 Pa.s
    sig_wax = 17.0          # [dynes/cm], 1 dyn/cm = 1e-7 N/m     
    rho_wax = 47.7          # [lbm/ft3], 1 lbm/ft3 = 16.0185 kg/m3
    
    if InjType == 'like':                
        K_prop = ((mu_prop*sig_prop/rho_prop)/(mu_wax*sig_wax/rho_wax))**0.25 # https://ntrs.nasa.gov/citations/19760023196     
        
        # https://ntrs.nasa.gov/citations/19760023196
        d_j = d_f = d_ox             #orifice diameter
        v_j = (v_iox + v_if)/2.0    
        D_f = D_ox = (1.6e5*(v_j*3.28084)**(-1.0)*(p_center/p_j)**(-0.1)*(d_j*39.3701)**0.57*K_prop)*0.0254e-6
        # print(D_f)
    else: #liquid - liquid   
        K_prop = ((mu_prop*sig_prop/rho_prop)/(mu_wax*sig_wax/rho_wax))**0.25 # https://ntrs.nasa.gov/citations/19760023196
                
        # unlike doublet https://ntrs.nasa.gov/citations/19720010642
        # tested for d_f <= d_ox
        P_D = (rho_f*v_if**2.0)/(rho_ox*v_iox**2.0)     #momentum ratio
        #droplet size
        D_f = (2.91e4*(v_if*3.28084)**(-0.76)*(d_f*39.3701)**(0.29)*(p_center/p_j)**(-0.65)*P_D**(0.165)*(d_ox/d_f)**(0.023)*K_prop)*1e-6
        D_ox = (2.72e4*(v_iox*3.28084)**(-0.57)*(d_ox*39.3701)**(0.65)*(p_center/p_j)**(-0.3)*P_D**(-0.25)*(d_ox/d_f)**(-0.17)*K_prop)*1e-6
    
    A_f =  n_f*((np.pi/4) * d_f**2)
    A_ox = n_ox*((np.pi/4) * d_ox**2)  
    print('A_f =', A_f, 'A_ox =', A_ox)
    return v_iox, v_if, D_f, D_ox, dp, eta_s

def injector2(default, propellant, v_iox, v_if, D_f, D_ox, p_inj, eta_s):
    '''
    Computes chamber pressure after injector.
    
    ''' 
    
    #Default variables
    C_d = default.Cd
    InjType = default.InjType
    
    #Propellant densities
    rho_ox = propellant.o_dens
    rho_f = propellant.f_dens_l
    
    p_c = p_inj / (1 + eta_s)
    
    zeta = 1.0/ C_d**2.0
    dp_ox = zeta * 0.5 * rho_ox*v_iox**2.0
    dp_f = zeta * 0.5 * rho_f*v_if**2.0
    
    if (dp_ox/p_c) < eta_s:
        print('dp_ox (', InjType,') <', eta_s,' p_c!')
    elif (dp_f/p_c) < eta_s:
        print('dp_f (', InjType,') <', eta_s,' p_c!')
        
    return p_c, dp_ox, dp_f

def validateInj():
    p_c = 180e5
    C_d = 0.7
    m = 160
    OF = 3.55
    rho_ox = 1141.0
    rho_f = 71.0
    mu_prop = 2.69e-3    
    sig_prop = 17.0   
    rho_prop = 47.7
    v_iox, v_if, D_f, D_ox, dp, eta_s, InjType = injector1(p_c, C_d, m, OF, rho_ox, rho_f, mu_prop, sig_prop, rho_prop)
    p_inj = 182.1e5
    # print(injector1(p_c, C_d, m, OF, rho_ox, rho_f, mu_prop, sig_prop, rho_prop))
    # print(injector2(v_iox, v_if, D_f, D_ox, p_inj, C_d, rho_ox, rho_f, eta_s, InjType))
    
    d_f = d_ox = (0.0135 +0.281)*0.0254/2.0 
    v_i = np.arange(0, 200.1,0.1)
    P_D = (rho_f*v_if**2.0)/(rho_ox*v_iox**2.0)
    
    plt.figure()
    ax = plt.axes()
    D_f = (2.91e4*(v_i*3.28084)**(-0.76)*(d_f*39.3701)**(0.29)*(1)**(-0.65)*P_D**(0.165)*(d_ox/d_f)**(0.023)*6)*1e-6
    D_ox = (2.72e4*(v_i*3.28084)**(-0.57)*(d_ox*39.3701)**(0.65)*(1)**(-0.3)*P_D**(-0.25)*(d_ox/d_f)**(-0.17)*6)*1e-6
    ax.plot(v_i, D_f*1e6, label = 'Droplet size fuel')
    ax.plot(v_i, D_ox*1e6, label = 'Droplet size ox')
    ax.set_xlabel('Injection velocity [m/s]')
    ax.set_ylabel('Droplet size [microns]')
    # plt.ylim((0,120))
    ax.legend()
    plt.show()
# validateInj()