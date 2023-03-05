## Notice:  - orifice areas, currently average from literature & A_iox=A_if
#           - C_d assumed to be 0.7


import numpy as np

def injector1(C_d, m, OF, rho_ox, rho_f, mu_prop, sig_prop, rho_prop):
    
    #PLACEHOLDERS (incl. userinput)
    # C_d = 0.7
    d_o = d_f = (0.0135 +0.281)*0.0254/2     # 3.74 mm # http://libgen.rs/book/index.php?md5=3D236B9BDD4070690CA83952058D9A1F p.113
    #D_i = (0.02 +0.08)*0.0254/2         # 1.27 mm
    p_center = 1 
    p_j = 1
    d_j = d_o    
    
    ##User Input: 
    # Unlike-doublet,-triplet,-quadletA,-quadletB,-pentad or like-doublet impingement? (pintle not yet)
    # orificie diameter, restrict to usual range : d_o, d_f # http://libgen.rs/book/index.php?md5=3D236B9BDD4070690CA83952058D9A1F p.113
    
    def Massflow(m, OF):
        m_ox = (OF/(OF+1))*m
        m_f = m - m_ox
        return m_ox, m_f
    
    #0.5-0.7. Changes with changing Re, is correct for only one set of conditions
    def Discharge(A_s, A_I):
        K_comp = 1
        zeta_c = K_comp*0.5*(1 - A_s/A_I)**(3/4)    # For sudden contraction
        zeta_e = (1 - A_s/A_I)**2                   # For sudden expansion
        zeta = zeta_c + zeta_e
        C_d = (1/zeta)**0.5
        return C_d, zeta
    
    t = 0
    InjTypes = ['like', 'unlike', 'pintle']
    InjType = InjTypes[t]
    
    m_ox, m_f = Massflow(m, OF)
        
    A_iox = A_if = 0.25*np.pi*d_o**2      ##!
    v_iox = m_ox / (rho_ox*A_iox)
    v_if = m_f / (rho_f*A_if)
        
    ### Droplet size (and ideally orifice ratio/ size in future)
    #liquid - liquid      
    if InjType == 'unlike':
        mu_wax = 2.69e-3    # [lbm/(ft-s)], 1 lbm/(ft-s) = 1.4881639 Pa.s
        sig_wax = 17        # [dynes/cm], 1 dyn/cm = 1e-7 N/m     
        rho_wax = 47.7          # [lbm/ft3], 1 lbm/ft3 = 16.0185 kg/m3
        K_prop = ((mu_prop*sig_prop/rho_prop)/(mu_wax*sig_wax/rho_wax))**0.25 # https://ntrs.nasa.gov/citations/19760023196
        k_doublet = k_quadletA = k_quadletB = 1
        k_triplet = 1.58
        k_pentad = 9.35
        k = [k_doublet, k_triplet, k_quadletA, k_quadletB, k_pentad] 
        b_doublet = b_quadletA = b_quadletB = 2/3
        b_triplet = 4/7
        b_pentad = 4/5
        b = [b_doublet, b_triplet, b_quadletA, b_quadletB, b_pentad]  #in order as above
        R = k*((rho_ox/rho_f)*(m_ox/m_f)**2)**b #orifice area ratio A_iox/A_if
        A_iox = R * A_if
        
        # unlike doublet https://ntrs.nasa.gov/citations/19720010642
        # tested for d_f <= d_o 
        P_D = (rho_f*v_if**2)/(rho_ox*v_iox**2)
        #droplet size
        D_f = 2.91e4*v_if**(-0.76)*d_f**(0.29)*(p_center/p_j)**(-0.65)*P_D**(0.165)*(d_o/d_f)**(0.023)*K_prop
        D_o = 2.72e4*v_iox**(-0.57)*d_o**(0.65)*(p_center/p_j)**(-0.3)*P_D**(-0.25)*(d_o/d_f)**(-0.17)*K_prop
        
    elif InjType == 'like':
        mu_wax = 2.69e-3    # [lbm/(ft-s)], 1 lbm/(ft-s) = 1.4881639 Pa.s
        sig_wax = 17        # [dynes/cm], 1 dyn/cm = 1e-7 N/m     
        rho_wax = 47.7          # [lbm/ft3], 1 lbm/ft3 = 16.0185 kg/m3
        K_prop = ((mu_prop*sig_prop/rho_prop)/(mu_wax*sig_wax/rho_wax))**0.25 # https://ntrs.nasa.gov/citations/19760023196
        ### D,D_f,D_o = droplet size
        # https://ntrs.nasa.gov/citations/19760023196 
        v_j = (v_iox + v_if)/2
        D_f = D_o = 1.6e5*v_j**(-1)*(p_center/p_j)**(-0.1)*d_j**0.57*K_prop
        
        ## like doublet : https://ntrs.nasa.gov/citations/19720010642
        #D_f = D_o = 4.85e4*v_j**(-0.75)*d_j**0.57*(p_center/p_j)**(-0.52)
        
            #single jet
        #D_f = D_o = 15.9e4*v_j**(-1)*d_j**0.57*(p_center/p_j)**(-0.1)

    return v_iox, v_if, D_f, D_o

def injector2(v_iox, v_if, D_f, D_o, p_inj, C_d, rho_ox, rho_f):
    
    t = 0
    InjTypes = ['like', 'unlike', 'pintle']
    InjType = InjTypes[t]
    
    ## For throttled 0.2, for Unthrottled engines 0.3
    if InjType == 'unlike':
        eta_dp = 0.1 #10-15% [Humble et al.]
    elif InjType == 'like':
        eta_dp = 0.2 #20-25% [Humble et al.]
    elif InjType == 'pintle':
        eta_dp = 0.1 #Oscar [DARE, Sparrow]
    
    zeta = 1/ C_d**2
    dp_ox = zeta * 0.5 * rho_ox*v_iox**2
    dp_f = zeta * 0.5 * rho_f*v_if**2
    
    p_c = p_inj - dp_ox
    p_c = p_inj - dp_f
    
    if (dp_ox/p_c) < eta_dp:
        print('dp_ox (', InjType,') <', eta_dp,' p_c!')
        
    return p_c, dp_ox, dp_f