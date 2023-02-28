def injector(m, OF, A_s, A_I, rho_ox, rho_f):
    
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
    InjTypes = ['like', 'unlike', 'non', 'hybrid']
    InjType = InjTypes[t]
    
    m_ox, m_f = Massflow(m, OF)
    #liquid - liquid      
    if InjType == 'unlike':
        k_doublet = k_quadletA = k_quadletB = 1
        k_triplet = 1.58
        k_pentad = 9.35
        k = [k_doublet, k_triplet, k_quadletA, k_quadletB, k_pentad] 
        b_doublet = b_quadletA = b_quadletB = 2/3
        b_triplet = 4/7
        b_pentad = 4/5
        b = [b_doublet, b_triplet, b_quadletA, b_quadletB, b_pentad]  #in order as above
        R = k*((rho_ox/rho_f)*(m_ox/m_f)**2)**b #orifice area ratio
    
    zeta = 0.7
    
    v_iox = m_ox / (rho_ox*A_i_ox)
    v_if = m_f / (rho_f*A_i_f)
    dp_ox = zeta * 0.5 * rho_ox*v_iox**2
    dp_f = zeta * 0.5 * rho_ox*v_iox**2
    p_c = dp_ox / 0.1 #pintle
    
    ## like doublet : https://ntrs.nasa.gov/citations/19720010642
    D = 4.85e4*v_i**(-0.75)*d_j**0.57*(p_c/p_j)**(-0.52)
    
        #single jet
    D = 15.9e4*v_i**(-1)*d_j**0.57*(p_c/p_j)**(-0.1)
    
    ## unlike doublet (same ref)
    # single jet
    d_f <= d_o #for both
    D_f = 2.91e4*v_f**(-0.76*d_f**0.29)*(p_c/p_j)**(-0.65)*P_d**(0.165*(d_o/d_f)**(0.023))
    D_o = 2.72e4*v_o**(-0.57*d_o**0.65)*(p_c/p_j)**(-0.3)*P_d**(-0.25*(d_o/d_f)**(-0.17))
    return v_i_ox, v_i_f, dp_ox, dp_f

# V_inj and droplet size lvl 0
# p_c at end