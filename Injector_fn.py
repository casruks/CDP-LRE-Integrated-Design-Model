import numpy as np
# injector(1712, 212, 1.2, 9.80665, 1137, 870, 1e-3, 1e-3, 5.52e5, 6.9e5, 1)

def injector(F, Isp, OF, g0, rho_ox, rho_f, D_i_ox, D_i_f, dp_ox, dp_f, K_comp):

    def Massflow(OF, F, Isp, n, g0):
        m = (F/(Isp*g0))/n
        m_ox = (OF/(OF+1))*m
        m_f = m - m_ox
        return [m_ox, m_f, m]

    #0.5-0.7. Changes with changing Re, is correct for only one set of conditions
    def Discharge(K_comp, A_s, A_I):
        
        zeta_c = K_comp*0.5*(1 - A_s/A_I)**(3/4)    # For sudden contraction
        zeta_e = (1 - A_s/A_I)**2                   # For sudden expansion
        zeta = zeta_c + zeta_e
        C_d = (1/zeta)**0.5
        return C_d, zeta

    def Velocity(C_d, dp, rho):
        return C_d * (2*dp/rho)**0.5

    m = Massflow(OF, F, Isp, 1, g0)
    A_i = [0.25*np.pi*D_i_ox**2, 0.25*np.pi*D_i_f**2]

    A_I_ox = 1
    A_I_f = 1
    v_i = [Velocity(Discharge(K_comp, A_i[0], A_I_ox)[0], dp_ox, rho_ox), \
           Velocity(Discharge(K_comp, A_i[1], A_I_f)[0], dp_f, rho_f)]
    
    n = [m[0]/(rho_ox*v_i[0]*A_i[0]), m[1]/(rho_f*v_i[1]*A_i[1])]
    
    m_i = [m[0]/n[0], m[1]/n[1]] 

    #Durting start transients approx:
    p_inj = [4*dp_ox, 4*dp_f]

    Q_i = [m_i[0]/ rho_ox, m_i[1] / rho_f]    

    return p_inj, v_i, m, m_i, Q_i, n
