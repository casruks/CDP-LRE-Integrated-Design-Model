import numpy as np
# injector(1712, 212, 1.2, 9.80665, 1137, 870, 1.93e-3, 1.77e-3, 5.52e5, 6.9e5, 1)

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
    #C_d = Discharge(K_comp, A_i[0], A_I_ox)[0]
    #C_d = Discharge(K_comp, A_i[1], A_I_f)[0], placeholders Q1 are below:
    v_i = [Velocity(0.87, dp_ox, rho_ox), \
           Velocity(0.91, dp_f, rho_f)]
    
    n = [m[0]/(rho_ox*v_i[0]*A_i[0]), m[1]/(rho_f*v_i[1]*A_i[1])]
    
    m_i = [m[0]/n[0], m[1]/n[1]] 

    #Durting start transients approx:
    p_inj = [4*dp_ox, 4*dp_f]

    Q_i = [m_i[0]/ rho_ox, m_i[1] / rho_f]    

    return p_inj, v_i, m, m_i, Q_i, n

## Reference values - LRE combustor Q1
C_dox = 0.87
C_df = 0.91
rho_ox = 1137
rho_f = 817
OF = 1.2
p_c = 20.68e5
p_if = 27.58e5
p_iox = 26.2e5
dp_ox = p_iox - p_c
dp_f = p_if - p_c
n = 5
F = 1712
Isp = 212

v_i_ox = C_dox * (2*dp_ox/rho_ox)**0.5
v_i_f = C_df * (2*dp_f/rho_f)**0.5

m = F/ (Isp*9.80665)
m_ox = (OF/(OF+1))*m
m_f = m - m_ox

m_i_ox = m_ox / n
m_i_f = m_f / n

Q_i_ox = m_i_ox / rho_ox 
Q_i_f = m_i_f / rho_f

A_i_ox = Q_i_ox / v_i_ox
A_i_f = Q_i_f / v_i_f

D_i_ox = (4*A_i_ox/np.pi)**0.5
D_i_f = (4*A_i_f/np.pi)**0.5
print(D_i_ox*1000, D_i_f*1000)
