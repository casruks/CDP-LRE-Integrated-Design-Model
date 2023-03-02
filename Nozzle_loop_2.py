def Nozzle_loop(Pc,F,Pamb,Tc,Propellant,Material,Nozzle_type,MR,eps,At,m_p,Dc,Default):
    ## INPUTS:
    # Default
    # Pc
    # Tc
    # Propellant
    # Material
    # Nozzle_type
    # MR
    # eps
    # At
    # m_p
    # Dc
    ## Outputs
    # t_noz Thickness of the nozzle
    # x_noz Points on the axis at which properties are evaluated
    # y_noz Geometry of the nozzle
    # Tw_ad_noz Adiabatic wall temperature in the nozzle
    # h_c_noz Convective heat coefficient in the nozzle
    # P_noz Pressure in the nozzle
    # T_noz Temperature in the nozzle

    from rocketcea.cea_obj import add_new_fuel, add_new_oxidizer
    from rocketcea.cea_obj_w_units import CEA_Obj
    from rocketcea import blend
    import math as mth
    import matplotlib.pyplot as plt
    import numpy as geek
    
    Ox=Propellant.Ox_name
    Fuel=Propellant.Fuel_name
    frozen_state=Propellant.Frozen_state
    
    ispObj = CEA_Obj( oxName=Ox, fuelName=Fuel,cstar_units='m/s',pressure_units='bar',temperature_units='K',isp_units='sec',density_units='kg/m^3',specific_heat_units='J/kg-K',viscosity_units='poise')

    T_w=Material.T_wall
    theta_con=Default.Theta_con
    Theta_conical=(Default.Theta_conical*mth.pi)/180
    Theta_bell=(Default.Theta_bell*mth.pi)/180
    TH_exit_bell=(Default.TH_exit_bell*mth.pi)/180

    # Definition of the geometry of the nozzle
    R_t=mth.sqrt(At/(mth.pi))
    R_u=R_t*Default.R_u_ratio
    con_ratio=Dc**2/((2*R_t)**2)

    if Nozzle_type==0:
        L_nozzle_div=((mth.sqrt(eps)-1)*R_t+R_u*(1/mth.cos(Theta_conical)-1))/mth.tan(Theta_conical)
        L_nozzle_con=((mth.sqrt(con_ratio)-1)*R_t+R_u*(1/mth.cos(theta_con)-1))/mth.tan(theta_con)
        L_tot=L_nozzle_con+L_nozzle_div

        xp=L_nozzle_con+R_u*mth.sin(Theta_conical)
        yp=R_t+(1-mth.cos(Theta_conical))*R_u

        xp_con=L_nozzle_con-R_u*mth.sin(theta_con)
        yp_con=R_t+(1-mth.cos(theta_con))*R_u

        x_con=geek.linspace(0,xp_con,num=1000)
        a_con=(yp_con-Dc/2)/(xp_con)
        y_con=Dc/2+a_con*x_con # With this part we have defined the coordinates for the convergent geometry

        x_throat1=geek.linspace(xp_con,L_nozzle_con,num=25)
        th_step=[]
        y_throat1=[]
        for i in x_throat1:
            th_step_cur=mth.asin((L_nozzle_con-i)/R_u)
            th_step.append(th_step_cur)
            y_cur=R_t+(1-mth.cos(th_step_cur))*R_u
            y_throat1.append(y_cur) # With this part we have defined the coordinates for the first part of the throat

        x_throat2=geek.linspace(L_nozzle_con,xp,num=25)#With this part we have defined the coordinates for the second part of the throat

        th_step=[]
        y_throat2=[]
        for i in x_throat2:
            th_step_cur=mth.asin((i-L_nozzle_con)/R_u)
            th_step.append(th_step_cur)
            y_cur=R_t+(1-mth.cos(th_step_cur))*R_u
            y_throat2.append(y_cur)

        x_div=geek.linspace(xp,L_tot,num=500)
        a_div=(mth.sqrt(eps)*R_t-yp)/(L_tot-xp)
        y_div=yp+a_div*(x_div-xp) # With this part we have defined the coordinates for the divergent part of the nozzle


        x_noz=geek.concatenate((x_con,x_throat1,x_throat2,x_div))
        y_noz=geek.concatenate((y_con,y_throat1,y_throat2,y_div)) 
    else:
        ye=mth.sqrt(eps)*R_t
        xp=0.382*R_t*mth.sin(Theta_bell)
        yp=1.382*R_t-0.382*R_t*mth.cos(Theta_bell)
        a=(mth.tan(mth.pi/2-TH_exit_bell)-mth.tan(mth.pi/2-Theta_bell))/(2*(ye-yp))
        b=mth.tan(mth.pi/2-Theta_bell)-2*(mth.tan(mth.pi/2-TH_exit_bell)-mth.tan(mth.pi/2-Theta_bell))/(2*(ye-yp))*yp
        c=xp-a*yp**2-b*yp
        L_nozzle_div=a*ye**2+b*ye+c
        L_nozzle_con=((mth.sqrt(con_ratio)-1)*R_t+R_u*(1/mth.cos(theta_con)-1))/mth.tan(theta_con)
        L_tot=L_nozzle_con+L_nozzle_div

        xp_con=L_nozzle_con-R_u*mth.sin(theta_con)
        yp_con=R_t+(1-mth.cos(theta_con))*R_u
        
        x_con=geek.linspace(0,xp_con,num=100)
        a_con=(yp_con-Dc/2)/(xp_con)
        y_con=Dc/2+a_con*x_con # With this we have defined the geometry of the convergent part

        y_throat1=[]
        th_step=[]
        x_throat1=geek.linspace(xp_con,L_nozzle_con,num=25)
        for i in x_throat1:
            th_step_cur=mth.asin((L_nozzle_con-i)/R_u)
            th_step.append(th_step_cur)
            y_cur=R_t+(1-mth.cos(th_step_cur))*R_u
            y_throat1.append(y_cur) # With this part we have defined the coordinates for the first part of the throat

        x_throat2=geek.linspace(L_nozzle_con,xp,num=25)
        th_step=[]
        y_throat2=[]
        for i in x_throat2:
            th_step_cur=mth.asin((i-L_nozzle_con)/(0.382*R_t))
            th_step.append(th_step_cur)
            y_cur=R_t+(1-mth.cos(th_step_cur))*R_t*0.382
            y_throat2.append(y_cur)


        x_div=geek.linspace(xp,L_tot,num=500)
        Delta=b**2-4*a*(c-x_div)
        y_div=(-b+mth.sqrt(Delta))/(2*a) # With this part we have defined the coordinates for the divergent part of the nozzle
        
        x_noz=geek.concatenate((x_con, x_throat1, x_throat2, x_div))
        y_noz=geek.concatenate((y_con, y_throat1, y_throat2, y_div));
    
    A_noz=mth.pi*y_noz**2

    # We are now going to determine the pressure in each section of the nozzle
    
    Ts=ispObj.get_Temperatures(Pc=Pc,MR=MR,eps=eps,frozen=frozen_state,frozenAtThroat=frozen_state)
    rhos=ispObj.get_Densities(Pc=Pc,MR=MR,eps=eps,frozen=frozen_state,frozenAtThroat=frozen_state)
    cps=ispObj.get_HeatCapacities(Pc=Pc,MR=MR,eps=eps,frozen=frozen_state,frozenAtThroat=frozen_state)
    gs=ispObj.get_Throat_MolWt_gamma(Pc=Pc,MR=MR,esp=eps,frozen=frozen_state)
    gsc=ispObj.get_Chamber_MolWt_gamma(Pc=Pc,MR=MR,esp=eps)

    Transp_c=ispObj.get_Chamber_Transport(Pc=Pc,MR=MR,eps=eps,frozen=frozen_state)
    Transp_t=ispObj.get_Throat_Transport(Pc=Pc,MR=MR,eps=eps,frozen=frozen_state)

    mu_c=Transp_c[1]
    mu_t=Transp_t[1]

    T_t=Ts[1]
    rho_t=rhos[1]
    rho_c=rhos[0]
    cp_c=cps[0]
    cp_t=cps[1]
    g_t=gs[1]
    g_c=gsc[1]
    R_t=cp_t*(g_t-1)/g_t
    R_c=cp_c*(g_c-1)/g_c

    P_t=rho_t*R_t*T_t # Pressure in the throat

    # For the convergent part we assume a linear decreas in pressure
    x_1=geek.concatenate((x_con,x_throat1))
    a_P_con=(P_t-Pc)/(x_1[-1]-x_1[0])
    P_con=Pc+a_P_con*x_1

    y_2=geek.concatenate((y_throat2,y_div))
    A_div=mth.pi*y_2**2 #Areas sampled at the various locations in the divergent

    P_div=[]
    for i in A_div:
        eps_it=i/At
        dp=ispObj.get_PcOvPe(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state,frozenAtThroat=frozen_state)
        Pe_it=Pc/dp
        P_div.append(Pe_it);
    
    P_noz=geek.concatenate((P_con,P_div)) # Pressure throughout the nozzle

    # We assume that properties change linearly in the convergent

    a_T_con=(T_t-Tc)/(x_1[-1]-x_1[0])
    T_con=Tc+a_T_con*x_1

    T_div=[]
    for i in A_div:
        eps_it=i/At
        T_s=ispObj.get_Temperatures(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state,frozenAtThroat=frozen_state)
        T_it=T_s[2]
        T_div.append(T_it);
    
    T_noz=geek.concatenate((T_con,T_div)) # Temperature throughout the nozzle

    a_R_con=(R_t-R_c)/(x_1[-1]-x_1[0])
    R_con=R_c+a_R_con*x_1 # R of gas in the convergent

    a_g_con=(g_t-g_c)/(x_1[-1]-x_1[0])
    g_con=g_c+a_g_con*x_1 # gamma of gas in the convergent

    a_mu_con=(mu_t-mu_c)/(x_1[-1]-x_1[0])
    mu_con=mu_c+a_mu_con*x_1  # viscosity of gases in the convergent

    a_cp_con=(cp_t-cp_c)/(x_1[-1]-x_1[0])
    cp_con=cp_c+a_cp_con*x_1 # heat capacity in the convergent

    rho_con=P_con/(R_con*T_con)
    u_con=mth.sqrt(T_con*g_con*R_con)
    y_1=geek.concatenate((y_con,y_throat1))
    y_2=geek.concatenate((y_throat2,y_div))
    A_con=mth.pi*y_1**2

    v_con=m_p/(rho_con*A_con)
    M_con=v_con/u_con

    Pr_con=4*g_con/(9*g_con-5) # These are all approximations

    r_con=Pr_con**(1/3)

    Tw_ad_con=T_con*(1+r_con*(g_con-1)/2*M_con**2)

    x_2=geek.concatenate((x_throat2,x_div))
    rho_div=[]
    cp_div=[]
    g_div=[]
    v_div=[]
    mu_div=[]
    Pr_div=[]
    k_div=[]
    Mach_div=[]

    for i in A_div:
        eps_it=i/At
        rhos_div=ispObj.get_Densities(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state,frozenAtThroat=frozen_state)
        cps_div=ispObj.get_HeatCapacities(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state,frozenAtThroat=frozen_state)
        gs_div=ispObj.get_exit_MolWt_gamma(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state,frozenAtThroat=frozen_state)
        u_sound=ispObj.get_SonicVelocities(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state,frozenAtThroat=frozen_state)
        Mach_n=ispObj.get_MachNumber(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state,frozenAtThroat=frozen_state)
        Transp=ispObj.get_Exit_Transport(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state)

        u_it=u_sound(2)
        v_it=u_it*Mach_n
        rho_div_it=rhos_div[2]
        cp_div_it=cps_div[2]
        g_div_it=gs_div[2]
        mu_div_it=Transp[1]
        Pr_div_it=Transp[3]
        k_div_it=Trans[2]
        rho_div.append(rho_div_it)
        cp_div.append(cp_div_it)
        g_div.append(g_div_it)
        v_div.append(v_it)
        mu_div.append(mu_div_it)
        Pr_div.append(Pr_div_it)
        Mach_div.append(Mach_n)
        k_div.append(k_div_it);
    
    r_div=Pr_div**(1/3)
    Tw_ad_div=T_div*(1+r_div*(g_div-1)/2*Mach_div**2)

    Tw_ad_noz=geek.concatenate((Tw_ad_con,Tw_ad_div))

    T_f_div=0.5*T_w+0.28*T_div+0.22*Tw_ad_div # Film temperature in the divergent calculated as in Ziebland
    T_f_con=0.5*T_w+0.28*T_con+0.22*Tw_ad_con # Film temperature in convergent

    a_b=0.026
    h_c_div=1.213*a_b*m_p**0.8*mu_div**0.2*cp_div*Pr_div**(-0.6)*(2*y_2)**(-1)*(Tc/T_f_div)**0.68
    h_c_con=1.213*a_b*m_p**0.8*mu_con**0.2*cp_con*Pr_con**(-0.6)*(2*y_1)**(-1)*(Tc/T_f_con)**0.68

    h_c_noz=geek.concatenate((h_c_con,h_c_div))

    sig=Material.sigma_y
    SF=Default.Safety_factor

    t_noz=SF*P_noz*y_noz/sig #Thickness of the wall in the nozzle

    return t_noz,x_noz,y_noz,Tw_ad_noz,h_c_noz,P_noz,T_noz;



