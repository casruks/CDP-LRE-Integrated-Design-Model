def Nozzle_loop_1(Pc,F,Pamb,Propellant,Default):
    ## INPUTS:
    # Pc= Chamber pressure in bar
    # F=Thrust in Newtons
    # Pamb= Ambient pressure in bar
    # Propellant.Ox_name= string with oxydizer name for RocketCEA
    # Propellant.Fuel_name= string with fuel name for Rocket CEA
    # Propellant.Frozen_state: 0= equilibrium everywhere, 1= equilibrium in the convergent and frozen in the divergent
    # Default.MR: 0=calculate for best c*, 'every_other_value'=Mixture ratio to use
    # Default.De_max= maximum exit diameter in meters
    # Default.toll_c_star= tollerance for the best MR calculation
    # Default.toll_F_obj= tollerance on the target force to be achieved during mass flow rate calculation
    # Default.Max_iterations_mass_flow= max iterations for the throat area and mass flow calculations
    # Default.toll_P_adapted= tollerance on difference between exit pressure and ambient pressure if nozzle can be adapted

    ## OUTPUTS:
    # m_p = mass flow rate kg/s
    # Tc= combustion chamber temperature K
    # MR= mixture ratio
    # At= throat area m^2
    # eps= expansion ratio
    # Isp= Specific impulse s

    from rocketcea.cea_obj_w_units import CEA_Obj
    import math as mth
    import matplotlib.pyplot as plt
    import numpy as geek
    
    Ox=Propellant.Ox_name
    Fuel=Propellant.Fuel_name
    frozen_state=Propellant.Frozen_state

    MR=Default.MR
    De_max=Default.De_max
    toll_c_star=Default.toll_c_star
    toll_F_obj=Default.toll_F_obj
    Max_iterations_mass_flow=Default.Max_iterations_mass_flow
    toll_P_adapted=Default.toll_P_adapted

    ispObj = CEA_Obj( oxName=Ox, fuelName=Fuel,cstar_units='m/s',pressure_units='bar',temperature_units='K',isp_units='sec',density_units='kg/m^3',specific_heat_units='J/kg-K',viscosity_units='poise',thermal_cond_units='W/cm-degC')

    if MR==0:
        MR_1=0.01
        MR_2=150

        diff=MR_1-MR_2
        MR_prev=MR_1;

        while abs(diff)>toll_c_star:
            
            MR_curr=(MR_2+MR_1)/2
            diff=MR_curr-MR_prev

            c_star_curr=ispObj.get_Cstar(Pc=Pc,MR=MR_curr)

            MR_curr_low=MR_curr-0.001
            MR_curr_high=MR_curr+0.001

            c_star_curr_low=ispObj.get_Cstar(Pc=Pc,MR=MR_curr_low)

            if c_star_curr>c_star_curr_low:
                MR_1=MR_curr;
            else:
                MR_2=MR_curr;
            MR_prev=MR_curr;
        MR=MR_curr #Mixture ratio
        c_star=c_star_curr;
    else:
        c_star=ispObj.get_Cstar(Pc=Pc,MR=MR);

    Tc=ispObj.get_Tcomb(Pc=Pc,MR=MR) # Function that returns the combustion chamber temperature

    v_eff=-1
    it=0

    At=0.001 # Initial value for throat area, to be better defined !!!!!!!!!!!!!!

   
    Ae_max=mth.pi*De_max**2/4 #Maximum exit area of the nozzle

 
    variation=10 # Value to enter the loop

    while v_eff<0 or variation>toll_F_obj:

        eps_max=Ae_max/At # Maximum expansion ratio with this throat area
        Pratio_max=ispObj.get_PcOvPe(Pc=Pc,MR=MR,eps=eps_max,frozen=frozen_state,frozenAtThroat=frozen_state)
        Pe_max=Pc/Pratio_max

        if Pe_max>Pamb:
            difference=0.001
            Pe=Pe_max
            Ae=Ae_max; #In this case the maximum exit area cannot reach adapted conditions and we will simply take that
        else:
            difference=0.1
            Ae_1=At+0.001; # In this case we enter the loop and start the iteration with Minimum exit area with this throat area
        
        while difference>toll_P_adapted:

            eps_1= Ae_1/At # First expansion ratio 
            Pratio_1=ispObj.get_PcOvPe(Pc=Pc,MR=MR,eps=eps_1,frozen=frozen_state,frozenAtThroat=frozen_state)
            Pe_1=Pc/Pratio_1 # Exit pressure 1
            
            difference=abs(Pe_1-Pamb)/Pamb 

            if difference>0.01:

                Ae_2=Ae_1*Pe_1/Pamb # Change value for next iteration
                Ae_1=Ae_2;
            else:
                Ae=Ae_1 # Finish iteration if we have reached the convergence
                Pe=Pe_1;
            
        eps_actual=Ae/At;
        
        Isp_it=ispObj.estimate_Ambient_Isp(Pc=Pc,MR=MR,eps=eps_actual,Pamb=Pamb,frozen=frozen_state,frozenAtThroat=frozen_state) # Calculates Isp for this iteration
        v_eff_it=Isp_it[0]*9.80665; # Calculates effective velocity for this iteration

        m_p_it=Pc*100000*At/c_star # Mass flow rate from continuity equation in the throat

        F_it=m_p_it*v_eff_it # Force computed at this iteration
    
        variation=abs(F_it-F)/F # Variation from target force

        if variation>toll_F_obj: # If we have to iterate again, we change the throat area
            mp2=m_p_it*(F/F_it);
            At=c_star*mp2/(Pc*100000);
    
        it=it+1;
        if it>Max_iterations_mass_flow:
            break;
        v_eff=v_eff_it;
    
    m_p=m_p_it # Propellant mass flow
    eps=Ae/At # Expansion ratio
    Isp=Isp_it # Isp value

    # Combustion chamber properties

    rhos=ispObj.get_Densities(Pc=Pc,MR=MR,eps=eps_actual,frozen=frozen_state,frozenAtThroat=frozen_state)
    Transp_c=ispObj.get_Chamber_Transport(Pc=Pc,MR=MR,eps=eps_actual,frozen=frozen_state)

    rho_c=rhos[0]
    cp_c=Transp_c[0]
    mu_c=Transp_c[1]
    k_c=Transp_c[2]/100
    Pr_c=Transp_c[3]

    return m_p,Tc,MR,At,eps,Isp[0],rho_c,cp_c,mu_c,k_c,Pr_c



