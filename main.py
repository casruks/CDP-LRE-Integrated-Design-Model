import Turbomachinery_code as Turbo
import CombustionChamber as Comb
import Injector as Inj
import Igniters as Ign
import Reliability as Rel
import Nozzle_loop_1 as Nz_1
import Nozzle_loop_2 as Nz_2
import Cooling
import Materials as Ms
import numpy as np
import Aux_classes as aux
import math
import GUI


#Main Function
def Main(d : list[aux.Data], default : aux.Default, prop : aux.Propellant, com : GUI.Communicate):
    p_old = 0.0
    p_new = 19000000/9999000*d[-1].Thrust+500000
    bool = 0 #this variable is used to show the combustor function we are in the first loop
    regCool=Cooling.RegenerativeCool(); #inicialise cooling

    it=0;
    #First loop, pressure convergence
    while abs(p_new-p_old)/p_new > default.pres_tol:
        #PROGRESS BAR
        prog = min(abs(p_new-p_old)/p_new, 0.99)
        com.progress.emit(int((1.0-prog)*100))

        p_old = p_new
        
        #Compute nozzle (1)
        #Inputs: 
        ## Chamber pressure in bars
        ## Thrust in Newton
        ## Ambient pressure in bars
        ## Propellants class
        ## Default class
        d[-1].m_nozz,d[-1].Tc,d[-1].O_F,d[-1].At,d[-1].Eps,d[-1].Isp,rho_c,cp_c,mu_c,k_c,Pr_c,d[-1].cstar,errors_nz1,warnings_nz1 = Nz_1.Nozzle_loop_1(p_new/100000.0,d[-1].Thrust,d[-1].Pa/100000.0,prop,default,default.Nozzle_type)
        # Outputs:
        ## mass flow rate in kg/s
        ## Chamber temperatures in K
        ## Oxidizer to fuel ratio
        ## Throat area in m^2
        ## Expansion ratio
        ## Specific impulse in s
        ## Density of gas in comustion chamber kg/m^3
        ## Specific heat at constant pressure in chamber in J/kgK
        ## Viscosity (mu) in combustion chamber in Pa*s
        ## Conduction constant in combustion chamber in W/(m*degC), it's in degC but as long as it is used with temperature differences it shouldn't make a difference
        ## Prandtl number in combustion chamber
        ## Error value: 1== there has been a fatal error
        #if(error==1): return False
        
        if errors_nz1!=0:
            return errors_nz1,0,0,0,0,0,0,0,0,warnings_nz1,0,0,0,0,0,0,0,0,0;
            
        #Compute injector (1)
        ## Added A_est, representing estimated total orifice area (A_ox+A_f) for sanity check with combustion chamber dimensions..
        d[-1].v_iox, d[-1].v_if, d[-1].D_f, d[-1].D_ox, d[-1].dp, d[-1].eta_s, d[-1].m_ox, d[-1].m_f, d[-1].n_ox, d[-1].n_f, A_est, er_inj1, wr_inj1 = Inj.injector1(default, prop, p_new, d[-1].m_nozz, d[-1].O_F)

        if er_inj1!=0:
            return errors_nz1,0,0,0,er_inj1,0,0,0,0,warnings_nz1,0,0,0,wr_inj1,0,0,0,0,0;
        #outputs
        #Compute Chamber
        #inputs
        ## Pressure at chamber in bars
        ## Throat area in m^2
        ## Material
        ## default
        ## fuel injection velocity
        ## oxidizer injection velocity
        ## Chamber temperature in K
        ## oxidizer to fuel ratio
        ## bool variable that specifies if inside loop or not
        ## density, cp, miu,k, prandlt - this all comes from nozzle just leave it like that
        ## IMPORTANT - this function is currently using hardcoded droplet diameter, bc droplet diameter coming from injector does not make sense.
        d[-1].h_comb, d[-1].Dc, d[-1].ThicknessChamber, d[-1].Chamber_L, d[-1].Re_c,wr_comb,er_comb= Comb.CombustionChamber(p_new, d[-1].At, prop, default.chamber_mat_select, default, d[-1].v_if, d[-1].v_iox, d[-1].D_f, d[-1].D_ox, bool,rho_c,cp_c,mu_c/10,k_c,Pr_c,A_est,d[-1].O_F)
        
        if er_comb!=0:
            return errors_nz1,0,er_comb,0,er_inj1,0,0,0,0,warnings_nz1,0,wr_comb,0,wr_inj1,0,0,0,0,0;
        #outputs
        ## conductive heat transfer coefficient
        ## chamber diameter in m
        ## chamber thickness in m
        ## chamber length in m
        ## reynolds number

        #COmpute nozzle (2)
        #Inputs
        ## Chamber pressure in bars
        ## Combustion chamber temperature in K
        ## Propellants class
        ## Material properties through material class
        ## Nozzle type (from default class)
        ## O_F ratio
        ## Expansion ratio
        ## Throat area in m^2
        ## Mass flow rate in kg/s
        ## Combustion chamber diameter in m
        ## Default class
        t_noz,x_noz,y_noz,Tw_ad_noz,h_c_noz,d[-1].Dt,d[-1].De,d[-1].L_con,d[-1].L_div,d[-1].L_total,x_noz_cool,y_noz_cool,errors_nz2,warnings_nz2=Nz_2.Nozzle_loop(p_new/100000.0, d[-1].Tc, prop, default.noz_mat_select, default.Nozzle_type, d[-1].O_F, d[-1].Eps, d[-1].At, d[-1].m_nozz, d[-1].Dc, default)
       
        #Outputs
        ## Array of thickness at the discretized points in the nozzle (corresponding to x_noz) in m
        ## Array of discretized positions in the nozzle where all variables and properties are calcualted (positions in m)
        ## Array of radius in the nozzle at the various discretization points (corresponding to x_noz) in m
        ## Array of adiabatic wall temperature at the various discretization points (corresopnding to x_noz) in K
        ## Array of coefficients of convective heat transfer at the various discretization points (corresponding to x_noz) in W/(m^2K)
        ## Diameter of the throat in m
        ## Exit diameter in m
        ## Length of nozzle convergent in m
        ## Length of nozzle divergent in m
        ## Total length of the nozzle in m
        ## Array of x coordinates used for the cooling properties
        ## Array of y coordinates used for the cooling properties
        d[-1].x_nozz=x_noz
        d[-1].y_nozz=y_noz
        if errors_nz2!=0:
            return errors_nz1,errors_nz2,er_comb,0,er_inj1,0,0,0,0,warnings_nz1,warnings_nz2,wr_comb,0,wr_inj1,0,0,0,0,0;
        #Compute regenerative
       
        #Tf_cool,dptcool=regCool.Run_for_Toperating1D(Tw_ad_noz, h_c_noz, t_noz,prop,Ms.Rhenium,default.A,default.T_fuel_tanks,d[-1].m_nozz/(1.0+d[-1].O_F)/default.n,x_noz_cool[-1],y_noz_cool)
        #Tf_cool,dptcool_c=regCool.Run_for_Toperating0D(d[-1].Tc, d[-1].h_comb, d[-1].ThicknessChamber,prop,Ms.Rhenium,d[-1].Chamber_L*default.Dr,Tf_cool,d[-1].m_nozz/(1.0+d[-1].O_F)/default.n,d[-1].Chamber_L)
        #dptcool=dptcool+dptcool_c
        Coolobj = Cooling.CoolingClass()
        Coolobj_c = Cooling.CoolingClass()
        Coolobj.Q = 0
        Coolobj_c.Q = 0
        nozzle_mass = Ms.Nozzle_mass(x_noz,y_noz,t_noz,default.noz_mat_select)
        chamber_mass= Ms.Chamber_mass(d[-1].Dc,d[-1].Chamber_L,d[-1].ThicknessChamber,default.chamber_mat_select)
        #Inputs 
        #Inicial temperature in K
        # heat capacitance of the material in J/K
        #operation time in seconds
        #mass of the nozzle in kg
        #emissivity
        #adiabatic wall temperature of the nozzle np.array
        #convection coefficient in W/(m2 K) np.array
        #Thickness of the nozzle in m np.array
        #coating thickness in m 
        #propellant
        #material of the nozzle
        #material of the coating
        #Hydralic diamter of the cooling system piping in m 
        #Segment area of contact in m2
        #inicial coolant temperature in K
        #coolant mass flow in kg/s
        #length of the nozzle in m
        #Radius of the nozzle in m
        # Option to overwrite area with input area bool
        # which regenerative function to call (default 0, do not change!)
        #


        alpha = 2 * math.pi * default.perimeter_percentage / default.n
        y_for_cooling_channel = np.amin(y_noz_cool)
        #A_nozzle = [x_noz_cool[-1] * y_for_cooling_channel * alpha]
        A_nozzle=Cooling.Nozzle_area_calculation(alpha,y_noz_cool,x_noz_cool)
        t_noz_cooling_reg=[default.cooling_thickness for i in range(len(Tw_ad_noz))]
        (
            Tf_cool,
            Tw_wall_nozzle_calculated,
            dptcool,
            _,
            type_variable_nozzle,
            T_outer_wall_nozzle,
            err_nozzle_cooling,
            warn_nozzle_cooling,
        ) = Coolobj.Run_cooling(
            default.T0,
            default.noz_mat_select.heat_cap,
            default.operationtime,
            nozzle_mass,
            default.eps,
            Tw_ad_noz,
            h_c_noz,
            t_noz,
            np.array([default.default_coating_thickness for i in range(len(t_noz))]),
            prop,
            default.noz_mat_select,
            default.default_coating,
            default.Dr,
            A_nozzle,
            default.T_fuel_tanks,
            d[-1].m_nozz / (1.0 + d[-1].O_F) / default.n,
            x_noz_cool[-1],
            y_noz_cool,
            True,
            default.regenerative_case,
            t_noz_cooling_reg
        )
        Coolobj.Q = Coolobj.Q * default.n

        #errors.append((err_nozzle_cooling))
        #warnings.append((warn_nozzle_cooling))

        #for i in range(len(errors)):
        #    if errors[i]!=0:
        #        return False, errors, warnings;
        #Outputs
        #end temperature of the coolant in K
        #wall temperature array
        #pressure loss
        #_
        #errors
        #warnings

        #Inputs
        #Inicial temperature in K
        # heat capacitance of the material in J/K
        #operation time in seconds
        #mass of the chamber in kg
        #emissivity
        #adiabatic wall temperature of the chamber np.array
        #convection coefficient in W/(m2 K) np.array
        #Thickness of the chamber in m np.array
        #default thickness of the coating to
        #propellant
        #material of the chamber
        #deafult material of the coating
        #Hydralic diamter of the cooling system piping in m 
        #Segment area of contact in m2
        #inicial coolant temperature in K
        #coolant mass flow in kg/s
        #length of the chamber in m
        #Radius of the chamber in m~
        # Option to overwrite area with input area bool
        # which regenerative function to call (default 0, do not change!)
        
        A_chamber = [d[-1].Chamber_L * y_for_cooling_channel * alpha]
        (
            Tf_cool,
            Tw_wall_chamber_calculated,
            dptcool_c,
            _,
            type_variable_chamber,
            T_outer_wall_chamber,
            err_chamber_cooling,
            warn_chamber_cooling,
        ) = Coolobj_c.Run_cooling(
            default.T0,
            default.chamber_mat_select.heat_cap,
            default.operationtime,
            chamber_mass,
            default.eps,
            np.array([d[-1].Tc]),
            np.array([d[-1].h_comb]),
            np.array([d[-1].ThicknessChamber]),
            np.array([default.default_coating_thickness]),
            prop,
            default.chamber_mat_select,
            default.default_coating,
            default.Dr,
            np.array(A_chamber),
            Tf_cool,
            d[-1].m_nozz / (1.0 + d[-1].O_F) / default.n,
            d[-1].Chamber_L,
            np.array([d[-1].Dc / 2]),
            True,
            default.regenerative_case,
            np.array([default.cooling_thickness]),
        )
        Coolobj_c.Q = Coolobj_c.Q * default.n
        #Outputs
        #end temperature of the coolant in K
        #wall temperature array
        #pressure loss
        #_
        #errors
        #warnings

        #errors.append((err_chamber_cooling))
        #warnings.append((warn_chamber_cooling))

        #for i in range(len(errors)):
        #    if errors[i]!=0:
        #        return False, errors, warnings;
        dptcool = dptcool + dptcool_c
        dptcool_cooling = dptcool
        error_out_cool=0
        err_cooling=0
        maximum_thermal_stress,safety_factor_cooling, max_temperature_inner,max_temperature_outer,error_out_cool=Cooling.outputs(T_outer_wall_chamber,T_outer_wall_nozzle,Tw_wall_chamber_calculated,Tw_wall_nozzle_calculated,default.noz_mat_select,default.chamber_mat_select, default.default_coating, d[-1].ThicknessChamber, t_noz, default.default_coating_thickness,y_noz,d[-1].Dc/2 ,0)
        err_cooling=err_cooling|error_out_cool
        err_cooling=err_cooling|err_chamber_cooling
        err_cooling=err_cooling|err_nozzle_cooling
        d[-1].type_variable_nozzle=type_variable_nozzle
        d[-1].type_variable_chamber=type_variable_chamber
        d[-1].T_after_cool=Tf_cool
        d[-1].maximum_thermal_stress=maximum_thermal_stress
        d[-1].safety_factor_cooling=safety_factor_cooling
        d[-1].max_temperature_inner=max_temperature_inner
        d[-1].max_temperature_outer=max_temperature_outer
        if type_variable_nozzle !=2 and type_variable_chamber!=2:
            d[-1].dptcool_cooling=0;
        elif isinstance(dptcool,(list,np.ndarray,tuple)):
            d[-1].dptcool_cooling=max(dptcool)
        else:
            d[-1].dptcool_cooling=dptcool
        if isinstance(Tf_cool,(np.ndarray,list,tuple)):
            Tf_cool=max(Tf_cool)
        #Compute Turbo
        #dp_cool=np.max(dptcool)
        #dp_cool=0.004
        #Tf_cool=400
        
    
        d[-1].ptinj, d[-1].W_Opump, d[-1].W_Fpump, d[-1].W_turb, d[-1].fuel_frac, d[-1].turbo_m, error_t = Turbo.TurboM(default, prop, d[-1].O_F, d[-1].Pa, Tf_cool, d[-1].dptcool_cooling, d[-1].m_nozz)
        
        if error_t!=0:
            return errors_nz1,errors_nz2,er_comb,error_t,er_inj1,0,0,0,0,warnings_nz1,warnings_nz2,wr_comb,0,wr_inj1,0,0,0,0,0;
        
        #Compute Injector (2)
        p_new, dp_ox, dp_f, er_inj2, wr_inj2 = Inj.injector2(default, prop, d[-1].v_iox, d[-1].v_if, d[-1].ptinj, d[-1].eta_s, d[-1].Pa)

        it=it+1
        if it>default.max_it_tool:
            break
        if p_new/d[-1].Pa<4:
            p_new=p_new+4*d[-1].Pa;
        if er_inj2!=0:
            return errors_nz1,errors_nz2,er_comb,error_t,er_inj2,0,0,0,0,warnings_nz1,warnings_nz2,wr_comb,0,wr_inj2,0,0,0,0,0;
        #print("P_new: " + str(p_new))
              
    d[-1].Pc = p_new
    bool = 1 #Shows the combustor it is out of the loop in order to compute mass!
    
    er_inj=er_inj1|er_inj2
    wr_inj=wr_inj1|wr_inj2

    # Plot nozzle contour
    #x_plot1=np.linspace(y_noz[0],y_noz[-1],len(y_noz))
    #y_plot=np.linspace(y_noz[0],y_noz[-1],len(y_noz))
    #x_plot2=np.linspace(y_noz[0],y_noz[-1],len(y_noz))
    
    #for i in range(len(x_noz)):
    #    x_plot1[i]=y_noz[i]
    #    y_plot[i]=x_noz[-1]-x_noz[i]
    #    x_plot2[i]=-x_plot1[i]
    
    #x_h1=np.linspace(x_plot2[0],x_plot1[0],100)
    #y_h1=np.linspace(x_plot2[0],x_plot1[0],100)
    #for i in range(len(x_h1)):
    #    y_h1[i]=y_plot[0]
    
    #x_h2=np.linspace(x_plot2[-1],x_plot1[-1],100)
    #y_h2=np.linspace(x_plot2[-1],x_plot1[-1],100)
    #for i in range(len(x_h2)):
    #    y_h2[i]=y_plot[-1]
    #plt.plot(x_plot1,y_plot,color='k')
    #plt.plot(x_plot2,y_plot,color='k')
    #plt.plot(x_h1,y_h1,color='k')
    #plt.plot(x_h2,y_h2,color='k')
    #plt.xlim(0, x_plot[-1]*1.1)
    #plt.ylim(0,y_plot[-1]*1.1)
    #ax = plt.gca()
    #ax.set_aspect('equal', adjustable='box')
    #plt.show()


    #Computes igniters mass
    #inputs
    ## massflow
    ## propellant
    ## default
    ## chamber temperature
    ## oxidizer to fuel ratio of main chamber
    Ign_propellant_mass, Ign_fuel_mass, Ign_ox_mass,wr_ign = Ign.Igniters(d[-1].m_nozz,prop,default,d[-1].Tc,d[-1].O_F,default.type)
    er_ign=0
    #Ign_propellant_mass=100
    #outputs
    ## mass used for igniter
    ## By this order: igniter propellant total mass, igniter fuel mass, igniter oxidizer mass, warnings

    #Compute reliability 
    Reliability, wr_rel = Rel.Reliability(default, d[-1].time, d[-1].Thrust, d[-1].Thrust, default.val)
    #d[-1].rel=max(Reliability)*100
    d[-1].rel=Reliability*100
    #Compute Mass:
    ChamberMass = Comb.CombustionChamber(p_new, d[-1].At, prop, default.chamber_mat_select, default, d[-1].v_if, d[-1].v_iox, d[-1].D_f, d[-1].D_ox, bool,rho_c,cp_c,mu_c/10,k_c,Pr_c,A_est,d[-1].O_F)[5]
    ##ChamberMass= Ms.Chamber_mass(d[-1].Dc,d[-1].Chamber_L,d[-1].ThicknessChamber,Ms.Rhenium)
    #ChamberMass = Comb.CombustionChamber(p_new, d[-1].At, prop, Ms.Rhenium, default, d[-1].v_if, d[-1].v_iox, d[-1].D_f, d[-1].D_ox, bool,rho_c,cp_c,mu_c/10,k_c,Pr_c,A_est,d[-1].O_F)[5]
    #Inputs for Ms.Mass
    ## chamber pressure
    ##nozzle material from aux
    ##valve material from aux
    ##Area ratio
    ##Throat Area
    ##Mass flow rate
    ##Safety Factor
    ##Cycle type (has to be changed in aux to affect Mass function)
    ##nozzle x coordinate
    ##nozzle y coordinate
    ##nozzle thickness
    ##Density of the oxidizer
    ##Density of the Fuel
    ##Oxidizer to Fuel Ratio
    
    Mass,t1,t2= ChamberMass + Ign_propellant_mass + Ms.Mass(p_new,aux.Default.noz_mat_select,aux.Default.valv_mat_select,d[-1].Eps,d[-1].At, d[-1].m_nozz,aux.Default.Safety_factor,aux.Default.cycle_type,x_noz,y_noz,t_noz,prop.o_dens,prop.f_dens_l,d[-1].O_F)
    M,error_mass,warning_mass=Ms.Mass(p_new,aux.Default.noz_mat_select,aux.Default.valv_mat_select,d[-1].Eps,d[-1].At, d[-1].m_nozz,aux.Default.Safety_factor,aux.Default.cycle_type,x_noz,y_noz,t_noz,prop.o_dens,prop.f_dens_l,d[-1].O_F)
    d[-1].Mtot=Mass
    d[-1].chambermass=ChamberMass
    d[-1].Igniter_compound=Ign_propellant_mass
    d[-1].Mnoz=Ms.Nozzle_mass(x_noz,y_noz,t_noz,aux.Default.noz_mat_select)

    if error_mass!=0:
            return errors_nz1,errors_nz2,er_comb,error_t,er_inj,er_ign,0,error_mass,0,warnings_nz1,warnings_nz2,wr_comb,0,wr_inj,wr_ign,0,warning_mass,0,0;
    #Output:
    ##Total Mass of the engine
    
    # Total propellant mass

    d[-1].Mprop=d[-1].turbo_m*d[-1].time
    #Compute Cost:
    #Inputs for Ms.Cost
    ##Mass of the engine
    ##Tech Readiness Factor
    ##Prior Experience Factor
    ##Reliability
    ##Learning Factor
    ##Cycle type
    Cost,error_cost,warning_cost= Ms.Cost(Mass,aux.Default.tech_ready,aux.Default.exp_factor, Reliability,aux.Default.learn_factor,aux.Default.cycle_type)
    if error_cost!=0:
            return errors_nz1,errors_nz2,er_comb,error_t,er_inj,er_ign,0,error_mass,error_cost,warnings_nz1,warnings_nz2,wr_comb,0,wr_inj,wr_ign,0,warning_mass,warning_cost,0;
    d[-1].cost=Cost 
    #Outputs:
    ##Cost of the engine
    
    ligament = aux.Default.cooling_thickness*(0.5)
    P_n = default.perimeter_percentage/default.n
    l = 2*np.pi*min(y_noz)*P_n
    w = 2*np.pi*min(y_noz)*(1-P_n)
    Life, Life_error, Life_warning = Ms.Life(aux.Default.noz_mat_select, maximum_thermal_stress, max_temperature_inner, max_temperature_outer, ligament, l, w, p_new )
    d[-1].Life=Life
    error_cost=error_cost|Life_error
    warning_cost=warning_cost|Life_warning
    #Reuseability Funtion:
    #inputs:
    ##Nozzle Material
    ##Maximum thermal stress
    ##Maximum temperature on the gas side 
    ##Maximum temperature for the coolant side
    ##Ligament thickness
    ##coolant channel width
    ##Ribs width
    ##Chamber Pressure
    Reuseability, Reuseability_error, Reuseability_warning = Ms.Reuseability(aux.Default.Reuses,aux.Data.time) 
    d[-1].Reusability=Reuseability
    #Missing Inputs: H,l,w (geometry of the cooling chanels written as 0.445, in the function above^)
    #They are currently inputs for the cooling function and are not returned anywhere. This is necessary before implementing into the fuction. Currently default values are used 
    #Output:
    #Low cycle fatigue life of the thrust chamber
   

    print("Calculations finished")
    return errors_nz1,errors_nz2,er_comb,error_t,er_inj,er_ign,0,error_mass,error_cost,warnings_nz1,warnings_nz2,wr_comb,0,wr_inj,wr_ign,0,warning_mass,warning_cost, wr_rel

if __name__ == '__main__':
    Main([aux.Data(15000,10,1000)], aux.Default(0), aux.Propellant(0))



