import os
import sys
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtWidgets import QDialog, QApplication, QWidget, QMainWindow, QLabel, QDialogButtonBox, QScrollArea, QVBoxLayout, QSizePolicy, QSpacerItem
from PyQt5.QtWidgets import QMessageBox
import main
import Aux_classes as aux
import csv
from datetime import datetime

set_images_path = "../Turbomachinery_code/"
class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        loadUi("GUI.ui",self)

        # Initialize Font
        screen = QApplication.primaryScreen()
        self.ScrSize = screen.size()
        fontSize = int(self.ScrSize.width()/100)
        self.setStyleSheet(''' font-size: ''' + str(fontSize) + '''px; ''')
        self.ensurePolished()
        self.label_inputs.setStyleSheet(''' font-size: ''' + str(fontSize+8) + '''px; ''')
        self.label_outputs.setStyleSheet(''' font-size: ''' + str(fontSize+8) + '''px; ''')

        # Initialize graphics


        #tabs
        self.tab_inputs.setCurrentIndex(0)
        self.tab_outputs.setCurrentIndex(0)

        #Elements
        self.line_O_F_2.hide()
        self.label_O_F_2.hide()
        self.label_O_F_no.hide()

        # Connects
        #Main
        self.OptOF(0)
        self.line_thrust.editingFinished.connect(self.checkThrust);  self.checkThrust()
        self.line_time.editingFinished.connect(self.checktime); self.checktime()
        self.line_pa.editingFinished.connect(self.checkPa); self.checkPa()
        self.line_O_F.editingFinished.connect(self.checkOF); self.checkOF()
        self.line_Otank_pres.editingFinished.connect(self.checkOtankPres); self.checkOtankPres()
        self.line_Ftank_pres.editingFinished.connect(self.checkFtankPres); self.checkFtankPres()

        #Nozzle
        self.nozzle_changed(0)
        self.line_De_max.editingFinished.connect(self.checkDeMax);  self.checkDeMax()
        self.line_eps_max.editingFinished.connect(self.checkEpsMax);  self.checkEpsMax()
        self.line_conv_ang.editingFinished.connect(self.checkConvAng);  self.checkConvAng()
        self.line_div_ang.editingFinished.connect(self.checkDivAng);  self.checkDivAng()
        self.line_par_0ang.editingFinished.connect(self.checkParAng0);  self.checkParAng0()
        self.line_par_fang.editingFinished.connect(self.checkParAngf);  self.checkParAngf()
        self.line_curv_conv.editingFinished.connect(self.checkCurvConv);  self.checkCurvConv()
        self.line_curv_div.editingFinished.connect(self.checkCurvDiv);  self.checkCurvDiv()

        #Combustion chamber
        self.line_SF.editingFinished.connect(self.checkSF);  self.checkSF()
        self.line_drop_ratio.editingFinished.connect(self.checkDropRat); self.checkDropRat()
        self.line_kloads.editingFinished.connect(self.checkKloads); self.checkKloads()

        #Turbo
        self.cycle_changed(0)
        self.line_Opump_eff.editingFinished.connect(self.checkOpumpEff);  self.checkOpumpEff()
        self.line_Fpump_eff.editingFinished.connect(self.checkFpumpEff);  self.checkFpumpEff()
        self.line_turb_eff.editingFinished.connect(self.checkTurbEff);  self.checkTurbEff()
        self.line_mech_eff.editingFinished.connect(self.checkMechEff);  self.checkMechEff()
        self.line_Ele_power.editingFinished.connect(self.checkElePow);  self.checkElePow()
        self.line_valve_losses.editingFinished.connect(self.checkValveLoss);  self.checkValveLoss()
        self.line_lGG.editingFinished.connect(self.checkLGG); self.checkLGG()

        #Cooling
        self.coating_changed(0)
        self.line_OxTank_temp.editingFinished.connect(self.checkOTankTemp);  self.checkOTankTemp()
        self.line_FuelTank_temp.editingFinished.connect(self.checkFTankTemp);  self.checkFTankTemp()
        self.line_T0_wall.editingFinished.connect(self.checkT0wall);  self.checkT0wall()
        self.line_coating_thick.editingFinished.connect(self.checkCoatThick);  self.checkCoatThick()
        self.line_DH.editingFinished.connect(self.checkDH);  self.checkDH()

        #Injectors
        self.Injector_changed(0); self.define_inj_dp(0)
        self.line_Cd.editingFinished.connect(self.checkCd);  self.checkCd()
        self.line_inj_dp.editingFinished.connect(self.checkInjDp); self.checkInjDp()

        #Ignitors
        self.igniter_changed(0)
        self.line_tign.editingFinished.connect(self.checktIgn);  self.checktIgn()
        self.line_O_F_ign.editingFinished.connect(self.checkOFIgn);  self.checkOFIgn()
        self.line_corr_ign.editingFinished.connect(self.checkCorrIgn);  self.checkCorrIgn()

        #Materials
        self.line_dens_nozz.editingFinished.connect(self.checkDensNozz); self.checkDensNozz()
        self.line_yield_nozz.editingFinished.connect(self.checkYieldNozz); self.checkYieldNozz()
        self.line_OpT_nozz.editingFinished.connect(self.checkOpTNozz); self.checkOpTNozz()
        self.line_E_nozz.editingFinished.connect(self.checkENozz); self.checkENozz()
        self.line_k_nozz.editingFinished.connect(self.checkKNozz); self.checkKNozz()
        self.line_cost_nozz.editingFinished.connect(self.checkCostNozz); self.checkCostNozz()
        self.line_dens_cham.editingFinished.connect(self.checkDensCham); self.checkDensCham()
        self.line_yield_cham.editingFinished.connect(self.checkYieldCham); self.checkYieldCham()
        self.line_OpT_cham.editingFinished.connect(self.checkOpTCham); self.checkOpTCham()
        self.line_E_cham.editingFinished.connect(self.checkECham); self.checkECham()
        self.line_k_cham.editingFinished.connect(self.checkKCham); self.checkKCham()
        self.line_cost_cham.editingFinished.connect(self.checkCostCham); self.checkCostCham()

    #Run
    def Run(self):
        #get propellant, cycle, nozzle, materials, etc...
        main.prop = aux.Propellant(self.combo_prop.currentIndex())
        main.default.cycle_type = self.combo_cycle.currentIndex()

        #run main
        err_nozz, err_chamber, err_turbo, err_inj, err_ign, err_cool, err_mass, warn_nozz, warn_chamber, warn_turbo, warn_inj, warn_ign, warn_cool, warn_mass = main.Main(main.dat)
        
        #Errors
        if(err_nozz or err_chamber or err_turbo or err_inj or err_ign or err_cool or err_mass): 
            msg = QMessageBox(self)
            msg.setIcon(QMessageBox.Critical)
            msg.setWindowTitle("Calculation error!")

            #Nozzle
            if (err_nozz & (1<<0)):
                msg.setText("Error in nozzle calculation - Chamber pressure is lower than ambient pressure")
            if (err_nozz & (1<<1)):
                msg.setText("Error in nozzle calculation - Mass flow rate is negative")
            if (err_nozz & (1<<2)):
                msg.setText("Error in nozzle calculation - Chamber temperature is negative")
            if (err_nozz & (1<<3)):
                msg.setText("Error in nozzle calculation - Throat area is negative")
            if (err_nozz & (1<<4)):
                msg.setText("Error in nozzle calculation - Expansion ratio is lower than 1")
            if (err_nozz & (1<<5)):
                msg.setText("Error in nozzle calculation - Isp is negative")
            if (err_nozz & (1<<6)):
                msg.setText("Error in nozzle calculation - Combustion chamber properties are negative")
            if (err_nozz & (1<<7)):
                msg.setText("Error in nozzle calculation - Divergent flow loss is higher than 1")
            if (err_nozz & (1<<8)):
                msg.setText("Error in nozzle calculation - Divergent nozzle length is negative")
            if (err_nozz & (1<<9)):
                msg.setText("Error in nozzle calculation - Mass flow rate in iteration becomes negative")
            if (err_nozz & (1<<10)):
                msg.setText("Error in nozzle calculation - Throat area in iteration is negative")
            
            #Combustion chamber
            if (err_chamber & (1<<0)):
                msg.setText("Error in combustion chamber calculation - Chamber length is higher than 2m")

            #Turbomachinery
            if (err_turbo & (1<<0)):
                msg.setText("Error in cycle calculation - Cycle not yet supported")
            if (err_turbo & (1<<1)):
                msg.setText("Error in cycle calculation - Optimization did not converge")
            if (err_turbo & (1<<2)):
                msg.setText("Error in cycle calculation - System could not be solved")
            if (err_turbo & (1<<3)):
                msg.setText("Error in cycle calculation - Potency equilibrium between pumps and turbine not satisfied")
            if (err_turbo & (1<<4)):
                msg.setText("Error in cycle calculation - Pump or turbine potency is negative")

            #Injector


            #Igniter


            #Cooling
            if (err_cool & (1<<0)):
                msg.setText("Error in cooling calculation - One of the inputs is negative")
            if (err_cool & (1<<1)):
                msg.setText("Error in cooling calculation - The temperature calculated for the heatsink is negative")
            if (err_cool & (1<<2)):
                msg.setText("Error in cooling calculation - The temperature calculated for radiation cooling is negative")
            if (err_cool & (1<<3)):
                msg.setText("Error in cooling calculation - One of the numerical outputs (coolant temperature, wall temperature or pressure loss) of the Regenerative Cooling explicit function Run1D() is negative")
            if (err_cool & (1<<4)):
                msg.setText("Error in cooling calculation - One of the numerical outputs (coolant temperature, wall temperature or pressure loss) of the Regenerative Cooling iterative function Run_for_Toperating1D is negative")
            if (err_cool & (1<<5)):
                msg.setText("Error in cooling calculation - One of the numerical outputs (coolant temperature, wall temperature or pressure loss) of the Regenerative Cooling iterative function Run1D_iterative_for_m is negative")
            if (err_cool & (1<<6)):
                msg.setText("Error in cooling calculation - The end coolant temperature is larger than the end wall temperature")
            if (err_cool & (1<<7)):
                msg.setText("Error in cooling calculation - End coolant temperature is negative or bigger than 1000 K")
            if (err_cool & (1<<8)):
                msg.setText("Error in cooling calculation - Wall temperature is negative or larger than operating temperature")
            if (err_cool & (1<<9)):
                msg.setText("Error in cooling calculation - Pressure loss is negative or larger 1 bar")
            if (err_cool & (1<<10)):
                msg.setText("Error in cooling calculation - mass flow is negative or larger 30 kg/s")
            if (err_cool & (1<<11)):
                msg.setText("Error in cooling calculation - Heat is negative")

            #Mass

            

            msg.exec_()
            return;

        #warnings
        if(warn_nozz or warn_chamber or warn_turbo or warn_inj or warn_ign or warn_cool or warn_mass):
            msg = QMessageBox(self)
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle("Warnings raised")
            msg.setText("Warnings raised during calculations:")
            scroll = QScrollArea(msg)
            scroll.setWidgetResizable(True)
            msg.content = QWidget()
            scroll.setWidget(msg.content)
            lay = QVBoxLayout(msg.content)
            msg.layout().addWidget(scroll, 1, 0, 1, 3)
            msg.layout().setColumnMinimumWidth(2,int(self.ScrSize.width()/3))
            msg.layout().setRowMinimumHeight(1,int(self.ScrSize.height()/4))

            #Nozzle
            if (warn_nozz & (1<<0)):
                lay.addWidget(QLabel("Warning in nozzle calculation - Chamber pressure is higher than 300 bar.", parent = msg))
            if (warn_nozz & (1<<1)):
                lay.addWidget(QLabel("Warning in nozzle calculation - Divergent angle is high.", parent = msg))
            if (warn_nozz & (1<<2)):
                lay.addWidget(QLabel("Warning in nozzle calculation - Divergent angle in bell nozzle is high.", parent = msg))
            if (warn_nozz & (1<<3)):
                lay.addWidget(QLabel("Warning in nozzle calculation - Exit angle in bell nozzle is high.", parent = msg))
            if (warn_nozz & (1<<4)):
                lay.addWidget(QLabel("Warning in nozzle calculation - Isp is higher than 550 s.", parent = msg))
            if (warn_nozz & (1<<5)):
                lay.addWidget(QLabel("Warning in nozzle calculation - Divergent flow loss is very high.", parent = msg))

            #Chamber
            if (warn_chamber & (1<<0)):
                lay.addWidget(QLabel("Warning in combustion chamber calculation - Droplet diameter at injection is smaller than 100 micrometers (this should not happen), chamber length is probably too small", parent = msg))
            if (warn_chamber & (1<<1)):
                lay.addWidget(QLabel("Warning in combustion chamber calculation - Droplet diameter at injection is higher than 250 micrometers, chamber will probably be too large", parent = msg))
            if (warn_chamber & (1<<2)):
                lay.addWidget(QLabel("Warning in combustion chamber calculation - Velocity of oxidizer at injection is smaller than 10 m/s, this is too small", parent = msg))
            if (warn_chamber & (1<<3)):
                lay.addWidget(QLabel("Warning in combustion chamber calculation - Velocity of oxidizer at injection is higher than 30 m/s, this is too high", parent = msg))
            if (warn_chamber & (1<<4)):
                lay.addWidget(QLabel("Warning in combustion chamber calculation - Velocity of fuel at injection is smaller than 10 m/s, this is too small.", parent = msg))
            if (warn_chamber & (1<<5)):
                lay.addWidget(QLabel("Warning in combustion chamber calculation - Velocity of fuel at injection is higher than 300 m/s, this is too high.", parent = msg))
            if (warn_chamber & (1<<6)):
                lay.addWidget(QLabel("Warning in combustion chamber calculation - Chamber length was increased to a value higher than what it takes to get full atomization of propellant, in order to make CC chamber fit the convergence ratio", parent = msg))
            if (warn_chamber & (1<<7)):
                lay.addWidget(QLabel("Warning in combustion chamber calculation - Chamber diameter was too small, all methods failed as such, Lstar was increased to make CC fulfill convergence ratio", parent = msg))
            if (warn_chamber & (1<<5)):
                lay.addWidget(QLabel("Warning in combustion chamber calculation - Chamber length is larger than 1m.", parent = msg))

            #Turbo

            #Injector

            #igniter
            if (warn_ign & (1<<0)):
                lay.addWidget(QLabel("Warning in igniter calculation - Chosen igniter type is not suitable for the chosen motor, or the correction factor is to small", parent = msg))
            if (warn_ign & (1<<1)):
                lay.addWidget(QLabel("Warning in igniter calculation - Mass is overestimated, the correction factor is to small", parent = msg))

            #Cooling

            #Mass

            lay.addItem(QSpacerItem(10,10,QSizePolicy.Expanding,QSizePolicy.Expanding))
            msg.exec_()

        #Out
        runNum = len(main.dat)
        self.combo_run.addItem("Run " + str(runNum))
        self.combo_run.setCurrentIndex(runNum-1)
        self.Output(runNum-1)

        #New data element
        main.dat.append(aux.Data(0,0,0))
        self.checkThrust()
        self.checktime()
        self.checkPa()
        self.OptOF(self.check_O_F.isChecked())
        

    # Output
    def Output(self, i : int):
        self.line_Isp.setText(str(main.dat[i].Isp))
        self.line_cstar.setText(str(main.dat[i].cstar))
        self.line_m.setText(str(main.dat[i].m))
        self.line_prop_mass.setText(str(main.dat[i].Mprop))
        self.line_mass.setText(str(main.dat[i].Mtot))
        self.line_cost.setText(str(main.dat[i].cost))
        self.line_reliability.setText(str(main.dat[i].rel))

        #Nozzle
        if(not main.dat[i].O_F): 
            self.line_O_F_2.setText(str(main.dat[i].O_F))
            self.line_O_F_2.show()
            self.label_O_F_2.show()
            self.label_O_F_no.show()
        else:
            self.line_O_F_2.hide()
            self.label_O_F_2.hide()
            self.label_O_F_no.hide()
        self.line_Dt.setText(str(main.dat[i].Dt))
        self.line_De.setText(str(main.dat[i].De))
        self.line_eps.setText(str(main.dat[i].Eps))
        self.line_len.setText(str(main.dat[i].L_total))
        self.line_len_conv.setText(str(main.dat[i].L_con))
        self.line_len_div.setText(str(main.dat[i].L_div))
        self.line_m_nozz.setText(str(main.dat[i].m_nozz))
        self.line_nozz_mass.setText(str(main.dat[i].Mnoz))

        #Combustion chamber
        self.line_Pc.setText(str(main.dat[i].Pc/1.0e5))
        self.line_Tc.setText(str(main.dat[i].Tc))
        self.line_Dc.setText(str(main.dat[i].Dc))
        self.line_Lc.setText(str(main.dat[i].Chamber_L))
        self.line_tc.setText(str(main.dat[i].ThicknessChamber))
        self.line_Mc.setText(str(main.dat[i].chambermass))

        #turbo
        self.line_Opump_power.setText(str(main.dat[i].W_Opump))
        self.line_Fpump_power.setText(str(main.dat[i].W_Fpump))
        self.line_Turb_power.setText(str(main.dat[i].W_turb))
        self.line_P_inj_in.setText(str(main.dat[i].ptinj))
        self.line_turbo_mass.setText(str(main.dat[i].turbo_m))

        #Cooling
        self.combo_cool_cham.setCurrentIndex(main.dat[i].type_variable_chamber)
        self.combo_cool_nozz.setCurrentIndex(main.dat[i].type_variable_nozzle)

        #Injectors
        self.line_n_oxinj.setText(str(main.dat[i].n_ox))
        self.line_n_finj.setText(str(main.dat[i].n_f))
        self.line_vel_ox.setText(str(main.dat[i].v_iox))
        self.line_vel_fuel.setText(str(main.dat[i].v_if))
        self.line_inj_pres.setText(str(main.dat[i].dp))

        #Ignitors
        self.line_Mign.setText(str(main.dat[i].Igniter_compound))

        #Materials


    #Export CSV
    def exportCSV(self):
        filename = "Export_" + datetime.now().strftime("%H_%M_%S") + ".csv"
        for data in main.dat:
            members = [attr for attr in dir(data) if not callable(getattr(data, attr)) and not attr.startswith("__")]
            with open(filename, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(members)
                writer.writerow([getattr(data, attr) for attr in members])


    #Blocks
    def OptOF(self, ischeck : bool):
        if ischeck:
            main.dat[-1].O_F = 0
            main.default.MR = 0
            self.line_O_F.setEnabled(False)
        else:
            self.checkOF()
            self.line_O_F.setEnabled(True)

    def define_inj_dp(self, ischeck : bool):
        main.default.dp_state = ischeck
        if ischeck:
            main.default.dp_user = float(self.line_inj_dp.text())
            self.line_inj_dp.setEnabled(True)
        else:
            self.line_inj_dp.setDisabled(True)

    def cycle_changed(self,i : int):
        if i == 2:
            self.checkLGG()
            self.label_lGG.show()
            self.line_lGG.show()
            self.label_no_lGG.show()
        else:
            self.label_lGG.hide()
            self.line_lGG.hide()
            self.label_no_lGG.hide()

        if i == 4:
            self.line_Ele_power.setEnabled(True)
            self.line_turb_eff.setEnabled(False)
        else:
            self.line_Ele_power.setEnabled(False)
            if i == 5:
                self.line_Opump_eff.setEnabled(False)
                self.line_Fpump_eff.setEnabled(False)
                self.line_turb_eff.setEnabled(False)
                self.line_mech_eff.setEnabled(False)
            else:
                self.line_Opump_eff.setEnabled(True)
                self.line_Fpump_eff.setEnabled(True)
                self.line_turb_eff.setEnabled(True)
                self.line_mech_eff.setEnabled(True)

    def nozzle_changed(self,i : int):
        if i==0:
            main.default.Nozzle_type = 0
            self.line_par_0ang.setEnabled(False)
            self.line_par_fang.setEnabled(False)
            self.line_div_ang.setEnabled(True)
        else:
            main.default.Nozzle_type = 1
            self.line_par_0ang.setEnabled(True)
            self.line_par_fang.setEnabled(True)
            self.line_div_ang.setEnabled(False)

    def Injector_changed(self, i : int):
        main.default.InjType = main.default.InjTypes[i]

    def igniter_changed(self, i : int):
        if(i==0):
            self.torch_changed(self.combo_torch.currentIndex())
            self.combo_hypergolic.hide()
            self.combo_torch.show()
        elif(i==2):
            self.hypergolic_changed(self.combo_hypergolic.currentIndex())
            self.combo_hypergolic.show()
            self.combo_torch.hide()
        else:
            main.default.type = "10"
            self.combo_hypergolic.hide()
            self.combo_torch.hide()

    def torch_changed(self, i : int):
        main.default.type = "0" + str(i)

    def hypergolic_changed(self, i : int):
        main.default.type = "2" + str(i)

    def nozz_mat_changed(self, i: int):
        main.default.noz_mat_select = main.default.material_list[i]
        if(i):
            self.label_dens_nozz.hide()
            self.label_yield_nozz.hide()
            self.label_E_nozz.hide()
            self.label_OpT_nozz.hide()
            self.label_k_nozz.hide()
            self.label_cost_nozz.hide()
            self.line_dens_nozz.hide()
            self.line_yield_nozz.hide()
            self.line_E_nozz.hide()
            self.line_OpT_nozz.hide()
            self.line_k_nozz.hide()
            self.line_cost_nozz.hide()
            self.label_kg_m3_nozz.hide()
            self.label_Pa_nozz.hide()
            self.label_Pa_nozz_2.hide()
            self.label_K_nozz.hide()
            self.label_W_mK_nozz.hide()
            self.label_eur_kg_nozz.hide()
        else:
            self.label_dens_nozz.show()
            self.label_yield_nozz.show()
            self.label_E_nozz.show()
            self.label_OpT_nozz.show()
            self.label_k_nozz.show()
            self.label_cost_nozz.show()
            self.line_dens_nozz.show()
            self.line_yield_nozz.show()
            self.line_E_nozz.show()
            self.line_OpT_nozz.show()
            self.line_k_nozz.show()
            self.line_cost_nozz.show()
            self.label_kg_m3_nozz.show()
            self.label_Pa_nozz.show()
            self.label_Pa_nozz_2.show()
            self.label_K_nozz.show()
            self.label_W_mK_nozz.show()
            self.label_eur_kg_nozz.show()

    def cham_mat_changed(self, i : int):
        main.default.noz_mat_select = main.default.material_list[i]
        if(i):
            self.label_dens_cham.hide()
            self.label_yield_cham.hide()
            self.label_E_cham.hide()
            self.label_OpT_cham.hide()
            self.label_k_cham.hide()
            self.label_cost_cham.hide()
            self.line_dens_cham.hide()
            self.line_yield_cham.hide()
            self.line_E_cham.hide()
            self.line_OpT_cham.hide()
            self.line_k_cham.hide()
            self.line_cost_cham.hide()
            self.label_kg_m3_cham.hide()
            self.label_Pa_cham.hide()
            self.label_Pa_cham_2.hide()
            self.label_K_cham.hide()
            self.label_W_mK_cham.hide()
            self.label_eur_kg_cham.hide()
        else:
            main.default.noz_mat_select = main.default.coating_list[0]
            self.label_dens_cham.show()
            self.label_yield_cham.show()
            self.label_E_cham.show()
            self.label_OpT_cham.show()
            self.label_k_cham.show()
            self.label_cost_cham.show()
            self.line_dens_cham.show()
            self.line_yield_cham.show()
            self.line_E_cham.show()
            self.line_OpT_cham.show()
            self.line_k_cham.show()
            self.line_cost_cham.show()
            self.label_kg_m3_cham.show()
            self.label_Pa_cham.show()
            self.label_Pa_cham_2.show()
            self.label_K_cham.show()
            self.label_W_mK_cham.show()
            self.label_eur_kg_cham.show()

    def coating_changed(self, i : int):
        main.default.default_coating = main.default.coating_list[i+1]

    #Checks
    def checkThrust(self):
        thrust = float(self.line_thrust.text())
        if thrust > 0.0 and thrust < 1.0e10:
            main.dat[-1].Thrust = thrust;
        else:
            self.line_thrust.setText("15000.0")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid thrust, try again.")
            msg.exec_()

    def checktime(self):
        var = float(self.line_time.text())
        if var > 0.0 and var < 1000:
            main.dat[-1].time = var;
        else:
            self.line_time.setText("10.0")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid thrust time, try again.")
            msg.exec_()

    def checkPa(self):
        var = float(self.line_pa.text())
        if var > 0.0 and var < 1.0e6:
            main.dat[-1].Pa = var;
        else:
            self.line_pa.setText("101325.0")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid ambient pressure, try again.")
            msg.exec_()

    def checkOF(self):
        var = float(self.line_O_F.text())
        if var > 0.0 and var < 100:
            main.dat[-1].O_F = var;
            main.default.MR = var;
        else:
            self.line_O_F.setText("5.0")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid O/F, try again.")
            msg.exec_()

    def checkOtankPres(self):
        var = float(self.line_Otank_pres.text())
        if var > 0.0 and var < 1.0e10:
            main.default.p_to = var;
        else:
            self.line_Otank_pres.setText("1.0e5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid oxidizer tank pressure, try again.")
            msg.exec_()

    def checkFtankPres(self):
        var = float(self.line_Ftank_pres.text())
        if var > 0.0 and var < 1.0e10:
            main.default.ptf = var;
        else:
            self.line_Ftank_pres.setText("1.0e5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid fuel tank pressure, try again.")
            msg.exec_()

    def checkOpumpEff(self):
        var = float(self.line_Opump_eff.text())
        if var > 0.0 and var <= 1.0:
            main.default.Eff_po = var;
        else:
            self.line_Opump_eff.setText("0.3")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid oxidizer pump efficiency, try again.")
            msg.exec_()

    def checkFpumpEff(self):
        var = float(self.line_Fpump_eff.text())
        if var > 0.0 and var <= 1.0:
            main.default.Eff_pf = var;
        else:
            self.line_Fpump_eff.setText("0.3")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid fuel pump efficiency, try again.")
            msg.exec_()

    def checkTurbEff(self):
        var = float(self.line_turb_eff.text())
        if var > 0.0 and var <= 1.0:
            main.default.Eff_t = var;
        else:
            self.line_turb_eff.setText("0.5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid turbine efficiency, try again.")
            msg.exec_()

    def checkMechEff(self):
        var = float(self.line_mech_eff.text())
        if var > 0.0 and var <= 1.0:
            main.default.Eff_m = var;
        else:
            self.line_mech_eff.setText("0.95")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid mechanical efficiency, try again.")
            msg.exec_()

    def checkElePow(self):
        var = float(self.line_Ele_power.text())
        if var > 0.0 and var <= 1.0e10:
            main.default.Wmotor = var;
        else:
            self.line_Ele_power.setText("1.0e6")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid electric motor power, try again.")
            msg.exec_()

    def checkValveLoss(self):
        var = float(self.line_valve_losses.text())
        if var > 0.0 and var <= 1.0e6:
            main.default.v_loss = var;
        else:
            self.line_valve_losses.setText("1000")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid valve losses, try again.")
            msg.exec_()

    def checkDeMax(self):
        var = float(self.line_De_max.text())
        if var > 0.0 and var <= 100:
            main.default.De_max = var;
        else:
            self.line_De_max.setText("5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid maximum exit diameter, try again.")
            msg.exec_()

    def checkEpsMax(self):
        var = float(self.line_eps_max.text())
        if var > 1.1 and var <= 1000:
            main.default.Eps_max = var;
        else:
            self.line_eps_max.setText("300")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid maximum expansion ratio, try again.")
            msg.exec_()

    def checkConvAng(self):
        var = float(self.line_conv_ang.text())
        if var > 0.0 and var < 90:
            main.default.Theta_con = var;
        else:
            self.line_conv_ang.setText("60")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid convergent angle, try again.")
            msg.exec_()

    def checkDivAng(self):
        var = float(self.line_div_ang.text())
        if var > 0.0 and var < 90:
            main.default.Theta_conical = var;
        else:
            self.line_div_ang.setText("15")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid divergent angle, try again.")
            msg.exec_()

    def checkParAng0(self):
        var = float(self.line_par_0ang.text())
        if var > float(self.line_par_fang.text()) and var < 90:
            main.default.Theta_bell = var;
        else:
            self.line_par_0ang.setText("55")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid parabola initial angle, try again.")
            msg.exec_()

    def checkParAngf(self):
        var = float(self.line_par_fang.text())
        if var > 0.0 and var < float(self.line_par_0ang.text()):
            main.default.TH_exit_bell = var;
        else:
            self.line_par_fang.setText("55")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid parabola final angle, try again.")
            msg.exec_()

    def checkCurvConv(self):
        var = float(self.line_curv_conv.text())
        if var > 0.0 and var < 4.0:
            main.default.R_u_ratio = var;
        else:
            self.line_curv_conv.setText("1.5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid throat curvature (convergent), try again.")
            msg.exec_()

    def checkCurvDiv(self):
        var = float(self.line_curv_div.text())
        if var > 0.0 and var < 2.0:
            main.default.R_u_bell = var;
        else:
            self.line_curv_div.setText("1.5")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid throat curvature (divergent), try again.")
            msg.exec_()

    def checkSF(self):
        var = float(self.line_SF.text())
        if var > 1.0 and var < 20.0:
            main.default.SF = var;
        else:
            self.line_SF.setText("1.2")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid safety factor, try again.")
            msg.exec_()

    def checkCd(self):
        var = float(self.line_Cd.text())
        if var > 0.0 and var < 2.0:
            main.default.Cd = var;
        else:
            self.line_Cd.setText("0.7")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid discharge coefficient, try again.")
            msg.exec_()

    def checktIgn(self):
        var = float(self.line_tign.text())
        if var > 0.0 and var < 6.0:
            main.default.ignburntime = var;
        else:
            self.line_tign.setText("4")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid igniter burn time, try again.")
            msg.exec_()

    def checkOFIgn(self):
        var = float(self.line_O_F_ign.text())
        if var > 0.0 and var < 6.0:
            main.default.ign_o_f = var;
        else:
            self.line_O_F_ign.setText("0.7")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid O/F of the igniter, try again.")
            msg.exec_()

    def checkCorrIgn(self):
        var = float(self.line_corr_ign.text())
        if var > 1.0 and var < 30.0:
            main.default.fudgefactor = var;
        else:
            self.line_corr_ign.setText("20")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid igniter correction factor, try again.")
            msg.exec_()

    def checkDropRat(self):
        var = float(self.line_drop_ratio.text())
        if var > 0.0 and var < 0.3:
            main.default.factor = var;
        else:
            self.line_drop_ratio.setText("0.2")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid final to initial droplet volume ratio, try again.")
            msg.exec_()

    def checkKloads(self):
        var = float(self.line_kloads.text())
        if var > 0.0:
            main.default.kloads = var;
        else:
            self.line_kloads.setText("1")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid correction factor for mass, try again.")
            msg.exec_()

    def checkOTankTemp(self):
        var = float(self.line_OxTank_temp.text())
        if var > 0.0 and var < 2000:
            main.default.T_ox_tanks = var;
        else:
            self.line_OxTank_temp.setText("60")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid temperature for oxidizer tank, try again.")
            msg.exec_()

    def checkFTankTemp(self):
        var = float(self.line_FuelTank_temp.text())
        if var > 0.0 and var < 2000:
            main.default.T_fuel_tanks = var;
        else:
            self.line_FuelTank_temp.setText("20")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid temperature for fuel tank, try again.")
            msg.exec_()

    def checkDensNozz(self):
        var = float(self.line_dens_nozz.text())
        if var > 0.0:
            main.default.material_list[0].density = var;
        else:
            self.line_dens_nozz.setText("21000")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid density for nozzle custom material, try again.")
            msg.exec_()

    def checkYieldNozz(self):
        var = float(self.line_yield_nozz.text())
        if var > 0.0:
            main.default.material_list[0].yieldstress_l = var;
        else:
            self.line_yield_nozz.setText("2300e6")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid ultimate yield stress for nozzle custom material, try again.")
            msg.exec_()

    def checkENozz(self):
        var = float(self.line_E_nozz.text())
        if var > 0.0:
            main.default.material_list[0].Emod = var;
        else:
            self.line_E_nozz.setText("471e9")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid elastic module for nozzle custom material, try again.")
            msg.exec_()

    def checkOpTNozz(self):
        var = float(self.line_OpT_nozz.text())
        if var > 0.0:
            main.default.material_list[0].OpTemp_u = var;
        else:
            self.line_OpT_nozz.setText("2200")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid operating temperature for nozzle custom material, try again.")
            msg.exec_()

    def checkKNozz(self):
        var = float(self.line_k_nozz.text())
        if var > 0.0:
            main.default.material_list[0].k = var;
        else:
            self.line_k_nozz.setText("48")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid conductivity for nozzle custom material, try again.")
            msg.exec_()
    
    def checkCostNozz(self):
        var = float(self.line_cost_nozz.text())
        if var > 0.0:
            main.default.material_list[0].cost = var;
        else:
            self.line_cost_nozz.setText("938")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid cost per kg for nozzle custom material, try again.")
            msg.exec_()

    def checkDensCham(self):
        var = float(self.line_dens_cham.text())
        if var > 0.0:
            main.default.coating_list[0].density = var;
        else:
            self.line_dens_cham.setText("21000")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid density for combustion chamber custom material, try again.")
            msg.exec_()

    def checkYieldCham(self):
        var = float(self.line_yield_cham.text())
        if var > 0.0:
            main.default.coating_list[0].yieldstress_l = var;
        else:
            self.line_yield_cham.setText("2300e6")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid ultimate yield stress for combustion chamber custom material, try again.")
            msg.exec_()

    def checkECham(self):
        var = float(self.line_E_cham.text())
        if var > 0.0:
            main.default.coating_list[0].Emod = var;
        else:
            self.line_E_cham.setText("471e9")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid elastic module for combustion chamber custom material, try again.")
            msg.exec_()

    def checkOpTCham(self):
        var = float(self.line_OpT_cham.text())
        if var > 0.0:
            main.default.coating_list[0].OpTemp_u = var;
        else:
            self.line_OpT_cham.setText("2200")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid operating temperature for combustion chamber custom material, try again.")
            msg.exec_()

    def checkKCham(self):
        var = float(self.line_k_cham.text())
        if var > 0.0:
            main.default.coating_list[0].k = var;
        else:
            self.line_k_cham.setText("48")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid conductivity for combustion chamber custom material, try again.")
            msg.exec_()
    
    def checkCostCham(self):
        var = float(self.line_cost_cham.text())
        if var > 0.0:
            main.default.coating_list[0].cost = var;
        else:
            self.line_cost_cham.setText("938")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid cost per kg for combustion chamber custom material, try again.")
            msg.exec_()

    def checkT0wall(self):
        var = float(self.line_T0_wall.text())
        if var > 0.0 and var < 1000:
            main.default.T0 = var;
        else:
            self.line_T0_wall.setText("298")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid initial wall temperature, try again.")
            msg.exec_()

    def checkCoatThick(self):
        var = float(self.line_coating_thick.text())
        if var >= 0.0 and var < 0.4:
            main.default.default_coating_thickness = var;
        else:
            self.line_coating_thick.setText("0")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid coating thickness, try again.")
            msg.exec_()

    def checkDH(self):
        var = float(self.line_DH.text())
        if var > 0.0 and var < 0.3:
            main.default.Dr = var;
        else:
            self.line_DH.setText("0.01")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid hydraulic diameter, try again.")
            msg.exec_()

    def checkInjDp(self):
        var = float(self.line_inj_dp.text())
        if var > 0.0 and var < 0.6:
            main.default.dp_user = var;
        else:
            self.line_inj_dp.setText("0.2")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid user defined pressure drop over injector, try again.")
            msg.exec_()
    def checkLGG(self):
        var = float(self.line_lGG.text())
        if var > 0.0 and var < 0.9:
            main.default.l_def = var;
        else:
            self.line_lGG.setText("0.1")
            msg = QMessageBox()
            msg.setWindowTitle("Input error!")
            msg.setText("Invalid percentage of total mass flow for GG, try again.")
            msg.exec_()

# Main
if __name__ == '__main__':
    #generate app
    app = QApplication(sys.argv)
    mainwindow = MainWindow()
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(mainwindow)
    widget.setWindowTitle("LRE Design Tool - 0.1")
    #widget.setWindowIcon(QtGui.QIcon("C:\Users\mail\OneDrive\EscritorioPortatil\Universidad\TU_Delft\CSDP\Turbomachinery_code\LRE_icon.jpg"))

    #Show
    widget.show()
    widget.showMaximized()

    #close
    try:
        sys.exit(app.exec_())
    except:
        print("Exiting")
