import os
import sys
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtWidgets import QDialog, QApplication, QWidget
from PyQt5.QtWidgets import QMessageBox
import main

set_images_path = "../Turbomachinery_code/"

class MainWindow(QDialog):
    CY = "EX"

    def __init__(self):
        super(MainWindow, self).__init__()
        loadUi("GUI.ui",self)

        # set initial image to corresponding initial cycle
        self.graphicsView.setStyleSheet("background-image: url("+set_images_path+"Turbo_cycle_types.jpeg);") 
        
        # Compile button
        self.pushButton_Rcompile.clicked.connect(self.compile)

        # Combo Box - Feed cycles
        listocc=["Expander Cycle", "Staged Combustion Cycle", "Coolant Bleed Cycle", "Gas Generator Cycle"]
        for feedcycle in listocc:
            self.comboBox.addItem(feedcycle)
        self.comboBox.currentIndexChanged.connect(self.combochanged)
        self.comboBox.setCurrentIndex(0)

    def combochanged(self): #Feed cycle selection
        self.label_currentStatus.setText("")
        if self.comboBox.currentText() == "Expander Cycle":
            self.graphicsView.setStyleSheet("background-image: url("+ set_images_path +"EX.png);") 
            self.CY = "EX"
        if self.comboBox.currentText() == "Staged Combustion Cycle":
            self.graphicsView.setStyleSheet("background-image: url("+ set_images_path +"SC.png);") 
            self.CY = "SC"
        if self.comboBox.currentText() == "Coolant Bleed Cycle":
            self.graphicsView.setStyleSheet("background-image: url("+ set_images_path +"CB.png);") 
            self.CY = "CB"
        if self.comboBox.currentText() == "Gas Generator Cycle":
            self.graphicsView.setStyleSheet("background-image: url("+ set_images_path +"GG.png);") 
            self.CY = "GG"
        return self.CY
    
    def compile(self): # action after clicking compile button
        msg = QMessageBox()
        msg.setWindowTitle("Input error!")
        x = 1
        if self.doubleSpinBox_T.value() <=0:
            x = 0
            msg.setText("Invalid thrust entry, try again.")
            msg.exec_()
        if self.doubleSpinBox_t.value() <=0:
            x = 0
            msg.setText("Invalid thrust time entry, try again.")
            msg.exec_()
        if self.doubleSpinBox_pa.value() <=0:
            x = 0
            msg.setText("Invalid ambient pressure entry, try again.")
            msg.exec_()
        if self.doubleSpinBox_ps.value() <=0:
            x = 0
            msg.setText("Invalid storage pressure entry, try again.")
            msg.exec_()
        if self.doubleSpinBox_OF.value() <=0:
            x = 0
            msg.setText("Invalid OF entry, try again.")
            msg.exec_()
        if x == 1:
            T = self.doubleSpinBox_T.value()*1000.0
            t = self.doubleSpinBox_t.value()
            p_a = self.doubleSpinBox_pa.value()*1.0e5
            p_s = self.doubleSpinBox_ps.value()
            OF = self.doubleSpinBox_OF.value()
            SafetyFactor = self.doubleSpinBox_SF.value() 
            No_reuses = self.spinBox_No.value() 
            self.label_currentStatus.setText("Compiling..")
        main.default.p_to = p_s*1.0e5
        main.default.ptf = p_s*1.0e5
        main.default.cycle_type = self.CY
        main.default.MR = OF
        main.default.Safety_factor = SafetyFactor
        #main.default. = No_reuses
        print("Running main...")
        Pc,Isp,m_p,M_p,Tc,Lc=main.Main(T, t, p_a)
        self.input_length.setText(str(Lc))
        self.line_pc.setText(str(Pc))
        self.line_Isp.setText(str(Isp))
        self.line_m.setText(str(m_p))
        self.line_mass.setText(str(M_p))
        self.line_T.setText(str(Tc))
# main
app = QApplication(sys.argv)
mainwindow = MainWindow()
widget = QtWidgets.QStackedWidget()
widget.addWidget(mainwindow)
#widget.setHeight(857)
#widget.setWidth(1024)
widget.setWindowTitle("LRE Design Tool - 0.0")
widget.show()
widget.showMaximized()
try:
    sys.exit(app.exec_())
except:
    print("Exiting")