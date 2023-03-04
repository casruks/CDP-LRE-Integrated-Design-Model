import os
import sys
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtWidgets import QDialog, QApplication, QWidget
from PyQt5.QtWidgets import QMessageBox

set_images_path = "C:/Users/casru/Documents/GitHub/CDP-LRE-Integrated-Design-Model/images/"

class MainWindow(QDialog):
    
    def __init__(self):
        super(MainWindow, self).__init__()
        loadUi("GUI.ui",self)
        # set initial image to corresponding initial cycle
        self.graphicsView.setStyleSheet("background-image: url("+set_images_path+"EX.png);") 
        
        # Compile button
        self.pushButton_Rcompile.clicked.connect(self.compile)

        # Combo Box - Feed cycles
        listocc=["Expander Cycle", "Staged Combustion Cycle", "Coolant Bleed Cycle", "Gas Generator Cycle"]
        for feedcycle in listocc:
            self.comboBox.addItem(feedcycle)
        self.comboBox.currentIndexChanged.connect(self.combochanged)

    def combochanged(self): #Feed cycle selection
        self.label_currentStatus.setText("")
        if self.comboBox.currentText() == "Expander Cycle":
            self.graphicsView.setStyleSheet("background-image: url("+ set_images_path +"EX.png);") 
            x = 10.
        if self.comboBox.currentText() == "Staged Combustion Cycle":
            self.graphicsView.setStyleSheet("background-image: url("+ set_images_path +"SC.png);") 
            x = 20.
        if self.comboBox.currentText() == "Coolant Bleed Cycle":
            self.graphicsView.setStyleSheet("background-image: url("+ set_images_path +"CB.png);") 
            x = 30.
        if self.comboBox.currentText() == "Gas Generator Cycle":
            self.graphicsView.setStyleSheet("background-image: url("+ set_images_path +"GG.png);") 
            x = 40.
        return print(x)
    
    def compile(self): # action after clicking compile button
        msg = QMessageBox()
        msg.setWindowTitle("Input error!")
        x = 1
        if self.spinBox_T.value() <=0:
            x = 0
            msg.setText("Invalid thrust entry, try again.")
            msg.exec_()
        if self.spinBox_t.value() <=0:
            x = 0
            msg.setText("Invalid thrust time entry, try again.")
            msg.exec_()
        if self.spinBox_pa.value() <=0:
            x = 0
            msg.setText("Invalid ambient pressure entry, try again.")
            msg.exec_()
        if self.spinBox_ps.value() <=0:
            x = 0
            msg.setText("Invalid storage pressure entry, try again.")
            msg.exec_()
        if self.spinBox_OF.value() <=0:
            x = 0
            msg.setText("Invalid OF entry, try again.")
            msg.exec_()
        if self.spinBox_No.value() <=0:
            x = 0
            msg.setText("Invalid No. of reuses entry, try again.")
            msg.exec_()
        if x == 1:
            T = self.spinBox_T.value()
            t = self.spinBox_t.value()
            p_a = self.spinBox_pa.value()
            p_s = self.spinBox_ps.value()
            OF = self.spinBox_OF.value()
            Feed_cycle = self.comboBox.currentText()
            SafetyFactor = self.spinBox_SF.value() 
            No_reuses = self.spinBox_No.value() 
            self.label_currentStatus.setText("Compiling..")
        return T, t, p_a, p_s, OF, Feed_cycle, SafetyFactor, No_reuses
        
# main
app = QApplication(sys.argv)
mainwindow = MainWindow()
widget = QtWidgets.QStackedWidget()
widget.addWidget(mainwindow)
widget.setFixedHeight(857)
widget.setFixedWidth(1024)
widget.show()
try:
    sys.exit(app.exec_())
except:
    print("Exiting")