# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'bino_form.ui'
#
# Created: Sat Sep  6 19:10:13 2014
#      by: PyQt4 UI code generator 4.9.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.setEnabled(True)
        Form.resize(810, 596)
        Form.setMouseTracking(True)
        Form.setAcceptDrops(True)
        Form.setAutoFillBackground(False)
        self.Close_Button = QtGui.QPushButton(Form)
        self.Close_Button.setGeometry(QtCore.QRect(380, 560, 75, 23))
        self.Close_Button.setObjectName(_fromUtf8("Close_Button"))
        self.inputfile = QtGui.QLineEdit(Form)
        self.inputfile.setGeometry(QtCore.QRect(149, 20, 241, 20))
        self.inputfile.setObjectName(_fromUtf8("inputfile"))
        self.Run_Button = QtGui.QPushButton(Form)
        self.Run_Button.setEnabled(False)
        self.Run_Button.setGeometry(QtCore.QRect(280, 560, 75, 23))
        self.Run_Button.setObjectName(_fromUtf8("Run_Button"))
        self.ciclos = QtGui.QLineEdit(Form)
        self.ciclos.setGeometry(QtCore.QRect(150, 50, 61, 20))
        self.ciclos.setObjectName(_fromUtf8("ciclos"))
        self.label_2 = QtGui.QLabel(Form)
        self.label_2.setGeometry(QtCore.QRect(20, 50, 111, 20))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.foodweb_checkbox = QtGui.QCheckBox(Form)
        self.foodweb_checkbox.setEnabled(True)
        self.foodweb_checkbox.setGeometry(QtCore.QRect(250, 102, 161, 17))
        self.foodweb_checkbox.setObjectName(_fromUtf8("foodweb_checkbox"))
        self.pl_ext_period = QtGui.QLineEdit(Form)
        self.pl_ext_period.setGeometry(QtCore.QRect(160, 400, 31, 20))
        self.pl_ext_period.setDragEnabled(True)
        self.pl_ext_period.setObjectName(_fromUtf8("pl_ext_period"))
        self.label_4 = QtGui.QLabel(Form)
        self.label_4.setGeometry(QtCore.QRect(110, 400, 51, 20))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.pl_ext_spike = QtGui.QLineEdit(Form)
        self.pl_ext_spike.setGeometry(QtCore.QRect(250, 400, 31, 20))
        self.pl_ext_spike.setDragEnabled(True)
        self.pl_ext_spike.setObjectName(_fromUtf8("pl_ext_spike"))
        self.label_5 = QtGui.QLabel(Form)
        self.label_5.setGeometry(QtCore.QRect(210, 400, 51, 20))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.label_6 = QtGui.QLabel(Form)
        self.label_6.setGeometry(QtCore.QRect(290, 400, 51, 20))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.pl_ext_numperiod = QtGui.QLineEdit(Form)
        self.pl_ext_numperiod.setGeometry(QtCore.QRect(340, 400, 31, 20))
        self.pl_ext_numperiod.setDragEnabled(True)
        self.pl_ext_numperiod.setObjectName(_fromUtf8("pl_ext_numperiod"))
        self.pl_ext_rate = QtGui.QLineEdit(Form)
        self.pl_ext_rate.setGeometry(QtCore.QRect(420, 400, 41, 20))
        self.pl_ext_rate.setDragEnabled(True)
        self.pl_ext_rate.setObjectName(_fromUtf8("pl_ext_rate"))
        self.label_7 = QtGui.QLabel(Form)
        self.label_7.setGeometry(QtCore.QRect(380, 400, 41, 20))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.pl_ext_start = QtGui.QLineEdit(Form)
        self.pl_ext_start.setGeometry(QtCore.QRect(520, 400, 41, 20))
        self.pl_ext_start.setDragEnabled(True)
        self.pl_ext_start.setObjectName(_fromUtf8("pl_ext_start"))
        self.label_8 = QtGui.QLabel(Form)
        self.label_8.setGeometry(QtCore.QRect(470, 400, 51, 20))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.pl_ext_species = QtGui.QLineEdit(Form)
        self.pl_ext_species.setGeometry(QtCore.QRect(630, 400, 61, 20))
        self.pl_ext_species.setDragEnabled(True)
        self.pl_ext_species.setObjectName(_fromUtf8("pl_ext_species"))
        self.label_9 = QtGui.QLabel(Form)
        self.label_9.setGeometry(QtCore.QRect(570, 400, 51, 20))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.pol_ext_period = QtGui.QLineEdit(Form)
        self.pol_ext_period.setGeometry(QtCore.QRect(160, 480, 31, 20))
        self.pol_ext_period.setDragEnabled(True)
        self.pol_ext_period.setObjectName(_fromUtf8("pol_ext_period"))
        self.pol_ext_numperiod = QtGui.QLineEdit(Form)
        self.pol_ext_numperiod.setGeometry(QtCore.QRect(340, 480, 31, 20))
        self.pol_ext_numperiod.setDragEnabled(True)
        self.pol_ext_numperiod.setObjectName(_fromUtf8("pol_ext_numperiod"))
        self.pol_ext_start = QtGui.QLineEdit(Form)
        self.pol_ext_start.setGeometry(QtCore.QRect(520, 480, 31, 20))
        self.pol_ext_start.setObjectName(_fromUtf8("pol_ext_start"))
        self.label_10 = QtGui.QLabel(Form)
        self.label_10.setGeometry(QtCore.QRect(470, 480, 41, 20))
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.label_11 = QtGui.QLabel(Form)
        self.label_11.setGeometry(QtCore.QRect(290, 480, 51, 20))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.pol_ext_rate = QtGui.QLineEdit(Form)
        self.pol_ext_rate.setGeometry(QtCore.QRect(420, 480, 41, 20))
        self.pol_ext_rate.setDragEnabled(True)
        self.pol_ext_rate.setObjectName(_fromUtf8("pol_ext_rate"))
        self.pol_ext_species = QtGui.QLineEdit(Form)
        self.pol_ext_species.setGeometry(QtCore.QRect(630, 480, 61, 20))
        self.pol_ext_species.setDragEnabled(True)
        self.pol_ext_species.setObjectName(_fromUtf8("pol_ext_species"))
        self.label_13 = QtGui.QLabel(Form)
        self.label_13.setGeometry(QtCore.QRect(210, 480, 51, 20))
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.label_14 = QtGui.QLabel(Form)
        self.label_14.setGeometry(QtCore.QRect(380, 480, 41, 20))
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.label_15 = QtGui.QLabel(Form)
        self.label_15.setGeometry(QtCore.QRect(570, 480, 51, 20))
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.label_16 = QtGui.QLabel(Form)
        self.label_16.setGeometry(QtCore.QRect(110, 480, 51, 20))
        self.label_16.setObjectName(_fromUtf8("label_16"))
        self.pol_ext_spike = QtGui.QLineEdit(Form)
        self.pol_ext_spike.setGeometry(QtCore.QRect(250, 480, 31, 20))
        self.pol_ext_spike.setDragEnabled(True)
        self.pol_ext_spike.setObjectName(_fromUtf8("pol_ext_spike"))
        self.label_12 = QtGui.QLabel(Form)
        self.label_12.setGeometry(QtCore.QRect(20, 370, 161, 20))
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.label_17 = QtGui.QLabel(Form)
        self.label_17.setGeometry(QtCore.QRect(20, 450, 171, 20))
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.ClearPlants = QtGui.QPushButton(Form)
        self.ClearPlants.setGeometry(QtCore.QRect(720, 400, 61, 21))
        self.ClearPlants.setObjectName(_fromUtf8("ClearPlants"))
        self.ClearPolin = QtGui.QPushButton(Form)
        self.ClearPolin.setGeometry(QtCore.QRect(720, 480, 61, 21))
        self.ClearPolin.setObjectName(_fromUtf8("ClearPolin"))
        self.output_suffix = QtGui.QLineEdit(Form)
        self.output_suffix.setGeometry(QtCore.QRect(650, 20, 131, 20))
        self.output_suffix.setObjectName(_fromUtf8("output_suffix"))
        self.label_3 = QtGui.QLabel(Form)
        self.label_3.setGeometry(QtCore.QRect(540, 20, 101, 20))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_18 = QtGui.QLabel(Form)
        self.label_18.setGeometry(QtCore.QRect(585, 64, 201, 20))
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.Comments_text = QtGui.QPlainTextEdit(Form)
        self.Comments_text.setGeometry(QtCore.QRect(480, 90, 301, 101))
        self.Comments_text.setObjectName(_fromUtf8("Comments_text"))
        self.label_report = QtGui.QLabel(Form)
        self.label_report.setGeometry(QtCore.QRect(250, 520, 111, 20))
        self.label_report.setText(_fromUtf8(""))
        self.label_report.setObjectName(_fromUtf8("label_report"))
        self.URL_report = QtGui.QLabel(Form)
        self.URL_report.setGeometry(QtCore.QRect(390, 510, 291, 40))
        self.URL_report.setText(_fromUtf8(""))
        self.URL_report.setOpenExternalLinks(True)
        self.URL_report.setObjectName(_fromUtf8("URL_report"))
        self.save_output_checkbox = QtGui.QCheckBox(Form)
        self.save_output_checkbox.setGeometry(QtCore.QRect(250, 120, 151, 17))
        self.save_output_checkbox.setObjectName(_fromUtf8("save_output_checkbox"))
        self.label_19 = QtGui.QLabel(Form)
        self.label_19.setGeometry(QtCore.QRect(20, 340, 161, 20))
        self.label_19.setObjectName(_fromUtf8("label_19"))
        self.random_removal = QtGui.QLineEdit(Form)
        self.random_removal.setGeometry(QtCore.QRect(200, 340, 41, 20))
        self.random_removal.setDragEnabled(True)
        self.random_removal.setObjectName(_fromUtf8("random_removal"))
        self.select_input_file = QtGui.QPushButton(Form)
        self.select_input_file.setGeometry(QtCore.QRect(20, 20, 101, 21))
        self.select_input_file.setObjectName(_fromUtf8("select_input_file"))
        self.Error_msg = QtGui.QLabel(Form)
        self.Error_msg.setGeometry(QtCore.QRect(240, 510, 391, 20))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 127, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 63, 63))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Midlight, brush)
        brush = QtGui.QBrush(QtGui.QColor(127, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Dark, brush)
        brush = QtGui.QBrush(QtGui.QColor(170, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Mid, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.BrightText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Shadow, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 127, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.AlternateBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 220))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.ToolTipBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.ToolTipText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 127, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 63, 63))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Midlight, brush)
        brush = QtGui.QBrush(QtGui.QColor(127, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Dark, brush)
        brush = QtGui.QBrush(QtGui.QColor(170, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Mid, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.BrightText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Shadow, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 127, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.AlternateBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 220))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.ToolTipBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.ToolTipText, brush)
        brush = QtGui.QBrush(QtGui.QColor(127, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 127, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 63, 63))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Midlight, brush)
        brush = QtGui.QBrush(QtGui.QColor(127, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Dark, brush)
        brush = QtGui.QBrush(QtGui.QColor(170, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Mid, brush)
        brush = QtGui.QBrush(QtGui.QColor(127, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.BrightText, brush)
        brush = QtGui.QBrush(QtGui.QColor(127, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Shadow, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.AlternateBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 220))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.ToolTipBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.ToolTipText, brush)
        self.Error_msg.setPalette(palette)
        self.Error_msg.setText(_fromUtf8(""))
        self.Error_msg.setObjectName(_fromUtf8("Error_msg"))
        self.URL_outputs = QtGui.QLabel(Form)
        self.URL_outputs.setGeometry(QtCore.QRect(480, 560, 131, 21))
        self.URL_outputs.setText(_fromUtf8(""))
        self.URL_outputs.setOpenExternalLinks(True)
        self.URL_outputs.setObjectName(_fromUtf8("URL_outputs"))
        self.TypeofMutualism = QtGui.QGroupBox(Form)
        self.TypeofMutualism.setGeometry(QtCore.QRect(20, 80, 191, 131))
        self.TypeofMutualism.setObjectName(_fromUtf8("TypeofMutualism"))
        self.TOM_LAp_abs_Button = QtGui.QRadioButton(self.TypeofMutualism)
        self.TOM_LAp_abs_Button.setEnabled(True)
        self.TOM_LAp_abs_Button.setGeometry(QtCore.QRect(10, 60, 161, 17))
        self.TOM_LAp_abs_Button.setChecked(False)
        self.TOM_LAp_abs_Button.setObjectName(_fromUtf8("TOM_LAp_abs_Button"))
        self.buttonGroup = QtGui.QButtonGroup(Form)
        self.buttonGroup.setObjectName(_fromUtf8("buttonGroup"))
        self.buttonGroup.addButton(self.TOM_LAp_abs_Button)
        self.TOM_May_Button = QtGui.QRadioButton(self.TypeofMutualism)
        self.TOM_May_Button.setGeometry(QtCore.QRect(10, 100, 161, 17))
        self.TOM_May_Button.setObjectName(_fromUtf8("TOM_May_Button"))
        self.buttonGroup.addButton(self.TOM_May_Button)
        self.TOM_LAp_u_Button = QtGui.QRadioButton(self.TypeofMutualism)
        self.TOM_LAp_u_Button.setEnabled(True)
        self.TOM_LAp_u_Button.setGeometry(QtCore.QRect(10, 80, 161, 17))
        self.TOM_LAp_u_Button.setObjectName(_fromUtf8("TOM_LAp_u_Button"))
        self.buttonGroup.addButton(self.TOM_LAp_u_Button)
        self.TOM_LAp_vh_Button = QtGui.QRadioButton(self.TypeofMutualism)
        self.TOM_LAp_vh_Button.setEnabled(True)
        self.TOM_LAp_vh_Button.setGeometry(QtCore.QRect(10, 40, 161, 17))
        self.TOM_LAp_vh_Button.setCheckable(True)
        self.TOM_LAp_vh_Button.setChecked(True)
        self.TOM_LAp_vh_Button.setObjectName(_fromUtf8("TOM_LAp_vh_Button"))
        self.buttonGroup.addButton(self.TOM_LAp_vh_Button)
        self.TOM_None_Button = QtGui.QRadioButton(self.TypeofMutualism)
        self.TOM_None_Button.setGeometry(QtCore.QRect(10, 20, 161, 17))
        self.TOM_None_Button.setObjectName(_fromUtf8("TOM_None_Button"))
        self.buttonGroup.addButton(self.TOM_None_Button)
        self.TypeofBlossom = QtGui.QGroupBox(Form)
        self.TypeofBlossom.setGeometry(QtCore.QRect(210, 220, 301, 41))
        self.TypeofBlossom.setTitle(_fromUtf8(""))
        self.TypeofBlossom.setObjectName(_fromUtf8("TypeofBlossom"))
        self.BType_Binary = QtGui.QRadioButton(self.TypeofBlossom)
        self.BType_Binary.setGeometry(QtCore.QRect(10, 10, 82, 17))
        self.BType_Binary.setChecked(True)
        self.BType_Binary.setObjectName(_fromUtf8("BType_Binary"))
        self.buttonGroup_2 = QtGui.QButtonGroup(Form)
        self.buttonGroup_2.setObjectName(_fromUtf8("buttonGroup_2"))
        self.buttonGroup_2.addButton(self.BType_Binary)
        self.BType_Gaussian = QtGui.QRadioButton(self.TypeofBlossom)
        self.BType_Gaussian.setGeometry(QtCore.QRect(100, 10, 121, 17))
        self.BType_Gaussian.setObjectName(_fromUtf8("BType_Gaussian"))
        self.buttonGroup_2.addButton(self.BType_Gaussian)
        self.plants_blossom_sd = QtGui.QLineEdit(self.TypeofBlossom)
        self.plants_blossom_sd.setEnabled(False)
        self.plants_blossom_sd.setGeometry(QtCore.QRect(240, 10, 41, 20))
        self.plants_blossom_sd.setDragEnabled(True)
        self.plants_blossom_sd.setObjectName(_fromUtf8("plants_blossom_sd"))
        self.label_21 = QtGui.QLabel(self.TypeofBlossom)
        self.label_21.setGeometry(QtCore.QRect(210, 10, 30, 20))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Symbol"))
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.label_21.setFont(font)
        self.label_21.setObjectName(_fromUtf8("label_21"))
        self.label_20 = QtGui.QLabel(Form)
        self.label_20.setGeometry(QtCore.QRect(20, 230, 151, 16))
        self.label_20.setObjectName(_fromUtf8("label_20"))
        self.plants_blossom = QtGui.QLineEdit(Form)
        self.plants_blossom.setGeometry(QtCore.QRect(170, 230, 31, 20))
        self.plants_blossom.setDragEnabled(True)
        self.plants_blossom.setObjectName(_fromUtf8("plants_blossom"))
        self.label_22 = QtGui.QLabel(Form)
        self.label_22.setGeometry(QtCore.QRect(570, 230, 51, 20))
        self.label_22.setObjectName(_fromUtf8("label_22"))
        self.blossom_pert_species = QtGui.QLineEdit(Form)
        self.blossom_pert_species.setGeometry(QtCore.QRect(630, 230, 61, 20))
        self.blossom_pert_species.setDragEnabled(True)
        self.blossom_pert_species.setObjectName(_fromUtf8("blossom_pert_species"))
        self.ClearBlossom = QtGui.QPushButton(Form)
        self.ClearBlossom.setGeometry(QtCore.QRect(720, 230, 61, 21))
        self.ClearBlossom.setObjectName(_fromUtf8("ClearBlossom"))
        self.label_23 = QtGui.QLabel(Form)
        self.label_23.setGeometry(QtCore.QRect(20, 280, 211, 16))
        self.label_23.setObjectName(_fromUtf8("label_23"))
        self.Bssvar_period = QtGui.QLineEdit(Form)
        self.Bssvar_period.setGeometry(QtCore.QRect(200, 300, 41, 20))
        self.Bssvar_period.setDragEnabled(True)
        self.Bssvar_period.setObjectName(_fromUtf8("Bssvar_period"))
        self.label_24 = QtGui.QLabel(Form)
        self.label_24.setGeometry(QtCore.QRect(140, 300, 51, 20))
        self.label_24.setObjectName(_fromUtf8("label_24"))
        self.label_25 = QtGui.QLabel(Form)
        self.label_25.setGeometry(QtCore.QRect(250, 300, 30, 20))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Symbol"))
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.label_25.setFont(font)
        self.label_25.setObjectName(_fromUtf8("label_25"))
        self.Bssvar_sd = QtGui.QLineEdit(Form)
        self.Bssvar_sd.setEnabled(True)
        self.Bssvar_sd.setGeometry(QtCore.QRect(270, 300, 41, 20))
        self.Bssvar_sd.setDragEnabled(True)
        self.Bssvar_sd.setObjectName(_fromUtf8("Bssvar_sd"))
        self.Bssvar_species = QtGui.QLineEdit(Form)
        self.Bssvar_species.setGeometry(QtCore.QRect(630, 300, 61, 20))
        self.Bssvar_species.setDragEnabled(True)
        self.Bssvar_species.setObjectName(_fromUtf8("Bssvar_species"))
        self.label_26 = QtGui.QLabel(Form)
        self.label_26.setGeometry(QtCore.QRect(570, 300, 51, 20))
        self.label_26.setObjectName(_fromUtf8("label_26"))
        self.ClearBssvar = QtGui.QPushButton(Form)
        self.ClearBssvar.setGeometry(QtCore.QRect(720, 300, 61, 21))
        self.ClearBssvar.setObjectName(_fromUtf8("ClearBssvar"))
        self.Modulation = QtGui.QGroupBox(Form)
        self.Modulation.setGeometry(QtCore.QRect(330, 270, 201, 111))
        self.Modulation.setObjectName(_fromUtf8("Modulation"))
        self.Bssvar_Type_none = QtGui.QRadioButton(self.Modulation)
        self.Bssvar_Type_none.setGeometry(QtCore.QRect(10, 20, 82, 17))
        self.Bssvar_Type_none.setChecked(True)
        self.Bssvar_Type_none.setObjectName(_fromUtf8("Bssvar_Type_none"))
        self.Bssvar_Type_sin_period = QtGui.QLineEdit(self.Modulation)
        self.Bssvar_Type_sin_period.setEnabled(False)
        self.Bssvar_Type_sin_period.setGeometry(QtCore.QRect(150, 80, 41, 20))
        self.Bssvar_Type_sin_period.setDragEnabled(True)
        self.Bssvar_Type_sin_period.setObjectName(_fromUtf8("Bssvar_Type_sin_period"))
        self.label_29 = QtGui.QLabel(self.Modulation)
        self.label_29.setGeometry(QtCore.QRect(90, 80, 51, 20))
        self.label_29.setObjectName(_fromUtf8("label_29"))
        self.Bssvar_Type_linear = QtGui.QRadioButton(self.Modulation)
        self.Bssvar_Type_linear.setGeometry(QtCore.QRect(10, 50, 91, 17))
        self.Bssvar_Type_linear.setObjectName(_fromUtf8("Bssvar_Type_linear"))
        self.Bssvar_Type_linear_slope = QtGui.QLineEdit(self.Modulation)
        self.Bssvar_Type_linear_slope.setEnabled(False)
        self.Bssvar_Type_linear_slope.setGeometry(QtCore.QRect(150, 50, 41, 20))
        self.Bssvar_Type_linear_slope.setDragEnabled(True)
        self.Bssvar_Type_linear_slope.setObjectName(_fromUtf8("Bssvar_Type_linear_slope"))
        self.label_30 = QtGui.QLabel(self.Modulation)
        self.label_30.setGeometry(QtCore.QRect(90, 50, 51, 20))
        self.label_30.setObjectName(_fromUtf8("label_30"))
        self.Bssvar_Type_sin = QtGui.QRadioButton(self.Modulation)
        self.Bssvar_Type_sin.setGeometry(QtCore.QRect(10, 80, 51, 17))
        self.Bssvar_Type_sin.setObjectName(_fromUtf8("Bssvar_Type_sin"))
        self.label_2.setBuddy(self.label_2)
        self.label_4.setBuddy(self.inputfile)
        self.label_5.setBuddy(self.inputfile)
        self.label_6.setBuddy(self.inputfile)
        self.label_7.setBuddy(self.inputfile)
        self.label_8.setBuddy(self.inputfile)
        self.label_9.setBuddy(self.inputfile)
        self.label_10.setBuddy(self.inputfile)
        self.label_11.setBuddy(self.inputfile)
        self.label_13.setBuddy(self.inputfile)
        self.label_14.setBuddy(self.inputfile)
        self.label_15.setBuddy(self.inputfile)
        self.label_16.setBuddy(self.inputfile)
        self.label_12.setBuddy(self.inputfile)
        self.label_17.setBuddy(self.inputfile)
        self.label_3.setBuddy(self.inputfile)
        self.label_18.setBuddy(self.inputfile)
        self.label_report.setBuddy(self.inputfile)
        self.label_19.setBuddy(self.inputfile)
        self.Error_msg.setBuddy(self.inputfile)
        self.label_21.setBuddy(self.inputfile)
        self.label_20.setBuddy(self.inputfile)
        self.label_22.setBuddy(self.inputfile)
        self.label_23.setBuddy(self.inputfile)
        self.label_24.setBuddy(self.inputfile)
        self.label_25.setBuddy(self.inputfile)
        self.label_26.setBuddy(self.inputfile)
        self.label_29.setBuddy(self.inputfile)
        self.label_30.setBuddy(self.inputfile)

        self.retranslateUi(Form)
        QtCore.QObject.connect(self.Close_Button, QtCore.SIGNAL(_fromUtf8("clicked()")), Form.close)
        QtCore.QObject.connect(self.Run_Button, QtCore.SIGNAL(_fromUtf8("clicked()")), self.ciclos.update)
        QtCore.QObject.connect(self.Run_Button, QtCore.SIGNAL(_fromUtf8("clicked()")), self.foodweb_checkbox.update)
        QtCore.QObject.connect(self.ClearPlants, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pl_ext_species.clear)
        QtCore.QObject.connect(self.ClearPlants, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pl_ext_rate.clear)
        QtCore.QObject.connect(self.ClearPlants, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pl_ext_numperiod.clear)
        QtCore.QObject.connect(self.ClearPlants, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pl_ext_spike.clear)
        QtCore.QObject.connect(self.ClearPlants, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pl_ext_period.clear)
        QtCore.QObject.connect(self.ClearPlants, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pl_ext_start.clear)
        QtCore.QObject.connect(self.ClearPolin, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pol_ext_species.clear)
        QtCore.QObject.connect(self.ClearPolin, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pol_ext_start.clear)
        QtCore.QObject.connect(self.ClearPolin, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pol_ext_rate.clear)
        QtCore.QObject.connect(self.ClearPolin, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pol_ext_numperiod.clear)
        QtCore.QObject.connect(self.ClearPolin, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pol_ext_spike.clear)
        QtCore.QObject.connect(self.ClearPolin, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pol_ext_period.clear)
        QtCore.QObject.connect(self.Run_Button, QtCore.SIGNAL(_fromUtf8("clicked()")), self.URL_report.repaint)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.Close_Button.setText(QtGui.QApplication.translate("Form", "Close", None, QtGui.QApplication.UnicodeUTF8))
        self.Run_Button.setText(QtGui.QApplication.translate("Form", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("Form", "Years", None, QtGui.QApplication.UnicodeUTF8))
        self.foodweb_checkbox.setText(QtGui.QApplication.translate("Form", "Foodweb", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("Form", "Years", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("Form", "Spike", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("Form", "Repeat", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setText(QtGui.QApplication.translate("Form", "Rate", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setText(QtGui.QApplication.translate("Form", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setText(QtGui.QApplication.translate("Form", "Species", None, QtGui.QApplication.UnicodeUTF8))
        self.label_10.setText(QtGui.QApplication.translate("Form", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setText(QtGui.QApplication.translate("Form", "Repeat", None, QtGui.QApplication.UnicodeUTF8))
        self.label_13.setText(QtGui.QApplication.translate("Form", "Spike", None, QtGui.QApplication.UnicodeUTF8))
        self.label_14.setText(QtGui.QApplication.translate("Form", "Rate", None, QtGui.QApplication.UnicodeUTF8))
        self.label_15.setText(QtGui.QApplication.translate("Form", "Species", None, QtGui.QApplication.UnicodeUTF8))
        self.label_16.setText(QtGui.QApplication.translate("Form", "Years", None, QtGui.QApplication.UnicodeUTF8))
        self.label_12.setText(QtGui.QApplication.translate("Form", "Plant perturbations", None, QtGui.QApplication.UnicodeUTF8))
        self.label_17.setText(QtGui.QApplication.translate("Form", "Pollinator perturbations", None, QtGui.QApplication.UnicodeUTF8))
        self.ClearPlants.setText(QtGui.QApplication.translate("Form", "Clear", None, QtGui.QApplication.UnicodeUTF8))
        self.ClearPolin.setText(QtGui.QApplication.translate("Form", "Clear", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("Form", "Output suffix", None, QtGui.QApplication.UnicodeUTF8))
        self.label_18.setText(QtGui.QApplication.translate("Form", "Comments for report", None, QtGui.QApplication.UnicodeUTF8))
        self.save_output_checkbox.setText(QtGui.QApplication.translate("Form", "Save output data", None, QtGui.QApplication.UnicodeUTF8))
        self.label_19.setText(QtGui.QApplication.translate("Form", "Random links removal", None, QtGui.QApplication.UnicodeUTF8))
        self.select_input_file.setText(QtGui.QApplication.translate("Form", "Input file", None, QtGui.QApplication.UnicodeUTF8))
        self.TypeofMutualism.setTitle(QtGui.QApplication.translate("Form", "Mutualism Model", None, QtGui.QApplication.UnicodeUTF8))
        self.TOM_LAp_abs_Button.setText(QtGui.QApplication.translate("Form", "Logistic |req|", None, QtGui.QApplication.UnicodeUTF8))
        self.TOM_May_Button.setText(QtGui.QApplication.translate("Form", "May", None, QtGui.QApplication.UnicodeUTF8))
        self.TOM_LAp_u_Button.setText(QtGui.QApplication.translate("Form", "Logistic u(req)", None, QtGui.QApplication.UnicodeUTF8))
        self.TOM_LAp_vh_Button.setText(QtGui.QApplication.translate("Form", "Verhulst", None, QtGui.QApplication.UnicodeUTF8))
        self.TOM_None_Button.setText(QtGui.QApplication.translate("Form", "None", None, QtGui.QApplication.UnicodeUTF8))
        self.BType_Binary.setText(QtGui.QApplication.translate("Form", "Binary", None, QtGui.QApplication.UnicodeUTF8))
        self.BType_Gaussian.setText(QtGui.QApplication.translate("Form", "Gaussian", None, QtGui.QApplication.UnicodeUTF8))
        self.label_21.setText(QtGui.QApplication.translate("Form", "s", None, QtGui.QApplication.UnicodeUTF8))
        self.label_20.setText(QtGui.QApplication.translate("Form", "Plants blossom prob.", None, QtGui.QApplication.UnicodeUTF8))
        self.label_22.setText(QtGui.QApplication.translate("Form", "Species", None, QtGui.QApplication.UnicodeUTF8))
        self.ClearBlossom.setText(QtGui.QApplication.translate("Form", "Clear", None, QtGui.QApplication.UnicodeUTF8))
        self.label_23.setText(QtGui.QApplication.translate("Form", "Plants blossom variability", None, QtGui.QApplication.UnicodeUTF8))
        self.label_24.setText(QtGui.QApplication.translate("Form", "Period", None, QtGui.QApplication.UnicodeUTF8))
        self.label_25.setText(QtGui.QApplication.translate("Form", "s", None, QtGui.QApplication.UnicodeUTF8))
        self.label_26.setText(QtGui.QApplication.translate("Form", "Species", None, QtGui.QApplication.UnicodeUTF8))
        self.ClearBssvar.setText(QtGui.QApplication.translate("Form", "Clear", None, QtGui.QApplication.UnicodeUTF8))
        self.Modulation.setTitle(QtGui.QApplication.translate("Form", "Modulation", None, QtGui.QApplication.UnicodeUTF8))
        self.Bssvar_Type_none.setText(QtGui.QApplication.translate("Form", "None", None, QtGui.QApplication.UnicodeUTF8))
        self.label_29.setText(QtGui.QApplication.translate("Form", "Period", None, QtGui.QApplication.UnicodeUTF8))
        self.Bssvar_Type_linear.setText(QtGui.QApplication.translate("Form", "linear", None, QtGui.QApplication.UnicodeUTF8))
        self.label_30.setText(QtGui.QApplication.translate("Form", "Slope", None, QtGui.QApplication.UnicodeUTF8))
        self.Bssvar_Type_sin.setText(QtGui.QApplication.translate("Form", "sin", None, QtGui.QApplication.UnicodeUTF8))

