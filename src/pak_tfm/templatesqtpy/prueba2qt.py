# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'prueba2.ui'
#
# Created: Mon Jun 18 16:00:08 2012
#      by: PyQt4 UI code generator 4.9.1
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
        Form.resize(301, 328)
        self.pushButton = QtGui.QPushButton(Form)
        self.pushButton.setGeometry(QtCore.QRect(200, 280, 75, 23))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.inputfile = QtGui.QLineEdit(Form)
        self.inputfile.setGeometry(QtCore.QRect(112, 20, 161, 20))
        self.inputfile.setObjectName(_fromUtf8("inputfile"))
        self.label = QtGui.QLabel(Form)
        self.label.setGeometry(QtCore.QRect(40, 21, 81, 20))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.Run_Button = QtGui.QPushButton(Form)
        self.Run_Button.setGeometry(QtCore.QRect(40, 280, 75, 23))
        self.Run_Button.setObjectName(_fromUtf8("Run_Button"))
        self.ClearAll_Button = QtGui.QPushButton(Form)
        self.ClearAll_Button.setGeometry(QtCore.QRect(120, 280, 75, 23))
        self.ClearAll_Button.setObjectName(_fromUtf8("ClearAll_Button"))
        self.ciclos = QtGui.QLineEdit(Form)
        self.ciclos.setGeometry(QtCore.QRect(160, 50, 113, 20))
        self.ciclos.setObjectName(_fromUtf8("ciclos"))
        self.label_2 = QtGui.QLabel(Form)
        self.label_2.setGeometry(QtCore.QRect(40, 50, 81, 20))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_2.setFont(font)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.Mutualism_checkbox = QtGui.QCheckBox(Form)
        self.Mutualism_checkbox.setGeometry(QtCore.QRect(100, 110, 131, 17))
        self.Mutualism_checkbox.setObjectName(_fromUtf8("Mutualism_checkbox"))
        self.foodweb_checkbox = QtGui.QCheckBox(Form)
        self.foodweb_checkbox.setGeometry(QtCore.QRect(100, 140, 91, 17))
        self.foodweb_checkbox.setObjectName(_fromUtf8("foodweb_checkbox"))
        self.label.setBuddy(self.inputfile)
        self.label_2.setBuddy(self.label_2)

        self.retranslateUi(Form)
        QtCore.QObject.connect(self.pushButton, QtCore.SIGNAL(_fromUtf8("clicked()")), Form.close)
        QtCore.QObject.connect(self.Run_Button, QtCore.SIGNAL(_fromUtf8("clicked()")), self.inputfile.update)
        QtCore.QObject.connect(self.ClearAll_Button, QtCore.SIGNAL(_fromUtf8("clicked()")), self.inputfile.clear)
        QtCore.QObject.connect(self.ciclos, QtCore.SIGNAL(_fromUtf8("returnPressed()")), self.ciclos.setFocus)
        QtCore.QObject.connect(self.ClearAll_Button, QtCore.SIGNAL(_fromUtf8("clicked()")), self.ciclos.clear)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setText(QtGui.QApplication.translate("Form", "Close", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("Form", "Input file", None, QtGui.QApplication.UnicodeUTF8))
        self.Run_Button.setText(QtGui.QApplication.translate("Form", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.ClearAll_Button.setText(QtGui.QApplication.translate("Form", "Clear All", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("Form", "Days", None, QtGui.QApplication.UnicodeUTF8))
        self.Mutualism_checkbox.setText(QtGui.QApplication.translate("Form", "Mutualism", None, QtGui.QApplication.UnicodeUTF8))
        self.foodweb_checkbox.setText(QtGui.QApplication.translate("Form", "Foodweb", None, QtGui.QApplication.UnicodeUTF8))

