#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Here we provide the necessary imports.
# The basic GUI widgets are located in QtGui module. 
import sys
import os
import b_sim
import sigmund_GLOBALS as sgGL
import sigmund_release as sgRL
import sigmund_graphs as sggraph
import sigmund_common as sgcom
# from PyQt4 import QtCore, QtGui
from bino_form import *
import re

class StartQT4(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.hay_error = False
        self.lista_err = []
        self.input_dir = './input'
        self.ciclos = 100
        self.release = sgRL.RELEASE_NUMBER
        self.daterelease = sgRL.RELEASE_STRING
        self.plants_blossom = 1.00
        self.plants_blossom_sd = 0.01
        self.typeofblossom = 'Binary'
        self.Bssvar_period = 0.1
        self.Bssvar_sd = 0.0
        self.Bssvar_Type_linear_slope = 1.0
        self.Bssvar_Type_sin_period = 50

        self.filename = ""
        self.input_file = ""
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.setWindowTitle('SIGMUND ' + "%0.2f " % (self.release / 100) +\
                             self.daterelease)
        self.setWindowIcon(QtGui.QIcon('upm.gif'))
        self.ui.ciclos.setText(str(self.ciclos))
        self.ui.plants_blossom.setText(str(self.plants_blossom))
        self.ui.plants_blossom_sd.setText(str(self.plants_blossom_sd))
        self.ui.Bssvar_period.setText(str(self.Bssvar_period))
        self.ui.Bssvar_sd.setText(str(self.Bssvar_sd))
        self.ui.Bssvar_Type_linear_slope.setText(str(self.Bssvar_Type_linear_slope))
        self.ui.Bssvar_Type_sin_period.setText(str(self.Bssvar_Type_sin_period))
        self.ui.blossom_pert_species.setText('ALL')
        self.haymut = False
        self.haypred = False
        self.May = False
        self.algorithm = 'Verhulst'
        self.typeofmodulation = 'None'
        self.Bssvar_modulationtype_list = []
        QtCore.QObject.connect(self.ui.inputfile, 
                               QtCore.SIGNAL("returnPressed()"), self.add_entry)
        QtCore.QObject.connect(self.ui.ciclos, 
                               QtCore.SIGNAL("returnPressed()"), self.add_entry)
        QtCore.QObject.connect(self.ui.foodweb_checkbox, 
                               QtCore.SIGNAL("stateChanged()"), self.add_entry)
        QtCore.QObject.connect(self.ui.Run_Button, QtCore.SIGNAL("clicked()"), 
                               self.add_entry)
        QtCore.QObject.connect(self.ui.select_input_file, 
                               QtCore.SIGNAL("clicked()"), self.select_file)
        self.ui.inputfile.setEnabled(False)
        self.ui.TOM_None_Button.clicked.connect(self.select_No_mutualism)
        self.ui.TOM_LAp_vh_Button.clicked.connect(self.select_LApproach_vh)
        self.ui.TOM_LAp_abs_Button.clicked.connect(self.select_LApproach_abs)
        self.ui.TOM_LAp_u_Button.clicked.connect(self.select_LApproach_u)
        self.ui.TOM_May_Button.clicked.connect(self.select_May)
        self.ui.BType_Binary.clicked.connect(self.select_BinaryBlossom)
        self.ui.BType_Gaussian.clicked.connect(self.select_GaussianBlossom)
        self.ui.Bssvar_Type_none.clicked.connect(self.select_NoneModulation)
        self.ui.Bssvar_Type_sin.clicked.connect(self.select_sinModulation)
        self.ui.Bssvar_Type_linear.clicked.connect(self.select_linearModulation)
        self.ui.ClearBlossom.clicked.connect(self.clear_blossom_pars)
        self.ui.ClearBssvar.clicked.connect(self.clear_Bssvar_pars)
           
    def select_LApproach_abs(self):
        self.algorithm = 'Logistic_abs'
        
    def select_LApproach_vh(self):
        self.algorithm = 'Verhulst'
        
    def select_LApproach_u(self):
        self.algorithm = 'Logistic_u'
        
    def select_No_mutualism(self):
        self.algorithm = 'NoMutualism'
        
    def select_May(self):
        self.algorithm = 'May'
        
    def select_BinaryBlossom(self):
        self.typeofblossom = 'Binary'
        self.ui.plants_blossom_sd.setEnabled(0)
    
    def select_GaussianBlossom(self):
        self.typeofblossom = 'Gaussian'
        self.ui.plants_blossom_sd.setEnabled(1)
        
    def select_NoneModulation(self):
        self.typeofmodulation = 'None'
        self.ui.Bssvar_Type_sin_period.setEnabled(0)
        self.ui.Bssvar_Type_linear_slope.setEnabled(0)
        
    def select_sinModulation(self):
        self.typeofmodulation = 'sin'
        self.ui.Bssvar_Type_sin_period.setEnabled(1)
        self.ui.Bssvar_Type_linear_slope.setEnabled(0)
        
    def select_linearModulation(self):
        self.typeofmodulation = 'linear'
        self.ui.Bssvar_Type_sin_period.setEnabled(0)
        self.ui.Bssvar_Type_linear_slope.setEnabled(1)
        
    def clear_blossom_pars(self):
        self.ui.BType_Binary.setChecked(True)    
        self.ui.BType_Binary.clicked.connect(self.select_BinaryBlossom)
        self.ui.blossom_pert_species.setText('ALL')
        self.plants_blossom = 1.00
        self.ui.plants_blossom.setText('1.0')
        self.plants_blossom_sd = 0.01
        self.ui.plants_blossom_sd.setText('0.01')
        
    def clear_Bssvar_pars(self):
        self.typeofmodulation = 'None'
        self.ui.Bssvar_Type_none.setChecked(True)  
        self.ui.Bssvar_Type_sin_period.setEnabled(0)
        self.ui.Bssvar_Type_linear_slope.setEnabled(0)
        self.Bssvar_period = 0.1
        self.Bssvar_sd = 0.0
        self.ui.Bssvar_sd.setText(str(self.Bssvar_sd))
        self.ui.Bssvar_period.setText(str(self.Bssvar_period))
        self.Bssvar_Type_linear_slope = 1.0
        self.Bssvar_Type_sin_period = 50
        self.ui.Bssvar_Type_linear_slope.setText(str(self.Bssvar_Type_linear_slope))
        self.ui.Bssvar_Type_sin_period.setText(str(self.Bssvar_Type_sin_period))
        self.ui.Bssvar_species.setText('')
        
    def error_exit(self):
        texto_aux = ""
        for i in self.lista_err:
            print (i)
            texto_aux += i
        self.ui.label_report.setText("")
        self.ui.URL_report.setText("")
        self.ui.Error_msg.setText("<html><head/><body><p align=center><span style=' font-size:12pt; font-weight:600; color:#ff0000;'>" +\
                                  texto_aux + "</span></p></body></html>")
        self.lista_err = []     
            
    def select_file(self):
        self.filename = QtGui.QFileDialog.getOpenFileName(self, 'Open File',
                                                          self.input_dir)
        a = self.filename.replace('\\', '/ ').split('/');
        self.input_file = a[-1]
        self.ui.inputfile.setText(self.input_file)
        self.ui.inputfile.setEnabled(False)
        self.ui.Run_Button.setEnabled(1)

    def create_list_species_affected(self, auxspec):
        auxlspec = []
        for numspe in auxspec:
            if len(numspe) == 1:
                auxlspec.append([int(numspe)])
            else:
                rango = numspe.split(':')
                auxlspec.append([j for j in range(int(rango[0]), 1 +\
                                                   int(rango[-1]))])
        listb = []
        for subla in auxlspec:
            for j in list(subla):
                listb.append(j)    
        listb = sorted(set(listb))   
        return listb
    
    def fetch_pert_data(self, xx_ext_period, xx_ext_spike, xx_ext_start, 
                        xx_ext_rate, xx_ext_numperiod, xx_ext_species):
        ret_extinction = {}
        print(xx_ext_period.text())
        if len(xx_ext_period.text())>0:
            try:
                ret_extinction['period'] = int( float(xx_ext_period.text())\
                                              * sgGL.DAYS_YEAR)
                ret_extinction['spike'] = float(xx_ext_spike.text())
                ret_extinction['start'] = int(xx_ext_start.text())
                ret_extinction['rate'] = float(xx_ext_rate.text())
                ret_extinction['numperiod'] = int(xx_ext_numperiod.text())
                
                print(ret_extinction['period'])
                print (''.join([str(item) for item in ret_extinction]))
                if (xx_ext_species.text().upper() == 'ALL'):
                    ret_extinction['species'] = ['ALL']
                else:
                    auxspec = xx_ext_species.text().split(',')
                    listb = self.create_list_species_affected(auxspec)
                    ret_extinction['species'] = listb
            except:         
                self.lista_err.append("ERROR: bad plant extinction format")
                self.error_exit()
                self.ui.Run_Button.setEnabled(1)
                self.ui.Close_Button.setEnabled(1)
                return
        return(ret_extinction)

    def trim_input_field_re(self):
        return(re.compile('\d+(\.\d+)?'))
    
    def txtfld2float(self,textfield):
        p = self.trim_input_field_re()
        if p.match(textfield) == None:
            el = 0
        else:
            el = float(textfield)
        return el

    def add_entry(self): 
        # Checking input file existence
        try:
            fh = open(self.filename, "r")
        except IOError:
            self.lista_err.append("ERROR: can\'t open file " + self.input_file)
            self.error_exit()
            return
        else:
            fh.close()
        # Testing that file exists       
        dirsal = 'output/'
        dirent = 'input/'
        dirs = os.path.dirname(dirsal)
        try:
            os.stat(dirs)
        except:
            os.makedirs(dirs)
        try:
            self.ciclos = int(self.ui.ciclos.text())
        except:            
            self.lista_err.append("ERROR: bad Days format ")
            self.error_exit()
            return    
        self.ui.label_report.setText("")
        self.ui.Run_Button.setEnabled(0)
        self.ui.Close_Button.setEnabled(0)
        self.ui.Error_msg.setText("")
        self.ui.URL_report.setText("Running")
        self.repaint()
        displayinic = 0
        dirsalida = 'output\\'
        input_fname = self.input_file.replace('_b.txt',
                                           '_a.txt').replace('_c.txt', '_a.txt')
        aux = input_fname.split('_a.txt')
        input_fname = aux[0]
        output_suffix = self.ui.output_suffix.text()
        comentario = self.ui.Comments_text.toPlainText()
        dirs = os.path.dirname(sys.argv[0])
        reportpath = os.path.join(dirs, dirsalida.replace('\\', '/'))
        fichr = reportpath + 'rep_' + input_fname + '_' + self.algorithm +\
                     '_' + output_suffix + '_' + str(int(self.ciclos)) + '.html'
        dispfichsal = fichr.split('/')
        linkname = '<a href=file:///' + fichr + '>' + dispfichsal[-1] + '</a>'
        print("%s" % comentario)
        outputdatasave = self.ui.save_output_checkbox.checkState() > 0
        self.haypred = self.ui.foodweb_checkbox.checkState() > 0
        haysup = 0
        el = self.txtfld2float(self.ui.random_removal.text())
        pb = self.txtfld2float(self.ui.plants_blossom.text())
        pb_sd = self.txtfld2float(self.ui.plants_blossom_sd.text())
        self.Bssvar_period = self.txtfld2float(self.ui.Bssvar_period.text())
        self.Bssvar_sd = self.txtfld2float(self.ui.Bssvar_sd.text())            
        self.Bssvar_modulationtype_list = []
        self.Bssvar_modulationtype_list.append(self.typeofmodulation)
        if (self.typeofmodulation == 'linear'):
            self.Bssvar_modulationtype_list.append(float(self.ui.Bssvar_Type_linear_slope.text()))
        else: 
            if (self.typeofmodulation == 'sin'):
                self.Bssvar_modulationtype_list.append(float(self.ui.Bssvar_Type_sin_period.text()))
        if (self.ui.Bssvar_species.text().upper() == 'ALL'):
            self.Bssvar_species = ['ALL']
        else:
            if (len(self.ui.Bssvar_species.text()) == 0):
                listb = []
            else:
                auxspec = self.ui.Bssvar_species.text().split(',')
                listb = self.create_list_species_affected(auxspec)
            self.Bssvar_species = listb
        # External perturbation data
        plants_extinction = self.fetch_pert_data(self.ui.pl_ext_period,
                                            self.ui.pl_ext_spike,
                                            self.ui.pl_ext_start,
                                            self.ui.pl_ext_rate,
                                            self.ui.pl_ext_numperiod,
                                            self.ui.pl_ext_species)      
        pols_extinction = self.fetch_pert_data(self.ui.pol_ext_period,
                                            self.ui.pol_ext_spike,
                                            self.ui.pol_ext_start,
                                            self.ui.pol_ext_rate,
                                            self.ui.pol_ext_numperiod,
                                            self.ui.pol_ext_species)
        blossom_perturbation = {}        
        if (self.ui.blossom_pert_species.text().upper() == 'ALL'):
            blossom_perturbation = ['ALL']
        else:
            auxspec = self.ui.blossom_pert_species.text().split(',')
            listb = self.create_list_species_affected(auxspec)
            blossom_perturbation = listb
        simulation_params = sgcom.SimulationConditions(filename = input_fname, 
                        year_periods = self.ciclos, 
                        hay_foodweb = self.haypred, hay_superpredadores = haysup,
                        data_save = outputdatasave, dirtrabajo = dirs, 
                        direntrada = dirent, dirsal = dirsalida,
                        eliminarenlaces = el, pl_ext = plants_extinction, 
                        pol_ext = pols_extinction, os = output_suffix, 
                        fichreport = fichr, com = comentario, 
                        algorithm = self.algorithm, plants_blossom_prob = pb, 
                        plants_blossom_sd = pb_sd, 
                        plants_blossom_type = self.typeofblossom, 
                        blossom_pert_list=blossom_perturbation,
                        release=self.release, Bssvar_period=self.Bssvar_period,
                        Bssvar_sd=self.Bssvar_sd,
                        Bssvar_modulationtype_list=self.Bssvar_modulationtype_list, \
                        Bssvar_species=self.Bssvar_species)  
        sim_ret_val = b_sim.bino_mutual (sim_cond = simulation_params)      
        self.ui.label_report.setText("Report file: ")
        self.ui.URL_report.setText(linkname)
        self.ui.Run_Button.setEnabled(1)
        self.ui.Close_Button.setEnabled(1)
        self.ui.URL_outputs.setText("<a href='file:///" +\
                                    reportpath.replace('\\', '/') +\
                                    "'>See all results</a>")       
        sggraph.mutual_render(simulation_params, sim_ret_val, displayinic,
                            self.ciclos * sgGL.DAYS_YEAR)
        if self.haypred:
            sggraph.food_render(simulation_params, sim_ret_val, displayinic, 
                                self.ciclos * sgGL.DAYS_YEAR)

if __name__ == "__main__": 
    app = QtGui.QApplication(sys.argv)
    myapp = StartQT4()
    myapp.show()
    sys.exit(app.exec_())
