#!/usr/bin/env python3
import sys
import os
import b_sim
import json
import pickle
import sigmund_GLOBALS as sgGL
import sigmund_release as sgRL
import sigmund_graphs as sggraph
import sigmund_common as sgcom
# from PyQt4 import QtCore, QtGui
from bino_form import *
import re

global simulation_params

class StartQT4(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.set_default_params()
        self.setWindowTitle('SIGMUND ' + "%0.2f " % (self.release / 100) +\
                             self.daterelease)
        self.setWindowIcon(QtGui.QIcon('upm.gif'))
        self.ui.TOM_None_Button.clicked.connect(self.select_No_mutualism)
        self.ui.TOM_LAp_vh_Button.clicked.connect(self.select_LApproach_vh)
        self.ui.TOM_LAp_abs_Button.clicked.connect(self.select_LApproach_abs)
        self.ui.TOM_LAp_u_Button.clicked.connect(self.select_LApproach_u)
        self.ui.TOM_May_Button.clicked.connect(self.select_May)
        self.ui.blossom_pert_species.setText('ALL')
        self.ui.BType_Binary.clicked.connect(self.select_BinaryBlossom)
        self.ui.BType_Gaussian.clicked.connect(self.select_GaussianBlossom)
        self.ui.Bssvar_Type_none.clicked.connect(self.select_NoneModulation)
        self.ui.Bssvar_Type_sin.clicked.connect(self.select_sinModulation)
        self.ui.Bssvar_Type_linear.clicked.connect(self.select_linearModulation)
        self.ui.ClearBlossom.clicked.connect(self.clear_blossom_pars)
        self.ui.ClearBssvar.clicked.connect(self.clear_Bssvar_pars)
        self.ui.clear_all_simparams.clicked.connect(self.clear_all_pars)
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
        QtCore.QObject.connect(self.ui.select_simulation, 
                               QtCore.SIGNAL("clicked()"), 
                               self.select_stored_simulation)
        QtCore.QObject.connect(self.ui.save_simulation, 
                               QtCore.SIGNAL("clicked()"), 
                               self.save_simulation_file)
        self.update_ui()

    def set_default_params(self):
        self.hay_error = False
        self.lista_err = []
        self.input_dir = './input'
        self.ciclos = 100
        self.release = sgRL.RELEASE_NUMBER
        self.daterelease = sgRL.RELEASE_STRING
        self.filename = ""
        self.input_file = ""
        self.input_fname_raw = ""
        self.haymut = False
        self.haypred = False
        self.May = False
        self.algorithm = 'Verhulst'
        self.typeofmodulation = 'None'
        self.Bssvar_modulationtype_list = []
        self.dirsal = sgGL.OUTPUTILES_PATH
        self.dirent = sgGL.INPUTFILES_PATH
        self.dirs = os.path.dirname(self.dirsal)
        self.clear_blossom_pars()
        self.clear_Bssvar_pars()
        
    def update_ui(self):
        self.ui.ciclos.setText(str(self.ciclos))
        self.ui.plants_blossom.setText(str(self.plants_blossom))
        self.ui.plants_blossom_sd.setText(str(self.plants_blossom_sd))
        self.ui.Bssvar_period.setText(str(self.Bssvar_period))
        self.ui.Bssvar_sd.setText(str(self.Bssvar_sd))
        self.ui.Bssvar_Type_linear_slope.setText(str(self.Bssvar_Type_linear_slope))
        self.ui.Bssvar_Type_sin_period.setText(str(self.Bssvar_Type_sin_period))
        self.ui.inputfile.setText(self.input_fname_raw)
        self.ui.inputfile.setEnabled(False)
        self.ui.simulation_file.setEnabled(False)
        self.repaint()

    def select_LApproach_abs(self):
        self.algorithm = 'Logistic_abs'
        self.ui.TOM_LAp_abs_Button.setChecked(True)
        
    def select_LApproach_vh(self):
        self.algorithm = 'Verhulst'
        self.ui.TOM_LAp_vh_Button.setChecked(True)
        
    def select_LApproach_u(self):
        self.algorithm = 'Logistic_u'
        self.ui.TOM_LAp_u_Button.setChecked(True)
        
    def select_No_mutualism(self):
        self.algorithm = 'NoMutualism'
        self.ui.TOM_None_Button.setChecked(True)
        
    def select_May(self):
        self.algorithm = 'May'
        self.ui.TOM_May_Button.setChecked(True)
        
    def select_BinaryBlossom(self):
        self.typeofblossom = 'Binary'
        self.ui.plants_blossom_sd.setEnabled(0)
        self.ui.BType_Binary.setChecked(True)
        self.ui.BType_Gaussian.setChecked(False)
    
    def select_GaussianBlossom(self):
        self.typeofblossom = 'Gaussian'
        self.ui.plants_blossom_sd.setEnabled(1)
        self.ui.BType_Gaussian.setChecked(True)
        self.ui.BType_Binary.setChecked(False)
        
    def select_NoneModulation(self):
        self.typeofmodulation = 'None'
        self.ui.Bssvar_Type_sin_period.setEnabled(0)
        self.ui.Bssvar_Type_linear_slope.setEnabled(0)
        self.ui.Bssvar_Type_none.setChecked(True)
        
    def select_sinModulation(self):
        self.typeofmodulation = 'sin'
        self.ui.Bssvar_Type_sin_period.setEnabled(1)
        self.ui.Bssvar_Type_linear_slope.setEnabled(0)
        self.ui.Bssvar_Type_sin.setChecked(True)
        
    def select_linearModulation(self):
        self.typeofmodulation = 'linear'
        self.ui.Bssvar_Type_sin_period.setEnabled(0)
        self.ui.Bssvar_Type_linear_slope.setEnabled(1)
        self.ui.Bssvar_Type_linear.setChecked(True)
        
    def clear_blossom_pars(self):
        self.ui.BType_Binary.setChecked(True)    
        self.ui.BType_Binary.clicked.connect(self.select_BinaryBlossom)
        self.typeofblossom = 'Binary'
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
        
    def clear_pert_data(self):
        self.ui.pl_ext_period.setText('')
        self.ui.pl_ext_spike.setText('')
        self.ui.pl_ext_numperiod.setText('')
        self.ui.pl_ext_rate.setText('')
        self.ui.pl_ext_start.setText('')
        self.ui.pl_ext_species.setText('')
        self.ui.pol_ext_period.setText('')
        self.ui.pol_ext_spike.setText('')
        self.ui.pol_ext_numperiod.setText('')
        self.ui.pol_ext_rate.setText('')
        self.ui.pol_ext_start.setText('')
        self.ui.pol_ext_species.setText('')
        
            
    def clear_all_pars(self):
        self.set_default_params()
        self.select_LApproach_vh()
        self.clear_pert_data()
        self.ui.random_removal.setText('')
        self.ui.output_suffix.setText('')
        self.ui.Comments_text.setPlainText('')
        self.ui.simulation_file.setText('')
        self.ui.inputfile.setText('')
        self.ui.foodweb_checkbox.setChecked(False)
        self.ui.save_output_checkbox.setChecked(False)
        self.ui.Run_Button.setEnabled(0)
        self.ui.clear_all_simparams.setEnabled(0)
        self.ui.save_simulation.setEnabled(0)
        self.update_ui()
        
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
  
    def listspecies2text(self,objlist):
        return(str(objlist).replace('[','').replace(']','').replace('\'',''))  
  
    def select_stored_simulation(self):
        self.simulation_file = QtGui.QFileDialog.getOpenFileName(self, 
                                                          'Open File',
                                        self.dirent+sgGL.SIMFILES_PATH,
                                                          filter='*.sim')
        if len(self.simulation_file) == 0:
            return
        a = self.simulation_file.replace('\\', '/ ').split('/');
        self.simulation_file = a[-1]
        self.ui.simulation_file.setText(self.simulation_file)
        self.ui.simulation_file.setEnabled(False)
        self.ui.Run_Button.setEnabled(1)
        try:
            fh = open(self.dirent+sgGL.SIMFILES_PATH+self.simulation_file, "rb")
        except IOError:
            self.lista_err.append("ERROR: can\'t open file " + \
                                  self.simulation_file)
            self.error_exit()
            return
        else:
            fh.close()
            with open(self.dirent+sgGL.SIMFILES_PATH+self.simulation_file, 
                      'rb') as f:
                storedsim = pickle.load(f)
                fentrada = storedsim.dirtrabajo + '\\' +storedsim.direntrada+\
                           storedsim.filename+'_a.txt'.replace('\\', '/ ')
                self.input_file = fentrada
                self.filename = storedsim.filename+'_a.txt'
                self.ui.inputfile.setText(self.filename)
                numspecies_a, numspecies_b = self.load_file_data(self.filename,
                                                          storedsim.direntrada)
                self.ciclos = storedsim.year_periods
                if (storedsim.eliminarenlaces==0):
                    self.ui.random_removal.setText('')
                else:
                    self.ui.random_removal.setText(str(storedsim.eliminarenlaces))
                if (storedsim.hay_foodweb):
                    self.ui.foodweb_checkbox.setChecked(True)
                if (storedsim.data_save):
                    self.ui.save_output_checkbox.setChecked(True)
                if len(storedsim.output_suff)>0:
                    self.ui.output_suffix.setText(storedsim.output_suff)
                if len(storedsim.com)>0:
                    self.ui.Comments_text.setPlainText(storedsim.com)            
                if (storedsim.algorithm == 'NoMutualism'):
                    self.select_No_mutualism()
                elif (storedsim.algorithm == 'Verhulst'):
                    self.select_LApproach_vh()
                elif (storedsim.algorithm == 'Logistic_abs'):
                    self.select_LApproach_abs()
                elif (storedsim.algorithm == 'Logistic_u'):
                    self.select_LApproach_u()
                elif (storedsim.algorithm == 'May'):
                    self.select_May()
                self.plants_blossom = storedsim.plants_blossom_prob
                self.plants_blossom_type = storedsim.plants_blossom_type
                if (self.plants_blossom_type == 'Gaussian'):
                    self.select_GaussianBlossom()
                if (self.plants_blossom_type == 'Binary'):
                    self.select_BinaryBlossom()
                self.ui.blossom_pert_species.setText(self.listspecies2text(storedsim.blossom_pert_list))
                self.plants_blossom_sd = storedsim.plants_blossom_sd
                self.typeofblossom = storedsim.plants_blossom_type
                self.Bssvar_period = storedsim.Bssvar_data.Bssvar_period
                self.Bssvar_sd = storedsim.Bssvar_data.Bssvar_sd
                self.typeofmodulation = 'None'
                self.Bssvar_modulationtype_list = storedsim.Bssvar_data.Bssvar_modulationtype_list
                self.ui.Bssvar_species.setText(self.listspecies2text(storedsim.Bssvar_data.Bssvar_species))
                if storedsim.Bssvar_data.Bssvar_modulationtype_list[0]=='linear':
                    self.Bssvar_Type_linear_slope = storedsim.Bssvar_data.Bssvar_modulationtype_list[1]
                    self.typeofmodulation = 'linear'
                    self.select_linearModulation()
                if storedsim.Bssvar_data.Bssvar_modulationtype_list[0]=='sin':
                    self.Bssvar_Type_sin_period = storedsim.Bssvar_data.Bssvar_modulationtype_list[1]
                    self.typeofmodulation = 'sin'
                    self.select_sinModulation()
                if storedsim.Bssvar_data.Bssvar_modulationtype_list[0]=='None':
                    self.typeofmodulation = 'None'
                    self.select_NoneModulation()
                if len(storedsim.pl_ext):
                    self.ui.pl_ext_species.setText(self.listspecies2text(storedsim.pl_ext['species']))
                    self.ui.pl_ext_period.setText(str(storedsim.pl_ext['period']//sgGL.DAYS_YEAR))
                    self.ui.pl_ext_spike.setText(str(storedsim.pl_ext['spike']))
                    self.ui.pl_ext_start.setText(str(storedsim.pl_ext['start']))
                    self.ui.pl_ext_rate.setText(str(storedsim.pl_ext['rate']))
                    self.ui.pl_ext_numperiod.setText(str(storedsim.pl_ext['numperiod']))
                if len(storedsim.pol_ext):
                    self.ui.pol_ext_species.setText(self.listspecies2text(storedsim.pol_ext['species']))   
                    self.ui.pol_ext_period.setText(str(storedsim.pol_ext['period']//sgGL.DAYS_YEAR))
                    self.ui.pol_ext_spike.setText(str(storedsim.pol_ext['spike']))
                    self.ui.pol_ext_start.setText(str(storedsim.pol_ext['start']))
                    self.ui.pol_ext_rate.setText(str(storedsim.pol_ext['rate']))
                    self.ui.pol_ext_numperiod.setText(str(storedsim.pol_ext['numperiod']))
                self.haypred = storedsim.hay_foodweb
                self.algorithm = storedsim.algorithm
                self.update_ui()

    def get_root_file_name(self, name):
        return(name.replace('_b.txt', '_a.txt').replace('_c.txt', '_a.txt'))
        
    def load_file_data(self,selected_filename,selected_path):    
        try:
            fh = open(selected_path+selected_filename, "r")
        except IOError:
            self.lista_err.append("ERROR: can\'t open file " + selected_filename)
            self.error_exit()
            return
        else:
            fh.close()
#             self.input_fname_raw = selected_filename.replace('_b.txt', '_a.txt').\
#                                                 replace('_c.txt', '_a.txt')
            self.input_fname_raw = self.get_root_file_name(selected_filename)
            l_minputchar_x = sgcom.dlmreadlike(self.input_fname_raw, selected_path)
            numspecies_a = len(l_minputchar_x) - sgGL.LINES_INFO_MATRIX
            numspecies_b = len(l_minputchar_x[0])
            self.ui.species_number.setText("Plant species: %d  Pollinator species: %d" %\
                             (numspecies_a, numspecies_b))
            self.repaint()
            return(numspecies_a, numspecies_b)
    
    def select_file(self):
        self.filename = QtGui.QFileDialog.getOpenFileName(self, 'Open File',
                                                          self.input_dir)
        if len(self.filename) == 0:
            return
        a = self.filename.replace('\\', '/ ').split('/');
        self.input_file = a[-1]
        self.ui.inputfile.setText(self.input_file)
        self.ui.inputfile.setEnabled(False)
        numspecies_a, numspecies_b = self.load_file_data(self.input_file,self.dirent) 
        self.ui.Run_Button.setEnabled(1)
    
    def save_simulation_file(self):
        """
        Create and show the Save FileDialog
        """
        a = self.ui.inputfile.text()
        a = self.get_root_file_name(a)
        a = a.split('_a.txt')
        output_suffix = self.ui.output_suffix.text()
        simfile_name =  self.input_dir+'/'+sgGL.SIMFILES_PATH + a[0] + '_' +\
             sgcom.create_file_suffix(self.algorithm,output_suffix,self.ciclos)+\
                        '.sim'
        simulation_selected_filename = QtGui.QFileDialog.getSaveFileName(self,
                "Save simulation parameters",
                simfile_name)
        if len(simulation_selected_filename)>0:
            simulation_params.write2file(simulation_selected_filename)
            
    def fetch_pert_data(self, xx_ext_period, xx_ext_spike, xx_ext_start, 
                        xx_ext_rate, xx_ext_numperiod, xx_ext_species):
        ret_extinction = {}

        if len(xx_ext_period.text())>0:
            try:
                ret_extinction['period'] = int( float(xx_ext_period.text())\
                                              * sgGL.DAYS_YEAR)
                ret_extinction['spike'] = float(xx_ext_spike.text())
                ret_extinction['start'] = float(xx_ext_start.text())
                ret_extinction['rate'] = float(xx_ext_rate.text())
                ret_extinction['numperiod'] = int(xx_ext_numperiod.text())
                ret_extinction['species']  = sgcom.create_list_species_affected(xx_ext_species.text())
            except:         
                self.lista_err.append("ERROR: Incorrect forced perturbation format")
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

#     def create_file_suffix(self,output_suffix):
#         return('_'.join([self.algorithm,output_suffix,str(int(self.ciclos))]) )
    
    def add_entry(self):
        global simulation_params
        try:
            os.stat(self.dirs)
        except:
            os.makedirs(self.dirs)
        try:
            self.ciclos = int(self.ui.ciclos.text())
        except:            
            self.lista_err.append("ERROR: Incorrect Years format ")
            self.error_exit()
            return    
        self.ui.label_report.setText("")
        self.ui.Run_Button.setEnabled(0)
        self.ui.Close_Button.setEnabled(0)
        self.ui.Error_msg.setText("")
        displayinic = 0
        self.input_fname_raw = self.input_file.replace('_b.txt',
                                           '_a.txt').replace('_c.txt', '_a.txt')
        aux = self.input_fname_raw.split('/')
        aux = aux[-1].split('_a.txt')
        input_fname = aux[0]
        output_suffix = self.ui.output_suffix.text()
        comentario = self.ui.Comments_text.toPlainText()
        #dirs = os.path.dirname(sys.argv[0])
        dirs = os.getcwd()
        reportpath = os.path.join(dirs, self.dirsal.replace('\\', '/'))
        file_suffix = sgcom.create_file_suffix(self.algorithm,output_suffix,self.ciclos)
        fichr = sgcom.create_fichreport_name(reportpath,input_fname,file_suffix)     
#         fichr = reportpath + 'rep_' + input_fname +'_'+\
#                 self.create_file_suffix(output_suffix)+ '.html'
        dispfichsal = fichr.split('/')
        linkname = '<a href=file:///' + fichr + '>' + dispfichsal[-1] + '</a>'
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
        if (len(self.ui.Bssvar_species.text()) == 0):
            self.Bssvar_species = []
        else:
            self.Bssvar_species = sgcom.create_list_species_affected(self.ui.Bssvar_species.text())
        Blossomvar_data = sgcom.BlossomVariability(self.Bssvar_period,
                                                    self.Bssvar_sd,
                                                self.Bssvar_modulationtype_list,
                                                self.Bssvar_species)
        # External perturbation datacom.
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

        blossom_perturbation = sgcom.create_list_species_affected(self.ui.blossom_pert_species.text())
        simulation_params = sgcom.SimulationConditions(filename = input_fname, 
                year_periods = self.ciclos, 
                hay_foodweb = self.haypred, hay_superpredadores = haysup,
                data_save = outputdatasave, dirtrabajo = dirs, 
                direntrada = self.dirent, dirsal = self.dirsal,
                eliminarenlaces = el, pl_ext = plants_extinction, 
                pol_ext = pols_extinction, output_suff = output_suffix, 
                fichreport = fichr, com = comentario, 
                algorithm = self.algorithm, plants_blossom_prob = pb, 
                plants_blossom_sd = pb_sd, 
                plants_blossom_type = self.typeofblossom, 
                blossom_pert_list=blossom_perturbation,
                release=self.release,
                Bssvar_data = Blossomvar_data )  
        
#         simulation_params.write2file()
        self.ui.Run_Button.setEnabled(1)
        self.ui.Close_Button.setEnabled(1)
        self.ui.URL_report.setText("Running")
        self.repaint()
        sim_ret_val = b_sim.bino_mutual (sim_cond = simulation_params)
        try:
            #sim_ret_val = b_sim.bino_mutual (sim_cond = simulation_params)
            pass
        except:
            self.lista_err.append("Simulation stopped, see details")
            self.error_exit()
            self.repaint() 
        else:   
            self.ui.label_report.setText("Report file: ")
            self.ui.URL_report.setText(linkname)
            self.ui.URL_outputs.setText("<a href='file:///" +\
                                        reportpath.replace('\\', '/') +\
                                        "'>See all results</a>")       
            sggraph.mutual_render(simulation_params, sim_ret_val, displayinic,
                                self.ciclos * sgGL.DAYS_YEAR)
            if self.haypred:
                sggraph.food_render(simulation_params, sim_ret_val, displayinic, 
                                    self.ciclos * sgGL.DAYS_YEAR)
            self.ui.save_simulation.setEnabled(True)
            self.ui.clear_all_simparams.setEnabled(True)

if __name__ == "__main__": 
    app = QtGui.QApplication(sys.argv)
    myapp = StartQT4()
    myapp.show()
    sys.exit(app.exec_())
