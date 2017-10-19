import sys
from PyQt4 import QtCore, QtGui, uic
import danfinitions as dan
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import curve_fit
import math
import pandas as pd

qtCreatorFile = "gui_cats_v1.ui" # Enter file here.
 
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
 
class MyApp(QtGui.QMainWindow, Ui_MainWindow):

    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        
        self.btnTopDownAll.clicked.connect(self.top_down_show_data)
        self.btnPlotFracSweep.clicked.connect(self.plot_errorbar_things)
        self.btnRunFracSweep.clicked.connect(self.calc_confidence_intervals)
        self.btnMake_DCM.clicked.connect(self.similarity_matrix)
        self.btnPlot_DCM.clicked.connect(self.show_similarity_matrix)
        self.btnCheckFiles.clicked.connect(self.get_preface)
        self.btnLoadFiles.clicked.connect(self.load_file_set)
        self.btnLoadFiles.setEnabled(False)
        self.btnLoadMom.clicked.connect(self.load_mother_file)
        self.btnLoadDad.clicked.connect(self.load_father_file)
        self.btnPlotThisRun.clicked.connect(self.plot_single_run)
        self.btnPlotParents.clicked.connect(self.plot_parents)
        self.btnCalcScores.clicked.connect(self.calculate_scores_better)
        self.btnPlotAllScores.clicked.connect(self.plot_all_scores)
        self.btnPlotParents.setEnabled(False)
        self.btnCalcScores.setEnabled(False)
        self.btnPlotAllScores.setEnabled(False)
        self.btnPlotThisRun.setEnabled(False)
        self.btnConstRrange.setEnabled(False)
        self.btnExtendRrange.setEnabled(False)
        self.btnPlotCompareRuns.setEnabled(False)
        
        self.btnPlotFTM.clicked.connect(self.show_feature_track_map)
        self.btnPlot_DCM.setEnabled(False)
        self.btnMake_DCM.setEnabled(False)
        self.radioButton_useGr.clicked.connect(self.update_axis_labels)
        self.radioButton_useId.clicked.connect(self.update_axis_labels)
        self.radioButton_useSq.clicked.connect(self.update_axis_labels)
        
        
        self.actionSave_All_Data.triggered.connect(self.save_alldata_one_file)
        self.actionSave_All_Data_in_Range_3D.triggered.connect(self.save_alldata_inRange_one_file)
        self.action2_Column_r_difference.triggered.connect(lambda: self.save_parent_diff(2))
        self.actionOne_file_r_mom_r_dad_r_difference.triggered.connect(lambda: self.save_parent_diff(4))
        self.actionAll_Fits_OneBigFile.triggered.connect(self.save_all_fits)
        self.actionMisfit_to_data_One_big_file_3D.triggered.connect(self.save_misfit_onebigfile)
        self.actionEach_fit_2Col.triggered.connect(lambda: self.save_each_fit(2))
        self.actionEach_fit_3Col.triggered.connect(lambda: self.save_each_fit(3))
        self.btnConstRrange.clicked.connect(lambda: self.rdep_fit_const_rrange(1))
        self.btnExtendRrange.clicked.connect(lambda: self.rdep_fit_const_rrange(2))
        self.actionSave_RtoR_inRng.triggered.connect(self.save_alldifference_inRange_one_file)
        self.actionReadTemperatureList.triggered.connect(self.read_in_run_variable)
        self.actionChoice_vs_Run_Number2col.triggered.connect(lambda: self.save_choice_vs_run_number(2))
        self.actionChoice_vs_Run_Number3col.triggered.connect(lambda: self.save_choice_vs_run_number(3))
        self.actionChoice_vs_Run_Number4col.triggered.connect(lambda: self.save_choice_vs_run_number(4))
        
 
        self.btnPlotCompareRuns.clicked.connect(self.PlotTwoRuns)
        #self.btnPlotNormError.clicked.connect(self.plot_normalized_error_for_run)
        
        self.labAlarm.setText("   ")
        self.labWarning.setText("")

        #point to checks on Box Car Mode checkbox
        self.checkBox_boxCar.stateChanged.connect(self.possible_warning)
        
        self.checkBox_defCMap.stateChanged.connect(self.cmap_settings_update)
        self.checkBox_ft_Int.pressed.connect(self.undo_ft_momanddad)
        
        self.checkBox_ft_upfMom.pressed.connect(self.undo_ft_int_check)
        self.checkBox_ft_upfDad.pressed.connect(self.undo_ft_int_check)

        #point to checks on spinboxes
        self.doubleSpinBox_Rmax.valueChanged.connect(self.possible_warning)
        self.doubleSpinBox_Rmin.valueChanged.connect(self.possible_warning)
        self.spinBox_RunMin.valueChanged.connect(self.possible_warning)
        self.spinBox_RunMax.valueChanged.connect(self.possible_warning)
        self.doubleSpinBox_deltaR.valueChanged.connect(self.possible_warning)
        self.doubleSpinBox_Rtop.valueChanged.connect(self.possible_warning)
        
        #disable all menu items until conditions are met
        self.actionSave_All_Data.setEnabled(False)
        self.actionSave_All_Data_in_Range_3D.setEnabled(False)
        self.menuSave_Parent_Differences.setEnabled(False)
        self.menuSave_Fits.setEnabled(False)
        self.menuSave_Fitted_Choice.setEnabled(False)
        self.actionSave_RtoR_inRng.setEnabled(False)
        
        #setin up variables like a BOSS
        self.file_preface = "GST"
        self.number_of_runs = 0
        self.rdata = []
        self.grdata = []
        self.fileList = []
        self.mom = []
        self.dad = []
        self.df_score = pd.DataFrame()
        self.run_variables = []
        self.df_all_data = []
        self.conf_int_stuff =[]
        
        self.recalc_params = np.zeros(5)
        
        self.use_xaxis_label = 'r ($\AA$)'
        self.use_yaxis_label = 'G(r)'
        
        self.rmax = 10000.0
        self.rmin = 0.0
        self.MomExists = False
        self.DadExists = False
        self.DataExists = False
        self.have_done_a_calc = False
        self.runMin = 0
        self.runMax = 0
        self.deltaR = 0
        self.rTop = 0
        self.deadZone = 0     

        self.datamin = 0.0
        self.datamax = 1.0
        
    def save_choice_vs_run_number(self,mode):
        self.labAlarm.setText(" ! ")
        self.notice_spinbox_values()

        if len(self.run_variables) > 0:
            x = self.run_variables[self.runMin:self.runMax+1]
        else:
            x = np.arange(self.runMin,self.runMax+1)# len(self.fit_dad_fraction))
            
        max_res = np.nanmax(self.residuals)
        norm_res = []
        for i in range(len(self.residuals)):
            norm_res.append(self.residuals[i]/max_res)
            
        #now for writing out
        default_name = "choice_vs_run_number_"+str(mode)+"col.dat"
        file = QtGui.QFileDialog.getSaveFileName(self, "Save File", default_name, ".dat")
        file_created_info = "# This file contains scores of runs based on linear combination of two parents.\n"
        
        if mode == 2:
            if len(self.run_variables) > 0:
                file_created_info += "# 2 Column format with mom-fraction vs. temperatures\n"
            else:
                file_created_info += "# 2 Column format with mom-fraction vs. run number\n"
        elif mode == 3:
            if len(self.run_variables) > 0:
                file_created_info += "# 3 Column format with mom-fraction & dad-fraction vs. temperatures\n"
            else:
                file_created_info += "# 3 Column format with mom-fraction & dad-fraction vs. run number\n"        
        elif mode == 4:
            if len(self.run_variables) > 0:
                file_created_info += "# 4 Column format with mom-fraction & dad-fraction & normalized residuals vs. temperatures\n"
                file_created_info += "# Residuals were normalized (divided) by a factor of "+str(max_res)+"\n"
            else:
                file_created_info += "# 4 Column format with mom-fraction & dad-fraction & normalized residuals vs. run number\n"   
                file_created_info += "# Residuals were normalized (divided) by a factor of "+str(max_res)+"\n"
                
        
        file_created_info += "# runs between "+str(self.runMin)+" and "+str(self.runMax)+"\n"
        file_created_info += "# plotting in R between "+str(self.rmin)+" and "+str(self.rmax)+"\n"
        
        file_created_info += "# Mother file (value 1): "+str(self.labMomFile.text())+"\n"
        file_created_info += "# Father file (value -1): "+str(self.labDadFile.text())+"\n"

        if self.recalc_params[4] == 0 : # calculated with different parameters
            file_created_info += "# Note: Fits were performed on a different range than shown in these plots.\n"
            if self.recalc_params[0] != self.rmin:
                file_created_info += "# Rmin used in fits was "+str(self.recalc_params[0])+"\n"
            if self.recalc_params[1] != self.rmax:
                file_created_info += "# Rmax used in fits was "+str(self.recalc_params[1])+"\n"
            if self.recalc_params[2] != self.runMin:
                file_created_info += "# Run_Min (lowest run) used in fitting was "+str(int(self.recalc_params[2]))+"\n"
            if self.recalc_params[3] != self.runMax:
                file_created_info += "# Run_Max (highest run) used in fitting was "+str(int(self.recalc_params[3]))+"\n"        

        date_is = time.strftime("%m/%d/%Y")
        time_is = time.strftime("%H:%M:%S")
        file_created_info += "# This data generated with CATS on "+str(date_is)+" at "+str(time_is)+".\n"
        file_created_info += "\n"
        outfile = open(file,'w')
        outfile.write(file_created_info)
        
        for runs in np.arange(self.runMin,self.runMax+1):
            this_mom_score = self.fit_mom_fraction[runs] 
            this_dad_score = self.fit_dad_fraction[runs]
            if mode == 2:
                outfile.write(str(x[runs])+" "+str(this_mom_score)+"\n")
            elif mode == 3:
                outfile.write(str(x[runs])+" "+str(this_mom_score)+" "+str(this_dad_score)+"\n")
            elif mode == 4:
                outfile.write(str(x[runs])+" "+str(this_mom_score)+" "+str(this_dad_score)+" "+str(norm_res[runs])+"\n")
        outfile.close()
        print "saved ok"

        self.labAlarm.setText("   ")   

    def cmap_settings_update(self):
        print "doing things!!!"
        if self.checkBox_defCMap.isChecked():
            print "restoring defaults"
            self.labRunNumber_4.setEnabled(False)
            self.PTE_cmap_choice.setEnabled(False)
            self.labRunNumber_5.setEnabled(False)
            self.doubleSpinBox_cmin.setEnabled(False)
            self.labRunNumber_6.setEnabled(False)
            self.doubleSpinBox_cmax.setEnabled(False)
            self.checkBox_showCbar.setEnabled(False)
            self.PTE_cmap_choice.setEnabled(False)
            self.PTE_cmap_choice.setPlainText('viridis')
            #self.doubleSpinBox_cmin.setValue(1.0)
            #self.doubleSpinBox_cmax.setValue(1.0)

        else:
            print "enabling controls"
            self.labRunNumber_4.setEnabled(True)
            self.PTE_cmap_choice.setEnabled(True)
            self.labRunNumber_5.setEnabled(True)
            self.doubleSpinBox_cmin.setEnabled(True)
            self.labRunNumber_6.setEnabled(True)
            self.doubleSpinBox_cmax.setEnabled(True)
            self.checkBox_showCbar.setEnabled(True)
            
            self.doubleSpinBox_cmin.setValue(self.datamin)
            self.doubleSpinBox_cmax.setValue(self.datamax)

        
    def possible_warning(self):
        #print "possible warning!!!"
        #for giving warning if parameters in spinboxes are different than those used to calculate fits recently
        self.notice_spinbox_values()
        still_ok = True
        if self.rmin != self.recalc_params[0] or self.rmax != self.recalc_params[1] or self.runMin != self.recalc_params[2] or self.runMax != self.recalc_params[3]:
            still_ok = False
            if self.have_done_a_calc:
                if still_ok == False:
                    self.labWarning.setText("You should probably recalculate ->")
                    self.recalc_params[4] = 0 # 0 shows a recalc has not been performed
        
        
        
        #for checking if boxcar mode is activated
        if self.checkBox_boxCar.isChecked():
            self.deltaR = self.rmax - self.rmin
            self.doubleSpinBox_deltaR.setValue(self.deltaR)
            self.doubleSpinBox_deltaR.setEnabled(False)
            self.doubleSpinBox_Rtop.setSingleStep(self.deltaR)
            self.doubleSpinBox_Rtop.setMinimum(self.rmax)
        else:
            self.deltaR = self.doubleSpinBox_deltaR.value()
            self.doubleSpinBox_deltaR.setEnabled(True)
            self.doubleSpinBox_Rtop.setMaximum(np.amax(self.rdata[0]))
            self.doubleSpinBox_Rtop.setSingleStep(1.0)
            self.doubleSpinBox_Rtop.setMinimum(self.rmax)
    
        self.deltaR = self.doubleSpinBox_deltaR.value()
        num_boxes = 1 + int(math.floor((self.doubleSpinBox_Rtop.value() - self.doubleSpinBox_Rmax.value()) / self.doubleSpinBox_deltaR.value()))
        that_val = str(self.use_xaxis_label)[0]
        self.labBoxWidth.setText("There will be "+str(num_boxes)+" boxes of width "+str(self.deltaR)+"\n each, up to maximum of "+that_val+" = "+str(num_boxes*self.deltaR))
    
    def undo_ft_momanddad(self):
        self.checkBox_ft_upfMom.setChecked(False)
        self.checkBox_ft_upfDad.setChecked(False)
            
    def undo_ft_int_check(self):    
        if self.checkBox_ft_Int.isChecked():
            self.checkBox_ft_Int.setChecked(False)
    
    def calculate_scores_better(self):
        self.notice_spinbox_values()
        self.labWarning.setText("")
        self.have_done_a_calc = True
        self.recalc_params[0] = self.rmin
        self.recalc_params[1] = self.rmax
        self.recalc_params[2] = self.runMin
        self.recalc_params[3] = self.runMax
        self.recalc_params[4] = 1 #set this flag to 1 to show recalc HAS been performed
        #make 'mini-mom' and 'mini-dad' for lsf procedure, which have dimensions rmin to rmax
        mini_mom = []
        mini_dad = []
        mini_r = []
        mini_data = []
        for i in range(len(self.rdata[0])):
            if self.rdata[0][i] >= self.rmin and self.rdata[0][i] <= self.rmax:
                mini_mom.append(self.mom[0][i])
                mini_dad.append(self.dad[0][i])
                mini_r.append(self.rdata[0][i])
        #have mini of parents and r, now go over each run to be used        
        for j in range(self.runMin,self.runMax+1):
            this_mini_data = []
            for i in range(len(self.rdata[0])):
                if self.rdata[0][i] >= self.rmin and self.rdata[0][i] <= self.rmax:
                    this_mini_data.append(self.grdata[j][i])
            mini_data.append(this_mini_data)
        #now that I have all mini-data, can do least_squares on each one
        
        #ok, now that we have mini-datasets, score them in 'mini-score'
        rng_was = self.runMax - self.runMin + 1
        
        def combo_data_simple(this_mini_r,p1):
            pmom = np.exp(-(p1**2))
            pdad = 1-pmom
            return pmom*this_mini_mom[:]+pdad*this_mini_dad[:]

        #loop over each dataset
		
        self.fit_mom_fraction = np.zeros(len(self.rdata))
        self.fit_dad_fraction = np.zeros(len(self.rdata))
        self.residuals = np.zeros(len(self.rdata))
        #offset = len(self.rdata) - rng_was
        is_on = 0
        #print "length of mini_data is: "+str(len(mini_data))
        #print "range to loop over up to "+str(len(self.rdata))
        #print "run min/max is: "+str(self.runMin)+" "+str(self.runMax)
        for this_run in range(len(self.rdata)):
#        for this_run in range(0,rng_was):
            #progress bar should probably go here
            completed_per = 100.0 * float(this_run)/float(rng_was-2)
            self.progressBar.setValue(completed_per)
            if this_run <= self.runMax and this_run >= self.runMin:
                lsf_results = np.zeros([2,len(mini_r)])
                lsf_residuals = []            
                
                this_mini_data = np.array(mini_data[is_on])
                this_mini_r = np.array(mini_r)
                this_mini_dad = np.array(mini_dad)
                this_mini_mom = np.array(mini_mom)
                is_on += 1
                try:
                    popt, pcov = curve_fit(combo_data_simple, this_mini_r, this_mini_data, p0 = (.5), maxfev=1200)
                    
                except:
                    print "Something went wrong, sorry about that."
                    popt = 0.0
                    break
                res = this_mini_data[:] - combo_data_simple(mini_r,popt[0])
                fres = sum(res**2)
                self.fit_mom_fraction[this_run]=(np.exp(-(popt[0]**2)))
                self.fit_dad_fraction[this_run]=(1.0 - (np.exp(-(popt[0]**2))))
                self.residuals[this_run]=(fres)


                self.btnPlotAllScores.setEnabled(True)
                self.checkBox_overlay_mom.setEnabled(True)
                self.checkBox_overlay_dad.setEnabled(True)
                self.checkBox_ft_upfMom.setEnabled(True)
                self.checkBox_ft_upfDad.setEnabled(True)
                
        self.menuSave_Fits.setEnabled(True)
        #print "scored!"
        self.menuSave_Fitted_Choice.setEnabled(True)
               
    def load_mother_file(self):
        #open file explorer
        name = QtGui.QFileDialog.getOpenFileName(self, "First Target Dataset (mother)")
        self.mom = []
        this_list = []
        this_list.append(name)
        momr,momgr = dan.read_in_all_runs_smart_list(this_list)
        self.mom.append(momgr[0])
        print "mom is loaded"
        self.labMomFile.setText(str(os.path.basename(str(name))))
        self.MomExists = True
        if self.DadExists:
            self.btnPlotParents.setEnabled(True)
            if self.DataExists:
                self.btnCalcScores.setEnabled(True)
                self.menuSave_Parent_Differences.setEnabled(True)
                self.btnConstRrange.setEnabled(True)
                self.btnExtendRrange.setEnabled(True)
                self.btnRunFracSweep.setEnabled(True)


        if self.DataExists == False:
            self.rdata.append(momr[0])

    def load_father_file(self):
        #open file explorer
        name = QtGui.QFileDialog.getOpenFileName(self, "Second Target Dataset (father)")
        self.dad = []
        this_list = []
        this_list.append(name)
        dadr,dadgr = dan.read_in_all_runs_smart_list(this_list)
        self.dad.append(dadgr[0])
        print "dad is loaded"
        self.labDadFile.setText(str(os.path.basename(str(name))))
        self.DadExists = True
        if self.MomExists:
            self.btnPlotParents.setEnabled(True)
            self.menuSave_Parent_Differences.setEnabled(True)
            if self.DataExists:
                self.btnCalcScores.setEnabled(True)
                self.btnConstRrange.setEnabled(True)
                self.btnExtendRrange.setEnabled(True)        
                self.btnRunFracSweep.setEnabled(True)

    def load_file_set(self):
        self.rdata = []
        self.labAlarm.setText(" ! ")    

#    self.rmax = float(self.doubleSpinBox_Rmax.value())
#    	print "loading files with rmax :"+str(self.rmax)
    	self.rdata, self.grdata = dan.read_in_all_runs_smart_list(self.fileList)
        #self.labSmile1.setText(" :)")
        self.labAlarm.setText("   ")    
        self.datamin = np.array(self.grdata).min().min()
        self.datamax = np.array(self.grdata).max().max()

        print "data min was "+str(self.datamin)
        
    	print "loaded "+str(len(self.grdata))
    	print "saw an rmax of "+str(np.amax(self.rdata[0]))
        self.doubleSpinBox_Rmax.setValue(np.amax(self.rdata[0]))    
        self.doubleSpinBox_Rmax.setMaximum(np.amax(self.rdata[0]))
        self.doubleSpinBox_Rmax.setMinimum(np.amin(self.rdata[0]))
        self.doubleSpinBox_Rmax.setSingleStep(np.amin(self.rdata[0][3]-self.rdata[0][2]))
    
        self.doubleSpinBox_Rmin.setValue(np.amax(self.rdata[0][0]))    
        self.doubleSpinBox_Rmin.setMaximum(np.amax(self.rdata[0]))
        self.doubleSpinBox_Rmin.setMinimum(np.amin(self.rdata[0]))
        self.doubleSpinBox_Rmin.setSingleStep(np.amin(self.rdata[0][3]-self.rdata[0][2]))
        
        self.spinBox_RunMin.setMinimum(0)
        self.spinBox_RunMax.setMinimum(0)
        self.spinBox_RunMax.setMaximum(len(self.rdata)-1)
        self.spinBox_RunMin.setMaximum(len(self.rdata)-1)

        self.spinBox_RunMin.setValue(0)
        self.spinBox_RunMax.setValue(len(self.grdata)-1)
        
        self.spinBox_RunNumber.setValue(0)
        self.spinBox_RunNumber.setMaximum(len(self.rdata)-1)
        self.spinBox_RunNumber.setMinimum(0)
 
        self.spinBox_RunNumber_2.setValue(1)
        self.spinBox_RunNumber_2.setMaximum(len(self.rdata)-1)
        self.spinBox_RunNumber_2.setMinimum(0) 
        
        self.doubleSpinBox_deltaR.setValue(1.0)
        self.doubleSpinBox_deltaR.setMinimum(np.amin(self.rdata[0]))
        self.doubleSpinBox_deltaR.setMaximum(np.amax(self.rdata[0]))
        
        #print "current value is : "+str(self.doubleSpinBox_Rtop.value())
        #print "want it to be : "+str(np.amax(self.rdata[0]))
        self.doubleSpinBox_Rtop.setMaximum(np.amax(self.rdata[0]))
        self.doubleSpinBox_Rtop.setValue(np.amax(self.rdata[0]))
        self.doubleSpinBox_Rtop.setMinimum(np.amin(self.rdata[0]))
        #print "but is now : "+str(self.doubleSpinBox_Rtop.value())
        
        self.DataExists = True
        if self.MomExists and self.DadExists:
            self.btnCalcScores.setEnabled(True)
            self.btnRunFracSweep.setEnabled(True)

        self.btnPlotCompareRuns.setEnabled(True)
        self.btnPlotThisRun.setEnabled(True)
        self.actionSave_All_Data.setEnabled(True)
        self.actionSave_All_Data_in_Range_3D.setEnabled(True)
        self.actionSave_RtoR_inRng.setEnabled(True)
        self.btnMake_DCM.setEnabled(True)
        self.checkBox_ft_Int.setEnabled(True)
        self.btnPlotFTM.setEnabled(True)
    
        self.labRmin_2.setEnabled(True)
        self.labRmax_2.setEnabled(True)
        self.doubleSpinBox_FTRmin.setEnabled(True)
        self.doubleSpinBox_FTRmax.setEnabled(True)

        self.btnTopDownAll.setEnabled(True)
        self.checkBox_useTDlimits.setEnabled(True)
        self.checkBox_overlayAvg.setEnabled(True)
        self.labRunNumber_7.setEnabled(True)
        self.PTE_overlay_color.setEnabled(True)
        self.doubleSpinBox_OLayOffset.setEnabled(True)
        self.doubleSpinBox_OLayScale.setEnabled(True)
        self.label_22.setEnabled(True)
        self.label_23.setEnabled(True)
        
        self.df_all_data = pd.DataFrame(index=self.rdata[0])
        for this_run in range(len(self.rdata)):
            self.df_all_data[this_run] = np.array(self.grdata[this_run])  
        self.doubleSpinBox_FTRmax.setValue(self.df_all_data.index[-1]) 
        self.doubleSpinBox_FTRmin.setValue(self.df_all_data.index[0])
        #self.doubleSpinBox_FTRmin.setMini(self.df_all_data.index[0])
           
        #set good initial guesses for overlay offsets and scales
        datamean = np.array(self.df_all_data.loc[:,:].mean(axis=1))
        yavg_max = max(datamean)
        yavg_min = min(datamean)
        num_data = (len(self.grdata))
        
        yspan = float(yavg_max - yavg_min)
        yspan = yspan/float(num_data)
        
        self.doubleSpinBox_OLayOffset.setValue(float(num_data/2))
        self.doubleSpinBox_OLayScale.setValue(yspan)
   
    def get_preface(self):
    	self.file_preface = str(self.plainTextEdit.text())
    	#now count how many files exist with that preface
    	myPath = os.getcwd()
    	self.fileList = glob.glob1(myPath,self.file_preface+"*")
    	self.number_of_runs = len(self.fileList)
    	self.labNumFiles.setText("I can see "+str(self.number_of_runs)+ " files with that preface")
        self.fileList.sort()	
	flist_view = open('file_list_used.txt','w')
	for i in range(len(self.fileList)):
		flist_view.write(str(self.fileList[i])+"\n")
	flist_view.close()

	if self.number_of_runs > 0:
    		self.btnLoadFiles.setEnabled(True)
    	else:
    		self.btnLoadFiles.setEnabled(False)
			
    def PlotTwoRuns(self):
        first_run = int(self.spinBox_RunNumber.value())
        second_run = int(self.spinBox_RunNumber_2.value())
        diff = np.subtract(self.grdata[second_run],self.grdata[first_run])
        plt.close("all")
        plt.figure()
        if self.checkBox_firstRun.isChecked():
            plt.plot(self.rdata[first_run],self.grdata[first_run])
        if self.checkBox_secondRun.isChecked():
            plt.plot(self.rdata[second_run],self.grdata[second_run])
        if self.checkBox_diff.isChecked():
            offset = self.doubleSpinBox_diffOffset.value()
            plt.plot(self.rdata[first_run],diff+offset)
        plt.xlabel(self.use_xaxis_label)
        plt.ylabel(self.use_yaxis_label)
        plt.show()
		
    def plot_normalized_error_for_run(self):
        run_index = int(self.spinBox_RunNumber.value())
        tot_error = 0.0
        norm_error = []
        for i in range(len(self.rdata[0])):
            this_diff = (self.fit_mom_fraction[run_index]*self.mom[0][i]+self.fit_dad_fraction[run_index]*self.dad[0][i]) - self.grdata[run_index][i]
            tot_error += (this_diff)**2
            norm_error.append(tot_error/float(i+1))
        plt.close("all")
        plt.figure()
        plt.xlim([self.doubleSpinBox_Rmin.value(),self.doubleSpinBox_Rmax.value()])
        plt.plot(self.rdata[run_index],norm_error)
        #plt.plot(self.rdata[run_index],self.grdata[run_index])
        if self.checkBox_compareNormError.isChecked():
            run_index2 = int(self.spinBox_RunNumber_2.value())
            tot_error = 0.0
            norm_error2 = []
            for i in range(len(self.rdata[0])):
                this_diff = (self.fit_mom_fraction[run_index2]*self.mom[0][i]+self.fit_dad_fraction[run_index2]*self.dad[0][i]) - self.grdata[run_index2][i]
                tot_error += (this_diff)**2
                norm_error2.append(tot_error/float(i+1))
            plt.plot(self.rdata[run_index2],norm_error2) 
        
        plt.xlabel(self.use_xaxis_label)
        plt.ylabel(self.use_yaxis_label)    
        plt.show()
           
    def plot_single_run(self):
        run_index = int(self.spinBox_RunNumber.value())
        plt.close("all")
        plt.figure()
        plt.plot(self.rdata[run_index],self.grdata[run_index])
        plt.xlabel(self.use_xaxis_label)
        plt.ylabel(self.use_yaxis_label)
        plt.show()
        
    def plot_all_scores(self):
        self.notice_spinbox_values()
        plt.close("all")
        plt.figure()
        if len(self.run_variables) > 0:
            x = self.run_variables[self.runMin:self.runMax+1]
        else:
            x = np.arange(self.runMin,self.runMax+1)# len(self.fit_dad_fraction))
        max_res = np.nanmax(self.residuals)
        norm_res = []
        for i in range(len(self.residuals)):
            norm_res.append(self.residuals[i]/max_res)
        plt.plot(x,self.fit_dad_fraction[self.runMin:self.runMax+1],x,self.fit_mom_fraction[self.runMin:self.runMax+1],x,norm_res[self.runMin:self.runMax+1])
        plt.xlabel(self.use_xaxis_label)
        plt.ylabel('$\phi$ / Residual')
        plt.show()

    def plot_parents(self):   
        plt.close("all")
        plt.figure()
        plt.plot(self.rdata[0],self.mom[0],self.rdata[0],self.dad[0])
        plt.xlabel(self.use_xaxis_label)
        plt.ylabel(self.use_yaxis_label)
        plt.show()
        
    def say_hi(self):
    	print "oh hey there"
    
    def rdep_fit_const_rrange(self,mode):
        # mode = 1 -> box-car style (though variable spacing)
        # mode = 2 -> draw out R-max each time
        self.deltaR = self.doubleSpinBox_deltaR.value()
        self.rTop = self.doubleSpinBox_Rtop.value()
        #self.deadZone = self.doubleSpinBox_deadWidth.value()
        self.notice_spinbox_values()
        if self.checkBox_boxCar.isChecked():
            self.deltaR = self.rmax - self.rmin
            self.doubleSpinBox_deltaR.setValue(self.deltaR)
            
        #figure out how many boxcars I need
        num_boxes = int(math.ceil((self.rTop - self.rmax) / self.deltaR))
        trueRtop = num_boxes * self.deltaR + self.rmax
        print "num_boxes start at "+str(num_boxes)+" "+str(trueRtop)
        while trueRtop > self.rTop:
            num_boxes -= 1
            trueRtop = num_boxes * self.deltaR + self.rmax
            print "dropped to "+str(num_boxes)+" "+str(trueRtop)
        print "setteled on this "+str(num_boxes)+" "+str(trueRtop)
        print "allocating a big'ol thingy"
        
        #lets keep track of each individual rmax
        setOfRmaxs = np.zeros(num_boxes+1)
        for i in range(len(setOfRmaxs)):
            setOfRmaxs[i] = i*self.deltaR + self.rmax
        print "first Rmax will be "+str(setOfRmaxs[0])
        print "last Rmax will be "+str(setOfRmaxs[len(setOfRmaxs)-1])
        
        fitsRdep_momfrac_constRrng = np.zeros((len(setOfRmaxs),len(self.grdata)))
        resRdep_constRrng = np.zeros((len(setOfRmaxs),len(self.grdata)))
        #now for the big calculation loop
        for i in range(len(setOfRmaxs)): #for a given rmax
            completed_per2 = 100.0 * float(i)/(float(len(setOfRmaxs)))
            self.progressBar_2.setValue(completed_per2)
            this_rTop = setOfRmaxs[i]
            this_rBottom = setOfRmaxs[i] - (self.rmax- self.rmin)
            if mode == 2: #extending
                this_rBottom = self.rmin
                
            print "going over "+str(this_rBottom)+" "+str(this_rTop)
#            for j in len(self.grdata): #for a specific dataset
            
            #make 'mini-mom' and 'mini-dad' for lsf procedure, which have dimensions rmin to rmax
            mini_mom = np.zeros(len(self.mom[0]))
            mini_dad = np.zeros(len(self.dad[0]))
            mini_r = np.zeros(len(self.rdata[0]))
            #mini_data = np.zeros(len(self.grdata[0]))
            mini_data = []
            print "making tiny parents"
            for r in range(len(self.rdata[0])):
                if self.rdata[0][r] >= this_rBottom and self.rdata[0][r] <= this_rTop:
                    mini_mom[r] = (self.mom[0][r])
                    mini_dad[r] = (self.dad[0][r])
                mini_r[r]=(self.rdata[0][r]) #can have every rvalue in the array
            #have mini of parents and r, now go over each run to be used        
            for run in range(self.runMin,self.runMax+1): #have each run, make minidata
                this_mini_data = np.zeros(len(self.grdata[0]))
                for r in range(len(self.rdata[0])):
                    if self.rdata[0][r] >= this_rBottom and self.rdata[0][r] <= this_rTop:
                        this_mini_data[r]=(self.grdata[run][r])
                mini_data.append(this_mini_data)
            #now that I have all mini-data, can do least_squares on each one           
            this_momFrac = np.zeros(len(self.rdata))
            this_dadFrac = np.zeros(len(self.rdata))
            
            def combo_data_simple(this_mini_r,p1):
                pmom = np.exp(-(p1**2))
                pdad = 1-pmom
                return pmom*this_mini_mom[:]+pdad*this_mini_dad[:]
                
            this_residuals = np.zeros(len(self.rdata))
            is_on = 0
            for this_run in range(len(self.rdata)):
                completed_per = 100.0 * float(this_run)/(float(self.runMax - self.runMin -1))
                self.progressBar.setValue(completed_per)
                #now decide if this run is between runMin and runMax
                if this_run <= self.runMax and this_run >= self.runMin:
                    lsf_results = np.zeros([2,len(mini_r)])
                    lsf_residuals = []            
                
                    this_mini_data = np.array(mini_data[is_on])
                    this_mini_r = np.array(mini_r)
                    this_mini_dad = np.array(mini_dad)
                    this_mini_mom = np.array(mini_mom)
                    is_on += 1
                    try:
                        popt, pcov = curve_fit(combo_data_simple, this_mini_r, this_mini_data, p0 = 0.5, maxfev = 1200)
                    except:
                        print "I'm just going to quit, this isn't working"
                        popt = 0
                        break
                    res = this_mini_data[:] - combo_data_simple(mini_r,popt[0])
                    fres = sum(res**2)
                    this_momFrac[this_run]=(np.exp(-(popt[0]**2.0))) #shorten to later line?
                    this_dadFrac[this_run]=(1.0 - (np.exp(-(popt[0]**2)))) #remove?
                    this_residuals[this_run]=(fres) #line maybe unnedded
                    
                    fitsRdep_momfrac_constRrng[i][this_run] = (np.exp(-(popt[0]**2.0))) #this_momFrac[this_run]
                    resRdep_constRrng[i][this_run] = (fres)
                    
        self.progressBar_2.setValue(100.0)

    	myPath = os.getcwd()    
        if mode == 1:
            default_name = "box_car_"
        elif mode == 2:
            default_name = "extending_rmax_"
        file_pref = QtGui.QFileDialog.getSaveFileName(self, "File preface?", default_name, "*")
        
        print "I'm just gonna print a couple of big'ol files for you"
        outfile1 = open(str(file_pref)+"choices.dat",'w')
        outfile2 = open(str(file_pref)+"residuals.dat",'w')
        for i in range(len(fitsRdep_momfrac_constRrng)): #Rmax choice
            for j in range(len(fitsRdep_momfrac_constRrng[0])): #run_number choice
                if j >= self.runMin and j <= self.runMax:
                    choice_was = ((1.0-fitsRdep_momfrac_constRrng[i][j]) - (fitsRdep_momfrac_constRrng[i][j]))                
                    if len(self.run_variables) > 0:
                        outfile1.write(str(setOfRmaxs[i])+" "+str(self.run_variables[j])+" "+str(choice_was)+"\n")
                        outfile2.write(str(setOfRmaxs[i])+" "+str(self.run_variables[j])+" "+str(resRdep_constRrng[i][j])+"\n")                    
                    else:
                        outfile1.write(str(setOfRmaxs[i])+" "+str(j)+" "+str(choice_was)+"\n")
                        outfile2.write(str(setOfRmaxs[i])+" "+str(j)+" "+str(resRdep_constRrng[i][j])+"\n")
            outfile1.write("\n")
            outfile2.write("\n")
        outfile1.close()
        outfile2.close()
        
        
        print "all done"
	
    def update_axis_labels(self):
        if self.radioButton_useGr.isChecked():
            print "G(r)-language"
            self.use_xaxis_label = 'r ($\AA$)'
            self.use_yaxis_label = 'G(r)'
            self.labRmin.setText('Rmin')
            self.labRmax.setText('Rmax')
            self.label_4.setText('deltaR')
            self.label_6.setText('Rtop')
            self.label_3.setText('R-range dependent Parent Fraction Series')
            self.label_5.setText('   Box-Car range : Starting with Rmin and Rmax on LHS,\nincrease Rmax and Rmin by deltaR up to the Rmax = Rtop')
            self.label_8.setText(' Extending R-range : Starting with Rmin/ max in Fit Parameters\n    increase Rmax by a value of Delta R each time up to Rtop')
            self.checkBox_useR_for_Fsweep.setText('Use Rmin/Rmax')
            self.labRmin_2.setText('R-local min')
            self.labRmax_2.setText('R-local max')
            self.btnConstRrange.setText('Constant R-range Width')
            self.btnExtendRrange.setText('Extending R-range Width')     
            self.checkBox_useTDlimits.setText('Use Limits of Rmin/Rmax')
          
            
        if self.radioButton_useId.isChecked():
            print "I(d)-language"
            self.use_xaxis_label = 'd ($\AA$)'
            self.use_yaxis_label = 'I(d)'
            self.labRmin.setText('dmin')
            self.labRmax.setText('dmax')
            self.label_4.setText('deltad')
            self.label_6.setText('dtop')
            self.label_3.setText('d-range dependent Parent Fraction Series')
            self.label_5.setText('   Box-Car range : Starting with dmin and dmax on LHS,\nincrease dmax and dmin by deltad up to the dmax = dtop')
            self.label_8.setText(' Extending d-range : Starting with dmin/ max in Fit Parameters\n    increase dmax by a value of Delta d each time up to dtop')
            self.checkBox_useR_for_Fsweep.setText('Use dmin/dmax')
            self.labRmin_2.setText('d-local min')
            self.labRmax_2.setText('d-local max')
            self.btnConstRrange.setText('Constant d-range Width')
            self.btnExtendRrange.setText('Extending d-range Width')            
            self.checkBox_useTDlimits.setText('Use Limits of dmin/dmax')
            
            
        if self.radioButton_useSq.isChecked():
            print "S(Q)-language"
            self.use_xaxis_label = 'Q ($\AA^{-1}$)'
            self.use_yaxis_label = 'S(Q)'
            self.labRmin.setText('Qmin')
            self.labRmax.setText('Qmax')
            self.label_4.setText('deltaQ')
            self.label_6.setText('Qtop')
            self.label_3.setText('Q-range dependent Parent Fraction Series')
            self.label_5.setText('   Box-Car range : Starting with Qmin and Qmax on LHS,\nincrease Qmax and Qmin by deltaQ up to the Qmax = Qtop')
            self.label_8.setText(' Extending Q-range : Starting with Qmin/ max in Fit Parameters\n    increase Qmax by a value of Delta Q each time up to Qtop')
            self.checkBox_useR_for_Fsweep.setText('Use Qmin/Qmax')
            self.labRmin_2.setText('Q-local min')
            self.labRmax_2.setText('Q-local max')
            self.btnConstRrange.setText('Constant Q-range Width')
            self.btnExtendRrange.setText('Extending Q-range Width')            
            self.checkBox_useTDlimits.setText('Use Limits of Qmin/Qmax')
        
        self.possible_warning()
        
    def notice_spinbox_values(self):
        #print "I see things!"
        self.rmax = float(self.doubleSpinBox_Rmax.value())
        self.rmin = float(self.doubleSpinBox_Rmin.value())
        self.runMin = int(self.spinBox_RunMin.value())#-1
        self.runMax = int(self.spinBox_RunMax.value())#-1
    
    
    def save_alldata_one_file(self):
        #save all data read in (regardless of R-range)
        self.labAlarm.setText(" ! ")
        default_name = "all_data_single_file.dat"
        file = QtGui.QFileDialog.getSaveFileName(self, "Save File", default_name, ".dat")
        file_created_info = "# This file contains all data files, for ease of 3D plotting.\n"
        file_created_info += "# Recommend plotting in gnuplot with commands: set pm3d map, splot 'this_file.dat'\n"
        date_is = time.strftime("%m/%d/%Y")
        time_is = time.strftime("%H:%M:%S")
        file_created_info += "# This data generated with CATS on "+str(date_is)+" at "+str(time_is)+".\n"
        outfile = open(file,'w')
        outfile.write(file_created_info)

        for i in range(len(self.rdata[0])):
            for j in range(len(self.rdata)):
                if len(self.run_variables) > 0:
                    outfile.write(str(self.rdata[j][i])+" "+str(self.run_variables[j])+" "+str(self.grdata[j][i])+"\n")
                else:
                    outfile.write(str(self.rdata[j][i])+" "+str(j)+" "+str(self.grdata[j][i])+"\n")
            outfile.write("\n")
        outfile.close()
        print "saved ok"
        self.labAlarm.setText("   ")
        
    def save_alldata_inRange_one_file(self):
        self.notice_spinbox_values()
        #save datafiles between Run_Min and Run_Max, Rmax to Rmin 
        self.labAlarm.setText(" ! ")
        default_name = "select_data_single_file.dat"
        file = QtGui.QFileDialog.getSaveFileName(self, "Save File", default_name, ".dat")
        file_created_info = "# This file contains selected data files, for ease of 3D plotting.\n"
        file_created_info += "# runs between "+str(self.runMin)+" and "+str(self.runMax)+"\n"
        file_created_info += "# plotting in R between "+str(self.rmin)+" and "+str(self.rmax)+"\n"
        file_created_info += "# Recommend plotting in gnuplot with commands: set pm3d map, splot '"+str(default_name)+"'\n"
        date_is = time.strftime("%m/%d/%Y")
        time_is = time.strftime("%H:%M:%S")
        file_created_info += "# This data generated with CATS on "+str(date_is)+" at "+str(time_is)+".\n"
        file_created_info += "\n"
        outfile = open(file,'w')
        outfile.write(file_created_info)

        for i in range(len(self.rdata[0])):
            if self.rdata[0][i] >= self.rmin and self.rdata[0][i] <= self.rmax:
                for j in range(len(self.rdata)):
                    if j >= self.runMin and j <= self.runMax:
                        if len(self.run_variables) > 0:
                            outfile.write(str(self.rdata[j][i])+" "+str(self.run_variables[j])+" "+str(self.grdata[j][i])+"\n")                        
                        else:
                            outfile.write(str(self.rdata[j][i])+" "+str(j)+" "+str(self.grdata[j][i])+"\n")
                outfile.write("\n")
        outfile.close()
        print "saved ok"
        self.labAlarm.setText("   ")    

    def save_alldifference_inRange_one_file(self):
        self.notice_spinbox_values()
        #save datafiles between Run_Min and Run_Max, Rmax to Rmin 
        self.labAlarm.setText(" ! ")   

        
        noise_lvl_string, ok = QtGui.QInputDialog.getText(self, 'Noise Level?', 'Enter threshold of change to report:')
        if len(noise_lvl_string) == 0:
            print "nothing entered, default is 0"
            noise_lvl = 0.0
        else:
            try:
                noise_lvl = float(noise_lvl_string)
            except:
                print "weird "+str(noise_lvl_string)
                noise_lvl = 0.0
        
        
        default_name = "runTOrun_differences_single_file.dat"
        file = QtGui.QFileDialog.getSaveFileName(self, "Save File", default_name, ".dat")
        file_created_info = "# This file contains selected run-to-run differences, for ease of 3D plotting.\n"
        file_created_info += "# runs between "+str(self.runMin)+" and "+str(self.runMax)+"\n"
        file_created_info += "# plotting in R between "+str(self.rmin)+" and "+str(self.rmax)+"\n"
        file_created_info += "# Recommend plotting in gnuplot with commands: set pm3d map, splot '"+str(default_name)+"'\n"
        date_is = time.strftime("%m/%d/%Y")
        time_is = time.strftime("%H:%M:%S")
        file_created_info += "# This data generated with CATS on "+str(date_is)+" at "+str(time_is)+".\n"
        file_created_info += "\n"
        outfile = open(file,'w')
        outfile.write(file_created_info)

        
        for i in range(len(self.rdata[0])):
            if self.rdata[0][i] >= self.rmin and self.rdata[0][i] <= self.rmax:
                for j in range(1,len(self.rdata)):
                    if j >= self.runMin and j <= self.runMax:
                        this_diff = self.grdata[j][i] - self.grdata[j-1][i]
                        if abs(this_diff) < noise_lvl:
                            this_diff = 0.0
                        if len(self.run_variables) > 0:
                            outfile.write(str(self.rdata[j][i])+" "+str(self.run_variables[j])+" "+str(this_diff)+"\n")
                        else:
                            outfile.write(str(self.rdata[j][i])+" "+str(j)+" "+str(this_diff)+"\n")
                outfile.write("\n")
        outfile.close()
        print "saved ok"
        self.labAlarm.setText("   ")   

        
    def save_parent_diff(self,choice):
        # save the difference between parents (dad - mom)
        self.labAlarm.setText(" ! ")
        default_name = "parental_differences_"+str(choice)+"col.dat"
        file = QtGui.QFileDialog.getSaveFileName(self, "Save File", default_name, ".dat")
        
        diff = []
        tot_diff = 0.0
        for i in range(len(self.rdata[0])):
            diff.append(self.dad[0][i] - self.mom[0][i])
            tot_diff += self.dad[0][i] - self.mom[0][i]
        if choice == 2:
            file_created_info = "# This file contains the difference between the two parent files (dad - mom).\n"
        elif choice == 4:
            file_created_info = "# This file contains the two parent files (mom, dad), and difference between them (dad - mom).\n"
        file_created_info += "# Mother file : "+str(self.labMomFile.text())+"\n"
        file_created_info += "# Father file : "+str(self.labDadFile.text())+"\n"
        file_created_info += "# Sum of difference was : "+str(tot_diff)+" or "+str(tot_diff/float(len(self.rdata[0])+1))+" per r-value.\n"
        date_is = time.strftime("%m/%d/%Y")
        time_is = time.strftime("%H:%M:%S")
        file_created_info += "# This data generated with CATS on "+str(date_is)+" at "+str(time_is)+".\n"
        file_created_info += "\n"        
        outfile = open(file,'w')
        outfile.write(file_created_info)
        
        for i in range(len(self.rdata[0])):
            if choice == 2:
                outfile.write(str(self.rdata[0][i])+" "+str(diff[i])+"\n")
            elif choice == 4:
                outfile.write(str(self.rdata[0][i])+" "+str(self.mom[0][i])+" "+str(self.dad[0][i])+" "+str(diff[i])+"\n")
        print "saved ok"
        self.labAlarm.setText("   ")
    
    def save_misfit_onebigfile(self):
        self.labAlarm.setText(" ! ")
        #update rmin/rmax and runMin runMax values (not sure if a good idea)
        self.notice_spinbox_values()

        default_name = "misfits_OneFile.dat"
        file = QtGui.QFileDialog.getSaveFileName(self, "Save File", default_name, ".dat")
        file_created_info = "# This file contains all the misfit of data for ease of 3D plotting.\n"
        file_created_info += "# runs between "+str(self.runMin)+" and "+str(self.runMax)+"\n"
        file_created_info += "# plotting in R between "+str(self.rmin)+" and "+str(self.rmax)+"\n"
        file_created_info += "# Mother file : "+str(self.labMomFile.text())+"\n"
        file_created_info += "# Father file : "+str(self.labDadFile.text())+"\n"
        if self.recalc_params[4] == 0 : # calculated with different parameters
            file_created_info += "# Note: Fits were performed on a different range than shown in these plots.\n"
            if self.recalc_params[0] != self.rmin:
                file_created_info += "# Rmin used in fits was "+str(self.recalc_params[0])+"\n"
            if self.recalc_params[1] != self.rmax:
                file_created_info += "# Rmax used in fits was "+str(self.recalc_params[1])+"\n"
            if self.recalc_params[2] != self.runMin:
                file_created_info += "# Run_Min (lowest run) used in fitting was "+str(int(self.recalc_params[2]))+"\n"
            if self.recalc_params[3] != self.runMax:
                file_created_info += "# Run_Max (highest run) used in fitting was "+str(int(self.recalc_params[3]))+"\n"
        file_created_info += "# Recommend plotting in gnuplot with commands: set pm3d map, splot 'this_file.dat'\n"
        date_is = time.strftime("%m/%d/%Y")
        time_is = time.strftime("%H:%M:%S")
        file_created_info += "# This data generated with CATS on "+str(date_is)+" at "+str(time_is)+".\n"
        file_created_info += "\n"     
        outfile = open(file,'w')
        outfile.write(file_created_info)
        #using current Rmin / Rmax values (and going over runs as defined by runMin/runMax)
        for i in range(len(self.rdata[0])):
            if self.rdata[0][i] >= self.rmin and self.rdata[0][i] <= self.rmax: # if in r-range we want
                for j in range(len(self.rdata)): # go over each run
                    if j>= self.runMin and j <= self.runMax:
                        # calculate fit
                        this_fit_val = self.fit_mom_fraction[j]*self.mom[0][i]+self.fit_dad_fraction[j]*self.dad[0][i]
                        this_diff = this_fit_val - self.grdata[j][i]
                        if len(self.run_variables) > 0:
                            outfile.write(str(self.rdata[0][i])+" "+str(self.run_variables[j])+" "+str(this_diff)+"\n")                        
                        else:
                            outfile.write(str(self.rdata[0][i])+" "+str(j)+" "+str(this_diff)+"\n")
                outfile.write("\n")
        outfile.close()
        print "saved ok"
        self.labAlarm.setText("   ")    
    
    def save_each_fit(self,choice):
        # save the difference between parents (dad - mom)
        self.labAlarm.setText(" ! ")
        self.notice_spinbox_values()

        default_name = "each_fit_"+str(choice)+"col_"
        preface = QtGui.QFileDialog.getSaveFileName(self, "Pick a prefix", default_name, "")
        print "default name "+str(preface)
        for each_fit in range(len(self.rdata)):
            if each_fit <= self.runMax and each_fit >= self.runMin:
                this_file_name = preface+str(each_fit)+".dat"
                outfile = open(this_file_name,'w')
                if choice == 2:
                    file_created_info = "# This file contains the fit using these two parents.\n"
                elif choice == 3:
                    file_created_info = "# This file contains the fit using these two parents, and associated mismatch with data.\n"
                file_created_info += "# Mother file : "+str(self.labMomFile.text())+"\n"
                file_created_info += "# Father file : "+str(self.labDadFile.text())+"\n"
                file_created_info += "# Fit was performed on dataset: "+str(self.fileList[each_fit])+"\n"
                if self.recalc_params[4] == 0 : # calculated with different parameters
                    file_created_info += "# Note: Fits were performed on a different range than shown in these plots.\n"
                    if self.recalc_params[0] != self.rmin:
                        file_created_info += "# Rmin used in fits was "+str(self.recalc_params[0])+"\n"
                    if self.recalc_params[1] != self.rmax:
                        file_created_info += "# Rmax used in fits was "+str(self.recalc_params[1])+"\n"
                    if self.recalc_params[2] != self.runMin:
                        file_created_info += "# Run_Min (lowest run) used in fitting was "+str(int(self.recalc_params[2]))+"\n"
                    if self.recalc_params[3] != self.runMax:
                        file_created_info += "# Run_Max (highest run) used in fitting was "+str(int(self.recalc_params[3]))+"\n"
                date_is = time.strftime("%m/%d/%Y")
                time_is = time.strftime("%H:%M:%S")
                file_created_info += "# This data generated with CATS on "+str(date_is)+" at "+str(time_is)+".\n"

                outfile.write(file_created_info)
                # metadata written into file, now write appropriate fits
                if choice == 2:
                    for r in range(len(self.rdata[0])):
                        if self.rdata[0][r] <= self.rmax and self.rdata[0][r] >= self.rmin:
                            this_fit_val = self.fit_mom_fraction[each_fit]*self.mom[0][r]+self.fit_dad_fraction[each_fit]*self.dad[0][r]
                            outfile.write(str(self.rdata[0][r])+" "+str(this_fit_val)+"\n")
                elif choice == 3:
                    for r in range(len(self.rdata[0])):
                        if self.rdata[0][r] <= self.rmax and self.rdata[0][r] >= self.rmin:
                            this_fit_val = self.fit_mom_fraction[each_fit]*self.mom[0][r]+self.fit_dad_fraction[each_fit]*self.dad[0][r]
                            this_diff = this_fit_val - self.grdata[each_fit][r]
                            outfile.write(str(self.rdata[0][r])+" "+str(this_fit_val)+" "+str(this_diff)+"\n")                            

                outfile.close()
        print "saved ok"
        self.labAlarm.setText("   ")
    
    def read_in_run_variable(self):
        
    	#now count how many files exist with that preface
    	myPath = os.getcwd()
    	#self.fileList = glob.glob1(myPath,self.file_preface+"*")
        
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file')
        with open(fname, "r") as data:
            read_file = data.readlines()
        
        print "file length was: "+str(len(read_file))
        file_names = []
        variables = []
        for i in range(len(read_file)):
            file_names.append(read_file[i].split()[0])
            self.run_variables.append(read_file[i].split()[1])
        #this_list = glob.glob1(myPath,fname)
        
        print "first file name is : "+str(file_names[0])
        print "first value is : "+str(self.run_variables[0])
        print "last value is: "+str(file_names[len(read_file)-1])+" "+str(self.run_variables[len(read_file)-1])
        #self.run_variables = variables        

    def similarity_matrix(self):
        #setup all data into dataframe
        df = pd.DataFrame(index=self.rdata[0])

        for this_run in range(len(self.rdata)):
            df[this_run] = np.array(self.grdata[this_run])
            
        score_matrix = np.zeros([len(df.columns),len(df.columns)])
        
        for y1i in range(len(df.columns)):
            col1 = df.columns[y1i]
            for y2i in range(len(df.columns)):
                col2 = df.columns[y2i]
                score_matrix[y1i,y2i] = dan.calc_rw(df.loc[self.rmin:self.rmax,col1],df.loc[self.rmin:self.rmax,col2])

        self.df_score = pd.DataFrame(data = score_matrix, index=df.columns, columns= df.columns)
        self.btnPlot_DCM.setEnabled(True)
        self.checkBox_overlay_else.setEnabled(True)
        self.spinBox_DCM_else_num.setEnabled(True)
        self.labDataset.setEnabled(True)
        self.spinBox_DCM_else_num.setMaximum(len(df.columns)-1)

    def calc_confidence_intervals(self):
        #self.df_all_data = pd.DataFrame(index=self.rdata[0])
        df_parents = pd.DataFrame(index=self.rdata[0])
        df_parents['mom'] = np.array(self.mom[0])
        df_parents['dad'] = np.array(self.dad[0])
        
        #for this_run in range(len(self.rdata)):
        #    self.df_all_data[this_run] = np.array(self.grdata[this_run])        
        
        use_conf = 1.0 - float(self.doubleSpinBox_confidence.value())
        use_num_combos = int(self.spinBox_NumCombosTry.value())
        
        if self.checkBox_useR_for_Fsweep.isChecked():
            srmin = float(self.doubleSpinBox_Rmin.value())
            srmax = float(self.doubleSpinBox_Rmax.value())
            best_phis, best_ress, wbp_highs, wbp_lows = dan.try_all_the_things(df_parents.loc[srmin:srmax,'mom'],df_parents.loc[srmin:srmax,'dad'],self.df_all_data.loc[srmin:srmax,:],confidence=use_conf,num_pts=use_num_combos)
        else:
            best_phis, best_ress, wbp_highs, wbp_lows = dan.try_all_the_things(df_parents.loc[:,'mom'],df_parents.loc[:,'dad'],self.df_all_data,confidence=use_conf,num_pts=use_num_combos)
        
        self.conf_int_stuff = [best_phis, best_ress, wbp_highs, wbp_lows]
        self.btnPlotFracSweep.setEnabled(True)
        
    def plot_errorbar_things(self):
    
        best_phis = self.conf_int_stuff[0]
        best_ress = self.conf_int_stuff[1]
        wbp_highs = self.conf_int_stuff[2]
        wbp_lows = self.conf_int_stuff[3]
        
        plt.figure()
        phi_error_high = (wbp_highs-best_phis)
        phi_error_low = (best_phis-wbp_lows)

        #phi_error_high = (wbp_highs-best_phis)
        #phi_error_low = (best_phis-wbp_lows)
        delta_conf = 100./float(self.spinBox_NumCombosTry.value())
        
        conf_title = str(100.*float(self.doubleSpinBox_confidence.value()))+' $\pm$ '+str(delta_conf)[:3]+'% Confidence Interval'
        plt.title(conf_title)
        
        plt.errorbar(self.df_all_data.columns,best_phis, yerr = [phi_error_low,phi_error_high],label='G(r)',linewidth=1.5,capsize=5,elinewidth=1.5)
        #plt.errorbar(df_all_SQ.columns,best_phis, yerr = [phi_error_low,phi_error_high],label='S(Q)',linewidth=1.5,capsize=5,elinewidth=1.5)
        plt.xlabel('runs')
        plt.ylabel('$\phi_{mom}$')
        plt.show()
    
    def top_down_show_data(self): 
        plt.figure()
        
        if self.checkBox_useTDlimits.isChecked():
            srmin = float(self.doubleSpinBox_Rmin.value())
            srmax = float(self.doubleSpinBox_Rmax.value())        
            runmin = int(self.spinBox_RunMin.value())
            runmax = int(self.spinBox_RunMax.value())
        else:
            srmin = self.df_all_data.index[0]
            srmax = self.df_all_data.index[-1]
            runmin = self.df_all_data.columns[0]
            runmax = self.df_all_data.columns[-1]        
        
        if self.checkBox_defCMap.isChecked():
            dan.top_down_plot(self.df_all_data.loc[srmin:srmax,runmin:runmax],cmap='viridis')
            plt.colorbar(label=self.use_yaxis_label)

        else:
            use_cmap = str(self.PTE_cmap_choice.text())
            use_cmin = float(self.doubleSpinBox_cmin.value())
            use_cmax = float(self.doubleSpinBox_cmax.value())
            dan.top_down_plot(self.df_all_data.loc[srmin:srmax,runmin:runmax],cmap=use_cmap,cmin=use_cmin,cmax=use_cmax)        
            if self.checkBox_showCbar.isChecked():
                plt.colorbar(label=self.use_yaxis_label)
        
        plt.xlabel(self.use_xaxis_label)
        plt.ylabel('runs')
        
        if self.checkBox_overlayAvg.isChecked():
            x = np.array(self.df_all_data.index)
            y = np.array(self.df_all_data.loc[:,:].mean(axis=1))
            usecolor = str(self.PTE_overlay_color.text())
            yoffset = float(self.doubleSpinBox_OLayOffset.value())
            yscale = float(self.doubleSpinBox_OLayScale.value())
            plt.plot(x,y*yscale + yoffset,color=usecolor,linewidth=1)
            
        
        plt.tight_layout()
        plt.show()
        
        
    def show_similarity_matrix(self):
        plt.figure()
        
        #main plot, use defaults or user provided
        if self.checkBox_defCMap.isChecked():
            dan.top_down_plot(self.df_score,cmap='viridis')
            plt.colorbar(label='Rw')
        else: #read from user provided
            use_cmap = str(self.PTE_cmap_choice.text())
            use_cmin = float(self.doubleSpinBox_cmin.value())
            use_cmax = float(self.doubleSpinBox_cmax.value())
            dan.top_down_plot(self.df_score,cmap=use_cmap,cmin=use_cmin,cmax=use_cmax)
       
            if self.checkBox_showCbar.isChecked():
                plt.colorbar(label='Rw')
                plt.show()
        
        plt.xlabel('Run #')
        plt.ylabel('Run #')
        
        if self.checkBox_overlay_else.isChecked():
            print ('putting something up ')
            overlay_num = int(self.spinBox_DCM_else_num.value())
            overx = np.array(self.df_score.index)
            overy = np.array(self.df_score[overlay_num])
            
            plt.plot(overx, overy*(len(self.df_score.columns)-1)/max(overy),'w')
            
            plt.show()
        
        if self.checkBox_overlay_mom.isChecked():
            print ('displaying mom-path')
            plt.plot(np.array(self.df_score.index),(1.-self.fit_mom_fraction)*(len(self.df_score.columns)-1),'r')
        if self.checkBox_overlay_dad.isChecked():
            print ('displaying dad-path')
            plt.plot(np.array(self.df_score.index),(1.-self.fit_dad_fraction)*(len(self.df_score.columns)-1),'b')
            
    
    def show_feature_track_map(self):
        print "well this is a fine mess"
        srmin = float(self.doubleSpinBox_FTRmin.value())
        srmax = float(self.doubleSpinBox_FTRmax.value())
        print "Display from rmin/rmax "+str(srmin)+ "to "+str(srmax)
        
        #setup things for first plot
        df_local = self.df_all_data.loc[srmin:srmax,:]
        this_yaxis= ''
        
        
        
        double_feature = False
        if self.checkBox_ft_Int.isChecked():
            this_yaxis = str('$\sum_{a}^{b}$')+str(self.use_yaxis_label)
            feature_y = np.array(df_local.sum(axis=0))
            feature_x = np.array(df_local.columns)
        if self.checkBox_ft_upfMom.isChecked():
            feature_x = np.array(df_local.columns)
            feature_y = self.fit_mom_fraction
            this_yaxis = str('$\phi_{mom}$')
        if self.checkBox_ft_upfDad.isChecked():
            feature_x = np.array(df_local.columns)
            if this_yaxis == '':
                this_yaxis = str('$\phi_{dad}$')
                feature_y = self.fit_dad_fraction
            else:
                this_yaxis += str('$ + \phi_{dad}$')
                other_feature_y = self.fit_dad_fraction
                double_feature = True
        
        try:
            feature_x
        except:
            print ("defaulting to intensity sum (you didn't select anything)")
            this_yaxis = str('$\sum_{a}^{b}$')+str(self.use_yaxis_label)
            feature_y = np.array(df_local.sum(axis=0))
            feature_x = np.array(df_local.columns) 
            self.checkBox_ft_Int.setChecked(True)
        
        plt.figure()
        plt.subplot(211)
        plt.plot(feature_x,feature_y,'k')
        if double_feature:
            plt.plot(feature_x, other_feature_y, 'r')
        plt.ylabel(this_yaxis)
        #plt.xlabel(self.use_yaxis_label)
        #plt.ylabel(self.use_xaxis_label)
        #plt.xlabel('run #')
        
        plt.subplot(212)
        dan.top_down_plot(self.df_all_data.loc[srmin:srmax,:].T)
        #plt.xlabel(self.use_yaxis_label)
        plt.ylabel(self.use_xaxis_label)
        plt.xlabel('run #')    
        plt.tight_layout()
        
    def save_all_fits(self):
        self.labAlarm.setText(" ! ")
        
        #update rmin/rmax and runMin runMax values (not sure if a good idea)
        self.notice_spinbox_values()

        default_name = "all_fits_OneFile.dat"
        file = QtGui.QFileDialog.getSaveFileName(self, "Save Menu", default_name, ".dat")
        file_created_info = "# This file contains all the fits for ease of 3D plotting.\n"
        file_created_info += "# runs between "+str(self.runMin)+" and "+str(self.runMax)+"\n"
        file_created_info += "# plotting in R between "+str(self.rmin)+" and "+str(self.rmax)+"\n"
        file_created_info += "# Mother file : "+str(self.labMomFile.text())+"\n"
        file_created_info += "# Father file : "+str(self.labDadFile.text())+"\n"        
        if self.recalc_params[4] == 0 : # calculated with different parameters
            file_created_info += "# Note: Fits were performed on a different range than shown in these plots.\n"
            if self.recalc_params[0] != self.rmin:
                file_created_info += "# Rmin used in fits was "+str(self.recalc_params[0])+"\n"
            if self.recalc_params[1] != self.rmax:
                file_created_info += "# Rmax used in fits was "+str(self.recalc_params[1])+"\n"
            if self.recalc_params[2] != self.runMin:
                file_created_info += "# Run_Min (lowest run) used in fitting was "+str(int(self.recalc_params[2]))+"\n"
            if self.recalc_params[3] != self.runMax:
                file_created_info += "# Run_Max (highest run) used in fitting was "+str(int(self.recalc_params[3]))+"\n"
        file_created_info += "# Recommend plotting in gnuplot with commands: set pm3d map, splot 'this_file.dat'\n"
        date_is = time.strftime("%m/%d/%Y")
        time_is = time.strftime("%H:%M:%S")
        file_created_info += "# This data generated with CATS on "+str(date_is)+" at "+str(time_is)+".\n"
        file_created_info += "\n"
        outfile = open(file,'w')
        outfile.write(file_created_info)
        #using current Rmin / Rmax values (and going over runs as defined by runMin/runMax)
        for i in range(len(self.rdata[0])):
            if self.rdata[0][i] >= self.rmin and self.rdata[0][i] <= self.rmax: # if in r-range we want
                for j in range(len(self.rdata)): # go over each run
                    if j>= self.runMin and j <= self.runMax:
                        # calculate fit
                        this_fit_val = self.fit_mom_fraction[j]*self.mom[0][i]+self.fit_dad_fraction[j]*self.dad[0][i]
                        this_diff = this_fit_val - self.grdata[j][i]
                        if len(self.run_variables) > 0:
                            outfile.write(str(self.rdata[0][i])+" "+str(self.run_variables[j])+" "+str(this_fit_val)+"\n")                        
                        else:
                            outfile.write(str(self.rdata[0][i])+" "+str(j)+" "+str(this_fit_val)+"\n")
                outfile.write("\n")
        outfile.close()
        print "saved ok"
        self.labAlarm.setText("   ")    
        
        
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MyApp()
    window.show()
    sys.exit(app.exec_())
