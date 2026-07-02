# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:30:20 2023

@author: timtic
"""

from    tkinter import *
import  numpy as np                         


class PolArStat_CV_Write_Read():
    Zero_IDX_E  = 1
    Zero_IDX_I  = 1
    ReadResist  = 1
    Exper_Notes = "No notes found"
    Output_File = None
    ###############################################################################
    # Define a function for writing the final CV-outputs into a .txt file, which 
    # has been specified before the measurement was initialized.
    ############################################################################### 
    def SAVE_DATA_AFTER_STOP(STOP_STATE, EXP_PARAMS, DATA):
        EXP_PAR_LABEL = np.array(["E_in_vs_RE_in_V","E_v1_vs_RE_in_V","E_v2_vs_RE_in_V","E_fi_vs_RE_in_V","n_cycles","Scanr_in_mV/s","Cond_t_in_s","Rread_in_Ohm"])
        #==================================================================
        #   Write output in the output-txt file
        #==================================================================
        PolArStat_CV_Write_Read.Output_File.write("Stop_State\t")
        PolArStat_CV_Write_Read.Output_File.write(STOP_STATE)
        PolArStat_CV_Write_Read.Output_File.write("\n")
        for i in range(len(EXP_PAR_LABEL)):
            PolArStat_CV_Write_Read.Output_File.write(EXP_PAR_LABEL[i])
            PolArStat_CV_Write_Read.Output_File.write("\t") 
            PolArStat_CV_Write_Read.Output_File.write(str(EXP_PARAMS[i]))
            PolArStat_CV_Write_Read.Output_File.write("\n")
        PolArStat_CV_Write_Read.Output_File.write("Notes:\t")
        PolArStat_CV_Write_Read.Output_File.write(PolArStat_CV_Write_Read.Exper_Notes)
        PolArStat_CV_Write_Read.Output_File.write("\n")
        PolArStat_CV_Write_Read.Output_File.write("======================================================\n\n")
        PolArStat_CV_Write_Read.Output_File.write("Ramp-Index\ttime in ms\tE_WE_vs_RE in V\tI in mA\tCyc.No.\n\n")
        for i in range(len(DATA[::,0])-1):
            PolArStat_CV_Write_Read.Output_File.write(str(DATA[i+1,0]))
            PolArStat_CV_Write_Read.Output_File.write("\t")
            PolArStat_CV_Write_Read.Output_File.write(str(0.001*DATA[i+1,1]))
            PolArStat_CV_Write_Read.Output_File.write("\t")
            PolArStat_CV_Write_Read.Output_File.write(str(-0.000249*(DATA[i+1,2]-PolArStat_CV_Write_Read.Zero_IDX_E)))
            PolArStat_CV_Write_Read.Output_File.write("\t")
            PolArStat_CV_Write_Read.Output_File.write(str(-0.12452*(DATA[i+1,3]-PolArStat_CV_Write_Read.Zero_IDX_I)/PolArStat_CV_Write_Read.ReadResist))
            PolArStat_CV_Write_Read.Output_File.write("\t")
            PolArStat_CV_Write_Read.Output_File.write(str(DATA[i+1,4]))
            PolArStat_CV_Write_Read.Output_File.write("\n")
        PolArStat_CV_Write_Read.Output_File.close()
        
    ###############################################################################
    # Define a function for writing a text output to a certain text-field
    # the INSERT is specified by the "from   tkinter import *" (see above)
    ############################################################################### 
    def WRITE_TEXT_OUTPUT(WHERE_TO_WRITE, TEXT_TO_WRITE): 
        WHERE_TO_WRITE.insert(INSERT, "%s \n" %TEXT_TO_WRITE )
        
    ###############################################################################
    # Define a function for writing a data output to a certain text-field
    # the INSERT is specified by the "from   tkinter import *" (see above)
    ###############################################################################   
    def WRITE_DATA_OUTPUT(WHERE_TO_WRITE, DATA_TO_WRITE):
        for i in range(len(DATA_TO_WRITE[::,0])):
            WHERE_TO_WRITE.insert(INSERT, "%.3f \t"  %((1e-6)*DATA_TO_WRITE[i,0])   )
            WHERE_TO_WRITE.insert(INSERT, "%.4f \t" %(-0.000249*(DATA_TO_WRITE[i,1]-PolArStat_CV_Write_Read.Zero_IDX_E))   )
            WHERE_TO_WRITE.insert(INSERT, "%.4f \t" %(-0.12452*(DATA_TO_WRITE[i,2]-PolArStat_CV_Write_Read.Zero_IDX_I)/PolArStat_CV_Write_Read.ReadResist)   )
            WHERE_TO_WRITE.insert(INSERT, "%.f \n"  %DATA_TO_WRITE[i,3]   )
        
        
        
        
        