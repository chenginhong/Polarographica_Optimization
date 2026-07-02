# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:30:20 2023

@author: timtic
"""

from    tkinter import *
import  numpy as np                         


class PolArStat_CA_Write_Read():
    Zero_IDX_E  = 1
    Zero_IDX_I  = 1
    ReadResist  = 1
    Exper_Notes = "No notes found"
    Output_File = None
    ###############################################################################
    # Define a function for writing the final CA-outputs into a .txt file, which 
    # has been specified before the measurement was initialized.
    ############################################################################### 
    def SAVE_DATA_AFTER_STOP(STOP_STATE, EXP_PARAMS, DATA):
        EXP_PAR_LABEL = np.array(["E1_vs_RE_in_V", "E2_vs_RE_in_V", "E3_vs_RE_in_V", "E4_vs_RE_in_V", "E5_vs_RE_in_V", "t1_in_s", "t2_in_s", "t3_in_s", "t4_in_s", "t5_in_s", "n_reps", "Rread_in_Ohm"])
        #==========================================================================
        #   Write output in the output-txt file
        #==========================================================================
        PolArStat_CA_Write_Read.Output_File.write("Stop_State\t")
        PolArStat_CA_Write_Read.Output_File.write(STOP_STATE)
        PolArStat_CA_Write_Read.Output_File.write("\n")
        for i in range(len(EXP_PAR_LABEL)):
            PolArStat_CA_Write_Read.Output_File.write(EXP_PAR_LABEL[i])
            PolArStat_CA_Write_Read.Output_File.write("\t") 
            PolArStat_CA_Write_Read.Output_File.write(str(EXP_PARAMS[i]))
            PolArStat_CA_Write_Read.Output_File.write("\n")
        PolArStat_CA_Write_Read.Output_File.write("Notes:\t")
        PolArStat_CA_Write_Read.Output_File.write(PolArStat_CA_Write_Read.Exper_Notes)
        PolArStat_CA_Write_Read.Output_File.write("\n")
        PolArStat_CA_Write_Read.Output_File.write("======================================================\n")
        PolArStat_CA_Write_Read.Output_File.write("Loop.No.\t Step.No.\t t/s\t E_vs_RE/V\t I/mA\n\n")
        for i in range(len(DATA[::,0])-1):
            PolArStat_CA_Write_Read.Output_File.write(str(DATA[i+1,0]))
            PolArStat_CA_Write_Read.Output_File.write("\t")
            PolArStat_CA_Write_Read.Output_File.write(str(DATA[i+1,1]))
            PolArStat_CA_Write_Read.Output_File.write("\t")
            PolArStat_CA_Write_Read.Output_File.write(str(0.001*DATA[i+1,2]))
            PolArStat_CA_Write_Read.Output_File.write("\t")
            PolArStat_CA_Write_Read.Output_File.write(str(-0.000249*(DATA[i+1,3]-PolArStat_CA_Write_Read.Zero_IDX_E)))
            PolArStat_CA_Write_Read.Output_File.write("\t")
            PolArStat_CA_Write_Read.Output_File.write(str(-0.12452*(DATA[i+1,4]-PolArStat_CA_Write_Read.Zero_IDX_I)/PolArStat_CA_Write_Read.ReadResist))
            PolArStat_CA_Write_Read.Output_File.write("\n")
        PolArStat_CA_Write_Read.Output_File.close()

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
            WHERE_TO_WRITE.insert(INSERT, "%.1f \t"  %(DATA_TO_WRITE[i,0])   )
            WHERE_TO_WRITE.insert(INSERT, "%.1f \t"  %(DATA_TO_WRITE[i,1])   )
            WHERE_TO_WRITE.insert(INSERT, "%.3f \t"  %(0.001*DATA_TO_WRITE[i,2])   )
            WHERE_TO_WRITE.insert(INSERT, "%.4f \t"  %(-0.000249*(DATA_TO_WRITE[i,3]-PolArStat_CA_Write_Read.Zero_IDX_E))   )
            WHERE_TO_WRITE.insert(INSERT, "%.4f \n"  %(-0.12452*(DATA_TO_WRITE[i,4]-PolArStat_CA_Write_Read.Zero_IDX_I)/PolArStat_CA_Write_Read.ReadResist)   )
        
    
