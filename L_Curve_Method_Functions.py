# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:50:26 2023

#Author:       Tim Tichter
#Years:        2022
#Function:     POLAROGRAPHICAS L-curve-Method for optimizing the regularitation parameter
#              of the DRT transformation. Likewise, all functions are based on the DRT script.
#              "DRT_Native_Gauss_ColeCole.py". No processing is done with the DRT-tools translation.
"""

##############################################################################################
#  Load all required modules for computing DRTs first
##############################################################################################
import numpy                           as np
from scipy.optimize                    import nnls
from math                              import ceil,floor

##############################################################################################
#  Re-import the DRT functions. This will also load the KK, Anti-KK and combi-KK
#  back-multiplication functions!
##############################################################################################
from DRT_Native_Gauss_ColeCole         import *


##############################################################################################
##############################################################################################
# Define the function for performing the L-curve Check with the native DRT
##############################################################################################
##############################################################################################

def L_Curve_Native(freq_raw, Z_raw, dampers_to_try, low_ext, high_ext, resol, DRT_Type, Reg_Type, BackCalc_Type):
    #------------------------------------------------------------
    # compute one DRT to find the length of the output-array
    #------------------------------------------------------------
    Initial_DRT        = DRT_Native(freq_raw, Z_raw, damp = dampers_to_try[0], low_ext, high_ext, resol, DRT_Type, Reg_Type)
    #------------------------------------------------------------
    # initialize empty arrays for storing a) the individual DRTs
    # b) the back-calculated EIS-data and c) the deviation of #
    # measured EIS-data and the fitting result.
    #------------------------------------------------------------
    Output_Taus_DRTs   = np.zeros((len(Initial_DRT[0]), 2*len(dampers_to_try)))
    Back_Calc_EIS      = np.zeros((len(Initial_DRT[0]), 3*len(dampers_to_try)))
    Deviations         = np.zeros(len(dampers_to_try))
    #------------------------------------------------------------    
    # iterate through all damper-values provided as input
    #------------------------------------------------------------
    for damper in range(len(dampers_to_try)):
        #--------------------------------------------------------
        # decide which Type of Back-calculation should be done
        # this will define which kind of DRT will be computed
        #--------------------------------------------------------
        if BackCalc_Type == "KK":
             KK_enable     = 1
             AntiKK_enable = 0
             Combi_enable  = 0 
        elif BackCalc_Type == "AntiKK":
            KK_enable     = 0
            AntiKK_enable = 1
            Combi_enable  = 0              
        elif BackCalc_Type == "CombiKK":
            KK_enable     = 0
            AntiKK_enable = 0
            Combi_enable  = 1 
            
            
        TempResult                     = DRT_Native(freq_raw, Z_raw, damp = dampers_to_try[0], low_ext, high_ext, resol, DRT_Type, Reg_Type)
        
        Output_Taus_DRTs[::,2*i+0]     = TempResult[0]
        Output_Taus_DRTs[::,2*i+1]     = TempResult[1]
        #-----------------------------------------------------------------------------------------
        # Back-calculation from DRTs to impedances in any permutation, KK, anti KK, combi
        #-----------------------------------------------------------------------------------------
        if KK_enable == 1:
            Z_Back           = KramersKronig_Multiplication(MATRIX = IM_Output[4], Real_DRT = RE_Output[1], Imag_DRT = IM_Output[1])
            Z_Back = Z_Back + Real_Offset     
        if AntiKK_enable == 1:
            Z_Back           = AntiKramersKronig_Multiplication(MATRIX = IM_Output[4], Real_DRT = RE_Output[1], Imag_DRT = IM_Output[1])
            if RemoveRealOffset == 1:
                Z_Back = Z_Back + Real_Offset
        if CombiKK_enable == 1:
            Z_Back           = CombiKramersKronig_Multiplication(MATRIX = CO_Output[4], Comb_DRT = CO_Output[1])
            if RemoveRealOffset == 1:
                Z_Back = Z_Back + Real_Offset














































