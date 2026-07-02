# -*- coding: utf-8 -*-
"""
Created on Tue May 11 11:55:09 2021

@author: gisel
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate                     import InterpolatedUnivariateSpline
from cmath                                 import *
from scipy.optimize                        import fsolve 
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
from tkinter                               import *
from tkinter.filedialog                    import askopenfilename
from tkinter.filedialog                    import asksaveasfilename



def Load_Computed_CA():
    Path_CA       = askopenfilename()
    data3D        = np.loadtxt(Path_CA, delimiter = '\t')
    Chrono_t3D    = data3D[::,0]
    Chrono_I3D    = data3D[::,1]
    x_Compare     = np.log10(np.logspace(-4, np.log10(Chrono_t3D[-1]), len(Chrono_t3D)))
    Cutoff_Definer(x_loaded = Chrono_t3D , y_loaded = Chrono_I3D)
    plt.figure(figsize = (5,5), dpi = 100)
    plt.plot(np.log10(Chrono_t3D), np.log10(Chrono_I3D), color = 'black', label = 'Loaded Data', linestyle = ':')
    plt.plot(x_Compare, np.log10(1/((10**x_Compare)*np.pi)**0.5), color = 'black', label = 'Analytical', linestyle = '-', linewidth = 0.5)
    plt.xlabel(r'log$_{10}$($t$ / s)', fontsize = 15)
    plt.ylabel(r'log$_{10}$(norm. flux)', fontsize = 15)
    plt.tight_layout()
    plt.show()
    

def Cutoff_Definer(x_loaded, y_loaded):
    Fenster = Toplevel()                                                         
    Fenster.title("Cutoff Definer")                         
    Fenster.geometry("300x200")                                            
    LowCutoff_label      = Label(Fenster,text= u'Low-Cutoff')
    LowCutoff_label.place(x = 20, y = 20)
    LowCutoff_Eingabe    = Entry(Fenster)  
    LowCutoff_Eingabe.insert(END, -0.5)                                             
    LowCutoff_Eingabe.place(x = 150, y = 20)
    HighCutoff_label       = Label(Fenster,text= r'High-Cutoff')
    HighCutoff_label.place(x =20, y = 45)
    HighCutoff_Eingabe     = Entry(Fenster)
    HighCutoff_Eingabe.insert(END, 0.5) 
    HighCutoff_Eingabe.place(x = 150, y = 45)
    Subdivisions_label    = Label(Fenster,text= r'Num. of Subdivisions')
    Subdivisions_label.place(x =20, y = 70)
    Subdivisions_Eingabe     = Entry(Fenster)
    Subdivisions_Eingabe.insert(END, 5) 
    Subdivisions_Eingabe.place(x = 150, y = 70)   
    NumOfTimeIncr_label    = Label(Fenster,text= r'Num. of dt in conv')
    NumOfTimeIncr_label.place(x =20, y = 95)
    NumOfTimeIncr_Eingabe     = Entry(Fenster)
    NumOfTimeIncr_Eingabe.insert(END, 100000) 
    NumOfTimeIncr_Eingabe.place(x = 150, y = 95)
    def Next():
        LowCutoff     = (float(LowCutoff_Eingabe.get()))
        HighCutoff    = (float(HighCutoff_Eingabe.get()))
        NumOfTimeIncr = (int(NumOfTimeIncr_Eingabe.get()))
        Subdivisons   = (int(Subdivisions_Eingabe.get()))
        Fenster.destroy()
        Conv_Func_Calcer(x_loaded, y_loaded, LowCutoff, HighCutoff, NumOfTimeIncr, Subdivisons)
        
    Next = Button(Fenster, text="Next",command=Next)
    Next.place(x = 150, y = 120, width = 60, height = 22)
    



def Conv_Func_Calcer(x_loaded, y_loaded, LowCutoff, HighCutoff, NumOfTimeIncr, Subdivisons):   
    #===================================================
    # Interpolating the current function
    #===================================================
    Chrono_t3D    = x_loaded
    Chrono_I3D    = y_loaded
    Inter_x       = np.log10(np.logspace(-4, np.log10(Chrono_t3D[-1]), len(Chrono_t3D)))
    slice_Point_A = ((np.log10(x_loaded)-LowCutoff)**2).argmin()
    slice_Point_B = ((np.log10(x_loaded)-HighCutoff)**2).argmin()
    Start_x       = Inter_x[:slice_Point_A]
    Start_1D      = np.log10(1/((10**Inter_x[:slice_Point_A])*np.pi)**0.5)
    End_x         = np.log10(Chrono_t3D[slice_Point_B:])
    End_3D        = np.log10(Chrono_I3D[slice_Point_B:])
    conc_x        = np.concatenate([Start_x , End_x ])
    conc_y        = np.concatenate([Start_1D, End_3D])
    Interpolation = InterpolatedUnivariateSpline(conc_x, conc_y, k=3)
    Inter_y       = Interpolation(Inter_x)
    Chrono_t      = 10**Inter_x
    Chrono_I      = 10**Inter_y
    #===================================================
    # Interpolating the current function -- completed
    #      Now, initialize M(t)-computation!
    #===================================================  
    Interp_log         = InterpolatedUnivariateSpline(np.log10(Chrono_t), Chrono_I, k=3)
    #---------------------------------------------------
    # Create the first set of dt to t = 1 s
    #---------------------------------------------------
    dt1                = 1.0/float(NumOfTimeIncr)
    t_Resol1           = np.arange(dt1, 1.0, dt1)
    #---------------------------------------------------
    # Interpolate on a logarithmic axis and create 
    # delta_I list.
    #---------------------------------------------------
    I_Resol1           = Interp_log(np.log10(t_Resol1))
    delta_I1           = I_Resol1[1::]-I_Resol1[0:-1:] 
    H1                 = np.zeros(len(I_Resol1))
    for i in range(len(I_Resol1)-1):  
        if i == 0:
            H1[i+1] = 2*(dt1/np.pi)**0.5
        if i != 0:
            H1[i+1] = ((1 - np.sum(H1[i::-1]*delta_I1[:i+1:]))/I_Resol1[0]) 
    H_old_Interpol   = InterpolatedUnivariateSpline(t_Resol1, H1, k = 3)
    #---------------------------------------------------
    # Create the nth-set function
    #---------------------------------------------------
    def Iterator(RecentUpperBound, PreviousUpperBound):
        dtnew              = RecentUpperBound/float(NumOfTimeIncr)
        t_Resolnew         = np.arange(dtnew, RecentUpperBound, dtnew)
        #---------------------------------------------------
        # Interpolate on a logarithmic axis and create 
        # delta_I list.
        #---------------------------------------------------
        I_Resolnew           = Interp_log(np.log10(t_Resolnew))
        delta_Inew           = I_Resolnew[1::]-I_Resolnew[0:-1:] 
        Hnew                 = np.zeros(len(I_Resolnew))
        for i in range(len(I_Resolnew)-1):  
            if i*dtnew < PreviousUpperBound-2*dtnew:
                Hnew[i] = H_old_Interpol(t_Resolnew[i])
            if i*dtnew >= PreviousUpperBound-2*dtnew:
                Hnew[i] = ((1 - np.sum(Hnew[i::-1]*delta_Inew[:i+1:]))/I_Resolnew[1])   
        Hnew_Interpol = InterpolatedUnivariateSpline(t_Resolnew, Hnew, k = 3)
        return Hnew_Interpol, t_Resolnew, Hnew
    #---------------------------------------------------
    # Create the nth-set 
    #---------------------------------------------------
    UpperBounds = np.logspace(1, np.log10(Chrono_t[-1]), Subdivisons)
    for bb in range(len(UpperBounds)):
        if bb == 0:
            Result = Iterator(UpperBounds[bb], 1)
        if bb > 0:
            Result = Iterator(UpperBounds[bb], UpperBounds[bb-1])
        H_old_Interpol = Result[0]
        t_Resol_temp   = Result[1]
        H_temp         = Result[2]
    
    #---------------------------------------------------
    # Do the final interpolation 
    #---------------------------------------------------
    Hfinal_Interpol   = InterpolatedUnivariateSpline(t_Resol_temp, H_temp, k = 3)
    t_Resol_final     = np.logspace(-6,np.log10(Chrono_t[-1]-5), 500)
    H_final           = Hfinal_Interpol(t_Resol_final)
    t_Resol_final[0]  = 0
    H_final[0]        = 0
    print ("Convolution Function got calculated")
    
    #================================================
    #                  Saving
    #================================================
    OUTPATH_MT  = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))
    file = open(OUTPATH_MT, "w")
    for jj in range(len(t_Resol_final[0:-1:])):
        file.write(str(t_Resol_final[jj]))
        file.write("\t")
        file.write(str(H_final[jj]))
        file.write("\n")
    file.close()
    
    
    #================================================
    #                  PLOTTING
    #================================================
    plt.figure(figsize = (10,5), dpi = 100)
    plt.subplot(121)
    plt.plot(np.log10(Chrono_t), np.log10(Chrono_I), color = 'black', label = 'Interpol', linestyle = '--')
    plt.plot(np.log10(Chrono_t3D), np.log10(Chrono_I3D), color = 'black', label = 'Numerical', linestyle = ':')
    plt.plot(np.log10(Chrono_t), np.log10(1/(Chrono_t*np.pi)**0.5), color = 'black', label = 'Analytical', linestyle = '-', linewidth = 0.5)
    plt.axvline(np.log10(Chrono_t3D[slice_Point_A]), color = 'black', linewidth = 0.5, linestyle = ':')
    plt.axvline(np.log10(Chrono_t3D[slice_Point_B]), color = 'black', linewidth = 0.5, linestyle = ':')
    plt.legend(frameon = False, fontsize = 15)
         
    plt.subplot(122)
    plt.plot((t_Resol_final), 2*(t_Resol_final/np.pi)**0.5, color = 'k', linewidth = 0.5, label = "1D semi-inf. ")
    plt.plot((t_Resol_final),H_final , linestyle=':', linewidth = 2, color = 'k', label = "3D network")
    plt.legend(frameon = False, fontsize = 15, loc='lower right', bbox_to_anchor=(1, 0.12))
    plt.xlabel('$t$' " / s", fontsize=15)
    plt.ylabel('$M(t)$'  "/ s" '$^{0.5}$', fontsize=15)
    plt.tick_params(direction = 'in', length=4, width=0.5, colors='k', labelsize = 15)
    plt.show()

        
        
        
    