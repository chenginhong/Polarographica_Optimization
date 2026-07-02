# coding: utf-8

#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================
#Author:       Tim Tichter
#Date:         2019.07.16
#Function:     POLAROGRAPHICAS Cyclic Voltammetry Functions

#=======================================================================================================================
#=======================================================================================================================
#importing all required modules from Python
#=======================================================================================================================
#=======================================================================================================================
from tkinter                           import *
from tkinter.filedialog                import askopenfilename
from tkinter.filedialog                import asksaveasfilename
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases          import key_press_handler
from matplotlib.figure                 import Figure
from scipy.interpolate                 import interp1d
from scipy.interpolate                 import InterpolatedUnivariateSpline
from scipy.optimize                    import curve_fit
from scipy.optimize                    import nnls
from scipy.linalg                      import toeplitz
from scipy.special                     import kv, iv, gamma
from scipy.integrate                   import quad
from scipy.fft                         import fft, ifft
from scipy.signal                      import hilbert
from cmath                             import *
from math                              import ceil,floor
import mpmath as mp
mp.mp.dps = 25
mp.mp.pretty = True

#=======================================================================================================================
#=======================================================================================================================
#define some global variables, functions and warning functions
#=======================================================================================================================
#=======================================================================================================================

global F
global R
global File_Was_Loaded
F = 96485.0
R = 8.314
File_Was_Loaded = 0
def cot(phi):
    return 1.0/tan(phi)
def csc(phi):
    return 1.0/sin(phi)
def coth(x):
    return 1/tanh(x)

#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================



# In[ ]:

def Open_ImpFit_File():
    root = Toplevel()
    root.title("Your Impedance Data")
    global File_Was_Loaded
    File_Was_Loaded = 1
    global ExpData   #needed later
    ExpData      = 1
    global data
    if Delimiter == 1:
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='\t', skip_header = FromRowxx, skip_footer = ToRowxx)
    if Delimiter == 0:     
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='', skip_header = FromRowxx, skip_footer = ToRowxx)
    #============================================================================
    #Plotting of loaded file
    #============================================================================ 
    f  = Figure(figsize=(9, 4), dpi=100)        
    b1 = f.add_subplot(121)
    b2 = f.add_subplot(122)   
    #============================================================================
    #define entry window
    #============================================================================
    #============================================================================ 
    global Z_REAL
    global Z_IMAG
    global Frequency_Array
    #==================================================================================
    Frequency         =  data[::Readeveryxx,0:1:1] 
    Z_REAL            =  data[::Readeveryxx,1:2:1] * UmrPot            
    Z_IMAG            =  data[::Readeveryxx,2:3:1] * UmRZ_IMAG
    #==================================================
    #turn frequencies alsways from high to low
    if Frequency[1] > Frequency[0]:            
        Frequency     = Frequency[::-1]
        Z_REAL     = Z_REAL[::-1]             
        Z_IMAG         = Z_IMAG[::-1]
    global Weite
    Weite = Z_REAL.shape[0]
    global NumMeas
    NumMeas = Z_REAL.shape[1]
    for i in range(Weite):
        Z_IMAGarrays            = np.squeeze(np.transpose(Z_IMAG[0::,i:i+1:]))                          
        Z_REALarrays            = np.squeeze(np.transpose(Z_REAL[0::,i:i+1:])) 
        LenPotArrays            = len(Z_REALarrays)
        Frequency_Array         = Frequency  
    #=========================================================
    #final plotting
    #=========================================================
    b1.set_title("Nyquist-Plot")
    b1.plot(Z_REAL, -Z_IMAG, linestyle='-',marker='',color='k')   
    b1.set_xlabel('$\mathfrak{RE}(Z(j\omega)/\Omega)$', fontsize = 15)
    b1.set_ylabel('$-\mathfrak{IM}(Z(j\omega)/\Omega)$', fontsize = 15)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b1.set_aspect('equal', adjustable = 'datalim')
    b2.set_title("Bode-Plot")    
    b2.plot(np.log10(Frequency), -Z_IMAG, linestyle='-',marker='',color='b', label = '-Z.imag')
    b2.plot(np.log10(Frequency), Z_REAL, linestyle='-',marker='',color='r',label = 'Z.real')
    b2.set_xlabel('log$_{10}\,(\\nu\,/\,\mathrm{s}^{-1})$', fontsize=15)
    b2.set_ylabel('$-\mathfrak{IM}(Z(j\omega)/\Omega)$ and $\mathfrak{RE}(Z(j\omega)/\Omega)$', fontsize = 15)
    b2.plot()
    b2.legend(frameon = False, fontsize = 15)
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)        
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)




# In[ ]:

def Get_ImpFit_Data():
    Fenster = Toplevel()                                                         
    Fenster.title("Get-Data")                         
    Fenster.geometry("400x270")                                            

    FromRowxxx_label      = Label(Fenster,text="Skip Header")
    FromRowxxx_label.grid(row=0, column=0)
    FromRowxxx_Eingabe    = Entry(Fenster)                                               
    FromRowxxx_Eingabe.grid(row=0, column=1)

    ToRowxxx_label       = Label(Fenster,text="Skip Footer")
    ToRowxxx_label.grid(row=1, column=0)
    ToRowxxx_Eingabe     = Entry(Fenster)
    ToRowxxx_Eingabe.grid(row=1, column=1)
    
    Readeveryxxx_label   = Label(Fenster,text="Read every")
    Readeveryxxx_label.grid(row=2, column=0)
    Readeveryxxx_Eingabe = Entry(Fenster)
    Readeveryxxx_Eingabe.grid(row=2, column=1)
    
    
    UmrPot_Label = Label(Fenster,text="Z_real Factor to be Ohm")
    UmrPot_Label.grid(row=3, column=0)
    UmrPot_Eingabe = Entry(Fenster)
    UmrPot_Eingabe.grid(row=3, column=1)
        
    UmRZ_IMAG_Label = Label(Fenster,text="Z_imag Factor to be Ohm")
    UmRZ_IMAG_Label.grid(row=4, column=0)
    UmRZ_IMAG_Eingabe = Entry(Fenster)
    UmRZ_IMAG_Eingabe.grid(row=4, column=1)
    
    HasToBe_Label  = Label(Fenster,text="Data Order")
    HasToBe_Label.grid(row=6, column=0)
    HasToBe_Label1  = Label(Fenster,text="Freq___Z.real___Z.Imag")
    HasToBe_Label1.grid(row=6, column=1)
        
    Delimiterxxx_Label  = Label(Fenster,text="Delimiter")
    Delimiterxxx_Label.grid(row=7, column=0)
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Tab", variable=var5).grid(row=7,column=1, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Space", variable=var6).grid(row=7,column=2, sticky=W)



    def Accept_ImpFit():
              
        global FromRowxx
        FromRowxx     = (int(FromRowxxx_Eingabe.get()))
        global ToRowxx
        ToRowxx       = (int(ToRowxxx_Eingabe.get()))
        global Readeveryxx
        Readeveryxx   = (int(Readeveryxxx_Eingabe.get()))
        global UmRZ_IMAG#
        UmRZ_IMAG      = (float(UmRZ_IMAG_Eingabe.get()))
        global UmrPot#
        UmrPot        = (float(UmrPot_Eingabe.get()))
        
        global DesiReactxx
        global FirstComesxx

        global Delimiter
        Delimiter         = var5.get()
    
    
    def Next_ImpFit():
        Open_ImpFit_File()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    Accept = Button(Fenster, text="Accept",command=Accept_ImpFit)
    Accept.grid(row=8, column=0) 
    
    Next = Button(Fenster, text="Next",command=Next_ImpFit)
    Next.grid(row=9, column=0)  


# In[ ]:
    
    
def No_Loaded_File_Warner():
    Fenster = Toplevel()                                                         
    Fenster.title("No loaded file found")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "Caution! You loaded no file that\nshould be analyzed. Go to\nopen file first and select\na file that should be analyzed.")


#---------------------------------------------------------------------------------------
#IMPEDANCE-STUFF
#---------------------------------------------------------------------------------------

def Semi_Inf_Plan_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 1
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()
def Fin_Plan_Trans_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 1
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()
def Fin_Plan_Ref_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 1
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()
def Semi_Inf_Zyl_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 1
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()
def Semi_Inf_Sphe_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 1
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()
def Finit_Int_Sphe_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 1
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()


#-----------------------------------------------------

def Semi_Inf_Plan_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 1
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()
def Fin_Plan_Trans_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 1
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()
def Fin_Plan_Ref_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 1
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()
def Semi_Inf_Zyl_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 1
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()
def Semi_Inf_Sphe_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 1
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()
def Finit_Int_Sphe_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 1
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()




# In[ ]:

def Imp_Entry_Window():
    
    Fenster = Toplevel()  
    if IMPFIT == 0:
        Fenster.geometry("400x600")
    if IMPFIT == 1:
        Fenster.geometry("400x600")
        
    if Imp1 == 1:
        Fenster.title("Semi-Infinite Planar Diffusion Impedance")
    if Imp2 == 1:
        Fenster.title("Finite Planar Transmissive Diffusion Impedance")
    if Imp3 == 1:
        Fenster.title("Finite Planar Refelctive Diffusion Impedance")
    if Imp4 == 1:
        Fenster.title("Semi-Infinite Cylindrical Diffusion Impedance")
    if Imp5 == 1:
        Fenster.title("Semi-Infinite Spherical Diffusion Impedance")
    if Imp6 == 1:
        Fenster.title("Finite Spherical internal Diffusion Impedance")
        
  
    
    n_Label = Label(Fenster, text="n*")
    n_Label.grid(row=0, column=0)
    n_Eingabe = Entry(Fenster)     
    n_Eingabe.insert(END, 1)                                           
    n_Eingabe.grid(row=0, column=1)
        
    
    alpha_Label = Label(Fenster, text="alpha*")
    alpha_Label.grid(row=1, column=0)
    alpha_Eingabe = Entry(Fenster)    
    alpha_Eingabe.insert(END, 0.5)                                            
    alpha_Eingabe.grid(row=1, column=1)
    
    
    kzero_Label = Label(Fenster,text="k_zero [cm/s]*")
    kzero_Label.grid(row=2, column=0)
    kzero_Eingabe = Entry(Fenster)
    kzero_Eingabe.insert(END, 0.001) 
    kzero_Eingabe.grid(row=2, column=1)
    
    
    Ezero_Label = Label(Fenster,text="E_zero vs Ref [V]*")
    Ezero_Label.grid(row=3, column=0)
    Ezero_Eingabe = Entry(Fenster)
    Ezero_Eingabe.insert(END, 0) 
    Ezero_Eingabe.grid(row=3, column=1)
    
    Eeq_Label = Label(Fenster,text="E_eq vs Ref [V]*")
    Eeq_Label.grid(row=4, column=0)
    Eeq_Eingabe = Entry(Fenster)
    Eeq_Eingabe.insert(END, 0) 
    Eeq_Eingabe.grid(row=4, column=1)
    
    T_Label = Label(Fenster,text="T [°C]*")
    T_Label.grid(row=5, column=0)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.insert(END, 25) 
    T_Eingabe.grid(row=5, column=1)
    
    if IMPFIT == 0:
    
        Freq_High_Label = Label(Fenster,text="Highest Frequency [Hz]*")
        Freq_High_Label.grid(row=6, column=0)
        Freq_High_Eingabe = Entry(Fenster)
        Freq_High_Eingabe.insert(END, 1000000) 
        Freq_High_Eingabe.grid(row=6, column=1)

        Freq_Low_Label = Label(Fenster,text="Lowest Frequency [Hz]*")
        Freq_Low_Label.grid(row=7, column=0)
        Freq_Low_Eingabe = Entry(Fenster)
        Freq_Low_Eingabe.insert(END, 0.1) 
        Freq_Low_Eingabe.grid(row=7, column=1)

        Freq_Num_Label = Label(Fenster,text="Number of Freq's*")
        Freq_Num_Label.grid(row=8, column=0)
        Freq_Num_Eingabe = Entry(Fenster)
        Freq_Num_Eingabe.insert(END, 100) 
        Freq_Num_Eingabe.grid(row=8, column=1)

    
    if Imp2 == 1:
        d_Label = Label(Fenster,text="distance 10^-6[m]*")
        d_Label.grid(row=9, column=0)
        d_Eingabe = Entry(Fenster)
        d_Eingabe.insert(END, 20) 
        d_Eingabe.grid(row=9, column=1)
    
    if Imp3 == 1:
        d_Label = Label(Fenster,text="distance 10^-6[m]*")
        d_Label.grid(row=9, column=0)
        d_Eingabe = Entry(Fenster)
        d_Eingabe.insert(END, 20) 
        d_Eingabe.grid(row=9, column=1)
    
    if Imp4 == 1:
        r_Label = Label(Fenster,text="radius 10^-6[m]*")
        r_Label.grid(row=9, column=0)
        r_Eingabe = Entry(Fenster)
        r_Eingabe.insert(END, 20) 
        r_Eingabe.grid(row=9, column=1)
    
    if Imp5 == 1:
        r_Label = Label(Fenster,text="radius 10^-6[m]*")
        r_Label.grid(row=9, column=0)
        r_Eingabe = Entry(Fenster)
        r_Eingabe.insert(END, 20) 
        r_Eingabe.grid(row=9, column=1)
        
    if Imp6 == 1:
        r_Label = Label(Fenster,text="radius 10^-6[m]*")
        r_Label.grid(row=9, column=0)
        r_Eingabe = Entry(Fenster)
        r_Eingabe.insert(END, 20) 
        r_Eingabe.grid(row=9, column=1)
        
    
   
    Ru_Label = Label(Fenster,text="Ru in Ohm*")
    Ru_Label.grid(row=10, column=0)
    Ru_Eingabe = Entry(Fenster)
    Ru_Eingabe.insert(END, 0.5) 
    Ru_Eingabe.grid(row=10, column=1)
    
    Capacity_Label = Label(Fenster,text="C in Microfarad*")
    Capacity_Label.grid(row=11, column=0)
    Capacity_Eingabe = Entry(Fenster)
    Capacity_Eingabe.insert(END, 50) 
    Capacity_Eingabe.grid(row=11, column=1)
    
    A_Label = Label(Fenster,text="A in cm^2*")
    A_Label.grid(row=12, column=0)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.insert(END, 1) 
    A_Eingabe.grid(row=12, column=1)
        
    c_Label = Label(Fenster,text="c_total in mol/L*")
    c_Label.grid(row=13, column=0)
    c_Eingabe = Entry(Fenster)
    c_Eingabe.insert(END, 0.1) 
    c_Eingabe.grid(row=13, column=1)
        
    D_Label = Label(Fenster,text="D 10^-6[cm^2/s]*")
    D_Label.grid(row=14, column=0)
    D_Eingabe = Entry(Fenster)
    D_Eingabe.insert(END, 1) 
    D_Eingabe.grid(row=14, column=1)
    
    var1 = IntVar()
    Checkbutton(Fenster, text="Unequal D", variable=var1).grid(row=15, column=0, sticky=W)
    
    Dred_Label = Label(Fenster,text="D_red 10^-6[cm^2/s]")
    Dred_Label.grid(row=16, column=0)
    Dred_Eingabe = Entry(Fenster)
    Dred_Eingabe.grid(row=16, column=1)
    
    Dox_Label = Label(Fenster,text="D_ox 10^-6[cm^2/s]")
    Dox_Label.grid(row=17, column=0)
    Dox_Eingabe = Entry(Fenster)
    Dox_Eingabe.grid(row=17, column=1)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Preceeding ch eq", variable=var2).grid(row=18, column=0, sticky=W)
    
    k_prec_f_Label = Label(Fenster,text="kf_chem_preceding [1/s]")
    k_prec_f_Label.grid(row=19, column=0)
    k_prec_f_Eingabe = Entry(Fenster)
    k_prec_f_Eingabe.grid(row=19, column=1)
    
    k_prec_b_Label = Label(Fenster,text="kb_chem_preceding [1/s]")
    k_prec_b_Label.grid(row=20, column=0)
    k_prec_b_Eingabe = Entry(Fenster)
    k_prec_b_Eingabe.grid(row=20, column=1)
    
    var3 = IntVar()
    Checkbutton(Fenster, text="Succeding ch eq", variable=var3).grid(row=21, column=0, sticky=W)
    
    k_succ_f_Label = Label(Fenster,text="kf_chem_succeding [1/s]")
    k_succ_f_Label.grid(row=22, column=0)
    k_succ_f_Eingabe = Entry(Fenster)
    k_succ_f_Eingabe.grid(row=22, column=1)
    
    k_succ_b_Label = Label(Fenster,text="kb_chem_succeding [1/s]")
    k_succ_b_Label.grid(row=23, column=0)
    k_succ_b_Eingabe = Entry(Fenster)
    k_succ_b_Eingabe.grid(row=23, column=1)
    
    
    var4 = IntVar()
    Checkbutton(Fenster, text="Save as txt", variable=var4).grid(row=24, column=0, sticky=W)
    
    
    #-------------------------------------------------------
    #Akzeptierfunktion schreiben
    #-------------------------------------------------------
    
    def AcceptParams():
        global n
        global alpha
        global kzero
        global Ezero
        global Eeq
        global T
        global distance
        global r
        global precedeq
        global succedeq
        global Unequal_D
        global AstxtSaver
        global A
        global c
        global Ru
        global FreqUp
        global FreqLow
        global FreqNum
        global Cap
    
        
    #------------------------------------------------------------
    #Felder lesen
    #------------------------------------------------------------
        
        n           = (float(n_Eingabe.get()))
        alpha       = (float(alpha_Eingabe.get()))
        kzero       = (float(kzero_Eingabe.get()))
        T           = (float(T_Eingabe.get())) + 273.15
        Cap         = (float(Capacity_Eingabe.get()))*0.000001
        
        if IMPFIT == 0:
            FreqUp      = np.log10(float(Freq_High_Eingabe.get()))
            FreqLow     = np.log10(float(Freq_Low_Eingabe.get()))
            FreqNum     = int(float(Freq_Num_Eingabe.get()))
        
    
        if Imp2 == 1:
            distance     = (float(d_Eingabe.get()))*0.0001
        if Imp3 == 1:
            distance     = (float(d_Eingabe.get()))*0.0001
        if Imp4 == 1:
            r            = (float(r_Eingabe.get()))*0.0001
        if Imp5 == 1:
            r            = (float(r_Eingabe.get()))*0.0001
        if Imp6 == 1:
            r            = (float(r_Eingabe.get()))*0.0001
    
         
        global Dox
        global Dred
        
        Unequal_D   = var1.get()
        if Unequal_D ==0:
            Dox       = (float(D_Eingabe.get()))*0.000001
            Dred      = (float(D_Eingabe.get()))*0.000001
        
        if Unequal_D ==1:
            Dox      = (float(Dox_Eingabe.get()))*0.000001
            Dred     = (float(Dred_Eingabe.get()))*0.000001
        
    
        Ru     = (float(Ru_Eingabe.get()))
        c      = (float(c_Eingabe.get()))*0.001
        A      = (float(A_Eingabe.get()))
        Ezero  = (float(Ezero_Eingabe.get()))
        Eeq    = (float(Eeq_Eingabe.get()))
        
        precedeq = var2.get()
        if precedeq == 1:
            global preced_kf
            global preced_kb
            preced_kf = (float(k_prec_f_Eingabe.get()))
            preced_kb = (float(k_prec_b_Eingabe.get()))
            
        succedeq = var3.get()
        if succedeq == 1:
            global succed_kf
            global succed_kb
            succed_kf = (float(k_succ_f_Eingabe.get()))  
            succed_kb = (float(k_succ_b_Eingabe.get())) 
            
        AstxtSaver  = var4.get()
        
        
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=25,column=1)
    button=Button(Fenster,text="Next",command=Imp_Type_Chooser).grid(row=26,column=1)
    
    
    


# In[ ]:

def Imp_Type_Chooser():
    if Imp1 == 1:
        if IMPFIT == 0:
            PlanarSemiinfImpCalculator()
        if IMPFIT == 1:
            ImpStartEndDefiner()
            PlanarSemiinfImpCalculator()
            
    if Imp2 == 1:
        if IMPFIT == 0:
            PlanarTransmImpCalculator()
        if IMPFIT == 1:
            ImpStartEndDefiner()
            PlanarTransmImpCalculator()
        
    if Imp3 == 1:
        if IMPFIT == 0:
            PlanarReflImpCalculator()
        if IMPFIT == 1:
            ImpStartEndDefiner()
            PlanarReflImpCalculator()
        
    if Imp4 == 1:
        if IMPFIT == 0:
            ZylSemiinfImpCalculator()        
        if IMPFIT == 1:
            ImpStartEndDefiner()
            ZylSemiinfImpCalculator()
        
    if Imp5 == 1:
        if IMPFIT == 0:
            SpherSemiinfImpCalculator()
        if IMPFIT == 1:
            ImpStartEndDefiner()
            SpherSemiinfImpCalculator()
        
    if Imp6 == 1:
        if IMPFIT == 0:
            SpherIntFinImpCalculator()
        if IMPFIT == 1:
            ImpStartEndDefiner()
            SpherIntFinImpCalculator()
    
    


# In[ ]:

def ImpStartEndDefiner():
    global Z_im_Meas
    global Z_re_Meas
    Z_im_Meas = np.squeeze(Z_IMAG)        #Nach dem Laden heißen die Arrays nur so wie früher
    Z_re_Meas = np.squeeze(Z_REAL)    #deshalb hier die Umbenennung
    


# In[ ]:

def PlanarSemiinfImpCalculator():
    
    root = Toplevel()
    root.title("Impedance plots")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]

    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp(-(1-alpha)*n*F*(Eeq-EZero)/(R*T))


    ZFar = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +kf_el/((1+Kp)*((Freq+lp)*Dred)**0.5)   +kf_el*Kp/((1+Kp)*((Freq)*Dred)**0.5)       +kb_el*Ks/((1+Ks)*((Freq+ls)*Dox)**0.5)   +kb_el/((1+Ks)*((Freq)*Dox)**0.5)      )                              

    #----------------------------------------------------------------------------------------------
   
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
        #==============================================
        #calculate Standard deviation
        #==============================================
        
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[-i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[-i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        
        global ReStDev
        global ImStDev
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        

        
        #==============================================
        
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b1.set_aspect('equal', adjustable = 'datalim')
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm  ', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    


# In[ ]:

def PlanarTransmImpCalculator():
    root = Toplevel()
    root.title("Impedance plots")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]
    
    
    
    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    global d
    
    d       = distance
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp(-(1-alpha)*n*F*(Eeq-EZero)/(R*T))


    ZFar = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +kf_el*np.tanh(d*((Freq+lp)/Dred)**0.5)/((1+Kp)*((Freq+lp)*Dred)**0.5)   +kf_el*Kp*np.tanh(d*(Freq/Dred)**0.5)/((1+Kp)*((Freq)*Dred)**0.5)       +kb_el*Ks*np.tanh(d*((Freq+ls)/Dox)**0.5)/((1+Ks)*((Freq+ls)*Dox)**0.5)   +kb_el*np.tanh(d*(Freq/Dox)**0.5)/((1+Ks)*((Freq)*Dox)**0.5)      )                              
    
    #----------------------------------------------------------------------------------------------
   
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
        
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        global ReStDev
        global ImStDev
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        
        
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b1.set_aspect('equal', adjustable = 'datalim')
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    
    


# In[ ]:

def PlanarReflImpCalculator():
    root = Toplevel()
    root.title("Impedance plots")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]
    
    
    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    global d
    
    d       = distance
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp(-(1-alpha)*n*F*(Eeq-EZero)/(R*T))


    def coth(x):
        return (np.tanh(x))**(-1)
    ZFar = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +kf_el*coth(d*((Freq+lp)/Dred)**0.5)/((1+Kp)*((Freq+lp)*Dred)**0.5)   +kf_el*Kp*coth(d*(Freq/Dred)**0.5)/((1+Kp)*((Freq)*Dred)**0.5)       +kb_el*Ks*coth(d*((Freq+ls)/Dox)**0.5)/((1+Ks)*((Freq+ls)*Dox)**0.5)   +kb_el*coth(d*(Freq/Dox)**0.5)/((1+Ks)*((Freq)*Dox)**0.5)      )                              

    #----------------------------------------------------------------------------------------------
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
        #==============================================
        #calculate Standard deviation
        #==============================================
        
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        
        global ReStDev
        global ImStDev
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        
          
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b1.set_aspect('equal', adjustable = 'datalim')
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    
    


# In[ ]:

def ZylSemiinfImpCalculator():
    root = Toplevel()
    root.title("Impedance plots")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]
    
    
    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp(-(1-alpha)*n*F*(Eeq-EZero)/(R*T))


    
    def Bessel_P_eq(x):
        return (kv(0,(r*((x+lp)/Dred)**0.5)))/(kv(1,(r*((x+lp)/Dred)**0.5))*((x+lp)*Dred)**0.5)
    def Bessel_P_n(x):
        return (kv(0,(r*((x)/Dred)**0.5)))/(kv(1,(r*((x)/Dred)**0.5))*((x)*Dred)**0.5)
    def Bessel_S_eq(x):
        return (kv(0,(r*((x+ls)/Dox)**0.5)))/(kv(1,(r*((x+ls)/Dox)**0.5))*((x+ls)*Dox)**0.5)
    def Bessel_S_n(x):
        return (kv(0,(r*((x)/Dox)**0.5)))/(kv(1,(r*((x)/Dox)**0.5))*((x)*Dox)**0.5)


    ZFar = np.empty(len(Freq),dtype = complex)
    for i in range(len(ZFar)):
        if np.abs(kv(1,(r*(Freq[i]/Dox)**0.5)).imag) != 0 :
            if np.abs(kv(1,(r*(Freq[i]/Dred)**0.5)).imag) != 0 :
                ZFar[i] = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +(kf_el/(1+Kp))*Bessel_P_eq(Freq[i])   +(kf_el*Kp/(1+Kp))*Bessel_P_n(Freq[i])  +(kb_el*Ks/(1+Ks))*Bessel_S_eq(Freq[i])   +(kb_el/(1+Ks))*Bessel_S_n(Freq[i]) )                             
        if np.abs(kv(1,(r*(Freq[i]/Dox)**0.5))) == 0: 
            ZFar[i] = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +(kf_el/((1+Kp)*((Freq[i]+lp)*Dred)**0.5)) +(kf_el*Kp/((1+Kp)*((Freq[i])*Dred)**0.5))  +(kb_el*Ks/((1+Ks)*((Freq[i]+lp)*Dox)**0.5))   +(kb_el/((1+Ks)*((Freq[i])*Dox)**0.5)) )                             
        if np.abs(kv(1,(r*(Freq[i]/Dred)**0.5))) == 0:  
            ZFar[i] = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +(kf_el/((1+Kp)*((Freq[i]+lp)*Dred)**0.5)) +(kf_el*Kp/((1+Kp)*((Freq[i])*Dred)**0.5))  +(kb_el*Ks/((1+Ks)*((Freq[i]+lp)*Dox)**0.5))   +(kb_el/((1+Ks)*((Freq[i])*Dox)**0.5)) )                             
    
    #----------------------------------------------------------------------------------------------
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
        
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        global ReStDev
        global ImStDev
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        
       
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b1.set_aspect('equal', adjustable = 'datalim')
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    


# In[ ]:

def SpherSemiinfImpCalculator():
    root = Toplevel()
    root.title("Impedance plots")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]
    
    
    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp(-(1-alpha)*n*F*(Eeq-EZero)/(R*T))


    ZFar = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +kf_el/((1+Kp)*(Dred/r + ((Freq+lp)*Dred)**0.5))   +kf_el*Kp/((1+Kp)*(Dred/r + ((Freq)*Dred)**0.5))       +kb_el*Ks/((1+Ks)*((Dox/r + ((Freq+ls)*Dox)**0.5) ))   +kb_el/((1+Ks)*((Dox/r + ((Freq)*Dox)**0.5) ))      )                              

    #----------------------------------------------------------------------------------------------
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
        
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        global ReStDev
        global ImStDev
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        
        
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b1.set_aspect('equal', adjustable = 'datalim')
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    


# In[ ]:

def SpherIntFinImpCalculator():
    root = Toplevel()
    root.title("Impedance plots")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]
            
            
    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp(-(1-alpha)*n*F*(Eeq-EZero)/(R*T))


    def coth(x):
        return (np.tanh(x))**(-1)

    ZFar = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +kf_el/(coth(r*((Freq+lp)/Dred)**0.5)*(1+Kp)*(Dred/r + ((Freq+lp)*Dred)**0.5))   +kf_el*Kp/(coth(r*(Freq/Dred)**0.5)*(1+Kp)*(Dred/r + ((Freq)*Dred)**0.5))       +kb_el*Ks/(coth(r*((Freq+lp)/Dox)**0.5)*(1+Ks)*((Dox/r + ((Freq+ls)*Dox)**0.5) ))   +kb_el/(coth(r*(Freq/Dox)**0.5)*(1+Ks)*((Dox/r + ((Freq)*Dox)**0.5) ))      )                              

    #----------------------------------------------------------------------------------------------
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
    
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        global ReStDev
        global ImStDev
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b1.set_aspect('equal', adjustable = 'datalim')
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    


# In[ ]:

def Simulated_IMPEDANZ_as_txt_saver():
    
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")
        
        
        if Imp1 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Planar Semi-infinit")
            fi.write("\n")
        if Imp2 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Planar Finit Transmissive")
            fi.write("\n")
        if Imp3 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Planar Finit Reflective")
            fi.write("\n")
        if Imp4 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Cylindrical semi infinte")
            fi.write("\n")
        if Imp5 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Spherical Semi-infinite")
            fi.write("\n")
        if Imp6 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Spherical Finite Internal")
            fi.write("\n")

        fi.write("n")
        fi.write("\t")
        fi.write(str(n))
        fi.write("\n")
        fi.write("T in [K]")
        fi.write("\t")
        fi.write(str(T))
        fi.write("\n")
        
        fi.write("D_red in [cm^2/s]")
        fi.write("\t")
        fi.write(str(Dred))
        fi.write("\n")
        fi.write("Dox in [cm^2/s]")
        fi.write("\t")
        fi.write(str(Dox))
        fi.write("\n")
        if precedeq == 1:
            #preceding
            fi.write("k_f_preceeding [1/s]")
            fi.write("\t")
            fi.write(str(kfp))
            fi.write("\n")
            fi.write("k_b_preceeding [1/s]")
            fi.write("\t")
            fi.write(str(kbp))
            fi.write("\n")
        if succedeq == 1:
            #succeding
            fi.write("k_f_succeeding [1/s]")
            fi.write("\t")
            fi.write(str(kfs))
            fi.write("\n")
            fi.write("k_b_succeeding [1/s]]")
            fi.write("\t")
            fi.write(str(kbs))
            fi.write("\n")      
        
        
        fi.write("alpha")
        fi.write("\t")
        fi.write(str(alpha))
        fi.write("\n")
        fi.write("k_zero in [cm/s]")
        fi.write("\t")
        fi.write(str(kzero))
        fi.write("\n")
            

        if Imp2 == 1:
            fi.write("d in [cm]")
            fi.write("\t")
            fi.write(str(distance))
            fi.write("\n")
        if Imp3 == 1:
            fi.write("d in [cm]")
            fi.write("\t")
            fi.write(str(distance))
            fi.write("\n")
        if Imp4 == 1:
            fi.write("radius in [cm]")
            fi.write("\t")
            fi.write(str(r))
            fi.write("\n")
        if Imp5 == 1:
            fi.write("radius in [cm]")
            fi.write("\t")
            fi.write(str(r))
            fi.write("\n")
        if Imp6 == 1:
            fi.write("radius in [cm]")
            fi.write("\t")
            fi.write(str(r))
            fi.write("\n")

        fi.write("E_zero vs. E_ref in [V]")
        fi.write("\t")
        fi.write(str(Ezero))
        fi.write("\n")
        fi.write("Eeq vs. E_ref in [V]")
        fi.write("\t")
        fi.write(str(Eeq))
        fi.write("\n")
        fi.write("Ru in [Ohm]")        
        fi.write("\t")
        fi.write(str(Ru))
        fi.write("\n")
        fi.write("A in [cm^2]")
        fi.write("\t")
        fi.write(str(A))
        fi.write("\n")
        fi.write("c in [mol/L]")
        fi.write("\t")
        fi.write(str(c*1000.0))
        
        
        for i in range(3):
            fi.write("\n")
        
        fi.write("Nyquist DATA")
        
        if IMPFIT == 0:
            fi.write("\n")
            fi.write("Frequency [Hz]")
            fi.write("\t")
            fi.write("Z_real_Calc [Ohm]")
            fi.write("\t")
            fi.write("Z_imag_Calc [Ohm]")
            fi.write("\t")
            fi.write("\n")
            for i in range(len(xxx)):
                fi.write(str(np.asscalar((Freq[i].imag)/(2*np.pi))))
                fi.write("\t")
                fi.write(str(np.asscalar(xxx[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(f[i])))
                fi.write("\t")
                fi.write("\n")

        
        if IMPFIT == 1:
            fi.write("\n")
            fi.write("\n")
            fi.write("\n")
            fi.write("Standard dev real part")
            fi.write("\t")
            fi.write(str(np.asscalar(ReStDev)))
            fi.write("\n")
            fi.write("Standard dev imag part")
            fi.write("\t")
            fi.write(str(np.asscalar(ImStDev)))
            fi.write("\n")
            fi.write("\n")
            fi.write("\n")
            
            fi.write("Frequency [Hz]")
            fi.write("\t")
            fi.write("Z_real_Calc [Ohm]")
            fi.write("\t")
            fi.write("Z_imag_Calc [Ohm]")
            fi.write("\t")
            fi.write("Z_real_Meas [Ohm]")
            fi.write("\t")
            fi.write("Z_imag_Meas [Ohm]")
            fi.write("\t")
            fi.write("\n")
            
            for i in range(len(Z_real_Meas)):
                fi.write(str(np.asscalar(Frequency_Array[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_real_Calc[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_imag_Calc[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_real_Meas[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_imag_Meas[i])))
                fi.write("\t")
                fi.write("\n")
    
    root.destroy()


# In[ ]:




# In[ ]:




# In[ ]:



