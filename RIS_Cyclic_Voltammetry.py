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


def No_Loaded_File_Warner():
    Fenster = Toplevel()                                                         
    Fenster.title("No loaded file found")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "Caution! You loaded no file that\nshould be fittet. Go to\nopen file first and select\na file that should be fitted.")

def Parameters_do_not_fit_warner():
    Fenster = Toplevel()                                                         
    Fenster.title("Entered Parameters are out of limits!")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "Caution! The entered Parameters\ncannot be used for the\ncalculation of a CV\nIt should be always\nD>10^-6 cm^2/s,  0.00001 cm < distance < 10 cm \nand r>0.1 micrometer!")

def Incorrect_Range_Warner():
    Fenster = Toplevel()                                                         
    Fenster.title("Warning! Incorrect numerical range.")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "Potential limits are not sufficiently\nfar away from E_zero. It should be\nat least\nXi_initial = -7 = nF(E_i-E_zero)/RT \nNevertheless, CV got calculated.")

def Counter_Error():
    Fenster = Toplevel()                                                         
    Fenster.title("Warning! Incorrect numerical range.")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "Counter Error! Your are out of the range of\nThe maximum numbe (50) or minimum number (1)\nof CVs. Go either to Clear- or Remove- option.")

def Unequal_length_warner():
    Fenster = Toplevel()                                                         
    Fenster.title("Warning! Incorrect numerical range.")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "Warning! Superposition-axis of data\nhave unequal length. This may lead\nTo unexpeted/unwanted behaviour!\nSuperposition may - or may not - be calculated.")

#=======================================================================================================================
#=======================================================================================================================
#define Open File Function for opening RIS_input file containing t and the random E(t)
#=======================================================================================================================
#=======================================================================================================================

def RIS_Definer():
    root = Toplevel()
    root.title("Your Data")
    global Input_Time
    global Input_Potential
    Path = askopenfilename()
    data = np.genfromtxt(Path, delimiter='\t')
    Input_Time        =  np.squeeze(data[::,0:1:1])
    Input_Potential   =  np.squeeze(data[::,1:2:1])           
    root.destroy()


#=======================================================================================================================
#=======================================================================================================================
#Superposition mode Function
#=======================================================================================================================
#=======================================================================================================================
def RISCV_Initializer():
    global Superposition_Mode
    global Counter
    Counter = 0
    Superposition_Mode = 1
    global Storage_x
    global Storage_y
    Storage_x = PotCalc
    Storage_y = np.zeros([len(PotCalc),50])
    
def RISCV_Append_Latest():
    global Counter
    global Storage_x
    global Storage_y
    Counter += 1
    if Counter >=1 and Counter <= 50:
        if len(Storage_x) != len(CurrCalc):
            Unequal_length_warner()
        Storage_y[::,(Counter-1)] = CurrCalc[0:len(Storage_x)]
    else:
        Counter_Error()
        
def RISCV_Remove_Latest():
    global Storage_y
    global Counter
    if Counter >=1 and Counter <= 50:
         Storage_y[::,(Counter-1)] = 0
    else:
        Counter_Error()
    Counter -= 1
    
def RISCV_E_vs_I_RS_Mode_Caller():
    global RS
    global ADD
    global t_vs_I
    t_vs_I = 0
    RS  = 1
    ADD = 0  
    RISCV_Superpos_Shower()

def RISCV_E_vs_I_Add_Mode_Caller():
    global RS
    global ADD
    global t_vs_I
    t_vs_I = 0
    RS  = 0
    ADD = 1  
    RISCV_Superpos_Shower()
    
def RISCV_E_vs_I_Double_Mode_Caller():
    global RS
    global ADD
    global t_vs_I
    t_vs_I = 0
    RS  = 1
    ADD = 1  
    RISCV_Superpos_Shower()
    
def RISCV_t_vs_I_RS_Mode_Caller():
    global RS
    global ADD
    global t_vs_I
    t_vs_I = 1
    RS  = 1
    ADD = 0  
    RISCV_Superpos_Shower()

def RISCV_t_vs_I_Add_Mode_Caller():
    global RS
    global ADD
    global t_vs_I
    t_vs_I = 1
    RS  = 0
    ADD = 1  
    RISCV_Superpos_Shower()
    
def RISCV_t_vs_I_Double_Mode_Caller():
    global RS
    global ADD
    global t_vs_I
    t_vs_I = 1
    RS  = 1
    ADD = 1  
    RISCV_Superpos_Shower()

def RISCV_Superpos_Shower():
    global Storage_x
    global Storage_y
    global yyy_RS
    global yyy_Add
    yyy_RS  = Storage_y[::,:Counter:]
    yyy_Add = np.squeeze(np.sum(yyy_RS, axis = 1)) 
    root = Toplevel()
    root.title("Simulated CV")
    if t_vs_I == 0:   
        Abbildung      = Figure(figsize=(5, 4), dpi=100)
        Superimpose    = Abbildung.add_subplot(111)
        if RS  == 1:
            Superimpose.plot(Storage_x[0:-1],yyy_RS[0:-1,::]) 
        if ADD == 1:
            Superimpose.plot(Storage_x[0:-1],yyy_Add[0:-1:], color = 'b')
        Superimpose.axvline(0,color='k')
        Superimpose.axhline(0,color='k')
    if t_vs_I == 1:    
        Abbildung      = Figure(figsize=(9, 4), dpi=100)
        Superimpose    = Abbildung.add_subplot(121)
        FFTSuperimpose = Abbildung.add_subplot(122)
        if RS  == 1:
            Superimpose.plot(Input_Time[0:-1],yyy_RS[0:-1,::]) 
        if ADD == 1:
            Superimpose.plot(Input_Time[0:-1],yyy_Add[0:-1:], color = 'b')
        FFTSuperimpose.plot(np.abs(fft(yyy_Add[0:-1:]))/np.max(np.abs(fft(yyy_Add[0:-1:]))))
        FFTSuperimpose.set_xlabel('Point no.', fontsize=12)
        FFTSuperimpose.set_ylabel('abs(fft(I / mA))/max(abs(fft(I / mA)))', fontsize=12)
        for axis in ['top','bottom','left','right']:
            FFTSuperimpose.spines[axis].set_linewidth(2)
            FFTSuperimpose.spines[axis].set_color('k')
        FFTSuperimpose.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        global Superposition_Mode
        Superposition_Mode =1
        IFFT_controller()
    
    
    Superimpose.set_xlabel('E vs. ' 'E$^0$' '(1$^{st}$ CV)' ' / V', fontsize=12)
    Superimpose.set_ylabel('I / mA', fontsize=12)
    for axis in ['top','bottom','left','right']:
        Superimpose.spines[axis].set_linewidth(2)
        Superimpose.spines[axis].set_color('k')
    Superimpose.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    Abbildung.tight_layout()
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)



def RISCV_Superpos_Saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))
    with open(root.filename,"w") as fi:
        fi.write("t in [s]")
        fi.write("\t")
        fi.write("Potential vs. E_zero in [V]")
        fi.write("\t")
        fi.write("Calculated Current in mA")
        if ADD == 1:
            fi.write("\n")
            fi.write("\n")
            for i in range(len(Storage_x)-1):
                fi.write(str(np.asscalar(Input_Time[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Storage_x[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(yyy_Add[i])))
                fi.write("\t")
                fi.write("\n")
        if RS == 1:
            fi.write("\n")
            fi.write("\n")
            fi.write("t in [s]")
            fi.write("\t")
            fi.write("Calculated Potential vs. E_zero in [V]")
            fi.write("\t")
            for j in range(len(yyy_RS[0,::])):
                fi.write("Indiv. Calc. Curr. in mA")
                fi.write("\t")
            fi.write("\n")
            for i in range(len(Storage_x)-1):
                fi.write(str(np.asscalar(Input_Time[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Storage_x[i])))
                fi.write("\t")
                for j in range(len(yyy_RS[0,::])):
                    fi.write(str(np.asscalar(yyy_RS[i,j])))
                    fi.write("\t")
                fi.write("\n")
                
    root.destroy()
                


    
#=======================================================================================================================
#=======================================================================================================================
#define CV Type Chooser and Model zeroer Functions. Model Zeroer is required that no interference between models will occur
#for RIS CV the same Mass Transfer Functions are required like for CV. So they are UNCHANGED
#The Difference is, that RIS_Eingabe instead of CV_Eingabe is called here and also the callinf functions are renamed
#=======================================================================================================================
#=======================================================================================================================


def Model_Zeroer():
    global model_a
    global model_b
    global model_b1
    global model_c
    global model_d
    global model_e
    global model_f
    global model_g
    global Statistical
    Statistical     = 0
    model_a         = 0
    model_b         = 0
    model_b1        = 0
    model_c         = 0
    model_d         = 0
    model_e         = 0
    model_f         = 0
    model_g         = 0
    global xxx
    global yyy
    global PsiArray
    global XiArray
    global CV_Interpolation_Hin
    global CV_Interpolation_Back


def RIS_Semi_Inf_Planar():
    Model_Zeroer()
    global model_a
    model_a = 1
    global FIT 
    FIT = 0
    Eingabe_RIS_CV_Simulator()
       
def RIS_Finit_Planar():
    Model_Zeroer()
    global model_b
    model_b = 1
    global FIT 
    FIT = 0
    Eingabe_RIS_CV_Simulator()
    
def RIS_Finit_Planar_Trans():
    Model_Zeroer()
    global model_b1
    model_b1 = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def RIS_Semi_Inf_Zyl_Ext():
    Model_Zeroer()
    global model_c
    model_c = 1
    global FIT 
    FIT = 0
    Eingabe_RIS_CV_Simulator()
    
def RIS_Finit_Zyl_Ext():
    Model_Zeroer()
    global model_d
    model_d = 1
    global FIT 
    FIT = 0
    Eingabe_RIS_CV_Simulator()
    
def RIS_Finit_Zyl_Int():
    Model_Zeroer()
    global model_e
    model_e = 1
    global FIT 
    FIT = 0
    Eingabe_RIS_CV_Simulator()
    
def RIS_Semi_Inf_Sphere_Ext():
    Model_Zeroer()
    global model_f
    model_f = 1
    global FIT 
    FIT = 0
    Eingabe_RIS_CV_Simulator()
    
def RIS_Finit_Sphere_Int():
    Model_Zeroer()
    global model_g
    model_g = 1
    global FIT 
    FIT = 0
    Eingabe_RIS_CV_Simulator()
    
    
    
#-------------------------------------------------------------------------------------
#Fitters
#-------------------------------------------------------------------------------------
    
    
def RIS_Semi_Inf_Planar_FITTER():
    Model_Zeroer()
    global model_a
    model_a = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def RIS_Finit_Planar_FITTER():
    Model_Zeroer()
    global model_b
    model_b = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
def RIS_Finit_Planar_Trans_FITTER():
    Model_Zeroer()
    global model_b1
    model_b1 = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def RIS_Semi_Inf_Zyl_Ext_FITTER():
    Model_Zeroer()
    global model_c
    model_c = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def RIS_Finit_Zyl_Ext_FITTER():
    Model_Zeroer()
    global model_d
    model_d = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
def RIS_Finit_Zyl_Int_FITTER():
    Model_Zeroer()
    global model_e
    model_e = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
def RIS_Semi_Inf_Sphere_Ext_FITTER():
    Model_Zeroer()
    global model_f
    model_f = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def RIS_Finit_Sphere_Int_FITTER():
    Model_Zeroer()
    global model_g
    model_g = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
        
        
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#Statistical Simulators
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

def RIS_Statistical_Finit_Planar():
    Model_Zeroer()
    global model_b
    global Statistical
    Statistical  = 1
    model_b = 1
    global FIT 
    FIT = 0
    Eingabe_RIS_CV_Simulator()
    
def RIS_Statistical_Finit_Zyl_Ext():
    Model_Zeroer()
    global model_d
    model_d = 1
    global Statistical
    Statistical  = 1
    global FIT 
    FIT = 0
    Eingabe_RIS_CV_Simulator()
    
def RIS_Statistical_Finit_Zyl_Int():
    Model_Zeroer()
    global model_e
    global Statistical
    Statistical  = 1
    model_e = 1
    global FIT 
    FIT = 0
    Eingabe_RIS_CV_Simulator()
    
def RIS_Statistical_Finit_Sphere_Int():
    Model_Zeroer()
    global model_g
    global Statistical
    Statistical  = 1
    model_g = 1
    global FIT 
    FIT = 0
    Eingabe_RIS_CV_Simulator()
    
    
    
#-------------------------------------------------------------------------------------
#Statistical Fitters
#-------------------------------------------------------------------------------------
    
    
def RIS_Statistical_Finit_Planar_FITTER():
    Model_Zeroer()
    global model_b
    global Statistical
    Statistical  = 1
    model_b = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
    
def RIS_Statistical_Finit_Zyl_Ext_FITTER():
    Model_Zeroer()
    global model_d
    global Statistical
    Statistical  = 1
    model_d = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
def RIS_Statistical_Finit_Zyl_Int_FITTER():
    Model_Zeroer()
    global model_e
    global Statistical
    Statistical  = 1
    model_e = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
    
def RIS_Statistical_Finit_Sphere_Int_FITTER():
    Model_Zeroer()
    global model_g
    global Statistical
    Statistical  = 1
    model_g = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_RIS_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()



# In[ ]:

def Eingabe_RIS_CV_Simulator():
    Fenster = Toplevel()  
    Fenster.geometry("850x630")
    Fenster.title("Set RIS-CV-Parameters")

    global Superposition_Mode
    Superposition_Mode = 0    
    
    colorbgr = Label(Fenster, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000)
    
    button=Button(Fenster,text="Klick to open Definition file with order\nt / s     tab.     (E vs. E°) / V    *",bg = '#80ffbf', command=RIS_Definer).place(x = 25, y = 10, width = 250, height = 45)
        
    SCStep_Label = Label(Fenster,text="Direction of forward scan or first Pot. Step*",bg = '#80ffbf')
    SCStep_Label.place(x = 25, y = 60)
    AnSweep = IntVar()
    Checkbutton(Fenster, text="Positive", bg = '#80ffbf', variable=AnSweep).place(x = 25, y = 85)
    CathSweep = IntVar()
    Checkbutton(Fenster, text="Negative", bg = '#80ffbf', variable=CathSweep).place(x = 100, y = 85)
    
    n_Label = Label(Fenster, text="n*", bg = '#80ffbf')
    n_Label.place(x = 300, y = 10)
    n_Eingabe = Entry(Fenster)                                               
    n_Eingabe.place(x = 425, y = 10)
    
    alpha_Label = Label(Fenster, text= u"\u03B1*",bg = '#80ffbf')
    alpha_Label.place(x = 300, y = 35)
    alpha_Eingabe = Entry(Fenster)                                               
    alpha_Eingabe.place(x = 425, y = 35)
    
    kzero_Label = Label(Fenster,text="k° [cm/s]*",bg = '#80ffbf')
    kzero_Label.place(x = 300, y = 60)
    kzero_Eingabe = Entry(Fenster)
    kzero_Eingabe.place(x = 425, y = 60)
    
    T_Label = Label(Fenster,text="T [°C]*",bg = '#80ffbf')
    T_Label.place(x = 300, y = 85)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.place(x = 425, y = 85)
    
    A_Label = Label(Fenster,text="A in cm^2*",bg = '#80ffbf')
    A_Label.place(x = 575, y = 10)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.place(x = 700, y = 10)
        
    c_Label = Label(Fenster,text="c in mol/L*",bg = '#80ffbf')
    c_Label.place(x = 575, y = 35)
    c_Eingabe = Entry(Fenster)
    c_Eingabe.place(x = 700, y = 35)    
    
    D_Label = Label(Fenster,text="D 10^-6[cm^2/s]*",bg = '#80ffbf')
    D_Label.place(x = 575, y = 60)
    D_Eingabe = Entry(Fenster)
    D_Eingabe.place(x = 700, y = 60)
    
    if FIT ==1:
        Ru_Label = Label(Fenster,text="Ru in Ohm*",bg = '#80ffbf')
        Ru_Label.place(x = 575, y = 85)
        Ru_Eingabe = Entry(Fenster)
        Ru_Eingabe.place(x = 700, y = 85)
    
    
    SepLine1 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine1.place(x = 0, y = 110, height = 10)
    
    
    colorbox1 = Label(Fenster,text="",bg = '#EFEFEF')
    colorbox1.place(x = 20, y = 130, width = 260, height = 85)
    var1 = IntVar()
    Checkbutton(Fenster, text="Unequal D", bg = '#EFEFEF', variable=var1).place(x = 25, y = 135)
    Df_Label = Label(Fenster,text="D_f 10^-6[cm^2/s]",bg = '#EFEFEF')
    Df_Label.place(x = 25, y = 160)
    Df_Eingabe = Entry(Fenster)
    Df_Eingabe.place(x = 150, y = 160)
    Db_Label = Label(Fenster,text="D_b 10^-6[cm^2/s]",bg = '#EFEFEF')
    Db_Label.place(x = 25, y = 185)
    Db_Eingabe = Entry(Fenster)
    Db_Eingabe.place(x = 150, y = 185)
    
    
    colorbox2 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox2.place(x = 295, y = 130, width = 260, height = 85)
    var2 = IntVar()
    Checkbutton(Fenster, text="Preceding ch. eq.",bg = '#EFEFEF', variable=var2).place(x = 300, y = 135)
    kp_Label = Label(Fenster,text="k_p [1/s]", bg = '#EFEFEF')
    kp_Label.place(x = 300, y = 160)
    kp_Eingabe = Entry(Fenster)
    kp_Eingabe.place(x = 425, y = 160)
    kmp_Label = Label(Fenster,text="k_-p [1/s]", bg = '#EFEFEF')
    kmp_Label.place(x = 300, y = 185)
    kmp_Eingabe = Entry(Fenster)
    kmp_Eingabe.place(x = 425, y = 185)
    
    
    colorbox3 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox3.place(x = 570, y = 130, width = 260, height = 85)
    var3 = IntVar()
    Checkbutton(Fenster, text="Following ch. eq.", bg = '#EFEFEF', variable=var3).place(x = 575, y = 135)
    kf_Label = Label(Fenster,text="k_f [1/s]", bg = '#EFEFEF')
    kf_Label.place(x = 575, y = 160)
    kf_Eingabe = Entry(Fenster)
    kf_Eingabe.place(x = 700, y = 160)
    kmf_Label = Label(Fenster,text="k_-f [1/s]", bg = '#EFEFEF')
    kmf_Label.place(x = 575, y = 185)
    kmf_Eingabe = Entry(Fenster)
    kmf_Eingabe.place(x = 700, y = 185)
    
    
    
    colorbox4 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox4.place(x = 20, y = 230, width = 260, height = 60)
    var15 = IntVar()
    Checkbutton(Fenster, text="Finite heterogeneous kinetics ",bg = '#EFEFEF', variable=var15).place(x = 25, y = 235)
    kmax_Label = Label(Fenster,text="fin. het. k [cm/s]",bg = '#EFEFEF')
    kmax_Label.place(x = 25, y = 260)
    kmax_Eingabe = Entry(Fenster)
    kmax_Eingabe.place(x = 150, y = 260)
    
    
    if Statistical  == 0:
        if model_b == 1 or model_b1 == 1:
            colorbox5 = Label(Fenster,text="", bg = '#EFEFEF')
            colorbox5.place(x = 295, y = 230, width = 260, height = 60)
            d_Label = Label(Fenster,text="distance 10^-6[m]*", bg = '#EFEFEF')
            d_Label.place(x = 300, y = 260)
            d_Eingabe = Entry(Fenster)
            d_Eingabe.place(x = 425, y = 260)
    
        if model_c == 1 or model_e == 1 or model_f == 1 or model_g == 1:
            colorbox5 = Label(Fenster,text="", bg = '#EFEFEF')
            colorbox5.place(x = 295, y = 230, width = 260, height = 60)
            r_Label = Label(Fenster,text="radius 10^-6[m]*", bg = '#EFEFEF')
            r_Label.place(x = 300, y = 260)
            r_Eingabe = Entry(Fenster)
            r_Eingabe.place(x = 425, y = 260)
        
        if model_d  == 1:
            colorbox5 = Label(Fenster,text="", bg = '#EFEFEF')
            colorbox5.place(x = 295, y = 230, width = 260, height = 60)
            d_Label = Label(Fenster,text="distance 10^-6[m]*", bg = '#EFEFEF')
            d_Label.place(x = 300, y = 235)
            d_Eingabe = Entry(Fenster)
            d_Eingabe.place(x = 425, y = 235)
            r_Label = Label(Fenster,text="radius 10^-6[m]*", bg = '#EFEFEF')
            r_Label.place(x = 300, y = 260)
            r_Eingabe = Entry(Fenster)
            r_Eingabe.place(x = 425, y = 260)
      
        
        
    if Statistical  == 1:
        colorbox5 = Label(Fenster,text="", bg = '#EFEFEF')
        colorbox5.place(x = 295, y = 230, width = 260, height = 60)
        DistFunc_Nodes_Label = Label(Fenster,text="Intervals of Dist.*", bg = '#EFEFEF')
        DistFunc_Nodes_Label.place(x = 300, y = 235)
        DistFunc_Nodes_Eingabe = Entry(Fenster)
        DistFunc_Nodes_Eingabe.place(x = 425, y = 235)
        
        if model_b == 1:
            Number_N_Label = Label(Fenster,text="Sheets per mm*", bg = '#EFEFEF')
        if model_d == 1 or model_e == 1:
            Number_N_Label = Label(Fenster,text="Cylinders per mm^2*", bg = '#EFEFEF')
        if model_g == 1:
            Number_N_Label = Label(Fenster,text="Spheres per mm^3*", bg = '#EFEFEF')
        Number_N_Label.place(x = 300, y = 260)
        Number_N_Eingabe = Entry(Fenster)
        Number_N_Eingabe.place(x = 425, y = 260)
        
        
        
        if model_d == 1:
            colorbox6 = Label(Fenster,text="", bg = '#EFEFEF')
            colorbox6.place(x = 570, y = 230, width = 260, height = 60)
            r_Label = Label(Fenster,text="Fiber radius 10^-6[m]*")
            r_Label.place(x = 575, y = 260)
            r_Eingabe = Entry(Fenster)
            r_Eingabe.place(x = 700, y = 260)
        
        if model_b == 1:
            colorbox6 = Label(Fenster,text="", bg = '#EFEFEF')
            colorbox6.place(x = 570, y = 230, width = 260, height = 60)
            sheet_thickness_Label = Label(Fenster,text="Sheet thick 10^-6[m]*")
            sheet_thickness_Label.place(x = 575, y = 260)
            sheet_thickness_Eingabe = Entry(Fenster)
            sheet_thickness_Eingabe.place(x = 700, y = 260)
    
    
    SepLine2 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine2.place(x = 0, y = 295, height = 10)
    
    
    colorbox7 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox7.place(x = 20, y = 310, width = 260, height = 60)
    ModDeltXi_Label = Label(Fenster,text="This field is inactive in RIS-CV mode.", bg = '#EFEFEF')
    ModDeltXi_Label.place(x = 25, y = 330)

    
    
    if FIT ==1:
        colorbox8 = Label(Fenster,text="", bg = '#EFEFEF')
        colorbox8.place(x = 295, y = 310, width = 535, height = 60)
        var7 = IntVar()
        Checkbutton(Fenster, text="Lin. Basecorr", bg = '#EFEFEF', variable=var7).place(x = 300, y = 315)
        LinStart_Label = Label(Fenster,text="Start vs. Ref. [V]", bg = '#EFEFEF')
        LinStart_Label.place(x = 300, y = 340)
        LinStart_Eingabe = Entry(Fenster)
        LinStart_Eingabe.place(x = 425, y = 340)
        LinEnd_Label = Label(Fenster,text="End vs. Ref. [V]", bg = '#EFEFEF')
        LinEnd_Label.place(x = 575, y = 340)
        LinEnd_Eingabe = Entry(Fenster)
        LinEnd_Eingabe.place(x = 700, y = 340)
    
    
    SepLine3 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine3.place(x = 0, y = 375, height = 10)
    
    
    if Statistical  == 1:
        colorbox9 = Label(Fenster,text="", bg = '#EFEFEF')
        colorbox9.place(x = 20, y = 390, width = 810, height = 60) 
        var16 = IntVar()
        Checkbutton(Fenster, text="Modify Distribution function", bg = '#EFEFEF', variable=var16).place(x = 25, y = 395)
        DFP_a_Label = Label(Fenster,text="Parameter a", bg = '#EFEFEF')
        DFP_a_Label.place(x = 25, y = 420)
        DFP_a_Eingabe = Entry(Fenster)
        DFP_a_Eingabe.place(x = 150, y = 420)
        DFP_b_Label = Label(Fenster,text="Parameter b", bg = '#EFEFEF')
        DFP_b_Label.place(x = 300, y = 420)
        DFP_b_Eingabe = Entry(Fenster)
        DFP_b_Eingabe.place(x = 425, y = 420)
        DFP_c_Label = Label(Fenster,text="Parameter c", bg = '#EFEFEF')
        DFP_c_Label.place(x = 575, y = 420)
        DFP_c_Eingabe = Entry(Fenster)
        DFP_c_Eingabe.place(x = 700, y = 420)
        
    if Statistical  == 0:
        NoEntryLabel = Label(Fenster,text="This field is inactive for your desired calculation.", bg = '#80ffbf')
        NoEntryLabel.place(x = 300, y = 410) 
    
    
    SepLine4 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine4.place(x = 0, y = 460, height = 10)
 
    
    
    colorbox10 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox10.place(x = 20, y = 470, width = 260, height = 120) 
    var6 = IntVar()
    Checkbutton(Fenster, text="Save as txt", bg = '#EFEFEF', variable=var6).place(x = 25, y = 475)
    if FIT == 0:
        if Statistical == 0:
            NoEntryLabel2 = Label(Fenster,text="No specified saving options provided.", bg = '#EFEFEF')
            NoEntryLabel2.place(x = 25, y = 500) 
        if Statistical  == 1:
            var12 = IntVar()
            Checkbutton(Fenster, text="Save Individuals", bg = '#EFEFEF',variable=var12).place(x = 25, y = 500)
            var14 = IntVar()
            Checkbutton(Fenster, text="Save Distribution function", bg = '#EFEFEF', variable=var14).place(x = 25, y = 525)
    if FIT == 1:
        if Statistical  == 0 :
            NoEntryLabel2 = Label(Fenster,text="No specified saving options provided.", bg = '#EFEFEF')
            NoEntryLabel2.place(x = 25, y = 500) 
        if Statistical  == 1:
            var12 = IntVar()
            Checkbutton(Fenster, text="Save Individuals", bg = '#EFEFEF',variable=var12).place(x = 25, y = 500)
            var14 = IntVar()
            Checkbutton(Fenster, text="Save Distribution function", bg = '#EFEFEF', variable=var14).place(x = 25, y = 525)
        
    
    colorbox11 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox11.place(x = 295, y = 470, width = 260, height = 120) 
    
    FFT_displayer_box = IntVar()
    Checkbutton(Fenster, text="Show abs(FFT(I))", bg = '#EFEFEF', variable=FFT_displayer_box).place(x = 300, y = 475)
    t_vs_I_displayer_box  = IntVar()
    Checkbutton(Fenster, text="Show t vs. I", bg = '#EFEFEF', variable=t_vs_I_displayer_box).place(x = 300, y = 500)
    
    if FIT == 0:
        if Statistical  == 1:
            var13 = IntVar()
            Checkbutton(Fenster, text="Show Distribution function", bg = '#EFEFEF', variable=var13).place(x = 300, y = 525)
            var11 = IntVar()
            Checkbutton(Fenster, text="Show Individuals", bg = '#EFEFEF', variable=var11).place(x = 300, y = 550)
    if FIT == 1:
        if Statistical  == 1:
            var13 = IntVar()
            Checkbutton(Fenster, text="Show Distribution function", bg = '#EFEFEF', variable=var13).place(x = 300, y = 475)
        if Statistical  == 0:
            NoEntryLabel3 = Label(Fenster,text="No specified displaying options provided.", bg = '#EFEFEF')
            NoEntryLabel3.place(x = 300, y = 475) 
        
        
    
    
    #-------------------------------------------------------
    #Akzeptierfunktion schreiben
    #-------------------------------------------------------
    
    def Accept_RISCV_Params():
        global n
        global alpha
        global kzero
        global kp
        global kmp
        global kf
        global kmf
        global T
        global distance
        global r
        global Precedeq
        global Succeedeq
        global Unequal_D
        global Cathsweeper
        global AstxtSaver
        global A
        global Area
        global c
        global Ru
        global LinBaseCorr
        global D_f
        global D_b
        global Show_Individuals
        global Save_Individuals
        global Show_Dist_Func
        global Save_Dist_Func
        global DistFunc_Nodes
        global DistFunc_Modifier
        global average_d_or_r
        global Number_N
        global DFP_a
        global DFP_b
        global DFP_c
        global sheet_thickness
        global FinKin
        global kfin
        global FFT_displayer
        global t_vs_I_displayer
        
#========================================================================
    #Read entries
#========================================================================
        
        n                = (float(n_Eingabe.get()))
        alpha            = (float(alpha_Eingabe.get()))
        kzero            = (float(kzero_Eingabe.get()))
        T                = (float(T_Eingabe.get())) + 273.15
        FinKin           = var15.get()
        kfin             = 0
        c                = (float(c_Eingabe.get()))*0.001
        A                = (float(A_Eingabe.get()))
        AstxtSaver       = var6.get()
        Cathsweeper      = CathSweep.get()
        FFT_displayer    = FFT_displayer_box.get()
        t_vs_I_displayer = t_vs_I_displayer_box.get()
            
        
        if FinKin == 1:
            kfin = kzero/(float(kmax_Eingabe.get()))
        
        
        
        if Statistical  == 0:
            if model_b == 1 or model_b1 == 1:
                distance     = (float(d_Eingabe.get()))*0.0001
            if model_c == 1  or model_f == 1 or model_e == 1  or model_g == 1:
                r            = (float(r_Eingabe.get()))*0.0001
            if model_d  == 1:
                distance     = (float(d_Eingabe.get()))*0.0001
                r            = (float(r_Eingabe.get()))*0.0001
            
        Unequal_D   = var1.get()
        if Unequal_D ==1:
            D_f     = (float(Df_Eingabe.get()))*0.000001
            D_b     = (float(Db_Eingabe.get()))*0.000001
        if Unequal_D ==0:
            global D
            D_f       = (float(D_Eingabe.get()))*0.000001
            D_b       = (float(D_Eingabe.get()))*0.000001
        
        Precedeq = var2.get()
        if Precedeq ==1:
            kp    = (float(kp_Eingabe.get()))
            kmp   = (float(kmp_Eingabe.get()))
           
        
        Succeedeq = var3.get()
        if Succeedeq ==1:
            kf   = (float(kf_Eingabe.get()))
            kmf  = (float(kmf_Eingabe.get()))
        
            
        
        if FIT == 1:  
            Ru     = (float(Ru_Eingabe.get()))
            Ezero  = (float(Ezero_Eingabe.get()))
            LinBaseCorr = var7.get()
            if LinBaseCorr == 1:
                global LineStart
                global LineEnd
                LineStart = (float(LinStart_Eingabe.get()))
                LineEnd   = (float(LinEnd_Eingabe.get()))
        
       
         
        
#========================================================================
#Define distribution function
#========================================================================
        if Statistical  == 1:
            
            Number_N = (float(Number_N_Eingabe.get()))
            #there are Number_N sheets per mm, cylinders per mm^2 or spheres per mm^3  
            #--> the following calculates r_average or d_average in µm
            if model_b ==1:
                average_d_or_r = (1000.0/(Number_N))    #distance between two sheet centers
            if model_d == 1 or model_e == 1:
                average_d_or_r = 2*(1000000.0/(np.pi*Number_N))**0.5
            if model_g == 1:
                average_d_or_r = 2 * (1000000000.0*3.0/(4.0*np.pi*Number_N))**0.333333333
                
            
            
            DistFunc_Modifier = var16.get()
            
            if DistFunc_Modifier == 0:
                if model_b ==1:
                    #1D dist-Function
                    DFP_a =  2.0
                    DFP_b =  2.0
                    DFP_c =  1.0
    
                if model_d ==1 or model_e == 1:
                    #2D dist-Function
                    DFP_a =  3.3095
                    DFP_b =  3.0328
                    DFP_c =  1.0787

                if model_g ==1:
                    #3D dist-Function
                    DFP_a =  4.8065
                    DFP_b =  4.06342
                    DFP_c =  1.16391
            
            if DistFunc_Modifier == 1:
                DFP_a =  (float(DFP_a_Eingabe.get()))
                DFP_b =  (float(DFP_b_Eingabe.get()))
                DFP_c =  (float(DFP_c_Eingabe.get()))
                
            if FIT == 0:
                Show_Individuals  = var11.get()
            Save_Individuals  = var12.get()
            Show_Dist_Func    = var13.get()
            Save_Dist_Func    = var14.get()
            DistFunc_Nodes    = 1+(int(DistFunc_Nodes_Eingabe.get()))
            
            
            if model_d  == 1:
                r               = (float(r_Eingabe.get()))*0.0001                 #Cylinder radius in cm (given entry in µm)
            if model_b == 1:
                sheet_thickness = (float(sheet_thickness_Eingabe.get()))*0.0001   #Sheet thickness in cm (given entry in µm)
                
#========================================================================
            

            
    button=Button(Fenster,text="Accept Parameters",bg = '#EFEFEF', command=Accept_RISCV_Params).place(x = 570, y = 470, width = 260, height = 55)
    button=Button(Fenster,text="Next",bg = '#EFEFEF', command=RISCV_Type_Chooser).place(x = 570, y =535, width = 260, height = 55)
    


# In[ ]:

def RISCV_Type_Chooser():
    if model_a == 1:
        RISCV_Start_End_Definer()
        if Statistical  == 0:
            RISCV_model_a_ConvolutionInverter()
        if Statistical  == 1:
            RISCV_DistFunc_Calcer()
    if model_b == 1:
        RISCV_Start_End_Definer()
        if Statistical  == 0:
            RISCV_model_b_ConvolutionInverter()
        if Statistical  == 1:
            RISCV_DistFunc_Calcer()
    if model_b1 == 1:
        RISCV_Start_End_Definer()
        if Statistical  == 0:
            RISCV_model_b1_ConvolutionInverter()
        if Statistical  == 1:
            RISCV_DistFunc_Calcer() 
    if model_c == 1:
        RISCV_Start_End_Definer()
        if Statistical  == 0:
            RISCV_model_c_ConvolutionInverter()
        if Statistical  == 1:
            RISCV_DistFunc_Calcer()
    if model_d == 1:
        RISCV_Start_End_Definer()
        if Statistical  == 0:
            RISCV_model_d_ConvolutionInverter()
        if Statistical  == 1:
            RISCV_DistFunc_Calcer()
    if model_e == 1:
        RISCV_Start_End_Definer()
        if Statistical  == 0:
            RISCV_model_e_ConvolutionInverter()
        if Statistical  == 1:
            RISCV_DistFunc_Calcer()
    if model_f == 1:
        RISCV_Start_End_Definer()
        if Statistical  == 0:
            RISCV_model_f_ConvolutionInverter()
        if Statistical  == 1:
            RISCV_DistFunc_Calcer()
    if model_g == 1:
        RISCV_Start_End_Definer()
        if Statistical  == 0:
            RISCV_model_g_ConvolutionInverter()
        if Statistical  == 1:
            RISCV_DistFunc_Calcer()
   


# In[ ]:

def RISCV_Start_End_Definer():    
    if FIT == 1: 
        global EEE
        global SSS
        global AsymmetrischesEnde
        StromarraysROH             = np.squeeze(np.transpose(Strom[0::,0:1:]))   
        UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,0:1:]))
        PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,0:1:]))
        global Stromarrays
        global Potenzialarrays
        Stromarrays = StromarraysROH
        UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::] 
        Potenzialarrays            = (PotenzialarraysROH[::]   -  Ru*Stromarrays/1000000) - Ezero   #referring all to E°
        if DesiReactxx ==0:
            Stromarrays        = -Stromarrays
            Potenzialarrays    = -Potenzialarrays
            EEE                = np.max(Potenzialarrays)    #Reduktion
            SSS                = Potenzialarrays[0]
            AsymmetrischesEnde = Potenzialarrays[-1]  
            Stromarrays        = -Stromarrays
            Potenzialarrays    = -Potenzialarrays
        
        if DesiReactxx ==1:
            EEE                = np.max(Potenzialarrays)    #Oxidation
            SSS                = Potenzialarrays[0]
            AsymmetrischesEnde = Potenzialarrays[-1]
        

        if LinBaseCorr == 1:
            global StartPotLinCorr
            global EndPotLinCorr
            global BaseSlope
            if DesiReactxx ==0:
                Potenzialarrays = Potenzialarrays[::-1]
                Stromarrays     = Stromarrays[::-1]
                HinScanEndindex = np.argmin(Potenzialarrays)
                LineStartIndex  = (np.abs(Potenzialarrays[-1:HinScanEndindex:-1]+Ezero-LineStart)).argmin()
                StartPotLinCorr = Potenzialarrays[LineStartIndex] #Startpot der Lin Korr vs. E_Zero
                LineEndIndex    = (np.abs(Potenzialarrays[-1:HinScanEndindex:-1]+Ezero-LineEnd)).argmin()
                EndPotLinCorr   = Potenzialarrays[LineEndIndex]   #Endpot der Lin Korr vs. E_Zero
                BaseSlope       = (Stromarrays[LineEndIndex]-Stromarrays[LineStartIndex])/(Potenzialarrays[LineEndIndex]-Potenzialarrays[LineStartIndex])
                Stromarrays     = (Stromarrays - Stromarrays[LineStartIndex])-BaseSlope*(Potenzialarrays-Potenzialarrays[LineStartIndex])
            if DesiReactxx ==1:
                HinScanEndindex = np.argmax(Potenzialarrays)
                LineStartIndex  = np.abs(Potenzialarrays[0:HinScanEndindex]+Ezero-LineStart).argmin()
                StartPotLinCorr = Potenzialarrays[LineStartIndex] #Startpot der Lin Korr vs. E_Zero
                LineEndIndex    = np.abs(Potenzialarrays[0:HinScanEndindex]+Ezero-LineEnd).argmin()
                EndPotLinCorr   = Potenzialarrays[LineEndIndex]   #Endpot der Lin Korr vs. E_Zero
                BaseSlope       = (Stromarrays[LineEndIndex]-Stromarrays[LineStartIndex])/(Potenzialarrays[LineEndIndex]-Potenzialarrays[LineStartIndex])
                Stromarrays     = (Stromarrays - Stromarrays[LineStartIndex])-BaseSlope*(Potenzialarrays-Potenzialarrays[LineStartIndex])
            
            
            
    


# In[ ]:

def RISCV_DistFunc_Calcer(): 

    global R_Array
    global P_Array
    
    y_Array = np.arange(0,100,0.01)
    P_Array = (DFP_c*DFP_b**(DFP_a/DFP_c) / gamma(DFP_a/DFP_c)) *(y_Array**(DFP_a-1) *np.exp(-DFP_b*y_Array**DFP_c))/average_d_or_r
    R_Array = average_d_or_r*y_Array
    
    IndexMax         =  np.asscalar(np.where(P_Array == np.max(P_Array))[0])
    Kurvenende       = (np.abs(P_Array[IndexMax::]-0.0005*P_Array[IndexMax])).argmin()   #Änderung auf 0.001
    Array_Ab_Maximum = P_Array[IndexMax::]
    Endwert          = Array_Ab_Maximum[Kurvenende]
    EndIndex         = np.asscalar(np.where(P_Array == Endwert)[0])

    global R_Array_Unchanged 
    global p_Array_Unchanged
    R_Array_Unchanged = R_Array       #new
    P_Array_Unchanged = P_Array       #new
    
    R_Array = R_Array[0:EndIndex] 
    P_Array = P_Array[0:EndIndex]
    
    global Probab_Array
    global R_Mean_Array  
    global RR_Mean_Array

    Probab_Array  = np.empty(DistFunc_Nodes) 
    R_Mean_Array  = np.empty(DistFunc_Nodes)
    RR_Mean_Array = np.empty(DistFunc_Nodes)
    
    DistFuncInterpolation = InterpolatedUnivariateSpline(R_Array_Unchanged,P_Array_Unchanged)
    R_Div = (R_Array[-1] - R_Array[0])/(1.0*float(DistFunc_Nodes))

    for i in range(len(Probab_Array)-1):
        R_Mean_Array[i] = R_Div*i  #Interval borders
        if (i+1) < (len(Probab_Array)-1):   #new
            Probab_Array[i] = quad(DistFuncInterpolation, R_Div*i,R_Div*(i+1) )[0]   #Wahscheinl. v. großem Segment
        if (i+1) == (len(Probab_Array)-1):
            Probab_Array[i] = quad(DistFuncInterpolation, R_Div*i, R_Array_Unchanged[-1] )[0]
            
                                                  
        if (i+1) < (len(Probab_Array)-1):
            R_Array_Resolved  = np.linspace(R_Div*i,R_Div*(i+1),50)
        if (i+1) == (len(Probab_Array)-1):
            R_Array_Resolved  = np.linspace(R_Div*i,R_Array_Unchanged[-1],5000)

        
        PR_Indiv_Resolved = np.zeros(len(R_Array_Resolved))         #Indiv. Probability of further resolved Integral-Part
        for j in range(len(R_Array_Resolved)-1):
            #berechnen von Mittlerem R von jedem einzelnen großen Segment.
            #Jedes große Segment wird in 50 Teile geiteilt und darüber dann gewichtet integriert um mittleres R
            #von großem Segment zu finden. Das wird dann mit der Gewichtung des großen Segments am Gesamtanteil verrechnet.
            PR_Indiv_Resolved[j] = 0.5*(R_Array_Resolved[j]+R_Array_Resolved[j+1])*quad(DistFuncInterpolation, R_Array_Resolved[j],R_Array_Resolved[j+1] )[0]
        
        RR_Mean_Array[i] = (1/Probab_Array[i])*np.sum(PR_Indiv_Resolved)    #gewichteter Flächenanteil 
                                                                            #Wahrscheinlichkeit, dass mittleres R aus großem
                                                                            #Segment auftaucht.
    Probab_Array = Probab_Array[0:-1]
    #Das sind eigentlich die Abstände von zwei Fasermitten, ALSO EIGENTLCIH KEINE R!!!
    R_Array      = R_Array
    R_Mean_Array = R_Mean_Array[0:-1]
   

    RISCV_DistFunc_User()
    
    if Show_Dist_Func == 1:
        RISCV_Dist_Func_Shower()

    


# In[ ]:

def RISCV_DistFunc_User():
    global distance 
    global a
    
    
    if FIT == 1:
        RISCV_Start_End_Definer()
        
    
    distance = 0.01    #to run calculation with default values before real calculation starts
    a        = 0.01    #to run calculation with default values before real calculation starts
    
    #HIER MÜSSEN ERST DIE CONVOLUTION INVERTERS GEZÜNDET WERDEN um len(yyy) zu bekommen
    if model_b == 1:
        RISCV_model_b_ConvolutionInverter()
    if model_d == 1:
        RISCV_model_d_ConvolutionInverter()
    if model_e == 1:
        RISCV_model_e_ConvolutionInverter()
    if model_g == 1:
        RISCV_model_g_ConvolutionInverter()
    
    
    global ChiSuperarray
    ChiSuperarray = np.zeros([len(yyy),len(R_Mean_Array)])
    
    
    #internal spherical finite diffusion --> calculation of area-weightning
    #----------------------------------------------------------------------------------------------------
    if model_g == 1:
        SUMMANDS    = np.zeros(len(Probab_Array))
        for i in range(len(Probab_Array)):
            SUMMANDS[i] = Probab_Array[i]*(RR_Mean_Array[i]/2.0)**2
        N_spheres = A/(4*np.pi*np.sum(SUMMANDS)) 
        
        INDIV_AREAS = np.zeros(len(Probab_Array))
        for i in range(len(Probab_Array)):
            INDIV_AREAS[i] = (N_spheres*Probab_Array[i])*(4*np.pi*(RR_Mean_Array[i]/2.0)**2)
    #----------------------------------------------------------------------------------------------------  
    
    #internal cylindrical finite diffusion --> calculation of area-weightning
    #----------------------------------------------------------------------------------------------------
    if model_e == 1:
        lunit = 0.001
        SUMMANDS    = np.zeros(len(Probab_Array))
        for i in range(len(Probab_Array)):
            SUMMANDS[i] = Probab_Array[i]*(RR_Mean_Array/2.0)[i]
        N_cyls = A/(2*np.pi*lunit*np.sum(SUMMANDS)) 
        
        INDIV_AREAS = np.zeros(len(Probab_Array))
        for i in range(len(Probab_Array)):
            INDIV_AREAS[i] = (N_cyls*Probab_Array[i])*(2*np.pi*lunit*(RR_Mean_Array/2.0)[i])
    #----------------------------------------------------------------------------------------------------  
     
        
    for i in range(len(Probab_Array)):   #is one shorter than RR_Mean_Array. Takes only the "useful" RR_Mean (not the last)
        
        #finite planar
        if model_b == 1:  
            distance = 0.5*RR_Mean_Array[i]*0.0001 - 0.5*sheet_thickness    #sheet distance/2 
        
        #internal cylindrical finite
        if model_e == 1:
            a = 0.5*RR_Mean_Array[i]*0.0001                                 #pore radii in µm
                
        #internal spherical finite
        if model_g == 1:
            a = 0.5*RR_Mean_Array[i]*0.0001                                 #pore radii in µm
        
        #external cylindrical finite
        if model_d == 1:
            #not necssary to define "a" here, because it will be re-defined in corresponding modell 
            distance = 0.5*RR_Mean_Array[i]*0.0001  - r                 #center to wall distances in µm

            
        INTERRUPT = 0
        
        if model_b == 1:
            if distance < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            RISCV_model_b_ConvolutionInverter()
        if model_d == 1:
            if distance < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            RISCV_model_d_ConvolutionInverter()
        if model_e == 1:
            if a < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            RISCV_model_e_ConvolutionInverter()
        if model_g == 1:
            if a < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            RISCV_model_g_ConvolutionInverter()
         
        for j in range(len(yyy)):
            if model_b == 1 or model_d == 1:
                ChiSuperarray[j,i] = yyy[j] *np.asscalar(Probab_Array[i])   #weightning by d-probability (planar fin and cylextfin)
            if model_e == 1 or model_g == 1:
                ChiSuperarray[j,i] = yyy[j] *INDIV_AREAS[i]/A       #weightning by area for internal cyl and internal sph

    
    global Chi_ges_Array
    Chi_ges_Array = np.empty(len(yyy))
    for i in range(len(yyy)):
        for j in range(len(Probab_Array)):
            Chi_ges_Array[i] = np.squeeze(np.sum(ChiSuperarray[i,::]))
            
    #====================================================================================================
    #Plotten von simuliertem
    #====================================================================================================
    if INTERRUPT == 0:
        global PotCalc     
        global CurrCalc 
        
        root = Toplevel()
        root.title("Simulated CV")
        
        Xi_i = -np.abs(xxx[0]*n*F/(R*T))

        if Xi_i > -6.9:
            Incorrect_Range_Warner()

        if FIT == 0:
            PotCalc  = xxx 
            CurrCalc = 1000*Chi_ges_Array
            Abbildung = Figure(figsize=(8, 4), dpi=100)
            b1 = Abbildung.add_subplot(121)
            b2 = Abbildung.add_subplot(122)
            if t_vs_I_displayer == 0:
                if Show_Individuals == 0:  
                    b1.plot(PotCalc[:-1:],CurrCalc[:-1:],color='r')
                if Show_Individuals == 1:
                    b1.plot(PotCalc[:-1:],1000*ChiSuperarray[:-1:,::])
                b1.set_xlabel('(E vs. E' '$^0$' ') / V', fontsize=12)
                b1.set_ylabel('I / mA', fontsize=12)
                b1.axvline(0,color='k')
                b1.axhline(0,color='k')
            if t_vs_I_displayer == 1:
                if Show_Individuals == 0:  
                    b1.plot(Input_Time[:-1:],CurrCalc[:-1:],color='r')
                if Show_Individuals == 1:
                    b1.plot(Input_Time[:-1:],1000*ChiSuperarray[:-1:,::])
                b1.set_xlabel('t / s', fontsize=12)
                b1.set_ylabel('I / mA', fontsize=12)   
            for axis in ['top','bottom','left','right']:
                b1.spines[axis].set_linewidth(2)
                b1.spines[axis].set_color('k')
            b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
            if FFT_displayer == 0:
                b2.plot(Input_Time[:-1:],Input_Potential[:-1:],color='r')
                b2.axvline(0,color='k')
                b2.axhline(0,color='k')
                b2.set_xlabel('t / s', fontsize=12)
                b2.set_ylabel('(E vs. ' '$E^0$' ') / V', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b2.spines[axis].set_linewidth(2)
                    b2.spines[axis].set_color('k')
                b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                if AstxtSaver == 1:
                    Simulated_RISCV_as_txt_saver()
                
            if FFT_displayer == 1:
                b2.plot(np.abs(fft(CurrCalc[:-1:]))/np.max(np.abs(fft(CurrCalc[:-1:]))),color='r')
                b2.set_xlabel('point no.', fontsize=12)
                b2.set_ylabel('abs(fft(I / mA))/max(abs(I / mA))', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b2.spines[axis].set_linewidth(2)
                    b2.spines[axis].set_color('k')
                b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                
                IFFT_controller()

            Abbildung.tight_layout()
            canvas = FigureCanvasTkAgg(Abbildung, master=root)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas, root)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

    


# In[ ]:

def RISCV_Dist_Func_Shower():   
    root = Toplevel()
    root.title("Distance distribution-function")
    
    Abbildung = Figure(figsize=(5, 4), dpi=100)
    b = Abbildung.add_subplot(111)
    for i in range(len(R_Mean_Array)):
        #mal 0.5. weil R_Mean_Array quasi durchmesser sind
        b.axvline(R_Mean_Array[i], color = 'k', linewidth = 0.5, linestyle ='-')   #das sind nicht die Mittelwerte, sondern die Stützstellen
        b.axvline(RR_Mean_Array[i], color = 'k', linewidth = 0.5, linestyle ='--') #Das sind die eigentlichen Mittelwerte
        
    #mal 0.5 weil R_Mean_Array quasi urchmesser sind und mal 2, damit Integral passt
    b.plot(R_Array,P_Array, color = 'b')
    b.set_xlabel('midpoint dist. x / micrometer', fontsize=12)
    b.set_ylabel('P(midpoint dist. x) / micrometer^-1', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b.spines[axis].set_linewidth(2)
        b.spines[axis].set_color('k')
    b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    Abbildung.tight_layout()
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    


#=================================================================================================================================
#=================================================================================================================================
#=================================================================================================================================
#Define Convolution functions in laplace Domain
#=================================================================================================================================
#=================================================================================================================================
#=================================================================================================================================


def PlanarSemiHin(s):
    return 1/s**0.5

def PlanarSemiBack(s):
    return 1/s**0.5

def PlanarFinitHin(s):
    return coth(distance*(s/D_f)**0.5)/s**0.5

def PlanarFinitBack(s):
    return coth(distance*(s/D_b)**0.5)/s**0.5

def PlanarFinitTransHin(s):
    return tanh(distance*(s/D_f)**0.5)/s**0.5

def PlanarFinitTransBack(s):
    return tanh(distance*(s/D_b)**0.5)/s**0.5

def ZylSemiHin(s):
    if np.isnan(kv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_f**0.5)*s**0.5) )) == False:
        return kv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_f**0.5)*s**0.5) )
    if np.isnan(kv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_f**0.5)*s**0.5) )) == True:
        return 1/s**0.5

def ZylSemiBack(s):
    if np.isnan(kv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_b**0.5)*s**0.5) )) == False:
        return kv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_b**0.5)*s**0.5) )
    if np.isnan(kv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_b**0.5)*s**0.5) )) == True:
        return 1/s**0.5

def ZylFinHin(s):
    return (kv(0,(a/D_f**0.5)*s**0.5)*iv(1,(d/D_f**0.5)*s**0.5) + iv(0,(a/D_f**0.5)*s**0.5)*kv(1,(d/D_f**0.5)*s**0.5))/(s**0.5 *(kv(1,(a/D_f**0.5)*s**0.5)*iv(1,(d/D_f**0.5)*s**0.5) - iv(1,(a/D_f**0.5)*s**0.5)*kv(1,(d/D_f**0.5)*s**0.5)))

def ZylFinBack(s):
    return (kv(0,(a/D_b**0.5)*s**0.5)*iv(1,(d/D_b**0.5)*s**0.5) + iv(0,(a/D_b**0.5)*s**0.5)*kv(1,(d/D_b**0.5)*s**0.5))/(s**0.5 *(kv(1,(a/D_b**0.5)*s**0.5)*iv(1,(d/D_b**0.5)*s**0.5) - iv(1,(a/D_b**0.5)*s**0.5)*kv(1,(d/D_b**0.5)*s**0.5)))

def ZylIntFinHin(s):
     if np.isnan(iv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_f**0.5)*s**0.5) )) == False:
        return iv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_f**0.5)*s**0.5) )
     if np.isnan(iv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_f**0.5)*s**0.5) )) == True:
        return 1/s**0.5

def ZylIntFinBack(s):
     if np.isnan(iv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_b**0.5)*s**0.5) )) == False:
        return iv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_b**0.5)*s**0.5) )
     if np.isnan(iv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_b**0.5)*s**0.5) )) == True:
        return 1/s**0.5

def SpherExtSemiHin(s):
    return 1/(s**0.5 + D_f**0.5/a)

def SpherExtSemiBack(s):
    return 1/(s**0.5 + D_b**0.5/a)

def SpherIntFinHin(s):
    return 1/(s**0.5 *coth(a*(s/D_f)**0.5) - D_f**0.5/a)

def SpherIntFinBack(s):
    return 1/(s**0.5 *coth(a*(s/D_b)**0.5) - D_b**0.5/a)



#=================================================================================================================================
#=================================================================================================================================
#=================================================================================================================================
#define all Talbot-Inverters
#=================================================================================================================================
#=================================================================================================================================
#=================================================================================================================================



def TalbotPlanarSemiHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanarSemiHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotPlanarSemiBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanarSemiBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotPlanarFinitHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanarFinitHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotPlanarFinitBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanarFinitBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real  

def TalbotPlanarFinitTransHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanarFinitTransHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotPlanarFinitTransBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanarFinitTransBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real  
 

def TalbotZylSemiHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylSemiHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotZylSemiBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylSemiBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   
 

def TalbotZylFinHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylFinHin(z)/z)*dz;
        
    return ((h/(2j*pi))*ans).real   

def TalbotZylFinBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylFinBack(z)/z)*dz;
        
    return ((h/(2j*pi))*ans).real   

def TalbotZylIntFinHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylIntFinHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotZylIntFinBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylIntFinBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotSpherExtSemiHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(SpherExtSemiHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotSpherExtSemiBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(SpherExtSemiBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   


def TalbotSpherIntFinHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(SpherIntFinHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotSpherIntFinBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print ("ERROR:   Inverse transform can not be calculated for t=0")
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(SpherIntFinBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

    
  
#=================================================================================================================================
#=================================================================================================================================
#=================================================================================================================================
#Invert the functions of Interest
#=================================================================================================================================
#=================================================================================================================================
#=================================================================================================================================



def RISCV_model_a_ConvolutionInverter():
    if D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    
    if D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array                      = np.logspace(-8,5,300)
        ft_Array_PlanarSemiHin       = np.empty(len(t_Array)) 
        ft_Array_PlanarSemiBack      = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_PlanarSemiHin[i]     = TalbotPlanarSemiHin(float(t_Array[i]),24)
            ft_Array_PlanarSemiBack[i]    = TalbotPlanarSemiBack(float(t_Array[i]),24)
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_PlanarSemiHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_PlanarSemiBack, k=3)
        
        RIS_CV_Calculator()
    

    
def RISCV_model_b_ConvolutionInverter():
    #also in statistical case distance will be passed as distance
    if distance < 0.00001 or distance > 10.0 or D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    if distance >= 0.00001 and distance <= 10.0 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array                       = np.logspace(-8,5,300)
        ft_Array_PlanarFinitHin       = np.empty(len(t_Array)) 
        ft_Array_PlanarFinitBack      = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_PlanarFinitHin[i]     = TalbotPlanarFinitHin(float(t_Array[i]),24)
            ft_Array_PlanarFinitBack[i]    = TalbotPlanarFinitBack(float(t_Array[i]),24)   
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_PlanarFinitHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_PlanarFinitBack, k=3)
        RIS_CV_Calculator()
    

def RISCV_model_b1_ConvolutionInverter():
    #also in statistical case distance will be passed as distance
    if distance < 0.00001 or distance > 10.0 or D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    if distance >= 0.00001 and distance <= 10.0 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array                                 = np.logspace(-8,5,300)
        ft_Array_PlanarFinitTransHin            = np.empty(len(t_Array)) 
        ft_Array_PlanarFinitTransBack           = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_PlanarFinitTransHin[i]     = TalbotPlanarFinitTransHin(float(t_Array[i]),24)
            ft_Array_PlanarFinitTransBack[i]    = TalbotPlanarFinitTransBack(float(t_Array[i]),24)   
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_PlanarFinitTransHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_PlanarFinitTransBack, k=3)
        RIS_CV_Calculator()
    
    
def RISCV_model_c_ConvolutionInverter():
    global a
    a = r
    
    if a>10 or a<0.00001 or D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
        
    if a<=10 and a>=0.00001 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array              = np.logspace(-8,5,300)
        ft_Array_ZylSemiHin  = np.empty(len(t_Array)) 
        ft_Array_ZylSemiBack = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_ZylSemiHin[i]     = TalbotZylSemiHin(float(t_Array[i]),24)
            ft_Array_ZylSemiBack[i]    = TalbotZylSemiBack(float(t_Array[i]),24)
        #----------------------------------------------------------------------------------------
        #Interpolationen finden
        #----------------------------------------------------------------------------------------
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_ZylSemiHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_ZylSemiBack, k=3)
        RIS_CV_Calculator()
    

    
def RISCV_model_d_ConvolutionInverter():
    global a
    a = r

    global d 
    d = a + distance
    
    if a>10 or a<0.00001 or distance>10 or distance<0.00001 or D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    if a<=10 and a>=0.00001 and distance<=10 and distance>=0.00001 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array              = np.logspace(-8,5,300)
        ft_Array_ZylSemiHin  = np.empty(len(t_Array)) 
        ft_Array_ZylFinHin   = np.empty(len(t_Array)) 
        ft_Array_PlanFinHin  = np.empty(len(t_Array)) 
        ft_Array_ZylSemiBack = np.empty(len(t_Array)) 
        ft_Array_ZylFinBack  = np.empty(len(t_Array)) 
        ft_Array_PlanFinBack = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_ZylSemiHin[i]     = TalbotZylSemiHin(float(t_Array[i]),24)
            ft_Array_PlanFinHin[i]     = TalbotPlanarFinitHin(float(t_Array[i]),24)
            ft_Array_ZylFinHin[i]      = TalbotZylFinHin(float(t_Array[i]),24)
            ft_Array_ZylSemiBack[i]    = TalbotZylSemiBack(float(t_Array[i]),24)
            ft_Array_PlanFinBack[i]    = TalbotPlanarFinitBack(float(t_Array[i]),24)
            ft_Array_ZylFinBack[i]     = TalbotZylFinBack(float(t_Array[i]),24)
        #--------------------------------------------------------------------------------------------    
        ft_Array_ZylSemi_revHin  = ft_Array_ZylSemiHin[::-1]
        ft_Array_ZylFin_revHin   = ft_Array_ZylFinHin[::-1] 
        ft_Array_PlanFin_revHin  = ft_Array_PlanFinHin[::-1]
        ft_Array_ZylSemi_revBack = ft_Array_ZylSemiBack[::-1]
        ft_Array_ZylFin_revBack  = ft_Array_ZylFinBack[::-1] 
        ft_Array_PlanFin_revBack = ft_Array_PlanFinBack[::-1]
        Schalter1Hin  = 0
        Schalter2Hin  = 0
        Schalter1Back = 0
        Schalter2Back = 0
        GesamtHin  = np.empty(len(t_Array))
        GesamtBack = np.empty(len(t_Array))
        #---------------------------------------------------------------------------------------
        #Hin
        #---------------------------------------------------------------------------------------
        for i in range(len(t_Array)):
            if Schalter1Hin == 0:
                if a <=0.01:
                    if distance >= 0.0001:
                        GesamtHin[i] = ft_Array_ZylFin_revHin[i]
                        if np.abs(ft_Array_PlanFin_revHin[i]- ft_Array_ZylFin_revHin[i]) <= 10**(-5):
                            Schalter1Hin = 1
                    if distance < 0.0001:
                        GesamtHin[i] = ft_Array_ZylFin_revHin[i]
                        if np.abs(ft_Array_PlanFin_revHin[i]- ft_Array_ZylFin_revHin[i]) <= 10**(-4):
                            Schalter1Hin = 1
                if a > 0.01:
                    GesamtHin[i] = ft_Array_ZylFin_revHin[i]
                    if np.abs(ft_Array_PlanFin_revHin[i]- ft_Array_ZylFin_revHin[i]) <= 10**(-3):
                        Schalter1Hin = 1
                if a > 0.1:
                    GesamtHin[i] = ft_Array_ZylFin_revHin[i]
                    if np.abs(ft_Array_PlanFin_revHin[i]- ft_Array_ZylFin_revHin[i]) <= 10**(-2):
                        Schalter1Hin = 1
                if a > 1:
                    GesamtHin[i] = ft_Array_ZylFin_revHin[i]
                    if np.abs(ft_Array_PlanFin_revHin[i]- ft_Array_ZylFin_revHin[i]) <= 10**(-1):
                        Schalter1Hin = 1
            if Schalter1Hin == 1:
                GesamtHin[i] = ft_Array_PlanFin_revHin[i]
        for i in range(len(t_Array)):
            if Schalter2Hin == 0:
                if np.abs(GesamtHin[i] - ft_Array_ZylSemi_revHin[i]) <= 10**(-5):
                    Schalter2Hin = 1
            if Schalter2Hin == 1:
                GesamtHin[i] = ft_Array_ZylSemi_revHin[i]  
        #---------------------------------------------------------------------------------------
        #Back
        #---------------------------------------------------------------------------------------       
            if Schalter1Back == 0:
                if a <=0.01:
                    if distance >= 0.0001:
                        GesamtBack[i] = ft_Array_ZylFin_revBack[i]
                        if np.abs(ft_Array_PlanFin_revBack[i]- ft_Array_ZylFin_revBack[i]) <= 10**(-5):
                            Schalter1Back = 1
                    if distance < 0.0001:
                        GesamtBack[i] = ft_Array_ZylFin_revBack[i]
                        if np.abs(ft_Array_PlanFin_revBack[i]- ft_Array_ZylFin_revBack[i]) <= 10**(-4):
                            Schalter1Back = 1
                if a > 0.01:
                    GesamtBack[i] = ft_Array_ZylFin_revBack[i]
                    if np.abs(ft_Array_PlanFin_revBack[i]- ft_Array_ZylFin_revBack[i]) <= 10**(-3):
                        Schalter1Back = 1
                if a > 0.1:
                    GesamtBack[i] = ft_Array_ZylFin_revBack[i]
                    if np.abs(ft_Array_PlanFin_revBack[i]- ft_Array_ZylFin_revBack[i]) <= 10**(-2):
                        Schalter1Back = 1
                if a > 1:
                    GesamtBack[i] = ft_Array_ZylFin_revBack[i]
                    if np.abs(ft_Array_PlanFin_revBack[i]- ft_Array_ZylFin_revBack[i]) <= 10**(-1):
                        Schalter1Back = 1
            if Schalter1Back == 1:
                GesamtBack[i] = ft_Array_PlanFin_revBack[i]
        for i in range(len(t_Array)):
            if Schalter2Back == 0:
                if np.abs(GesamtBack[i] - ft_Array_ZylSemi_revBack[i]) <= 10**(-5):
                    Schalter2Back = 1
            if Schalter2Back == 1:
                GesamtBack[i] = ft_Array_ZylSemi_revBack[i]     
        #----------------------------------------------------------------------------------------
        #Interpolationen finden
        #----------------------------------------------------------------------------------------
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,GesamtHin[::-1], k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,GesamtBack[::-1], k=3)
        RIS_CV_Calculator()
    
    
    
def RISCV_model_e_ConvolutionInverter():
    global a
    if Statistical == 0:
        a = r
    
    if a>10 or a<0.00001 or D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    if a<=10 and a>=0.00001 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array                = np.logspace(-8,5,300)
        ft_Array_ZylIntFinHin  = np.empty(len(t_Array)) 
        ft_Array_ZylIntFinBack = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_ZylIntFinHin[i]     = TalbotZylIntFinHin(float(t_Array[i]),24)
            ft_Array_ZylIntFinBack[i]    = TalbotZylIntFinBack(float(t_Array[i]),24)
        #----------------------------------------------------------------------------------------
        #Interpolationen finden
        #----------------------------------------------------------------------------------------
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_ZylIntFinHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_ZylIntFinBack, k=3)
        RIS_CV_Calculator()
    


def RISCV_model_f_ConvolutionInverter():
    global a
    a = r
    
    if a>10 and a<0.00001 and D_f <0.000001 and D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    if a<=10 and a>=0.00001 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array                      = np.logspace(-8,5,300)
        ft_Array_SphereExtSemiHin    = np.empty(len(t_Array)) 
        ft_Array_SphereExtSemiBack   = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_SphereExtSemiHin[i]     = TalbotSpherExtSemiHin(float(t_Array[i]),24)
            ft_Array_SphereExtSemiBack[i]    = TalbotSpherExtSemiBack(float(t_Array[i]),24)          
        #----------------------------------------------------------------------------------------
        #Interpolationen finden
        #----------------------------------------------------------------------------------------
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_SphereExtSemiHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_SphereExtSemiBack, k=3)
        RIS_CV_Calculator()
    
    

def RISCV_model_g_ConvolutionInverter():
    global a
    if Statistical == 0:
        a = r
    

    if a>10 or a<0.00001 or D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    if a<=10 and a>=0.00001 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array                     = np.logspace(-8,5,300)
        ft_Array_SphereIntFinHin    = np.empty(len(t_Array)) 
        ft_Array_SphereIntFinBack   = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_SphereIntFinHin[i]     = TalbotSpherIntFinHin(float(t_Array[i]),24)
            ft_Array_SphereIntFinBack[i]    = TalbotSpherIntFinBack(float(t_Array[i]),24)  
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_SphereIntFinHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_SphereIntFinBack, k=3)
        RIS_CV_Calculator()
    
    


# In[ ]:

def RIS_CV_Calculator():
    
    if Statistical == 0:
        root = Toplevel()
        root.title("Simulated CV")
    
    global Kp
    global Kf
    global p
    global f
    p               = 0.0
    Kp              = 1000000000000.0
    Kf              = 0.00000000001
    f               = 0.0
    if Precedeq == 1:
        p               = kp + kmp
        Kp              = kp/kmp
    if Succeedeq == 1:
        Kf              = kf/kmf
        f               = kf + kmf
    global t_Array_cont
    global Xi_i
    if FIT == 0:
        if Cathsweeper == 0:
            Xi_Array        = (n*F/(R*T))*Input_Potential  
        if Cathsweeper == 1:
            Xi_Array        = -(n*F/(R*T))*Input_Potential
        t_Array_cont    = Input_Time
        Exp_Array       = np.exp(-Xi_Array)
        Preced_Array    = np.exp(-p*t_Array_cont)
        Follow_Array    = np.exp(-f*t_Array_cont)
        Kinetik_Array   = D_f**0.5/(kzero*Exp_Array**(-alpha))
        
        
    
    Fin_Kin_p_Array = 1/(1 + kfin*Exp_Array**(-alpha))      #hier mit minus, weil Exp-Array schon umgedreht ist!
    Fin_Kin_s_Array = 1/(1 + kfin*Exp_Array**((1-alpha)))
     
    
    

    #-------------------------------------------------------------------------------------------------- 
    #Jetzt CV berechnen
    #--------------------------------------------------------------------------------------------------
    
    global WuFuInt_p
    WuFuInt_p        = np.empty(len(Xi_Array))
    Delta_WuFuInt_p  = np.empty(len(Xi_Array))
    for i in range(len(t_Array_cont)-1):
        WuFuInt_p[i]        = CV_InterpolationHin(t_Array_cont[i])
        Delta_WuFuInt_p[i]  = (CV_InterpolationHin(t_Array_cont[i+1]) - CV_InterpolationHin(t_Array_cont[i]))

    global WuFuInt_f
    WuFuInt_f        = np.empty(len(Xi_Array))
    Delta_WuFuInt_f  = np.empty(len(Xi_Array))

    for i in range(len(t_Array_cont)-1):
        WuFuInt_f[i]       = CV_InterpolationBack(t_Array_cont[i])
        Delta_WuFuInt_f[i] = (CV_InterpolationBack(t_Array_cont[i+1]) - CV_InterpolationBack(t_Array_cont[i]))
    
    
    I_Array  = np.zeros(len(Xi_Array))

    for i in range(len(Xi_Array)-1):
    
        #MUCH MUCH MUCH FASTER!!!!!!!!!!!!!!!!!!!!! 
        Summandenarray_p_GG   = Delta_WuFuInt_p[i::-1]*I_Array[:i+1:]*Preced_Array[i::-1]
        Summandenarray_f_GG   = Delta_WuFuInt_f[i::-1]*I_Array[:i+1:]*Follow_Array[i::-1]
        Summandenarray_p      = Delta_WuFuInt_p[i::-1]*I_Array[:i+1:]
        Summandenarray_f      = Delta_WuFuInt_f[i::-1]*I_Array[:i+1:]

        I_Array[i] = (Fin_Kin_p_Array[i]*(1/(1+Kp))*(Kp*c*n*F*A*(D_f**0.5) - Kp*np.sum(Summandenarray_p) - np.sum(Summandenarray_p_GG)) -(Fin_Kin_s_Array[i]*(Exp_Array[i]*(1/(1+Kf))*(D_f/D_b)**0.5)*(np.sum(Summandenarray_f) + Kf*np.sum(Summandenarray_f_GG))  ))/(  Kinetik_Array[i] +  Fin_Kin_p_Array[i]*WuFuInt_p[1]+ (Fin_Kin_s_Array[i]*(D_f/D_b)**0.5)*WuFuInt_f[1]*Exp_Array[i])

    #-----------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------
    global xxx
    global yyy
    xxx   = (Input_Potential)
    yyy   = I_Array
    
    if FIT ==0:
        if Cathsweeper == 1:
            yyy = -yyy
            
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON , nur wenn nicht statistisch--> der hat seinen eigenen Plotter
    #-------------------------------------------------------------------------------------------------------------------------
    
    if Statistical == 0:
        global PotCalc
        global CurrCalc
        
        Xi_i = Xi_Array[0]
        
        if Xi_i > -6.9:
            Incorrect_Range_Warner()
            
        if FIT == 0:
            PotCalc  = xxx 
            CurrCalc = 1000*yyy
            Abbildung = Figure(figsize=(8, 4), dpi=100)
            b1 = Abbildung.add_subplot(121)
            b2 = Abbildung.add_subplot(122)
            
            if t_vs_I_displayer == 0:
                b1.plot(PotCalc[:-1:],CurrCalc[:-1:],color='r')
                b1.set_xlabel('(E vs. E' '$^0$' ') / V', fontsize=12)
                b1.set_ylabel('I / mA', fontsize=12)
                b1.axvline(0,color='k')
                b1.axhline(0,color='k')
            if t_vs_I_displayer == 1:
                b1.plot(Input_Time[:-1:],CurrCalc[:-1:],color='r')
                b1.set_xlabel('t / s', fontsize=12)
                b1.set_ylabel('I / mA', fontsize=12)
            for axis in ['top','bottom','left','right']:
                b1.spines[axis].set_linewidth(2)
                b1.spines[axis].set_color('k')
            b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
            if FFT_displayer == 0:
                b2.plot(Input_Time[:-1:],Input_Potential[:-1:],color='r')
                b2.axvline(0,color='k')
                b2.axhline(0,color='k')
                b2.set_xlabel('t / s', fontsize=12)
                b2.set_ylabel('(E vs. ' '$E^0$' ') / V', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b2.spines[axis].set_linewidth(2)
                    b2.spines[axis].set_color('k')
                b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                if AstxtSaver == 1:
                    if Statistical == 0:
                        global FFTAstxtSaver
                        FFTAstxtSaver = 0
                        Simulated_RISCV_as_txt_saver()
                
            if FFT_displayer == 1:
                b2.plot(np.abs(fft(CurrCalc[:-1:]))/np.max(np.abs(fft(CurrCalc[:-1:]))),color='r')
                b2.set_xlabel('point no.', fontsize=12)
                b2.set_ylabel('abs(fft(I / mA))/max(abs(I / mA))', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b2.spines[axis].set_linewidth(2)
                    b2.spines[axis].set_color('k')
                b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                
                IFFT_controller()

            Abbildung.tight_layout()
            canvas = FigureCanvasTkAgg(Abbildung, master=root)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas, root)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)


def IFFT_controller():
    Fenster = Toplevel()                                                         
    Fenster.title("Analyze n-th harmonics control")                         
    Fenster.geometry("1000x300")   

    var1 = IntVar()
    Checkbutton(Fenster, text="Take DC-Base", variable=var1).grid(row=0, column=0, sticky=W)
    BaseFrom_Label = Label(Fenster, text="From low point")
    BaseFrom_Label.grid(row=0, column=1)
    BaseFrom_Eingabe = Entry(Fenster)                                               
    BaseFrom_Eingabe.grid(row=0, column=2)
    BaseTo_Label = Label(Fenster, text="To low point")
    BaseTo_Label.grid(row=0, column=3)
    BaseTo_Eingabe = Entry(Fenster)                                               
    BaseTo_Eingabe.grid(row=0, column=4)
    BaseUpFrom_Label = Label(Fenster, text="From up point")
    BaseUpFrom_Label.grid(row=0, column=5)
    BaseUpFrom_Eingabe = Entry(Fenster)                                               
    BaseUpFrom_Eingabe.grid(row=0, column=6)
    BaseUpTo_Label = Label(Fenster, text="To up point")
    BaseUpTo_Label.grid(row=0, column=7)
    BaseUpTo_Eingabe = Entry(Fenster)                                               
    BaseUpTo_Eingabe.grid(row=0, column=8)
    
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Take fundamental", variable=var2).grid(row=1, column=0, sticky=W)
    FundFrom_Label = Label(Fenster, text="From low point")
    FundFrom_Label.grid(row=1, column=1)
    FundFrom_Eingabe = Entry(Fenster)                                               
    FundFrom_Eingabe.grid(row=1, column=2)
    FundTo_Label = Label(Fenster, text="To low point")
    FundTo_Label.grid(row=1, column=3)
    FundTo_Eingabe = Entry(Fenster)                                               
    FundTo_Eingabe.grid(row=1, column=4)
    FundUpFrom_Label = Label(Fenster, text="From up point")
    FundUpFrom_Label.grid(row=1, column=5)
    FundUpFrom_Eingabe = Entry(Fenster)                                               
    FundUpFrom_Eingabe.grid(row=1, column=6)
    FundUpTo_Label = Label(Fenster, text="To up point")
    FundUpTo_Label.grid(row=1, column=7)
    FundUpTo_Eingabe = Entry(Fenster)                                               
    FundUpTo_Eingabe.grid(row=1, column=8)
    
    
    var3 = IntVar()
    Checkbutton(Fenster, text="Take 1-st harmonic", variable=var3).grid(row=2, column=0, sticky=W)
    FirstFrom_Label = Label(Fenster, text="From low point")
    FirstFrom_Label.grid(row=2, column=1)
    FirstFrom_Eingabe = Entry(Fenster)                                               
    FirstFrom_Eingabe.grid(row=2, column=2)
    FirstTo_Label = Label(Fenster, text="To low point")
    FirstTo_Label.grid(row=2, column=3)
    FirstTo_Eingabe = Entry(Fenster)                                               
    FirstTo_Eingabe.grid(row=2, column=4)
    FirstUpFrom_Label = Label(Fenster, text="From up point")
    FirstUpFrom_Label.grid(row=2, column=5)
    FirstUpFrom_Eingabe = Entry(Fenster)                                               
    FirstUpFrom_Eingabe.grid(row=2, column=6)
    FirstUpTo_Label = Label(Fenster, text="To up point")
    FirstUpTo_Label.grid(row=2, column=7)
    FirstUpTo_Eingabe = Entry(Fenster)                                               
    FirstUpTo_Eingabe.grid(row=2, column=8)
    
    
    var4 = IntVar()
    Checkbutton(Fenster, text="Take 2-nd harmonic", variable=var4).grid(row=3, column=0, sticky=W)
    SecondFrom_Label = Label(Fenster, text="From low point")
    SecondFrom_Label.grid(row=3, column=1)
    SecondFrom_Eingabe = Entry(Fenster)                                               
    SecondFrom_Eingabe.grid(row=3, column=2)
    SecondTo_Label = Label(Fenster, text="To low point")
    SecondTo_Label.grid(row=3, column=3)
    SecondTo_Eingabe = Entry(Fenster)                                               
    SecondTo_Eingabe.grid(row=3, column=4)
    SecondUpFrom_Label = Label(Fenster, text="From up point")
    SecondUpFrom_Label.grid(row=3, column=5)
    SecondUpFrom_Eingabe = Entry(Fenster)                                               
    SecondUpFrom_Eingabe.grid(row=3, column=6)
    SecondUpTo_Label = Label(Fenster, text="To up point")
    SecondUpTo_Label.grid(row=3, column=7)
    SecondUpTo_Eingabe = Entry(Fenster)                                               
    SecondUpTo_Eingabe.grid(row=3, column=8)
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Take 3-rd harmonic", variable=var5).grid(row=4, column=0, sticky=W)
    ThirdFrom_Label = Label(Fenster, text="From low point")
    ThirdFrom_Label.grid(row=4, column=1)
    ThirdFrom_Eingabe = Entry(Fenster)                                               
    ThirdFrom_Eingabe.grid(row=4, column=2)
    ThirdTo_Label = Label(Fenster, text="To low point")
    ThirdTo_Label.grid(row=4, column=3)
    ThirdTo_Eingabe = Entry(Fenster)                                               
    ThirdTo_Eingabe.grid(row=4, column=4)
    ThirdUpFrom_Label = Label(Fenster, text="From up point")
    ThirdUpFrom_Label.grid(row=4, column=5)
    ThirdUpFrom_Eingabe = Entry(Fenster)                                               
    ThirdUpFrom_Eingabe.grid(row=4, column=6)
    ThirdUpTo_Label = Label(Fenster, text="To up point")
    ThirdUpTo_Label.grid(row=4, column=7)
    ThirdUpTo_Eingabe = Entry(Fenster)                                               
    ThirdUpTo_Eingabe.grid(row=4, column=8)
    
    
    var6 = IntVar()
    Checkbutton(Fenster, text="Take 4-th harmonic", variable=var6).grid(row=5, column=0, sticky=W)
    FourthFrom_Label = Label(Fenster, text="From low point")
    FourthFrom_Label.grid(row=5, column=1)
    FourthFrom_Eingabe = Entry(Fenster)                                               
    FourthFrom_Eingabe.grid(row=5, column=2)
    FourthTo_Label = Label(Fenster, text="To low point")
    FourthTo_Label.grid(row=5, column=3)
    FourthTo_Eingabe = Entry(Fenster)                                               
    FourthTo_Eingabe.grid(row=5, column=4)
    FourthUpFrom_Label = Label(Fenster, text="From up point")
    FourthUpFrom_Label.grid(row=5, column=5)
    FourthUpFrom_Eingabe = Entry(Fenster)                                               
    FourthUpFrom_Eingabe.grid(row=5, column=6)
    FourthUpTo_Label = Label(Fenster, text="To up point")
    FourthUpTo_Label.grid(row=5, column=7)
    FourthUpTo_Eingabe = Entry(Fenster)                                               
    FourthUpTo_Eingabe.grid(row=5, column=8)
    
    
    var7 = IntVar()
    Checkbutton(Fenster, text="Take 5-th harmonic", variable=var7).grid(row=6, column=0, sticky=W)
    FifthFrom_Label = Label(Fenster, text="From low point")
    FifthFrom_Label.grid(row=6, column=1)
    FifthFrom_Eingabe = Entry(Fenster)                                               
    FifthFrom_Eingabe.grid(row=6, column=2)
    FifthTo_Label = Label(Fenster, text="To low point")
    FifthTo_Label.grid(row=6, column=3)
    FifthTo_Eingabe = Entry(Fenster)                                               
    FifthTo_Eingabe.grid(row=6, column=4)
    FifthUpFrom_Label = Label(Fenster, text="From up point")
    FifthUpFrom_Label.grid(row=6, column=5)
    FifthUpFrom_Eingabe = Entry(Fenster)                                               
    FifthUpFrom_Eingabe.grid(row=6, column=6)
    FifthUpTo_Label = Label(Fenster, text="To up point")
    FifthUpTo_Label.grid(row=6, column=7)
    FifthUpTo_Eingabe = Entry(Fenster)                                               
    FifthUpTo_Eingabe.grid(row=6, column=8)
    
    
    var8 = IntVar()
    Checkbutton(Fenster, text="Take 6-th harmonic", variable=var8).grid(row=7, column=0, sticky=W)
    SixthFrom_Label = Label(Fenster, text="From low point")
    SixthFrom_Label.grid(row=7, column=1)
    SixthFrom_Eingabe = Entry(Fenster)                                               
    SixthFrom_Eingabe.grid(row=7, column=2)
    SixthTo_Label = Label(Fenster, text="To low point")
    SixthTo_Label.grid(row=7, column=3)
    SixthTo_Eingabe = Entry(Fenster)                                               
    SixthTo_Eingabe.grid(row=7, column=4)
    SixthUpFrom_Label = Label(Fenster, text="From up point")
    SixthUpFrom_Label.grid(row=7, column=5)
    SixthUpFrom_Eingabe = Entry(Fenster)                                               
    SixthUpFrom_Eingabe.grid(row=7, column=6)
    SixthUpTo_Label = Label(Fenster, text="To up point")
    SixthUpTo_Label.grid(row=7, column=7)
    SixthUpTo_Eingabe = Entry(Fenster)                                               
    SixthUpTo_Eingabe.grid(row=7, column=8)
    
    
    var9 = IntVar()
    Checkbutton(Fenster, text="Save IFFT CVs", variable=var9).grid(row=8, column=0, sticky=W)
    
    
    def IFFTAccepter():
        global BaseInversion
        global FundInversion
        global v1stInversion
        global v2ndInversion
        global v3rdInversion
        global v4thInversion
        global v5thInversion
        global v6thInversion
        global FFTAstxtSaver
        global AstxtSaver
    
        BaseInversion  = var1.get()
        FundInversion  = var2.get()
        v1stInversion  = var3.get()
        v2ndInversion  = var4.get()
        v3rdInversion  = var5.get()
        v4thInversion  = var6.get()
        v5thInversion  = var7.get()
        v6thInversion  = var8.get()
        FFTAstxtSaver  = var9.get()
        
        
        if FFTAstxtSaver == 1:
            if Superposition_Mode == 1:
                AstxtSaver = 1              #required redefinition od AstxtSaver
    
        if BaseInversion == 1:
            global BaseInvStartPoint
            global BaseInvEndPoint
            BaseInvStartPoint  = int(BaseFrom_Eingabe.get())
            BaseInvEndPoint    = int(BaseTo_Eingabe.get())
            global BaseUpInvStartPoint
            global BaseUpInvEndPoint
            BaseUpInvStartPoint  = int(BaseUpFrom_Eingabe.get())
            BaseUpInvEndPoint    = int(BaseUpTo_Eingabe.get())
        
        
        if FundInversion == 1:
            global FundInvStartPoint
            global FundInvEndPoint
            FundInvStartPoint  = int(FundFrom_Eingabe.get())
            FundInvEndPoint    = int(FundTo_Eingabe.get())
            global FundUpInvStartPoint
            global FundUpInvEndPoint
            FundUpInvStartPoint  = int(FundUpFrom_Eingabe.get())
            FundUpInvEndPoint    = int(FundUpTo_Eingabe.get())
        
        if v1stInversion == 1:
            global v1stInvStartPoint
            global v1stInvEndPoint
            v1stInvStartPoint  = int(FirstFrom_Eingabe.get())
            v1stInvEndPoint    = int(FirstTo_Eingabe.get())
            global v1stUpInvStartPoint
            global v1stUpInvEndPoint
            v1stUpInvStartPoint  = int(FirstUpFrom_Eingabe.get())
            v1stUpInvEndPoint    = int(FirstUpTo_Eingabe.get())
        
        if v2ndInversion == 1:
            global v2ndInvStartPoint
            global v2ndInvEndPoint
            v2ndInvStartPoint  = int(SecondFrom_Eingabe.get())
            v2ndInvEndPoint    = int(SecondTo_Eingabe.get())
            global v2ndUpInvStartPoint
            global v2ndUpInvEndPoint
            v2ndUpInvStartPoint  = int(SecondUpFrom_Eingabe.get())
            v2ndUpInvEndPoint    = int(SecondUpTo_Eingabe.get())
        
        if v3rdInversion == 1:
            global v3rdInvStartPoint
            global v3rdInvEndPoint
            v3rdInvStartPoint  = int(ThirdFrom_Eingabe.get())
            v3rdInvEndPoint    = int(ThirdTo_Eingabe.get())
            global v3rdUpInvStartPoint
            global v3rdUpInvEndPoint
            v3rdUpInvStartPoint  = int(ThirdUpFrom_Eingabe.get())
            v3rdUpInvEndPoint    = int(ThirdUpTo_Eingabe.get())
        
        if v4thInversion == 1:
            global v4thInvStartPoint
            global v4thInvEndPoint
            v4thInvStartPoint  = int(FourthFrom_Eingabe.get())
            v4thInvEndPoint    = int(FourthTo_Eingabe.get())
            global v4thUpInvStartPoint
            global v4thUpInvEndPoint
            v4thUpInvStartPoint  = int(FourthUpFrom_Eingabe.get())
            v4thUpInvEndPoint    = int(FourthUpTo_Eingabe.get())
        
        if v5thInversion == 1:
            global v5thInvStartPoint
            global v5thInvEndPoint
            v5thInvStartPoint  = int(FifthFrom_Eingabe.get())
            v5thInvEndPoint    = int(FifthTo_Eingabe.get())
            global v5thUpInvStartPoint
            global v5thUpInvEndPoint
            v5thUpInvStartPoint  = int(FifthUpFrom_Eingabe.get())
            v5thUpInvEndPoint    = int(FifthUpTo_Eingabe.get())
        
        if v6thInversion == 1:
            global v6thInvStartPoint
            global v6thInvEndPoint
            v6thInvStartPoint  = int(SixthFrom_Eingabe.get())
            v6thInvEndPoint    = int(SixthTo_Eingabe.get())
            global v6thUpInvStartPoint
            global v6thUpInvEndPoint
            v6thUpInvStartPoint  = int(SixthUpFrom_Eingabe.get())
            v6thUpInvEndPoint    = int(SixthUpTo_Eingabe.get())
        
    
    button=Button(Fenster,text="Accept",command=IFFTAccepter).grid(row=9,column=0)                                        
    button=Button(Fenster,text="Next",command=IFFT_Calcer).grid(row=10,column=0)
        
    
def IFFT_Calcer():
    global BaseCV
    global Fundamental
    global FirstHarm
    global SecondHarm
    global ThirdHarm
    global FourthHarm
    global FifthHarm
    global SixthHarm
    
    global BaseToInvert
    global FundToInvert
    global FirstToInvert
    global SecondToInvert
    global ThirdToInvert
    global FourthToInvert
    global FifthToInvert
    global SixthToInvert
    
   
    if Superposition_Mode == 0:
        BaseToInvert   = fft(CurrCalc[:-1:])
        FundToInvert   = fft(CurrCalc[:-1:])
        FirstToInvert  = fft(CurrCalc[:-1:])
        SecondToInvert = fft(CurrCalc[:-1:])
        ThirdToInvert  = fft(CurrCalc[:-1:])
        FourthToInvert = fft(CurrCalc[:-1:])
        FifthToInvert  = fft(CurrCalc[:-1:])
        SixthToInvert  = fft(CurrCalc[:-1:])
    
    if Superposition_Mode == 1:
        BaseToInvert   = fft(yyy_Add[:-1:])
        FundToInvert   = fft(yyy_Add[:-1:])
        FirstToInvert  = fft(yyy_Add[:-1:])
        SecondToInvert = fft(yyy_Add[:-1:])
        ThirdToInvert  = fft(yyy_Add[:-1:])
        FourthToInvert = fft(yyy_Add[:-1:])
        FifthToInvert  = fft(yyy_Add[:-1:])
        SixthToInvert  = fft(yyy_Add[:-1:])
    
    
    #=======================================================================================
    #=======================================================================================
    #Time vs. IFFT BaseCV Current
    #=======================================================================================
    #=======================================================================================
    
    if BaseInversion == 1:
        root = Toplevel()                                                         
        root.title("Base-CV")                         
        root.geometry("500x500")   
        for i in range(len(BaseToInvert)):
            if i < BaseInvStartPoint:
                BaseToInvert[i] = 0 +0j
            if i > BaseUpInvEndPoint :
                BaseToInvert[i] = 0 +0j
            if i > BaseInvEndPoint and i< BaseUpInvStartPoint:
                BaseToInvert[i] = 0 +0j
        BaseCV = ifft(BaseToInvert)
        Abbildung = Figure(figsize=(4, 4), dpi=100)
        b1 = Abbildung.add_subplot(111)
        if Superposition_Mode == 0:
            b1.plot(t_Array_cont[:-1:],BaseCV[::],color='r')
        if Superposition_Mode == 1:
            b1.plot(Input_Time[:-1:],BaseCV[::],color='r')
        b1.set_xlabel('t / s', fontsize=12)
        b1.set_ylabel('I / mA', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        Abbildung.tight_layout()
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
    #=======================================================================================
    #=======================================================================================
    #Time vs. IFFT FundamentalCV Current
    #=======================================================================================
    #=======================================================================================
    
    if FundInversion == 1:
        root = Toplevel()                                                         
        root.title("Fundamental-inversion-CV")                         
        root.geometry("500x500")   
        for i in range(len(FundToInvert)):
            if i < FundInvStartPoint:
                FundToInvert[i] = 0+0j
            if i > FundUpInvEndPoint :
                FundToInvert[i] = 0+0j
            if i > FundInvEndPoint and i < FundUpInvStartPoint:
                FundToInvert[i] = 0+0j
        Fundamental = ifft(FundToInvert)
        Abbildung = Figure(figsize=(4, 4), dpi=100)
        b1 = Abbildung.add_subplot(111)
        if Superposition_Mode == 0:
            b1.plot(t_Array_cont[:-1:],Fundamental[::],color='r')
        if Superposition_Mode == 1:
            b1.plot(Input_Time[:-1:],Fundamental[::],color='r')
        b1.set_xlabel('t / s', fontsize=12)
        b1.set_ylabel('I / mA', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        Abbildung.tight_layout()
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
    #=======================================================================================
    #=======================================================================================
    #Time vs. IFFT First HarminicCV Current
    #=======================================================================================
    #=======================================================================================
    
    if v1stInversion == 1:
        root = Toplevel()                                                         
        root.title("First harmonic-inversion-CV")                         
        root.geometry("500x500")   
        for i in range(len(FirstToInvert)):
            if i < v1stInvStartPoint:
                FirstToInvert[i] = 0+0j
            if i > v1stUpInvEndPoint:
                FirstToInvert[i] = 0+0j
            if i > v1stInvEndPoint and i < v1stUpInvStartPoint:
                FirstToInvert[i] = 0+0j
        FirstHarm = ifft(FirstToInvert)
        Abbildung = Figure(figsize=(4, 4), dpi=100)
        b1 = Abbildung.add_subplot(111)
        if Superposition_Mode == 0:
            b1.plot(t_Array_cont[:-1:],FirstHarm[::],color='r')
        if Superposition_Mode == 1:
            b1.plot(Input_Time[:-1:],FirstHarm[::],color='r')
        b1.set_xlabel('t / s', fontsize=12)
        b1.set_ylabel('I / mA', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        Abbildung.tight_layout()
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
    #=======================================================================================
    #=======================================================================================
    #Time vs. IFFT Second HarminicCV Current
    #=======================================================================================
    #=======================================================================================
    
    if v2ndInversion == 1:
        root = Toplevel()                                                         
        root.title("Second harmonic-inversion-CV")                         
        root.geometry("500x500")   
        for i in range(len(SecondToInvert)):
            if i < v2ndInvStartPoint:
                SecondToInvert[i] = 0+0j
            if i > v2ndUpInvEndPoint:
                SecondToInvert[i] = 0+0j
            if i > v2ndInvEndPoint and i< v2ndUpInvStartPoint:
                SecondToInvert[i] = 0+0j
        SecondHarm = ifft(SecondToInvert)
        Abbildung = Figure(figsize=(4, 4), dpi=100)
        b1 = Abbildung.add_subplot(111)
        if Superposition_Mode == 0:
            b1.plot(t_Array_cont[:-1:],SecondHarm[::],color='r')
        if Superposition_Mode == 1:
            b1.plot(Input_Time[:-1:],SecondHarm[::],color='r')
        b1.set_xlabel('t / s', fontsize=12)
        b1.set_ylabel('I / mA', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        Abbildung.tight_layout()
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
    #=======================================================================================
    #=======================================================================================
    #Time vs. IFFT Third HarminicCV Current
    #=======================================================================================
    #=======================================================================================
    
    if v3rdInversion == 1:
        root = Toplevel()                                                         
        root.title("Third harmonic-inversion-CV")                         
        root.geometry("500x500")   
        for i in range(len(ThirdToInvert)):
            if i < v3rdInvStartPoint:
                ThirdToInvert[i] = 0+0j
            if i > v3rdUpInvEndPoint:
                ThirdToInvert[i] = 0+0j
            if i > v3rdInvEndPoint and i < v3rdUpInvStartPoint:
                ThirdToInvert[i] = 0+0j
        ThirdHarm = ifft(ThirdToInvert)
        Abbildung = Figure(figsize=(4, 4), dpi=100)
        b1 = Abbildung.add_subplot(111)
        if Superposition_Mode == 0:
            b1.plot(t_Array_cont[:-1:],ThirdHarm[::],color='r')
        if Superposition_Mode == 1:
            b1.plot(Input_Time[:-1:],ThirdHarm[::],color='r')
        b1.set_xlabel('t / s', fontsize=12)
        b1.set_ylabel('I / mA', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        Abbildung.tight_layout()
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
    #=======================================================================================
    #=======================================================================================
    #Time vs. IFFT Fourth HarminicCV Current
    #=======================================================================================
    #=======================================================================================
    
    if v4thInversion == 1:
        root = Toplevel()                                                         
        root.title("Fourth harmonic-inversion-CV")                         
        root.geometry("500x500")   
        for i in range(len(FourthToInvert)):
            if i < v4thInvStartPoint:
                FourthToInvert[i] = 0+0j
            if i > v4thUpInvEndPoint:
                FourthToInvert[i] = 0+0j
            if i > v4thInvEndPoint and i< v4thUpInvStartPoint:
                FourthToInvert[i] = 0+0j
        FourthHarm = ifft(FourthToInvert)
        Abbildung = Figure(figsize=(4, 4), dpi=100)
        b1 = Abbildung.add_subplot(111)
        if Superposition_Mode == 0:
            b1.plot(t_Array_cont[:-1:],FourthHarm[::],color='r')
        if Superposition_Mode == 1:
            b1.plot(Input_Time[:-1:],FourthHarm[::],color='r')
        b1.set_xlabel('t / s', fontsize=12)
        b1.set_ylabel('I / mA', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        Abbildung.tight_layout()
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
    #=======================================================================================
    #=======================================================================================
    #Time vs. IFFT Fifth HarminicCV Current
    #=======================================================================================
    #=======================================================================================
    
    if v5thInversion == 1:
        root = Toplevel()                                                         
        root.title("Fifth harmonic-inversion-CV")                         
        root.geometry("500x500")   
        for i in range(len(FifthToInvert)):
            if i < v5thInvStartPoint:
                FifthToInvert[i] = 0+0j
            if i > v5thUpInvEndPoint:
                FifthToInvert[i] = 0+0j
            if i > v5thInvEndPoint and i< v5thUpInvStartPoint:
                FifthToInvert[i] = 0+0j
        FifthHarm = ifft(FifthToInvert)
        Abbildung = Figure(figsize=(4, 4), dpi=100)
        b1 = Abbildung.add_subplot(111)
        if Superposition_Mode == 0:
            b1.plot(t_Array_cont[:-1:],FifthHarm[::],color='r')
        if Superposition_Mode == 1:
            b1.plot(Input_Time[:-1:],FifthHarm[::],color='r')
        b1.set_xlabel('t / s', fontsize=12)
        b1.set_ylabel('I / mA', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        Abbildung.tight_layout()
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
    #=======================================================================================
    #=======================================================================================
    #Time vs. IFFT Sixth HarminicCV Current
    #=======================================================================================
    #=======================================================================================
    
    if v6thInversion == 1:
        root = Toplevel()                                                         
        root.title("Sixth harmonic-inversion-CV")                         
        root.geometry("500x500")   
        for i in range(len(SixthToInvert)):
            if i < v6thInvStartPoint:
                SixthToInvert[i] = 0+0j
            if i > v6thUpInvEndPoint:
                SixthToInvert[i] = 0+0j
            if i > v6thInvEndPoint and i< v6thUpInvStartPoint:
                SixthToInvert[i] = 0+0j
        SixthHarm = ifft(SixthToInvert)
        Abbildung = Figure(figsize=(4, 4), dpi=100)
        b1 = Abbildung.add_subplot(111)
        if Superposition_Mode == 0:
            b1.plot(t_Array_cont[:-1:],SixthHarm[::],color='r')
        if Superposition_Mode == 1:
            b1.plot(Input_Time[:-1:],SixthHarm[::],color='r')
        b1.set_xlabel('t / s', fontsize=12)
        b1.set_ylabel('I / mA', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        Abbildung.tight_layout()
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    if AstxtSaver == 1:
        Simulated_RISCV_as_txt_saver()

      


def Simulated_RISCV_as_txt_saver():
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")

        fi.write("n")
        fi.write("\t")
        fi.write(str(n))
        fi.write("\n")
        fi.write("T in [K]")
        fi.write("\t")
        fi.write(str(T))
        fi.write("\n")
        if Unequal_D ==0:
            fi.write("Common D in [cm^2/s]")
            fi.write("\t")
            fi.write(str(D_f))
            fi.write("\n")
        if Unequal_D ==1:
            fi.write("D_foreward in [cm^2/s]")
            fi.write("\t")
            fi.write(str(D_f))
            fi.write("\n")
            fi.write("D_backward in [cm^2/s]")
            fi.write("\t")
            fi.write(str(D_b))
            fi.write("\n")
        fi.write("alpha_forward")
        fi.write("\t")
        fi.write(str(alpha))
        fi.write("\n")
        fi.write("k_zero in [cm/s]")
        fi.write("\t")
        fi.write(str(kzero))
        fi.write("\n")
        if FinKin ==1:
            fi.write("k_het. max. in [cm/s]")
            fi.write("\t")
            fi.write(str(kzero/kfin))
            fi.write("\n")
        
        fi.write("A in [cm^2]")
        fi.write("\t")
        fi.write(str(A))
        fi.write("\n")
        fi.write("c in [mol/cm^3]")
        fi.write("\t")
        fi.write(str(c))
        fi.write("\n")
        
        
       
        if Precedeq == 1:
            fi.write("k_p [1/s]")
            fi.write("\t")
            fi.write(str(kp))
            fi.write("\n")
            fi.write("k_-p [1/s]")
            fi.write("\t")
            fi.write(str(kmp))
            fi.write("\n")
            fi.write("Kp")
            fi.write("\t")
            fi.write(str(Kp))
            fi.write("\n")
            fi.write("p")
            fi.write("\t")
            fi.write(str(p))
            fi.write("\n")
            
        if Succeedeq == 1:
            fi.write("k_f [1/s]")
            fi.write("\t")
            fi.write(str(kf))
            fi.write("\n")
            fi.write("k_-f [1/s]")
            fi.write("\t")
            fi.write(str(kmf))
            fi.write("\n")
            fi.write("Kf")
            fi.write("\t")
            fi.write(str(Kf))
            fi.write("\n")
            fi.write("f")
            fi.write("\t")
            fi.write(str(f))
            fi.write("\n")
        
        if Statistical  == 0:
    
            if model_b == 1 or model_b1 == 1 :
                fi.write("d in [cm]")
                fi.write("\t")
                fi.write(str(distance))
                fi.write("\n")

            if model_c == 1 or model_e == 1 or model_f == 1 or model_g == 1:
                fi.write("radius in [cm]")
                fi.write("\t")
                fi.write(str(r))
                fi.write("\n")

            if model_d   == 1:
                fi.write("radius in [cm]")
                fi.write("\t")
                fi.write(str(r))
                fi.write("\n")
                fi.write("d in [cm]")
                fi.write("\t")
                fi.write(str(distance))
                fi.write("\n")
        
        if Statistical  == 1:
            if DistFunc_Modifier == 1:
                fi.write("Dist. Func. was modified")
                fi.write("\n")
            if DistFunc_Modifier == 0:
                fi.write("Default Dist. Func. was used")
                fi.write("\n")
            
            fi.write("Dist. Func. param. a")
            fi.write("\t")
            fi.write(str(DFP_a))
            fi.write("\n")
            
            fi.write("Dist. Func. param. b")
            fi.write("\t")
            fi.write(str(DFP_b))
            fi.write("\n")
            
            fi.write("Dist. Func. param. c")
            fi.write("\t")
            fi.write(str(DFP_c))
            fi.write("\n")
            
            if model_b == 1:
                fi.write("sheets per mm")   
            if model_d == 1 or model_e == 1:
                fi.write("cylinders per mm^2")
            if model_g == 1:
                fi.write("spheres per mm^3")
            fi.write("\t")
            fi.write(str(Number_N))
            fi.write("\n")
                
            
            if model_d   == 1:
                fi.write("radius of Fiber in 10^-6[m]")
                fi.write("\t")
                fi.write(str(r*10000))
                fi.write("\n")
            if model_b   == 1:
                fi.write("sheet thickness 10^-6[m]")
                fi.write("\t")
                fi.write(str(sheet_thickness*10000))
                fi.write("\n")
            fi.write("average midpoint distance 10^-6[m]")
            fi.write("\t")
            fi.write(str(average_d_or_r))
            fi.write("\n")
                
        
        for i in range(3):
            fi.write("\n")
        
        fi.write("CURVE DATA")
        
        
#---------------------------------------------------------------------------------------------------------------------------        
#Wenn der Fit null ist... also wenn nur berechnet wird...  

        if FIT == 0:  
            
            fi.write("\n")
            fi.write("t [s]")
            fi.write("\t")
            fi.write("E_vs_E_zero [V]")
            fi.write("\t")
            fi.write("I[mA]")
            fi.write("\t")
            if FFTAstxtSaver == 1:
                if Superposition_Mode == 0:
                    FTCurrCalc = fft(CurrCalc[:-1:])
                if Superposition_Mode == 1:
                    FTyyy_Add = fft(yyy_Add[:-1:])
                fi.write("Point No")
                fi.write("\t")
                fi.write("RE(fft(I / mA)")
                fi.write("\t")
                fi.write("IM(fft(I / mA)")
                fi.write("\t")
                fi.write("ABS(fft(I / mA)/max(ABS(fft(I / mA))")
                fi.write("\t")
                if BaseInversion == 1:
                    fi.write("Base CV I[mA]")
                    fi.write("\t")
                if FundInversion == 1:
                    fi.write("Fundumantal CV I[mA]")
                    fi.write("\t")
                if v1stInversion == 1:
                    fi.write("First harmonic CV I[mA]")
                    fi.write("\t")
                if v2ndInversion == 1:
                    fi.write("Second harmonic CV I[mA]")
                    fi.write("\t")
                if v3rdInversion == 1:
                    fi.write("Third harmonic CV I[mA]")
                    fi.write("\t")
                if v4thInversion == 1:
                    fi.write("Fourth harmonic CV I[mA]")
                    fi.write("\t")
                if v5thInversion == 1:
                    fi.write("Fifth harmonic CV I[mA]")
                    fi.write("\t")
                if v6thInversion == 1:
                    fi.write("Sixth harmonic CV I[mA]")
                    fi.write("\t")        
            
            fi.write("\n")
            for i in range(len(xxx)-1):
                fi.write(str(np.asscalar(Input_Time[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Input_Potential[i])))
                fi.write("\t")
                if Superposition_Mode == 0:
                    fi.write(str(np.asscalar(CurrCalc[i])))
                if Superposition_Mode == 1:
                    fi.write(str(np.asscalar(yyy_Add[i])))
                
                fi.write("\t")
                if FFTAstxtSaver == 1:
                    fi.write(str(i))
                    fi.write("\t")
                    if Superposition_Mode == 0:
                        fi.write(str(np.asscalar(FTCurrCalc[i].real)))
                        fi.write("\t")
                        fi.write(str(np.asscalar(FTCurrCalc[i].imag)))
                        fi.write("\t")
                        fi.write(str(np.asscalar(np.abs(FTCurrCalc[i]))/np.max(np.abs(FTCurrCalc))))
                        fi.write("\t")
                    if Superposition_Mode == 1:
                        fi.write(str(np.asscalar(FTyyy_Add[i].real)))
                        fi.write("\t")
                        fi.write(str(np.asscalar(FTyyy_Add[i].imag)))
                        fi.write("\t")
                        fi.write(str(np.asscalar(np.abs(FTyyy_Add[i]))/np.max(np.abs(FTyyy_Add))))
                        fi.write("\t")
                
                    if BaseInversion == 1:
                        fi.write(str(np.asscalar(BaseCV[i].real)))
                        fi.write("\t")
                    if FundInversion == 1:
                        fi.write(str(np.asscalar(Fundamental[i].real)))
                        fi.write("\t")
                    if v1stInversion == 1:
                        fi.write(str(np.asscalar(FirstHarm[i].real)))
                        fi.write("\t")
                    if v2ndInversion == 1:
                        fi.write(str(np.asscalar(SecondHarm[i].real)))
                        fi.write("\t")
                    if v3rdInversion == 1:
                        fi.write(str(np.asscalar(ThirdHarm[i].real)))
                        fi.write("\t")
                    if v4thInversion == 1:
                        fi.write(str(np.asscalar(FourthHarm[i].real)))
                        fi.write("\t")
                    if v5thInversion == 1:
                        fi.write(str(np.asscalar(FifthHarm[i].real)))
                        fi.write("\t")
                    if v6thInversion == 1:
                        fi.write(str(np.asscalar(SixthHarm[i].real)))
                        fi.write("\t")                
                
                fi.write("\n")
            fi.write("\n")
            fi.write("\n")
            
            
        if FIT == 0:  
            if Statistical  == 1:
                if Save_Individuals == 1:
                    fi.write("Mean center to wall in 10^-6[m]")
                    fi.write("\t")
                    for i in range(len(R_Mean_Array)):
                        fi.write(str(np.asscalar(0.5*RR_Mean_Array[i])))                        
                        fi.write("\t")
                    fi.write("\n")
                    
                    fi.write("Center to wall Integ-Limits in 10^-6[m]")
                    fi.write("\t")
                    for i in range(len(R_Mean_Array)):
                        fi.write(str(np.asscalar(0.5*R_Mean_Array[i])))                        
                        fi.write("\t")
                    fi.write("\n")
                    fi.write("\n")
                    fi.write("t [s]")
                    fi.write("\t")
                    fi.write("E_vs_E_zero [V]")
                    fi.write("\t")
                    for i in range(len(R_Mean_Array)):
                        fi.write("weighted I[mA]")
                        fi.write("\t")
                    fi.write("\n")
                    for i in range(len(Input_Potential)-1):
                        fi.write(str(np.asscalar(Input_Time[i])))
                        fi.write("\t")
                        fi.write(str(np.asscalar(Input_Potential[i])))
                        fi.write("\t")
                        for j in range(len(R_Mean_Array)):
                            fi.write(str(np.asscalar(1000*ChiSuperarray[i,j])))
                            fi.write("\t")       
                        fi.write("\n")
                    fi.write("\n")
                    fi.write("\n")
                
                    
        if FIT == 0:  
            if Statistical  == 1:
                if Save_Dist_Func == 1:
                    fi.write("center to wall / 10^-6 m")
                    fi.write("\t")
                    fi.write("P(center to wall) /(10^-6 m)^-1")
                    fi.write("\t")
                    fi.write("\n")
                    for i in range(len(R_Array_Unchanged)):
                        #weil es ja eigentlich Abstände von Fasermittelpunkten sind!!! durch 2
                        fi.write(str(np.asscalar(0.5*R_Array_Unchanged[i])))
                        fi.write("\t")
                        #Mal 2 damit Normierung stimmt!
                        fi.write(str(np.asscalar(2*P_Array_Unchanged[i])))
                        fi.write("\t")
                        fi.write("\n")
                    fi.write("\n")
                    fi.write("\n")
                    fi.write("\n")                       
                
    root.destroy()


