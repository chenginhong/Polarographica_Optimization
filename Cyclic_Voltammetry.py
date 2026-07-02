
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
from cmath                             import *
from math                              import ceil,floor
from scipy.optimize                    import fsolve

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
#define Open File Functions required for reading experimental data for CV fitting
#=======================================================================================================================
#=======================================================================================================================


def Open_CV_File():
    root = Toplevel()
    root.title("Your Data")
    root.geometry("500x500")
    global File_Was_Loaded
    File_Was_Loaded = 1
    global ExpData   #needed later to differentiate experimental and simulated data of ACCV
    ExpData      = 1
    Delim        = Delimiter
    global data
    if Delim == 1:
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='\t', skip_header = FromRowxx, skip_footer = ToRowxx)
    if Delim == 0:
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='', skip_header = FromRowxx, skip_footer = ToRowxx)
    #============================================================================
    #Plotting of loaded file
    #============================================================================
    f = Figure(figsize=(5, 5), tight_layout = True)
    b = f.add_subplot(111)
    #============================================================================
    global Potenzial
    global Strom
    if FirstComesxx  == 1:
        Potenzial     =  data[::Readeveryxx,0::2] * UmrPot            
        Strom         =  data[::Readeveryxx,1::2] * UmRStrom
    if FirstComesxx  == 0:
        Potenzial     =  data[::Readeveryxx,1::2] * UmrPot
        Strom         =  data[::Readeveryxx,0::2] * UmRStrom
    global Weite
    Weite = Potenzial.shape[0]
    global NumMeas
    NumMeas = Potenzial.shape[1]
    for i in range(Weite):
        Stromarrays             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))                          
        Potenzialarrays         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:])) 
        LenPotArrays            = len(Potenzialarrays)
        b.plot(Potenzial, 0.001*Strom, linestyle='-',marker='',color='k')
        b.set_xlabel('E vs. Ref. / V', fontsize=12)
        b.set_ylabel('I / mA', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b.spines[axis].set_linewidth(2)
            b.spines[axis].set_color('k')
        b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    
    
    
def Open_CV_File_PolArStat():
    root = Toplevel()
    root.title("Your Data")
    root.geometry("500x500")
    def Arraycondenser(ARRAY):
        NLB   = 0
        Out   = np.zeros(len(ARRAY[0,1::]))
        for i in range(len(ARRAY[::,0])-1):
            if ARRAY[i+1,0] != ARRAY[i,0]:
                B   = np.mean(ARRAY[NLB:i+1,1::], axis = 0)
                B[0] = ARRAY[i,1]
                Out = np.vstack([Out,B])        
                NLB = i+1
        return Out[1::,::]
    
    Path           = askopenfilename()
    data           = np.genfromtxt(Path, skip_header = 14, skip_footer = 1)
    data_condensed = Arraycondenser(data[1::,::])
    
    #============================================================================
    #Plotting of loaded file
    #============================================================================
    f = Figure(figsize=(5, 5), tight_layout = True)
    b = f.add_subplot(111)
    #============================================================================
    global Potenzial
    global Strom
    times_of_exp  = data_condensed[data_condensed[::,3] == Cycle_to_analyze, 0] 
    Potenzial     = data_condensed[data_condensed[::,3] == Cycle_to_analyze, 1]      
    Strom         = data_condensed[data_condensed[::,3] == Cycle_to_analyze, 2]
    
    PotenzialRAW     = data[data[::,4] == Cycle_to_analyze, 2]      
    StromRAW         = data[data[::,4] == Cycle_to_analyze, 3]
    
    b.plot(PotenzialRAW, StromRAW, linestyle='-',marker='',color='lightgrey', label = "raw data")
    b.plot(Potenzial, Strom, linestyle='-',marker='',color='k', label = "step-average")
    b.legend(frameon = False, fontsize = 12)
    b.set_xlabel('E vs. Ref. / V', fontsize=12)
    b.set_ylabel('I / mA', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b.spines[axis].set_linewidth(2)
        b.spines[axis].set_color('k')
    b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    txt_filename = asksaveasfilename(defaultextension='.txt')
    outfile_CV_Data = open(txt_filename, 'w')
    outfile_CV_Data.write("E_vs_Ref_in_V")
    outfile_CV_Data.write("\t")
    outfile_CV_Data.write("I_in_microampere")
    outfile_CV_Data.write("\n")
    for i in range(len(Potenzial)):
        outfile_CV_Data.write(str(Potenzial[i]))
        outfile_CV_Data.write("\t")
        outfile_CV_Data.write(str(1000*Strom[i]))
        outfile_CV_Data.write("\n")
    outfile_CV_Data.close()
    



def Get_CV_Data():
    Fenster = Toplevel()                                                         
    Fenster.title("Get CV Data")                         
    Fenster.geometry("400x270")
    
    def Classical_CV_Call():
        Get_CV_Data_Classical()
        def quit():
            Fenster.destroy()
        quit()
        
    def PolArStat_CV_Call():
        Get_CV_Data_PolArStat()
        def quit():
            Fenster.destroy()
        quit()
    
    
    Classical_CV = Button(Fenster, text="Open condensed CV data \nfrom PolArStat or\ndata from external device.",command=Classical_CV_Call, width = 20, height = 5)
    Classical_CV.place(x = 35, y = 25)
    PolArStatCV = Button(Fenster, text="Condense CV data from \nPolArStat",command=PolArStat_CV_Call, width = 20, height = 5)
    PolArStatCV.place(x = 210, y = 25)
    



def Get_CV_Data_PolArStat():
    Fenster = Toplevel()                                                         
    Fenster.title("Get Data")                         
    Fenster.geometry("400x270")                                            
       
    Cycle_to_analyze_label      = Label(Fenster,text="Cycle to export")
    Cycle_to_analyze_label.place(x = 35, y = 80)
    Cycle_to_analyze_Eingabe    = Entry(Fenster)   
    Cycle_to_analyze_Eingabe.insert(END, 1)                                             
    Cycle_to_analyze_Eingabe.place(x = 155, y = 80, width = 175, height = 22)
    

    def NextCVPolArStat():
        global Cycle_to_analyze
        Cycle_to_analyze   = (int(Cycle_to_analyze_Eingabe.get()))
        Open_CV_File_PolArStat()
        def quit():
            Fenster.destroy()
        quit()
    Next = Button(Fenster, text="Next",command=NextCVPolArStat)
    Next.place(x = 155, y = 120, width = 175, height = 35)
    

    

def Get_CV_Data_Classical():
    Fenster = Toplevel()                                                         
    Fenster.title("Get Data")                         
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
    UmRStrom_Label = Label(Fenster,text="Current Factor to be Microampere")
    UmRStrom_Label.grid(row=3, column=0)
    UmRStrom_Eingabe = Entry(Fenster)
    UmRStrom_Eingabe.grid(row=3, column=1)
    UmrPot_Label = Label(Fenster,text="Potential Factor to be Volt")
    UmrPot_Label.grid(row=4, column=0)
    UmrPot_Eingabe = Entry(Fenster)
    UmrPot_Eingabe.grid(row=4, column=1)
    DesiReactxxx_Label  = Label(Fenster,text="Desired Reaction")
    DesiReactxxx_Label.grid(row=5, column=0)
    var1 = IntVar()
    Checkbutton(Fenster, text="Oxidation", variable=var1).grid(row=5,column=1, sticky=W)
    var2 = IntVar()
    Checkbutton(Fenster, text="Reduction", variable=var2).grid(row=5,column=2, sticky=W)
    FirstComesxxx_Label  = Label(Fenster,text="First comes")
    FirstComesxxx_Label.grid(row=6, column=0)
    var3 = IntVar()
    Checkbutton(Fenster, text="Potential", variable=var3).grid(row=6,column=1, sticky=W)
    var4 = IntVar()
    Checkbutton(Fenster, text="Current", variable=var4).grid(row=6,column=2, sticky=W)
    Delimiterxxx_Label  = Label(Fenster,text="Delimiter")
    Delimiterxxx_Label.grid(row=7, column=0)
    var5 = IntVar()
    Checkbutton(Fenster, text="Tab", variable=var5).grid(row=7,column=1, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Space", variable=var6).grid(row=7,column=2, sticky=W)
    def AcceptCV():
        global FromRowxx
        FromRowxx     = (int(FromRowxxx_Eingabe.get()))
        global ToRowxx
        ToRowxx       = (int(ToRowxxx_Eingabe.get()))
        global Readeveryxx
        Readeveryxx   = (int(Readeveryxxx_Eingabe.get()))
        global UmRStrom#
        UmRStrom      = (float(UmRStrom_Eingabe.get()))
        global UmrPot#
        UmrPot        = (float(UmrPot_Eingabe.get()))
        global DesiReactxx
        global FirstComesxx
        DesiReactxx   = var1.get()
        FirstComesxx  = var3.get()
        global Delimiter
        Delimiter         = var5.get()
    def NextCV():
        Open_CV_File()
        def quit():
            Fenster.destroy()
        quit()
    Accept = Button(Fenster, text="Accept",command=AcceptCV)
    Accept.grid(row=8, column=0) 
    Next = Button(Fenster, text="Next",command=NextCV)
    Next.grid(row=9, column=0)  
    
    

    

#=======================================================================================================================
#=======================================================================================================================
#Superposition mode Function
#=======================================================================================================================
#=======================================================================================================================
def Initializer():
    global Superposition_Mode
    global Counter
    Counter = 0
    Superposition_Mode = 1
    global Storage_x
    global Storage_y
    Storage_x = PotCalc
    if FIT == 1:
        global Storage_x_DataMeas
        Storage_x_DataMeas = Pot
        
   
        
    Storage_y = np.zeros([len(PotCalc),50])
    
def Append_Latest():
    global Counter
    global Storage_x
    global Storage_y
    Counter += 1
    if Counter >=1 and Counter <= 50:
        if FIT == 0:
            if Statistical == 0:
                if len(Storage_x) != len(CurrCalc):
                    Unequal_length_warner()
                Storage_y[::,(Counter-1)] = CurrCalc[0:len(Storage_x)]
            if Statistical == 1:
                if len(Storage_x) != len(np.squeeze(np.sum(CurrCalc, axis = 1))):
                    Unequal_length_warner()
                Storage_y[::,(Counter-1)] = np.squeeze(np.sum(CurrCalc, axis = 1))[0:len(Storage_x)]
        if FIT == 1:
            Storage_y[::,(Counter-1)] = CurrCalc[0:len(Storage_x)]
    else:
        Counter_Error()
        
def Remove_Latest():
    global Storage_y
    global Counter
    if Counter >=1 and Counter <= 50:
         Storage_y[::,(Counter-1)] = 0
    else:
        Counter_Error()
    Counter -= 1
    
def RS_Mode_Caller():
    global RS
    global ADD
    RS  = 1
    ADD = 0  
    Superpos_Shower()

def Add_Mode_Caller():
    global RS
    global ADD
    RS  = 0
    ADD = 1  
    Superpos_Shower()
    
def Double_Mode_Caller():
    global RS
    global ADD
    RS  = 1
    ADD = 1  
    Superpos_Shower()

def Superpos_Shower():
    global Storage_x
    global Storage_y
    global yyy_RS
    global yyy_Add
    yyy_RS  = Storage_y[::,:Counter:]
    yyy_Add = np.squeeze(np.sum(yyy_RS, axis = 1)) 
    root = Toplevel()
    root.title("Simulated CV")
    Abbildung      = Figure(figsize=(5, 4), tight_layout=True)
    Superimpose    = Abbildung.add_subplot(111)
    if RS  == 1:
        Superimpose.plot(Storage_x[0:-1],yyy_RS[0:-1,::], color = 'r') 
    if ADD == 1:
        Superimpose.plot(Storage_x[0:-1],yyy_Add[0:-1:], color = 'b')
    if FIT == 1:
        NormLengthPot    = (1.0/float(len(Pot[0:-1])))
        NormLengthxxx    = (1.0/float(len(Storage_x[0:-1])))
        NormArrayPot     = np.arange(0,1,NormLengthPot)
        NormArrayxxx     = np.arange(0,1,NormLengthxxx)
        Fitinterpolation = InterpolatedUnivariateSpline(NormArrayxxx,yyy_Add[0:-1:],k=3)
        InterpolatedCurr = np.empty(len(NormArrayPot))
        for i in range(len(NormArrayPot)):
            InterpolatedCurr[i] = Fitinterpolation(NormArrayPot[i])
        Abweichungsarray = np.empty(len(Curr)-1)
        for i in range(len(Curr)-1):
            Abweichungsarray[i] = (((0.001*Curr[i] - InterpolatedCurr[i])/np.max(0.001*np.abs(Curr)))**2)/(float(len(Curr))-1)
        Varianz            = np.sum(Abweichungsarray)
        global Standardabweichung
        Standardabweichung = Varianz**0.5
        Superimpose.plot(Storage_x[0:-1],yyy_Add[0:-1:], color = 'lightblue', marker ='.', linestyle ='') #plotting again but only as markers
        Superimpose.plot(Storage_x_DataMeas[:-1:],0.001*Curr[:-1:],color='k')        
        Superimpose.annotate(u'\u03C3' '=%5.3f'  % Standardabweichung , xy=(0.05, 0.83), xycoords='axes fraction')
    Superimpose.axvline(0,color='k')
    Superimpose.axhline(0,color='k')
    Superimpose.set_xlabel('E vs. ' 'E$^0$' '(1$^{st}$ CV)' ' / V', fontsize=12)
    Superimpose.set_ylabel('I / mA', fontsize=12)
    for axis in ['top','bottom','left','right']:
        Superimpose.spines[axis].set_linewidth(2)
        Superimpose.spines[axis].set_color('k')
    Superimpose.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)



def Superpos_Saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))
    with open(root.filename,"w") as fi:
        if FIT == 1:
            fi.write("Standard deviation")
            fi.write("\t")
            fi.write(str(np.asscalar(Standardabweichung)))        
            fi.write("\n")
            fi.write("\n")
        fi.write("Additive Mode")
        fi.write("\n")
        fi.write("\n")
        fi.write("t in [s]")
        fi.write("\t")
        fi.write("Calculated Potential vs. E_zero in [V]")
        fi.write("\t")
        fi.write("Calculated Current in mA")
        fi.write("\n")
        fi.write("\n")
        for i in range(len(Storage_x)-1):
            fi.write(str(np.asscalar(t_Array_cont[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Storage_x[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(yyy_Add[i])))
            fi.write("\t")
            fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Randles-Sevcik Mode")
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
            fi.write(str(np.asscalar(t_Array_cont[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Storage_x[i])))
            fi.write("\t")
            for j in range(len(yyy_RS[0,::])):
                fi.write(str(np.asscalar(yyy_RS[i,j])))
                fi.write("\t")
            fi.write("\n")
        if FIT == 1:
            fi.write("\n")
            fi.write("\n")
            fi.write("Measured and Ru corr. Pot vs. E_zero in [V]")
            fi.write("\t")
            fi.write("Measured Current in mA")
            fi.write("\t")
            fi.write("\n")
            for i in range(len(Pot)):
                fi.write(str(np.asscalar(Pot[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(0.001*Curr[i])))
                fi.write("\t")
                fi.write("\n")            
    root.destroy()
                


    
#=======================================================================================================================
#=======================================================================================================================
#define CV Type Chooser and Model zeroer Functions. Model Zeroer is required that no interference between models will occur
#=======================================================================================================================
#=======================================================================================================================


def Model_Zeroer():
    global model_3D
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
    model_3D        = 0
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

    
    
def External_3DMT():
    Model_Zeroer()
    global model_3D
    model_3D = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()    

def Semi_Inf_Planar():
    Model_Zeroer()
    global model_a
    model_a = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
       
def Finit_Planar():
    Model_Zeroer()
    global model_b
    model_b = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Finit_Planar_Trans():
    Model_Zeroer()
    global model_b1
    model_b1 = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Semi_Inf_Zyl_Ext():
    Model_Zeroer()
    global model_c
    model_c = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Finit_Zyl_Ext():
    Model_Zeroer()
    global model_d
    model_d = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Finit_Zyl_Int():
    Model_Zeroer()
    global model_e
    model_e = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Semi_Inf_Sphere_Ext():
    Model_Zeroer()
    global model_f
    model_f = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Finit_Sphere_Int():
    Model_Zeroer()
    global model_g
    model_g = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
    
    
#-------------------------------------------------------------------------------------
#Fitters
#-------------------------------------------------------------------------------------
    

def External_3DMT_FITTER():
    Model_Zeroer()
    global model_3D
    model_3D = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def Semi_Inf_Planar_FITTER():
    Model_Zeroer()
    global model_a
    model_a = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def Finit_Planar_FITTER():
    Model_Zeroer()
    global model_b
    model_b = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
def Finit_Planar_Trans_FITTER():
    Model_Zeroer()
    global model_b1
    model_b1 = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def Semi_Inf_Zyl_Ext_FITTER():
    Model_Zeroer()
    global model_c
    model_c = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def Finit_Zyl_Ext_FITTER():
    Model_Zeroer()
    global model_d
    model_d = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
def Finit_Zyl_Int_FITTER():
    Model_Zeroer()
    global model_e
    model_e = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
def Semi_Inf_Sphere_Ext_FITTER():
    Model_Zeroer()
    global model_f
    model_f = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def Finit_Sphere_Int_FITTER():
    Model_Zeroer()
    global model_g
    model_g = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
        
        
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#Statistical Simulators
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

def Statistical_Finit_Planar():
    Model_Zeroer()
    global model_b
    global Statistical
    Statistical  = 1
    model_b = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Statistical_Finit_Zyl_Ext():
    Model_Zeroer()
    global model_d
    model_d = 1
    global Statistical
    Statistical  = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Statistical_Finit_Zyl_Int():
    Model_Zeroer()
    global model_e
    global Statistical
    Statistical  = 1
    model_e = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Statistical_Finit_Sphere_Int():
    Model_Zeroer()
    global model_g
    global Statistical
    Statistical  = 1
    model_g = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
    
    
#-------------------------------------------------------------------------------------
#Statistical Fitters
#-------------------------------------------------------------------------------------
    
    
def Statistical_Finit_Planar_FITTER():
    Model_Zeroer()
    global model_b
    global Statistical
    Statistical  = 1
    model_b = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
    
def Statistical_Finit_Zyl_Ext_FITTER():
    Model_Zeroer()
    global model_d
    global Statistical
    Statistical  = 1
    model_d = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
def Statistical_Finit_Zyl_Int_FITTER():
    Model_Zeroer()
    global model_e
    global Statistical
    Statistical  = 1
    model_e = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
    
def Statistical_Finit_Sphere_Int_FITTER():
    Model_Zeroer()
    global model_g
    global Statistical
    Statistical  = 1
    model_g = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()



# In[ ]:

def Eingabe_CV_Simulator():
    Fenster = Toplevel()  
    Fenster.geometry("850x670")
    Fenster.title("Set CV-Parameters")    
    
    colorbgr = Label(Fenster, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000)
    
    n_Label = Label(Fenster, text="n*", bg = '#80ffbf')
    n_Label.place(x = 25, y = 10)
    n_Eingabe = Entry(Fenster) 
    n_Eingabe.insert(END, 1)                                               
    n_Eingabe.place(x = 150, y = 10)
    
    
    alpha_Label = Label(Fenster, text= u"\u03B1*",bg = '#80ffbf')
    alpha_Label.place(x = 25, y = 35)
    alpha_Eingabe = Entry(Fenster)       
    alpha_Eingabe.insert(END, 0.5)                                         
    alpha_Eingabe.place(x = 150, y = 35)
    
    kzero_Label = Label(Fenster,text="k° [cm/s]*",bg = '#80ffbf')
    kzero_Label.place(x = 25, y = 60)
    kzero_Eingabe = Entry(Fenster)
    kzero_Eingabe.insert(END, 1000) 
    kzero_Eingabe.place(x = 150, y = 60)
    
    T_Label = Label(Fenster,text="T [°C]*",bg = '#80ffbf')
    T_Label.place(x = 25, y = 85)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.insert(END, 25) 
    T_Eingabe.place(x = 150, y = 85)
    
    Scanrate_Label = Label(Fenster,text="Scanrate [mV/s]*",bg = '#80ffbf')
    Scanrate_Label.place(x = 300, y = 10)
    Scanrate_Eingabe = Entry(Fenster)
    Scanrate_Eingabe.insert(END, 20) 
    Scanrate_Eingabe.place(x = 425, y = 10)
    
    A_Label = Label(Fenster,text="A in cm^2*",bg = '#80ffbf')
    A_Label.place(x = 300, y = 35)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.insert(END, 20) 
    A_Eingabe.place(x = 425, y = 35)
        
    c_Label = Label(Fenster,text="c in mol/L*",bg = '#80ffbf')
    c_Label.place(x = 300, y = 60)
    c_Eingabe = Entry(Fenster)
    c_Eingabe.insert(END, 0.1) 
    c_Eingabe.place(x = 425, y = 60)    
    
    D_Label = Label(Fenster,text="D 10^-6[cm^2/s]*",bg = '#80ffbf')
    D_Label.place(x = 300, y = 85)
    D_Eingabe = Entry(Fenster)
    D_Eingabe.insert(END, 1) 
    D_Eingabe.place(x = 425, y = 85)
    
    
    if FIT ==1:
        Ezero_Label = Label(Fenster,text="E_zero vs Ref [V]*",bg = '#80ffbf')
        Ezero_Label.place(x = 575, y = 10)
        Ezero_Eingabe = Entry(Fenster)
        Ezero_Eingabe.insert(END, 0) 
        Ezero_Eingabe.place(x = 700, y = 10)
        
        Ru_Label = Label(Fenster,text="Ru in Ohm*",bg = '#80ffbf')
        Ru_Label.place(x = 575, y = 35)
        Ru_Eingabe = Entry(Fenster)
        Ru_Eingabe.insert(END, 0) 
        Ru_Eingabe.place(x = 700, y = 35)
        
    if FIT ==0:
        
        Start_Label = Label(Fenster,text="E_start vs. Ezero [V]*",bg = '#80ffbf')
        Start_Label.place(x = 575, y = 10)
        Start_Eingabe = Entry(Fenster)
        Start_Eingabe.insert(END, -0.4) 
        Start_Eingabe.place(x = 700, y = 10)
        
        Switch_Label = Label(Fenster,text="E_switch vs. Ezero [V]*",bg = '#80ffbf')
        Switch_Label.place(x = 575, y = 35)
        Switch_Eingabe = Entry(Fenster)
        Switch_Eingabe.insert(END, 0.4) 
        Switch_Eingabe.place(x = 700, y = 35)
        
        End_Label = Label(Fenster,text="E_end vs. Ezero [V]*",bg = '#80ffbf')
        End_Label.place(x = 575, y = 60)
        End_Eingabe = Entry(Fenster)
        End_Eingabe.insert(END, -0.4) 
        End_Eingabe.place(x = 700, y = 60)
    
    SepLine0 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine0.place(x = 0, y = 110, height = 10)
    
    
    colorbox0 = Label(Fenster,text="",bg = '#EFEFEF')
    colorbox0.place(x = 20, y = 125, width = 810, height = 55)
    
    var00 = IntVar()
    Checkbutton(Fenster, text="Convolute Ohmic resistance and (non/ideal) double layer capacities", bg = '#EFEFEF', variable=var00).place(x = 25, y = 130)
    
    ROhm_Label = Label(Fenster,text="R in Ohm",bg = '#EFEFEF')
    ROhm_Label.place(x = 25, y = 155)
    ROhm_Eingabe = Entry(Fenster)
    ROhm_Eingabe.insert(END, 1e-3) 
    ROhm_Eingabe.place(x = 150, y = 155)
    
    Cdl_Label = Label(Fenster,text=r"C in Fs^{1-gamma}",bg = '#EFEFEF')
    Cdl_Label.place(x = 300, y = 155)
    Cdl_Eingabe = Entry(Fenster)
    Cdl_Eingabe.insert(END, 1e-6) 
    Cdl_Eingabe.place(x = 425, y =155)
    
    gammaCPE_Label = Label(Fenster,text="gamma CPE",bg = '#EFEFEF')
    gammaCPE_Label.place(x = 575, y =155)
    gammaCPE_Eingabe = Entry(Fenster)
    gammaCPE_Eingabe.insert(END, 1) 
    gammaCPE_Eingabe.place(x = 700, y = 155)
                
            
    
    SepLine1 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine1.place(x = 0, y = 185, height = 10)
    
    
    colorbox1 = Label(Fenster,text="",bg = '#EFEFEF')
    colorbox1.place(x = 20, y = 200, width = 260, height = 85)
    var1 = IntVar()
    Checkbutton(Fenster, text="Unequal D", bg = '#EFEFEF', variable=var1).place(x = 25, y = 205)
    Df_Label = Label(Fenster,text="D_f 10^-6[cm^2/s]",bg = '#EFEFEF')
    Df_Label.place(x = 25, y = 230)
    Df_Eingabe = Entry(Fenster)
    Df_Eingabe.place(x = 150, y = 230)
    Db_Label = Label(Fenster,text="D_b 10^-6[cm^2/s]",bg = '#EFEFEF')
    Db_Label.place(x = 25, y = 255)
    Db_Eingabe = Entry(Fenster)
    Db_Eingabe.place(x = 150, y = 255)
    
    
    colorbox2 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox2.place(x = 295, y = 200, width = 260, height = 85)
    var2 = IntVar()
    Checkbutton(Fenster, text="Preceding ch. eq.",bg = '#EFEFEF', variable=var2).place(x = 300, y = 205)
    kp_Label = Label(Fenster,text="k_p [1/s]", bg = '#EFEFEF')
    kp_Label.place(x = 300, y = 230)
    kp_Eingabe = Entry(Fenster)
    kp_Eingabe.place(x = 425, y = 230)
    kmp_Label = Label(Fenster,text="k_-p [1/s]", bg = '#EFEFEF')
    kmp_Label.place(x = 300, y = 255)
    kmp_Eingabe = Entry(Fenster)
    kmp_Eingabe.place(x = 425, y = 255)
    
    
    colorbox3 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox3.place(x = 570, y = 200, width = 260, height = 85)
    var3 = IntVar()
    Checkbutton(Fenster, text="Following ch. eq.", bg = '#EFEFEF', variable=var3).place(x = 575, y = 205)
    kf_Label = Label(Fenster,text="k_f [1/s]", bg = '#EFEFEF')
    kf_Label.place(x = 575, y = 230)
    kf_Eingabe = Entry(Fenster)
    kf_Eingabe.place(x = 700, y = 230)
    kmf_Label = Label(Fenster,text="k_-f [1/s]", bg = '#EFEFEF')
    kmf_Label.place(x = 575, y = 255)
    kmf_Eingabe = Entry(Fenster)
    kmf_Eingabe.place(x = 700, y = 255)
    
    
    
    colorbox4 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox4.place(x = 20, y = 300, width = 260, height = 60)
    var15 = IntVar()
    Checkbutton(Fenster, text="Finite heterogeneous kinetics ",bg = '#EFEFEF', variable=var15).place(x = 25, y = 305)
    kmax_Label = Label(Fenster,text="fin. het. k [cm/s]",bg = '#EFEFEF')
    kmax_Label.place(x = 25, y = 330)
    kmax_Eingabe = Entry(Fenster)
    kmax_Eingabe.place(x = 150, y = 330)
    
    
    if Statistical  == 0:
        if model_b == 1 or model_b1 == 1:
            colorbox5 = Label(Fenster,text="", bg = '#EFEFEF')
            colorbox5.place(x = 295, y = 300, width = 260, height = 60)
            d_Label = Label(Fenster,text="distance 10^-6[m]*", bg = '#EFEFEF')
            d_Label.place(x = 300, y = 330)
            d_Eingabe = Entry(Fenster)
            d_Eingabe.insert(END, 20) 
            d_Eingabe.place(x = 425, y = 330)
    
        if model_c == 1 or model_e == 1 or model_f == 1 or model_g == 1:
            colorbox5 = Label(Fenster,text="", bg = '#EFEFEF')
            colorbox5.place(x = 295, y = 300, width = 260, height = 60)
            r_Label = Label(Fenster,text="radius 10^-6[m]*", bg = '#EFEFEF')
            r_Label.place(x = 300, y = 300)
            r_Eingabe = Entry(Fenster)
            r_Eingabe.insert(END, 20) 
            r_Eingabe.place(x = 425, y = 300)
        
        if model_d  == 1:
            colorbox5 = Label(Fenster,text="", bg = '#EFEFEF')
            colorbox5.place(x = 295, y = 300, width = 260, height = 60)
            d_Label = Label(Fenster,text="distance 10^-6[m]*", bg = '#EFEFEF')
            d_Label.place(x = 300, y = 305)
            d_Eingabe = Entry(Fenster)
            d_Eingabe.insert(END, 50) 
            d_Eingabe.place(x = 425, y = 305)
            r_Label = Label(Fenster,text="radius 10^-6[m]*", bg = '#EFEFEF')
            r_Label.place(x = 300, y = 330)
            r_Eingabe = Entry(Fenster)
            r_Eingabe.insert(END, 20) 
            r_Eingabe.place(x = 425, y = 330)
      
        
        
    if Statistical  == 1:
        colorbox5 = Label(Fenster,text="", bg = '#EFEFEF')
        colorbox5.place(x = 295, y = 300, width = 260, height = 60)
        DistFunc_Nodes_Label = Label(Fenster,text="Intervals of Dist.*", bg = '#EFEFEF')
        DistFunc_Nodes_Label.place(x = 300, y = 305)
        DistFunc_Nodes_Eingabe = Entry(Fenster)
        DistFunc_Nodes_Eingabe.insert(END, 10) 
        DistFunc_Nodes_Eingabe.place(x = 425, y = 305)
        
        if model_b == 1:
            Number_N_Label = Label(Fenster,text="Sheets per mm*", bg = '#EFEFEF')
        if model_d == 1 or model_e == 1:
            Number_N_Label = Label(Fenster,text="Cylinders per mm^2*", bg = '#EFEFEF')
        if model_g == 1:
            Number_N_Label = Label(Fenster,text="Spheres per mm^3*", bg = '#EFEFEF')
        Number_N_Label.place(x = 300, y = 330)
        Number_N_Eingabe = Entry(Fenster)
        Number_N_Eingabe.insert(END, 100)
        Number_N_Eingabe.place(x = 425, y = 330)
        
        
        
        if model_d == 1:
            colorbox6 = Label(Fenster,text="", bg = '#EFEFEF')
            colorbox6.place(x = 570, y = 300, width = 260, height = 60)
            r_Label = Label(Fenster,text="Fiber radius 10^-6[m]*")
            r_Label.place(x = 575, y = 330)
            r_Eingabe = Entry(Fenster)
            r_Eingabe.insert(END, 5) 
            r_Eingabe.place(x = 700, y = 330)
        
        if model_b == 1:
            colorbox6 = Label(Fenster,text="", bg = '#EFEFEF')
            colorbox6.place(x = 570, y = 300, width = 260, height = 60)
            sheet_thickness_Label = Label(Fenster,text="Sheet thick 10^-6[m]*")
            sheet_thickness_Label.place(x = 575, y = 330)
            sheet_thickness_Eingabe = Entry(Fenster)
            sheet_thickness_Eingabe.insert(END, 2) 
            sheet_thickness_Eingabe.place(x = 700, y = 330)
    
    
    SepLine2 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine2.place(x = 0, y = 365, height = 10)
    
    
    colorbox7 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox7.place(x = 20, y = 380, width = 260, height = 60)
    var17 = IntVar()
    Checkbutton(Fenster, text="Mod. Xi-Resol.", bg = '#EFEFEF', variable=var17).place(x = 25, y = 385)
    ModDeltXi_Label = Label(Fenster,text="delta Xi", bg = '#EFEFEF')
    ModDeltXi_Label.place(x = 25, y = 410)
    ModDeltXi_Eingabe = Entry(Fenster)
    ModDeltXi_Eingabe.place(x = 150, y = 410)
    
    
    
    if FIT ==1:
        colorbox8 = Label(Fenster,text="", bg = '#EFEFEF')
        colorbox8.place(x = 295, y = 380, width = 535, height = 60)
        var7 = IntVar()
        Checkbutton(Fenster, text="Lin. Basecorr", bg = '#EFEFEF', variable=var7).place(x = 300, y = 385)
        LinStart_Label = Label(Fenster,text="Start vs. Ref. [V]", bg = '#EFEFEF')
        LinStart_Label.place(x = 300, y = 410)
        LinStart_Eingabe = Entry(Fenster)
        LinStart_Eingabe.place(x = 425, y = 410)
        LinEnd_Label = Label(Fenster,text="End vs. Ref. [V]", bg = '#EFEFEF')
        LinEnd_Label.place(x = 575, y = 410)
        LinEnd_Eingabe = Entry(Fenster)
        LinEnd_Eingabe.place(x = 700, y = 410)
    
    
    SepLine3 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine3.place(x = 0, y = 445, height = 10)
    
    
    if Statistical  == 1:
        colorbox9 = Label(Fenster,text="", bg = '#EFEFEF')
        colorbox9.place(x = 20, y = 460, width = 810, height = 60) 
        var16 = IntVar()
        Checkbutton(Fenster, text="Modify Distribution function", bg = '#EFEFEF', variable=var16).place(x = 25, y = 465)
        DFP_a_Label = Label(Fenster,text="Parameter a", bg = '#EFEFEF')
        DFP_a_Label.place(x = 25, y = 490)
        DFP_a_Eingabe = Entry(Fenster)
        DFP_a_Eingabe.place(x = 150, y = 490)
        DFP_b_Label = Label(Fenster,text="Parameter b", bg = '#EFEFEF')
        DFP_b_Label.place(x = 300, y = 490)
        DFP_b_Eingabe = Entry(Fenster)
        DFP_b_Eingabe.place(x = 425, y = 490)
        DFP_c_Label = Label(Fenster,text="Parameter c", bg = '#EFEFEF')
        DFP_c_Label.place(x = 575, y = 490)
        DFP_c_Eingabe = Entry(Fenster)
        DFP_c_Eingabe.place(x = 700, y = 490)
        
    if Statistical  == 0:
        NoEntryLabel = Label(Fenster,text="This field is inactive for your desired calculation.", bg = '#80ffbf')
        NoEntryLabel.place(x = 300, y = 480) 
    
    
    SepLine4 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine4.place(x = 0, y = 530, height = 10)
 
    
    
    colorbox10 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox10.place(x = 20, y = 540, width = 260, height = 120) 
    var6 = IntVar()
    Checkbutton(Fenster, text="Save as txt", bg = '#EFEFEF', variable=var6).place(x = 25, y = 545)
    if FIT == 0:
        var9 = IntVar()
        Checkbutton(Fenster, text="Save Chi vs. Xi", bg = '#EFEFEF', variable=var9).place(x = 25, y = 570)
        if Statistical  == 0:
            var10 = IntVar()
            Checkbutton(Fenster, text="Save convolution functions", bg = '#EFEFEF', variable=var10).place(x = 25, y = 595)
        if Statistical  == 1:
            var12 = IntVar()
            Checkbutton(Fenster, text="Save Individuals", bg = '#EFEFEF',variable=var12).place(x = 25, y = 595)
            var14 = IntVar()
            Checkbutton(Fenster, text="Save Distribution function", bg = '#EFEFEF', variable=var14).place(x = 25, y = 620)
    if FIT == 1:
        if Statistical  == 0 :
            NoEntryLabel2 = Label(Fenster,text="No specified saving options provided.", bg = '#EFEFEF')
            NoEntryLabel2.place(x = 25, y = 570) 
        if Statistical  == 1:
            var12 = IntVar()
            Checkbutton(Fenster, text="Save Individuals", bg = '#EFEFEF',variable=var12).place(x = 25, y = 570)
            var14 = IntVar()
            Checkbutton(Fenster, text="Save Distribution function", bg = '#EFEFEF', variable=var14).place(x = 25, y = 595)
        
    
    colorbox11 = Label(Fenster,text="", bg = '#EFEFEF')
    colorbox11.place(x = 295, y = 540, width = 260, height = 120) 
    if FIT == 0:
        var8 = IntVar()
        Checkbutton(Fenster, text="Show Xi vs. Chi", bg = '#EFEFEF', variable=var8).place(x = 300, y = 545)
        if Statistical  == 1:
            var13 = IntVar()
            Checkbutton(Fenster, text="Show Distribution function", bg = '#EFEFEF', variable=var13).place(x = 300, y = 570)
            var11 = IntVar()
            Checkbutton(Fenster, text="Show Individuals", bg = '#EFEFEF', variable=var11).place(x = 300, y = 595)
    if FIT == 1:
        if Statistical  == 1:
            var13 = IntVar()
            Checkbutton(Fenster, text="Show Distribution function", bg = '#EFEFEF', variable=var13).place(x = 300, y = 545)
        if Statistical  == 0:
            NoEntryLabel3 = Label(Fenster,text="No specified displaying options provided.", bg = '#EFEFEF')
            NoEntryLabel3.place(x = 300, y = 545) 
        
        
    
    

    
    #-------------------------------------------------------
    #Akzeptierfunktion schreiben
    #-------------------------------------------------------
    
    def AcceptParams():
        global n
        global alpha
        global kzero
        global kp
        global kmp
        global kf
        global kmf
        global Ezero
        global T
        global distance
        global r
        global Scanrate
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
        global DeltaXiModder
        global ConvoluteRC
        global ROhm
        global CCPE
        global gammaCPE
        global kmax
        
        
#========================================================================
    #Read entries
#========================================================================
        
        n             = (float(n_Eingabe.get()))
        alpha         = (float(alpha_Eingabe.get()))
        kzero         = (float(kzero_Eingabe.get()))
        T             = (float(T_Eingabe.get())) + 273.15
        Scanrate      = (float(Scanrate_Eingabe.get()))*0.001
        FinKin        = var15.get()
        kmax          = 1e23
        kfin          = 0
        c             = (float(c_Eingabe.get()))*0.001
        A             = (float(A_Eingabe.get()))
        AstxtSaver    = var6.get()
        DeltaXiModder = var17.get()
        ConvoluteRC   = var00.get()
        if ConvoluteRC == 1:
            ROhm      = (float(ROhm_Eingabe.get()))
            CCPE      = (float(Cdl_Eingabe.get()))
            gammaCPE  = (float(gammaCPE_Eingabe.get()))
        
        if DeltaXiModder == 1:
            global ModDeltXi
            ModDeltXi = (float(ModDeltXi_Eingabe.get()))
            
        
        if FinKin == 1:
            kmax = float(kmax_Eingabe.get())
            kfin = kzero/kmax
        
        
        
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
        
        
        if FIT ==0:
            global EEE
            global SSS
            global SWT
            EEE         = (float(End_Eingabe.get()))
            SSS         = (float(Start_Eingabe.get()))
            SWT         = (float(Switch_Eingabe.get()))
            if SSS < SWT:
                Cathsweeper = 0
            if SSS > SWT:
                Cathsweeper = 1
                
            
            
            global Just_Current_Template
            Just_Current_Template = var8.get()
            global Chi_vs_Xi_Saver
            Chi_vs_Xi_Saver = var9.get()
            
            if Statistical  == 0:
                global Convolution_Saver
                Convolution_Saver = var10.get()
        
        
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
            

            
    button=Button(Fenster,text="Accept Parameters",bg = '#EFEFEF', command=AcceptParams).place(x = 570, y = 540, width = 260, height = 55)
    button=Button(Fenster,text="Next",bg = '#EFEFEF', command=CV_Type_Chooser).place(x = 570, y =605, width = 260, height = 55)
    


# In[ ]:

def CV_Type_Chooser():
    if model_3D == 1:
        Start_End_Definer()
        External_MT_Loader()
    if model_a == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_a_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_b == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_b_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_b1 == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_b1_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer() 
    if model_c == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_c_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_d == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_d_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_e == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_e_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_f == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_f_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_g == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_g_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
   


# In[ ]:

def Start_End_Definer():    
    if FIT == 1: 
        global EEE
        global SSS
        global SWT
        global Stromarrays
        global Potenzialarrays
        
        StromarraysROH             = np.squeeze(np.transpose(Strom[0::,0:1:]))   
        UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,0:1:]))
        PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,0:1:]))
        Stromarrays                = StromarraysROH
        UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::] 
        Potenzialarrays            = (PotenzialarraysROH[::]   -  Ru*Stromarrays/1000000) - Ezero   #referring all to E°
        
        if DesiReactxx ==0:
            Stromarrays        = -Stromarrays
            Potenzialarrays    = -Potenzialarrays
            SWT                = np.max(Potenzialarrays)    #Reduktion
            SSS                = Potenzialarrays[0]
            EEE                = Potenzialarrays[-1]  
            Stromarrays        = -Stromarrays
            Potenzialarrays    = -Potenzialarrays
        
        if DesiReactxx ==1:
            SWT                = np.max(Potenzialarrays)    #Oxidation
            SSS                = Potenzialarrays[0]
            EEE                = Potenzialarrays[-1] 


        if LinBaseCorr == 1:
            global StartPotLinCorr
            global EndPotLinCorr
            global BaseSlope
            if DesiReactxx ==0:
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

def DistFunc_Calcer(): 

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
    global P_Array_Unchanged
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
   

    DistFunc_User()
    
    if Show_Dist_Func == 1:
        Dist_Func_Shower()

    


# In[ ]:

def DistFunc_User():
    global distance 
    global a
    
    
    if FIT == 1:
        Start_End_Definer()
        
    
    distance = 0.01    #to run calculation with default values before real calculation starts
    a        = 0.01    #to run calculation with default values before real calculation starts
    
    #HIER MÜSSEN ERST DIE CONVOLUTION INVERTERS GEZÜNDET WERDEN um len(yyy) zu bekommen
    if model_b == 1:
        model_b_ConvolutionInverter()
    if model_d == 1:
        model_d_ConvolutionInverter()
    if model_e == 1:
        model_e_ConvolutionInverter()
    if model_g == 1:
        model_g_ConvolutionInverter()
    
    
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
            distance = 0.5*RR_Mean_Array[i]*0.0001  - r   #center to wall distances in µm

            
        INTERRUPT = 0
        
        if model_b == 1:
            if distance < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            model_b_ConvolutionInverter()
        if model_d == 1:
            if distance < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            model_d_ConvolutionInverter()
        if model_e == 1:
            if a < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            model_e_ConvolutionInverter()
        if model_g == 1:
            if a < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            model_g_ConvolutionInverter()
         
        for j in range(len(yyy)):
            if model_b == 1 or model_d == 1:
                ChiSuperarray[j,i] = yyy[j] *np.asscalar(Probab_Array[i])   #weightning by d-probability (planar fin and cylextfin)
            if model_e == 1 or model_g == 1:
                ChiSuperarray[j,i] = yyy[j] *INDIV_AREAS[i]/A       #weightning by area for internal cyl and internal sph

    
    global Chi_ges_Array
    Chi_ges_Array = np.empty(len(yyy))
    for i in range(len(yyy)):
        for j in range(len(Probab_Array)):
            Chi_ges_Array[i] = np.sum(ChiSuperarray[i,::])
            
    #====================================================================================================
    #Plotten von simuliertem
    #====================================================================================================
    if INTERRUPT == 0:
        global PotCalc     
        global CurrCalc 
        
        root = Toplevel()
        root.title("Simulated CV")
        root.geometry("500x500")

        if Xi_i > -6.9:
            Incorrect_Range_Warner()

        if FIT == 0:
            PotCalc  = R*T/(n*F)*xxx 
            CurrCalc = (n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000*ChiSuperarray

            Abbildung = Figure(figsize=(5, 4), tight_layout = True)
            b = Abbildung.add_subplot(111)

            if Just_Current_Template == 0:
                if Show_Individuals == 0:
                    b.plot(PotCalc[0:-1],Chi_ges_Array[0:-1]*(n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000, color = 'r')     
                if Show_Individuals == 1:
                    b.plot(PotCalc[0:-1],CurrCalc[0:-1,::])
                b.axvline(0,color='k')
                b.axhline(0,color='k')
                to = 0
                if np.abs(EEE) > np.abs(SWT):  
                    to = np.abs(EEE)
                if np.abs(EEE) <= np.abs(SWT):  
                    to = np.abs(SWT)
                b.set_xlim([-(to+0.05*to),(to+0.05*to)])    
                b.set_xlabel('E vs. E_zero. / V', fontsize=12)
                b.set_ylabel('I / mA', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b.spines[axis].set_linewidth(2)
                    b.spines[axis].set_color('k')
                b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)

            if Just_Current_Template == 1:
                if Show_Individuals == 0:
                    b.plot(xxx[0:-1],Chi_ges_Array[0:-1], color = 'r') 
                    if xxx[0]<xxx[1]:
                        b.annotate(u'\u03A8' '=%5.3f'  % max(Chi_ges_Array[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                    if xxx[0]>xxx[1]:
                        b.annotate(u'\u03A8' '=%5.3f'  % min(Chi_ges_Array[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                if Show_Individuals == 1:
                    b.plot(xxx[0:-1],ChiSuperarray[0:-1,::])
                b.axvline(0,color='k')
                b.axhline(0,color='k')
                to = 0
                if np.abs(EEE) > np.abs(SWT):  
                    to = n*F*np.abs(EEE)/(R*T)
                if np.abs(EEE) <= np.abs(SWT):  
                    to = n*F*np.abs(SWT)/(R*T)
                b.set_xlim([-(to+0.1*to),(to+0.1*to)])    
                b.set_xlabel(u'\u03BE', fontsize=12)
                b.set_ylabel(u'\u03C7', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b.spines[axis].set_linewidth(2)
                    b.spines[axis].set_color('k')
                b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)

            canvas = FigureCanvasTkAgg(Abbildung, master=root)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas, root)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)


        if FIT == 1:
            PotCorr  = R*T/(n*F) 
            CurrCorr = (n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000000
            global Pot
            global Curr
            Pot              = Potenzialarrays
            Curr             = Stromarrays
            PotCalc          = xxx*PotCorr
            Chi_ges_Array    = CurrCorr*Chi_ges_Array       #Unterschied zu Vorher :)
            CurrCalc         = 0.001*Chi_ges_Array
            NormLengthPot    = (1.0/float(len(Pot)))
            NormLengthxxx    = (1.0/float(len(PotCalc)))
            NormArrayPot     = np.arange(0,1,NormLengthPot)
            NormArrayxxx     = np.arange(0,1,NormLengthxxx)
            Fitinterpolation = InterpolatedUnivariateSpline(NormArrayxxx,Chi_ges_Array,k=3)
            InterpolatedCurr = np.empty(len(NormArrayPot))
            for i in range(len(NormArrayPot)):
                InterpolatedCurr[i] = Fitinterpolation(NormArrayPot[i])
            Abweichungsarray = np.empty(len(Curr)-1)
            for i in range(len(Curr)-1):
                Abweichungsarray[i] = (((Curr[i] - InterpolatedCurr[i])/np.max(np.abs(Curr)))**2)/(float(len(Curr))-1)
            Varianz            = np.sum(Abweichungsarray)
            global Standardabweichung
            Standardabweichung = Varianz**0.5
            Abbildung = Figure(figsize=(5, 4))
            b = Abbildung.add_subplot(111)
            b.plot(PotCalc[1:-1:5],0.001*Chi_ges_Array[1:-1:5],color='r', marker = '.', linestyle = '')
            b.plot(PotCalc[:-1:5],0.001*Chi_ges_Array[:-1:5],color='r')
            b.plot(Pot[:-1:],0.001*Curr[:-1:],color='k')       
            b.annotate(u'\u03C3' '=%5.3f'  % Standardabweichung , xy=(0.05, 0.83), xycoords='axes fraction')
            b.axvline(0,color='k')
            b.axhline(0,color='k')
            b.set_xlabel('E vs. E_zero. / V', fontsize=12)
            b.set_ylabel('I / mA', fontsize=12)
            for axis in ['top','bottom','left','right']:
                b.spines[axis].set_linewidth(2)
                b.spines[axis].set_color('k')
            b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            canvas = FigureCanvasTkAgg(Abbildung, master=root)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas, root)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)



        if AstxtSaver == 1:
            Simulated_CV_as_txt_saver()

    


# In[ ]:

def Dist_Func_Shower():   
    root = Toplevel()
    root.title("Distance distribution-function")
    root.geometry("500x500")
    
    Abbildung = Figure(figsize=(5, 4), tight_layout=True)
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
#Invert the functions of Interest or load MT func externally
#=================================================================================================================================
#=================================================================================================================================
#=================================================================================================================================
def External_MT_Loader():
    Path_MT       = askopenfilename()
    data_MT       = np.loadtxt(Path_MT, delimiter = '\t')
    MT_Time      = data_MT[::,0]
    MT3D         = data_MT[::,1]
    global CV_InterpolationHin
    CV_InterpolationHin = InterpolatedUnivariateSpline(MT_Time,MT3D, k=3)
    global CV_InterpolationBack
    CV_InterpolationBack = InterpolatedUnivariateSpline(MT_Time,MT3D, k=3)
    CV_Calculator()



def model_a_ConvolutionInverter():
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
        
        CV_Calculator()
    

    
def model_b_ConvolutionInverter():
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
        CV_Calculator()
    

def model_b1_ConvolutionInverter():
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
        CV_Calculator()
    
    
def model_c_ConvolutionInverter():
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
        CV_Calculator()
    

    
def model_d_ConvolutionInverter():
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
        CV_Calculator()
    
    
    
def model_e_ConvolutionInverter():
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
        CV_Calculator()
    


def model_f_ConvolutionInverter():
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
        CV_Calculator()
    
    

def model_g_ConvolutionInverter():
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
        CV_Calculator()
    
 
    
 



# In[ ]:

def CV_Calculator():
    
    if Statistical == 0:
        root = Toplevel()
        root.title("Simulated CV")
        root.geometry("500x500")
    
    global Kp
    global Kf
    global p
    global f
    p               = 0.0
    Kp              = 1e9
    Kf              = 1e-9
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
        ScaRa           = Scanrate
        delta_Xi        = 0.01
        if DeltaXiModder == 1:
            delta_Xi = ModDeltXi
        E_switch        = SWT
        delta_E         = delta_Xi*R*T/(n*F)
        Num_of_E_forw   = int((np.abs(SSS)+np.abs(SWT))/delta_E)
        E_i             = np.abs(SWT) - (Num_of_E_forw+1)*delta_E
        sigma           = n*F*ScaRa/(R*T)
        Xi_i            = -n*F*np.abs(E_i)/(R*T)   #initial Xi
        Xi_f            = n*F*(np.abs(SWT))/(R*T)  #final Xi
        if Cathsweeper == 0:
            Xi_Asym     = n*F*(EEE)/(R*T)
        if Cathsweeper == 1:
            Xi_Asym     = n*F*(-EEE)/(R*T)
        Xi_Array_hin    = np.arange(Xi_i,Xi_f,delta_Xi) 
        Xi_Array_back   = np.arange(Xi_f,Xi_Asym-delta_Xi,-delta_Xi) 
        Xi_Array        = np.concatenate([Xi_Array_hin,Xi_Array_back])
        delta_t         = delta_Xi*R*T/(n*F*ScaRa)
        t_Array_cont    = np.zeros(len(Xi_Array))
        for i in range(len(Xi_Array )):
            t_Array_cont[i] =  i*delta_t
        Exp_Array       = np.exp(-Xi_Array)
        Preced_Array    = np.exp(-p*t_Array_cont)
        Follow_Array    = np.exp(-f*t_Array_cont)
        Kinetik_Array   = D_f**0.5/(kzero*Exp_Array**(-alpha))
        
    if FIT == 1:
        ScaRa           = Scanrate
        delta_Xi        = 0.01
        if DeltaXiModder == 1:
            delta_Xi = ModDeltXi
            
        delta_E         = delta_Xi*R*T/(n*F)
        Num_of_E_forw   = int((SWT-SSS)/delta_E)
        Num_of_E_backw  = int((SWT-EEE)/delta_E)
       
        
        
        sigma           = n*F*ScaRa/(R*T)
        Xi_i            = n*F*SSS/(R*T)                   #initial Xi
        Xi_Asym         = n*F*(SWT)/(R*T)  #Asymmetrisch endendes Xi
        Xi_f            = n*F*(EEE)/(R*T)                 #final Xi
        Xi_Array_hin    = np.arange(Xi_Asym, Xi_i, -delta_Xi) 
        Xi_Array_hin    = Xi_Array_hin[::-1]
        Xi_Array_back   = np.arange(Xi_Asym+delta_Xi, Xi_f-delta_Xi,-delta_Xi) 
        Xi_Array        = np.concatenate([Xi_Array_hin,Xi_Array_back])
        
        delta_t         = delta_Xi*R*T/(n*F*ScaRa)
        t_Array_cont    = np.zeros(len(Xi_Array))
        for i in range(len(Xi_Array )):
            t_Array_cont[i] =  i*delta_t
        Exp_Array       = np.exp(-Xi_Array)
        Preced_Array    = np.exp(-p*t_Array_cont)
        Follow_Array    = np.exp(-f*t_Array_cont)
        Kinetik_Array   = D_f**0.5/(kzero*Exp_Array**(-alpha))
        
    
    Fin_Kin_p_Array = 1/(1 + kfin*Exp_Array**(-alpha))      #hier mit minus, weil Exp-Array schon umgedreht ist!
    Fin_Kin_s_Array = 1/(1 + kfin*Exp_Array**((1-alpha)))
     
    
    #(kmax/(kmax + (kzero*Kn[m]**(alpha))*np.exp(-alpha*n*F*R_ohm*x/(R*T))))

    #-------------------------------------------------------------------------------------------------- 
    #Jetzt CV berechnen
    #--------------------------------------------------------------------------------------------------
    
    global WuFuInt_p
    WuFuInt_p        = np.empty(len(Xi_Array))
    Delta_WuFuInt_p  = np.empty(len(Xi_Array))
    for i in range(len(t_Array_cont)-1):
        WuFuInt_p[i]        = CV_InterpolationHin(t_Array_cont[i])
    for i in range(len(t_Array_cont)-2):
        Delta_WuFuInt_p[i]  = WuFuInt_p[i+1] - WuFuInt_p[i]
        
    global WuFuInt_f
    WuFuInt_f        = np.empty(len(Xi_Array))
    Delta_WuFuInt_f  = np.empty(len(Xi_Array))

    for i in range(len(t_Array_cont)-1):
        WuFuInt_f[i]       = CV_InterpolationBack(t_Array_cont[i])
    for i in range(len(t_Array_cont)-2):
        Delta_WuFuInt_f[i]  = WuFuInt_f[i+1] - WuFuInt_f[i]
    
    
    ############################################################################################
    # This function is used, no Ohmic resistances and double layer capacities are involved
    ############################################################################################
    if ConvoluteRC == 0:
        Chi_Array  = np.zeros(len(Xi_Array))
        
        for i in range(len(Xi_Array)-1):
        
            #MUCH MUCH MUCH FASTER!!!!!!!!!!!!!!!!!!!!! 
            Summandenarray_p_GG   = Delta_WuFuInt_p[i::-1]*Chi_Array[:i+1:]*Preced_Array[i::-1]
            Summandenarray_f_GG   = Delta_WuFuInt_f[i::-1]*Chi_Array[:i+1:]*Follow_Array[i::-1]
            Summandenarray_p      = Delta_WuFuInt_p[i::-1]*Chi_Array[:i+1:]
            Summandenarray_f      = Delta_WuFuInt_f[i::-1]*Chi_Array[:i+1:]
            Chi_Array[i]          = (Fin_Kin_p_Array[i]*(1/(1+Kp))*(Kp - Kp*np.sum(Summandenarray_p) - np.sum(Summandenarray_p_GG)) -(Fin_Kin_s_Array[i]*(D_f/D_b)**0.5)*(Exp_Array[i]/(1+Kf))*(np.sum(Summandenarray_f) + Kf*np.sum(Summandenarray_f_GG))  )/(  Kinetik_Array[i] +  Fin_Kin_p_Array[i]*WuFuInt_p[1]+ (Fin_Kin_s_Array[i]*(D_f/D_b)**0.5)*WuFuInt_f[1]*Exp_Array[i])
    
    ############################################################################################
    # This function is used, when R and C should be convoluted to give I
    ############################################################################################
    if ConvoluteRC == 1:
        def ContourIntegral_ConvRC(t, N, LapFunc):
            h = 2*np.pi/N;
            shift = 0.0;
            ans   = 0.0;
            if t == 0:
                print ("ERROR:   Inverse transform can not be calculated for t=0")
                return ("Error");
            for k in range(0,N):
                theta  = -np.pi + (k+1./2)*h;
                z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
                dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
                ans    = ans + exp(z*t)*(LapFunc(z))*dz; 
            return ((h/(2j*np.pi))*ans).real   
        
        def Talbot_ConvRC(LapFunc):
            t_Array             = np.logspace(-6,6,400)    #np.logspace(-8,5,300)
            ft_Array            = np.empty(len(t_Array)) 
            for i in range(len(t_Array)):
                ft_Array[i]     = ContourIntegral_ConvRC(float(t_Array[i]), 24, LapFunc) 
            TimeDomainFunc      = InterpolatedUnivariateSpline(t_Array, ft_Array, k=3)
            return TimeDomainFunc
        
        #=================================================================================
        #=================================================================================
        #Adjust Parameters for CV Simulation for CPE-CV function
        #=================================================================================
        #=================================================================================
        nu    =  Scanrate           ;    E_0   =  0           ;    E_Span = (SWT - SSS) + (SWT - EEE)
        R_ohm =  ROhm               ;    K_dl  =  CCPE        ;    gamma  = gammaCPE    
        #=================================================================================
        # derived units
        #=================================================================================
        Tau       = R_ohm*K_dl
        E_Span    = (SWT - SSS) + (SWT - EEE) 
        t         = t_Array_cont[1::]
        steps1    = int(round(len(t)*(SWT - SSS)/E_Span))
        steps2    = int(round(len(t)*(SWT - EEE)/E_Span))

        E         = (R*T/(n*F))*Xi_Array
        Kn        = np.exp(n*F*(E-E_0)/(R*T))
        timescale = np.logspace(-20,6,10000)   
        Preced    = np.exp(-p*t_Array_cont[1::])
        Follow    = np.exp(-f*t_Array_cont[1::])
        #=================================================================================
        # Laplace-domain solution of the CPE-current funtion
        #=================================================================================
        def Term_One_Lap(s):
            return nu*K_dl/((s**(-gamma) + Tau)*s**2)
        #=================================================================================
        # prepare the lists for the pure CPE-current without Faradaic current
        #=================================================================================
        Pure_Idl_t1                  =    np.zeros(len(t_Array_cont[1::]))
        Pure_Idl_t2                  =    np.zeros(len(t_Array_cont[1::]))
        Pure_Idl_t1                  =    Talbot_ConvRC(Term_One_Lap)(t_Array_cont[1::])
        Pure_Idl_t2[steps1::]        = -2*Talbot_ConvRC(Term_One_Lap)(t_Array_cont[steps1:-1:] - t_Array_cont[steps1]) 
        Pure_Idl_t                   =    Pure_Idl_t1 + Pure_Idl_t2 
        def h_Lap(s):
            return 1/(s*(1 + Tau*s**(gamma)))
        H_timescale     = Talbot_ConvRC(h_Lap)(timescale)
        H_interpol      = InterpolatedUnivariateSpline(timescale, H_timescale, k = 3)
        #================================================================================
        # Mass transfer function
        #================================================================================
              
        M_t_p  = WuFuInt_p
        M_t_f  = WuFuInt_f
        
        
        H_t   = H_interpol(t_Array_cont[1::])
        #=================================================================================
        #=================================================================================
        #Solve Convolution Integrals
        #=================================================================================
        #=================================================================================
        I_f    = np.zeros(len(E))
        I_ges  = np.zeros(len(E))
        I_dl   = np.zeros(len(E))
        Iconvm = 0
        for m in range(len(E)-1):
            if m > 0: 
                #Vectorized for loop for acceleration the convolution of lists-finally works :)
                IconvF_p       = np.sum(I_f[0:m-1:]*(M_t_p[m:1:-1]  - M_t_p[m-1:0:-1]))
                IconvF_p_kin   = np.sum(I_f[0:m-1:]*(M_t_p[m:1:-1]  - M_t_p[m-1:0:-1])*Preced[m-1:0:-1])
                IconvF_f       = np.sum(I_f[0:m-1:]*(M_t_f[m:1:-1]  - M_t_f[m-1:0:-1]))
                IconvF_f_kin   = np.sum(I_f[0:m-1:]*(M_t_f[m:1:-1]  - M_t_f[m-1:0:-1])*Follow[m-1:0:-1])
                if m > 1:
                    Iconvm    = np.sum(H_t[0:m-2:]*(I_f[m-1:1:-1] - I_f[m-2:0:-1]))
                if m == 1:      #Function will be defined only once INSIDE the for loop
                    def Rekur(x): 
                        fan = (kmax/(kmax + (kzero*Kn[m]**(alpha))*np.exp(-alpha*n*F*R_ohm*x/(R*T))))
                        fca = (kmax/(kmax + kzero*Kn[m]**(-(1-alpha))*np.exp((1-alpha)*n*F*R_ohm*x/(R*T))))
                        
                        return (  (  (  (  fan*(1/(1+Kp))*(n*F*A*c*Kp*D_f**0.5 - Kp*IconvF_p -  IconvF_p_kin ) 
                                    - (fca)*((1/(1+Kf))*(1/(Kn[m]*np.exp(-n*F*R_ohm*x/(R*T)))))*((D_f/D_b)**0.5)*(IconvF_f + Kf*IconvF_f_kin )  )/   
                                    ((D_f**0.5/(kzero*(Kn[m]**alpha)*np.exp(-alpha*n*F*R_ohm*x/(R*T))))
                                    + M_t_p[1]*fan+ M_t_f[1]*fca*((D_f/D_b)**0.5)*(1/(Kn[m]*np.exp(-n*F*R_ohm*x/(R*T))))   )   ) 
                                    - I_f[m-1])*H_t[1] + Pure_Idl_t[m] + Iconvm - x  )
                I_ges[m]  = fsolve(Rekur, I_ges[m-1])
                I_f[m]    = (I_ges[m] - Pure_Idl_t[m] - Iconvm)/H_t[1] + I_f[m-1]
                I_dl[m]   = I_ges[m] - I_f[m]
        
        
        Chi_Array = I_ges/(n*F*A*c*D_f**0.5)       # This is essentially NOT chi, but a pre version. It is I/(nFAc(D**0.5)). Chis is, however,
                                                 # I/(nFAc(D**0.5)\sigma) with \sigma = (nF\nu/RT)**0.5). This is done in the lines below
                                                 # So basicall yyy is the true Chi
    
    
    #-----------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------
    global xxx
    global yyy
    xxx = (Xi_Array)
    yyy = (Chi_Array[::])/sigma**0.5
    
    if len(xxx) > len(yyy):
        deltaLen = int(len(xxx)-len(yyy))
        xxx      = xxx[deltaLen::]
    if len(yyy) > len(xxx):
        deltaLen = int(len(yyy)-len(xxx))
        yyy      = yyy[deltaLen::]
     
    
    
    if FIT == 0:
        if Cathsweeper == 1:
            xxx = -(Xi_Array)
            yyy = -yyy
    if FIT ==1:
        if DesiReactxx == 0:
            xxx = -(Xi_Array)
            yyy = -yyy
            
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON , nur wenn nicht statistisch--> der hat seinen eigenen Plotter
    #-------------------------------------------------------------------------------------------------------------------------
    
    if Statistical == 0:
        global PotCalc
        global CurrCalc
        
        if Xi_i > -6.9:
            Incorrect_Range_Warner()
            
        if FIT == 0:
            PotCalc  = R*T/(n*F)*xxx 
            CurrCalc = (n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000*yyy
            Abbildung = Figure(figsize=(5, 4), tight_layout=True)
            b = Abbildung.add_subplot(111)

            if Just_Current_Template == 0:
                b.plot(PotCalc[:-1:],CurrCalc[:-1:],color='r')
                if Cathsweeper == 0:
                    b.annotate('I_p / mA' '=%5.3f'  % max(CurrCalc[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                if Cathsweeper == 1:
                    b.annotate('I_p / mA' '=%5.3f'  % min(CurrCalc[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                b.axvline(0,color='k')
                b.axhline(0,color='k')
                to = 0
                if np.abs(EEE) > np.abs(SWT):
                    to = np.abs(EEE)
                if np.abs(EEE) <= np.abs(SWT):
                    to = np.abs(SWT)
                    
                b.set_xlim([-(to+0.05*to),(to+0.05*to)])    
                b.set_xlabel('E vs. E_zero. / V', fontsize=12)
                b.set_ylabel('I / mA', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b.spines[axis].set_linewidth(2)
                    b.spines[axis].set_color('k')
                b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)

            if Just_Current_Template == 1:
                b.plot(xxx[:-1:],yyy[:-1:],color='r')
                if Cathsweeper == 0:
                    b.annotate(u'\u03A8' '=%5.3f'  % max(yyy[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                if Cathsweeper == 1:
                    b.annotate(u'\u03A8' '=%5.3f'  % min(yyy[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                b.axvline(0,color='k')
                b.axhline(0,color='k')
                to = 0
                if np.abs(EEE) > np.abs(SWT):
                    to = n*F*np.abs(EEE)/(R*T)
                if np.abs(EEE) <= np.abs(SWT):
                    to = n*F*np.abs(SWT)/(R*T)
                b.set_xlim([-(to+0.1*to),(to+0.1*to)])    
                b.set_xlabel(u'\u03BE', fontsize=12)
                b.set_ylabel(u'\u03C7', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b.spines[axis].set_linewidth(2)
                    b.spines[axis].set_color('k')
                b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            canvas = FigureCanvasTkAgg(Abbildung, master=root)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas, root)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)


        if FIT == 1:
            global Pot
            global Curr
            PotCorr          = R*T/(n*F) 
            CurrCorr         = (n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000000
            Pot              = Potenzialarrays
            Curr             = Stromarrays
            xxx              = xxx*PotCorr
            yyy              = yyy*CurrCorr
            PotCalc          = xxx
            CurrCalc         = 0.001*yyy
            NormLengthPot    = (1.0/len(Pot))
            NormLengthxxx    = (1.0/len(xxx))
            NormArrayPot     = np.arange(0,1,NormLengthPot)
            NormArrayxxx     = np.arange(0,1,NormLengthxxx)
            Fitinterpolation = InterpolatedUnivariateSpline(NormArrayxxx,yyy,k=3)
            InterpolatedCurr = np.empty(len(NormArrayPot))
            for i in range(len(NormArrayPot)):
                InterpolatedCurr[i] = Fitinterpolation(NormArrayPot[i])
            Abweichungsarray = np.empty(len(Curr)-1)
            for i in range(len(Curr)-1):
                Abweichungsarray[i] = (((Curr[i] - InterpolatedCurr[i])/np.max(np.abs(Curr)))**2)/(float(len(Curr))-1)
            Varianz            = np.sum(Abweichungsarray)
            global Standardabweichung
            Standardabweichung = Varianz**0.5
            Abbildung = Figure(figsize=(5, 4))
            b = Abbildung.add_subplot(111)
            b.plot(xxx[1:-1:5],0.001*yyy[1:-1:5],color='r', marker ='.', linestyle = '')   
            b.plot(xxx[:-1:],0.001*yyy[:-1:],color='r', label='Fit')   
            b.plot(Pot[:-1:],0.001*Curr[:-1:],color='k', label='Data')  
            b.annotate(u'\u03C3' '=%5.3f'  % Standardabweichung , xy=(0.05, 0.83), xycoords='axes fraction')
            b.axvline(0,color='k')
            b.axhline(0,color='k')
            b.set_xlabel('E vs. E_zero. / V', fontsize=12)
            b.set_ylabel('I / mA', fontsize=12)
            b.legend(frameon = False, fontsize = 12)
            for axis in ['top','bottom','left','right']:
                b.spines[axis].set_linewidth(2)
                b.spines[axis].set_color('k')
            b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            canvas = FigureCanvasTkAgg(Abbildung, master=root)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas, root)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        if AstxtSaver == 1:
            Simulated_CV_as_txt_saver()



# In[ ]:

def Simulated_CV_as_txt_saver():
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
        fi.write("Scanrate in [V/s]")
        fi.write("\t")
        fi.write(str(Scanrate))
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
        
        if ConvoluteRC == 1:
            fi.write("RC got convoluted")
            fi.write("\n")
            fi.write("R in Ohm")
            fi.write("\t")
            fi.write(str(ROhm))
            fi.write("\n")
            fi.write("Cdl in F*s^(1-gamma)")
            fi.write("\t")
            fi.write(str(CCPE))
            fi.write("\n")
            fi.write("gamma CPE")
            fi.write("\t")
            fi.write(str(gammaCPE))
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
        
        if DeltaXiModder == 1:
            fi.write("Xi-Resol. got modified")
            fi.write("\t")
            fi.write(str(ModDeltXi))
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
            
                
        
            
        if FIT == 1:
            
            fi.write("E_zero vs. E_ref in [V]")
            fi.write("\t")
            fi.write(str(Ezero))
            fi.write("\n")
            fi.write("R_u in [Ohm]")
            fi.write("\t")
            fi.write(str(Ru))
            fi.write("\n")
            fi.write("Standard deviation")
            fi.write("\t")
            fi.write(str(Standardabweichung))
            fi.write("\n")
            if LinBaseCorr == 1:
                fi.write("Base was corrected")
                fi.write("\n")
                fi.write("Slope Start Pot vs. E_zero [V]")
                fi.write("\t")
                fi.write(str(StartPotLinCorr))
                fi.write("\n")
                fi.write("Slope End Pot vs. E_zero [V]")
                fi.write("\t")
                fi.write(str(EndPotLinCorr))
                fi.write("\n")
                fi.write("Slope of Baseline [micA/V]")
                fi.write("\t")
                fi.write(str(BaseSlope))
                fi.write("\n")
            
        
        for i in range(3):
            fi.write("\n")
        
        fi.write("CURVE DATA")
        
        
#---------------------------------------------------------------------------------------------------------------------------        
#Wenn der Fit null ist... also wenn nur berechnet wird...  

        if FIT == 0:  
            fi.write("\n")
            fi.write("E_vs_E_zero [V]")
            fi.write("\t")
            fi.write("I[mA]")
            fi.write("\t")
            if Chi_vs_Xi_Saver == 1:
                fi.write("Xi")
                fi.write("\t")
                fi.write("Chi")
                fi.write("\t")
            if Statistical  == 0:
                if Convolution_Saver == 1:
                    fi.write("t")
                    fi.write("\t")
                    fi.write("Convolution_integral_forward")
                    fi.write("\t")
                    fi.write("Convolution_integral_backward")
                    fi.write("\t")              
            fi.write("\n")

        if FIT == 0:
            for i in range(len(xxx)-1):
                fi.write(str(np.asscalar(PotCalc[i])))
                fi.write("\t")
                if Statistical  == 0:
                    fi.write(str(np.asscalar(CurrCalc[i])))
                    fi.write("\t")
                if Statistical  == 1:
                    fi.write(str(np.asscalar(Chi_ges_Array[i]*(n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000)))
                    fi.write("\t")
                if Statistical  == 0:
                    if Chi_vs_Xi_Saver == 1:
                        fi.write(str(np.asscalar(xxx[i])))
                        fi.write("\t")
                        fi.write(str(np.asscalar(yyy[i])))
                        fi.write("\t")
                    if Convolution_Saver == 1:
                        fi.write(str(np.asscalar(t_Array_cont[i])))
                        fi.write("\t")
                        fi.write(str(np.asscalar(WuFuInt_p[i])))
                        fi.write("\t")
                        fi.write(str(np.asscalar(WuFuInt_f[i])))
                        fi.write("\t")
                if Statistical  == 1:
                    if Chi_vs_Xi_Saver == 1:
                        fi.write(str(np.asscalar(xxx[i])))
                        fi.write("\t")
                        fi.write(str(np.asscalar(Chi_ges_Array[i])))
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
                    fi.write("E_vs_E_zero [V]")
                    fi.write("\t")
                    for i in range(len(R_Mean_Array)):
                        fi.write("weighted I[mA]")
                        fi.write("\t")
                    fi.write("\n")
                    for i in range(len(PotCalc)-1):
                        fi.write(str(np.asscalar(PotCalc[i])))
                        fi.write("\t")
                        for j in range(len(R_Mean_Array)):
                            fi.write(str(np.asscalar(CurrCalc[i,j])))
                            fi.write("\t")       
                        fi.write("\n")
                    fi.write("\n")
                    fi.write("\n")
                
                    if Chi_vs_Xi_Saver == 1: 
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
                        fi.write("Xi")
                        fi.write("\t")
                        for i in range(len(R_Mean_Array)):
                            fi.write("weighted Chi")
                            fi.write("\t")
                        fi.write("\n")
                        for i in range(len(PotCalc)-1):
                            fi.write(str(np.asscalar(xxx[i])))
                            fi.write("\t")
                            for j in range(len(R_Mean_Array)):
                                fi.write(str(np.asscalar(CurrCalc[i,j]/((n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000))))
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

#---------------------------------------------------------------------------------------------------------------------------        
#Wenn der Fit gleich eins ist, also gefittet wird
#---------------------------------------------------------------------------------------------------------------------------        

        if FIT == 1:
            fi.write("\n")
            fi.write("Calculated Potential vs. E_zero in [V]")
            fi.write("\t")
            fi.write("Calculated Current in microampere")
            fi.write("\t")
            fi.write("\n")
            
            if Statistical  == 0:
                for i in range(len(xxx)-1):
                    fi.write(str(np.asscalar(xxx[i])))
                    fi.write("\t")
                    fi.write(str(np.asscalar(yyy[i])))
                    fi.write("\t")
                    fi.write("\n")
                fi.write("\n")
                fi.write("\n")
                fi.write("\n")
                
            if Statistical  == 1:
                for i in range(len(xxx)-1):
                    fi.write(str(np.asscalar(PotCalc[i])))
                    fi.write("\t")
                    fi.write(str(np.asscalar(Chi_ges_Array[i])))
                    fi.write("\t")
                    fi.write("\n")
                fi.write("\n")
                fi.write("\n")
                fi.write("\n")
                
                
            if Statistical  == 1:
                if Save_Individuals == 1:
                    fi.write("r to wall in 10^-6[m]")
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
                    fi.write("E_vs_E_zero [V]")
                    fi.write("\t")
                    for i in range(len(R_Mean_Array)):
                        fi.write("weighted I[mA]")
                        fi.write("\t")
                    fi.write("\n")
                    for i in range(len(PotCalc)-1):
                        fi.write(str(np.asscalar(PotCalc[i])))
                        fi.write("\t")
                        for j in range(len(R_Mean_Array)): #soll in milliamps sein, deshalb mal 1000
                            fi.write(str(np.asscalar(1000*ChiSuperarray[i,j]*(n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5))))
                            fi.write("\t")       
                        fi.write("\n")
                    fi.write("\n")
                    fi.write("\n")
                
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


#---------------------------------------------------------------------------------------------------------------------------        
     
            fi.write("R_Corrected measured Potential vs. E_zero in [V]")
            fi.write("\t")
            fi.write("Measured Current in microampere")
            fi.write("\t")
            fi.write("\n")
            
            for i in range(len(Pot)):
                fi.write(str(np.asscalar(Pot[i]-Ezero)))
                fi.write("\t")
                fi.write(str(np.asscalar(Curr[i])))
                fi.write("\t")
                fi.write("\n")            
                        
                
    root.destroy()


