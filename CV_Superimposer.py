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
from scipy.fft                          import fft, ifft
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



F = 96485.0
R = 8.314

def cot(phi):
    return 1.0/tan(phi)

def csc(phi):
    return 1.0/sin(phi)

def coth(x):
    return 1/tanh(x)


#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================
#Re-import written modules for Cyclic_Voltammetry

from Cyclic_Voltammetry import *
#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================

def Wrong_Number_of_CV_Warner():
    Fenster = Toplevel()                                                         
    Fenster.title("Wrong Number of CVs")                         
    Fenster.geometry("400x200")
    colorbgr = Label(Fenster, text= "", bg = '#FF0000')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000) 
    Te = Text(Fenster, height=3, width=30)
    Te.place(x = 75, y = 75)
    Te.insert(END, "Caution! Wrong number of\nCVs to superimpose. Number has\nto be min. 2 and max. 8.")
    

def No_Mode_Warner():
    Fenster = Toplevel()                                                         
    Fenster.title("Wrong Number of CVs")                         
    Fenster.geometry("400x200")
    colorbgr = Label(Fenster, text= "", bg = '#FF0000')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000) 
    Te = Text(Fenster, height=3, width=30)
    Te.place(x = 75, y = 75)
    Te.insert(END, "Caution! you have not selected a\nsuperposition mode. \nCheck RS-Mode OR Add-Mode")
    


#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================


def CVSIM_Superimposer():
    def Runner():
        global choose
        Add_Mode  = var1.get()
        RS_Mode   = var2.get()
        choose    = int(En1.get())
        if choose == 1:
            Semi_Inf_Planar()
        if choose == 2:
            Finit_Planar()
        if choose == 3:
            Finit_Planar_Trans()
        if choose == 4:
            Semi_Inf_Zyl_Ext()
        if choose == 5:
            Finit_Zyl_Ext()
        if choose == 6:
            Finit_Zyl_Int()
        if choose == 7:
            Semi_Inf_Sphere_Ext()
        if choose == 8:
            Finit_Sphere_Int()
        if choose == 9:
            Statistical_Finit_Planar()
        if choose == 10:
            Statistical_Finit_Zyl_Ext()
        if choose == 11:
            Statistical_Finit_Zyl_Int()
        if choose == 12:
            Statistical_Finit_Sphere_Int()
    Fenster = Toplevel()                                                         
    Fenster.title("Superposition mode") 
    Fenster.geometry("650x400")
    colorbgr = Label(Fenster, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000)  
    CVType_Label  = Label(Fenster,text="Superimpose up to 50 CVs by choosing the CV-Type number 1 to 12", bg = '#80ffbf', font = ('Arial',12,'bold'))
    CVType_Label.place(x = 25, y = 15)    
    Lab1 = Label(Fenster,text=" CV-Type Number", bg = '#80ffbf')
    Lab1.place(x = 25, y = 50 )
    En1  = Entry(Fenster)
    En1.place(x =200, y = 50 , width =90, height = 20)
    button0 =Button (Fenster,text="run",bg = '#EFEFEF', command=Runner).place(x = 200, y = 85, width = 90, height = 50)
    var1 = IntVar()
    Checkbutton(Fenster, text="Additive mode ",bg = '#80ffbf', variable=var1).place(x = 25, y = 85)
    var2 = IntVar()
    Checkbutton(Fenster, text="Randles-Sevcik mode",bg = '#80ffbf', variable=var2).place(x = 25, y = 110)
    Te = Text(Fenster, height=4, width=25, bg = '#80ffbf')
    Te.place(x = 350, y = 50, width = 250, height = 85)
    Te.insert(END, "Additive mode will add\nCVs to one single CV. \nRandles-Sevcik mode\nwill overlay them.")
    SepLine1 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine1.place(x = 0, y = 150, height = 20)
    Label_1 = Label(Fenster,text="1)   Semi-inf. planar", bg = '#80ffbf')
    Label_1.place(x = 25, y = 175 )
    Label_2 = Label(Fenster,text="2)   Fin. reflective planar", bg = '#80ffbf')
    Label_2.place(x = 25, y = 200 )
    Label_3 = Label(Fenster,text="3)   Fin. transmissive planar", bg = '#80ffbf')
    Label_3.place(x = 25, y = 225 )
    Label_4 = Label(Fenster,text="4)   Semi-inf. ext. cylindrical.", bg = '#80ffbf')
    Label_4.place(x = 25, y = 250 )
    Label_5 = Label(Fenster,text="5)   Fin. reflective ext. cylindrical", bg = '#80ffbf')
    Label_5.place(x = 200, y = 175 )
    Label_6 = Label(Fenster,text="6)   Fin. reflective int. cylindrical", bg = '#80ffbf')
    Label_6.place(x = 200, y = 200 )
    Label_7 = Label(Fenster,text="7)   Semi-inf. ext. spherical", bg = '#80ffbf')
    Label_7.place(x = 200, y = 225 )
    Label_8 = Label(Fenster,text="8)   Fin. reflective int. spherical", bg = '#80ffbf')
    Label_8.place(x = 200, y = 250 )
    Label_9 = Label(Fenster,text=" 9)   Statistical fin. reflective planar", bg = '#80ffbf')
    Label_9.place(x = 400, y = 175 )
    Label_10 = Label(Fenster,text="10)   Statistical fin. reflective ext. cylindrical", bg = '#80ffbf')
    Label_10.place(x = 400, y = 200 )
    Label_11 = Label(Fenster,text="11)   Statistical fin. reflective int. cylindrical", bg = '#80ffbf')
    Label_11.place(x = 400, y = 225 )
    Label_12 = Label(Fenster,text="12)   Statistical fin. reflective int. spherical", bg = '#80ffbf')
    Label_12.place(x = 400, y = 250 )
    SepLine2 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine2.place(x = 0, y = 275, height = 20) 
    def Displayer():
        if RS_Mode == 1 and Add_Mode == 1:
            Double_Mode_Caller()
        if RS_Mode == 0 and Add_Mode == 0:
            No_Mode_Warner()
        if RS_Mode == 1 and Add_Mode == 0:
            RS_Mode_Caller()
        if Add_Mode == 1 and RS_Mode == 0:
            Add_Mode_Caller()
    def Mode_Decider():
        global RS_Mode
        global Add_Mode
        Add_Mode  = var1.get()
        RS_Mode   = var2.get()
        Displayer()
    button1=Button(Fenster,text="Begin/Clear",bg = '#EFEFEF',        command=Initializer    ).place(x = 25, y =300 , width = 100, height = 85)
    button2=Button(Fenster,text="Append latest",bg = '#EFEFEF',      command=Append_Latest  ).place(x = 150, y =300 , width = 225, height = 35)
    button3=Button(Fenster,text="Remove latest",bg = '#EFEFEF',      command=Remove_Latest  ).place(x = 150, y =350 , width = 225, height = 35)
    button4=Button(Fenster,text="Show Superposition",bg = '#EFEFEF', command=Mode_Decider).place(x = 400, y =300 , width = 225, height = 35)
    button5=Button(Fenster,text="Save Superposition",bg = '#EFEFEF', command=Superpos_Saver ).place(x = 400, y =350 , width = 225, height = 35)
    
    
def CVFIT_Superimposer():
    def Runner():
        global choose
        choose    = int(En1.get())
        if choose == 1:
            Semi_Inf_Planar_FITTER()
        if choose == 2:
            Finit_Planar_FITTER()
        if choose == 3:
            Finit_Planar_Trans_FITTER()
        if choose == 4:
            Semi_Inf_Zyl_Ext_FITTER()
        if choose == 5:
            Finit_Zyl_Ext_FITTER()
        if choose == 6:
            Finit_Zyl_Int_FITTER()
        if choose == 7:
            Semi_Inf_Sphere_Ext_FITTER()
        if choose == 8:
            Finit_Sphere_Int_FITTER()
        if choose == 9:
            Statistical_Finit_Planar_FITTER()
        if choose == 10:
            Statistical_Finit_Zyl_Ext_FITTER()
        if choose == 11:
            Statistical_Finit_Zyl_Int_FITTER()
        if choose == 12:
            Statistical_Finit_Sphere_Int_FITTER()
    Fenster = Toplevel()                                                         
    Fenster.title("Superposition mode") 
    Fenster.geometry("650x400")
    colorbgr = Label(Fenster, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000)  
    CVType_Label  = Label(Fenster,text="Superimpose up to 50 CVs by choosing the CV-Type number 1 to 12", bg = '#80ffbf', font = ('Arial',12,'bold'))
    CVType_Label.place(x = 25, y = 15)    
    Lab1 = Label(Fenster,text=" CV-Type Number", bg = '#80ffbf')
    Lab1.place(x = 25, y = 50 )
    En1  = Entry(Fenster)
    En1.place(x =200, y = 50 , width =90, height = 20)
    button0 =Button (Fenster,text="run",bg = '#EFEFEF', command=Runner).place(x = 200, y = 85, width = 90, height = 50)
    var1 = IntVar()
    Checkbutton(Fenster, text="Additive mode ",bg = '#80ffbf', variable=var1).place(x = 25, y = 85)
    var2 = IntVar()
    Checkbutton(Fenster, text="Randles-Sevcik mode",bg = '#80ffbf', variable=var2).place(x = 25, y = 110)
    Te = Text(Fenster, height=4, width=25, bg = '#80ffbf')
    Te.place(x = 350, y = 50, width = 250, height = 85)
    Te.insert(END, "Additive mode will add\nCVs to one single CV. \nRandles-Sevcik mode\nwill overlay them.")
    SepLine1 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine1.place(x = 0, y = 150, height = 20)
    Label_1 = Label(Fenster,text="1)   Semi-inf. planar", bg = '#80ffbf')
    Label_1.place(x = 25, y = 175 )
    Label_2 = Label(Fenster,text="2)   Fin. reflective planar", bg = '#80ffbf')
    Label_2.place(x = 25, y = 200 )
    Label_3 = Label(Fenster,text="3)   Fin. transmissive planar", bg = '#80ffbf')
    Label_3.place(x = 25, y = 225 )
    Label_4 = Label(Fenster,text="4)   Semi-inf. ext. cylindrical.", bg = '#80ffbf')
    Label_4.place(x = 25, y = 250 )
    Label_5 = Label(Fenster,text="5)   Fin. reflective ext. cylindrical", bg = '#80ffbf')
    Label_5.place(x = 200, y = 175 )
    Label_6 = Label(Fenster,text="6)   Fin. reflective int. cylindrical", bg = '#80ffbf')
    Label_6.place(x = 200, y = 200 )
    Label_7 = Label(Fenster,text="7)   Semi-inf. ext. spherical", bg = '#80ffbf')
    Label_7.place(x = 200, y = 225 )
    Label_8 = Label(Fenster,text="8)   Fin. reflective int. spherical", bg = '#80ffbf')
    Label_8.place(x = 200, y = 250 )
    Label_9 = Label(Fenster,text=" 9)   Statistical fin. reflective planar", bg = '#80ffbf')
    Label_9.place(x = 400, y = 175 )
    Label_10 = Label(Fenster,text="10)   Statistical fin. reflective ext. cylindrical", bg = '#80ffbf')
    Label_10.place(x = 400, y = 200 )
    Label_11 = Label(Fenster,text="11)   Statistical fin. reflective int. cylindrical", bg = '#80ffbf')
    Label_11.place(x = 400, y = 225 )
    Label_12 = Label(Fenster,text="12)   Statistical fin. reflective int. spherical", bg = '#80ffbf')
    Label_12.place(x = 400, y = 250 )
    SepLine2 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine2.place(x = 0, y = 275, height = 20) 
    def Displayer():
        if RS_Mode == 1 and Add_Mode == 1:
            Double_Mode_Caller()
        if RS_Mode == 0 and Add_Mode == 0:
            No_Mode_Warner()
        if RS_Mode == 1 and Add_Mode == 0:
            RS_Mode_Caller()
        if Add_Mode == 1 and RS_Mode == 0:
            Add_Mode_Caller()
    def Mode_Decider():
        global RS_Mode
        global Add_Mode
        Add_Mode  = var1.get()
        RS_Mode   = var2.get()
        Displayer()
    button1=Button(Fenster,text="Begin/Clear",bg = '#EFEFEF',        command=Initializer    ).place(x = 25, y =300 , width = 100, height = 85)
    button2=Button(Fenster,text="Append latest",bg = '#EFEFEF',      command=Append_Latest  ).place(x = 150, y =300 , width = 225, height = 35)
    button3=Button(Fenster,text="Remove latest",bg = '#EFEFEF',      command=Remove_Latest  ).place(x = 150, y =350 , width = 225, height = 35)
    button4=Button(Fenster,text="Show Superposition",bg = '#EFEFEF', command=Mode_Decider).place(x = 400, y =300 , width = 225, height = 35)
    button5=Button(Fenster,text="Save Superposition",bg = '#EFEFEF', command=Superpos_Saver ).place(x = 400, y =350 , width = 225, height = 35)
    
    



    
    

          