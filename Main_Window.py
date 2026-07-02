
# coding: utf-8

#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================
#importing all required modules fr66om Python
#=======================================================================================================================
#=======================================================================================================================
import tkinter as tk                 
import numpy as np
import matplotlib.pyplot
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases          import key_press_handler
from matplotlib.figure                 import Figure
import warnings
# to filter the error in the Talbot inversion, which gets replaced by other functions if it occurs :)
warnings.filterwarnings("ignore", category=RuntimeWarning) 

#=======================================================================================================================
#=======================================================================================================================
## define fonts
#=======================================================================================================================
#=======================================================================================================================
font = {'family': 'Times New Roman', 'color':  'black','weight': 'normal','size': 15,}
matplotlib.pyplot.rcParams['mathtext.fontset'] = 'dejavuserif'
matplotlib.pyplot.rcParams['font.sans-serif'] = ['Times new Roman']
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
#===========================================================================================
#===========================================================================================
#===========================================================================================
#IMPORT ALL WRITTEN MODULES
#===========================================================================================
#===========================================================================================
#===========================================================================================
from Koutecky_Levich              import *
from Tafel                        import *
from Randles_Rev                  import *
from Randles_Irr                  import *
from Cottrell                     import *
from Impedance                    import *
from AC_Cyclic_Voltammetry        import *
from Cyclic_Voltammetry           import *
from RIS_Cyclic_Voltammetry       import *
from ACCV_Superimposer            import *
from CV_Superimposer              import *
from RISCV_Superimposer           import *
from LSA_Cyclic_Voltammetry       import *
from LSACV_Superimposer           import *
from STC_Cyclic_Voltammetry       import *
from STCCV_Superimposer           import *
from NST_Chronoamperometry        import *
from NSTCA_Superimposer           import *
from ChargeDischarge              import *
from CT_DATA_LOADER               import *
from DOUGLAS_GUNN                 import *
from Mass_Transfer_Calculator     import *
from DRT_Interface                import *
from Kramers_Kronig_Check         import *

#------------------------------------------------------------------------------
# time for timed events, serial for serial
# communication of Arduino with the PC. The 
# serial requires the installation of pyserial
# on the operating PC.
#------------------------------------------------------------------------------"""   
import time
import serial
#------------------------------------------------------------------------------
###############################################################################
#    Re-Import the modules which have be split to other files
###############################################################################
from PolArStat_CV_Script_GUI       import*
from PolArStat_CA_Script_GUI       import*
from PolArStat_Calibration_Script  import*

#==========================================================================================
#Define Data selection-functions
#==========================================================================================

def CV_Data():
    Get_CV_Data()
def CA_Data():
    Get_NSTCA_Data() 
def ACCV_Data():
    Get_ACCV_Data() 
def LSACV_Data():
    Get_LSACV_Data()
def STCCV_Data():
    Get_STCCV_Data()
def ImpFit_Data():
    Get_ImpFit_Data()
def ImpDRT_Data():
    Get_ImpDRT_Data()
    
def NotNow():
    Fenster = tk.Toplevel()                                                         
    Fenster.title("Function will follow soon")                         
   
    HammerDuck_image       = tk.PhotoImage(file = r"IMAGES\Hammerduck.png")       #Image as background in toplevel
    background_label       = tk.Label(Fenster, image=HammerDuck_image)            #Image as background in toplevel
    background_label.image = HammerDuck_image                                     #Image as background in toplevel
    w                      = HammerDuck_image.width()                             #Image as background in toplevel
    h                      = HammerDuck_image.height()                            #Image as background in toplevel
    xcoodinate             = (Fenster.winfo_screenwidth()  -  w)/2                #Pop up window in center of screen
    ycoodinate             = (Fenster.winfo_screenheight()  - h)/2                #Pop up window in center of screen
    Fenster.geometry('%dx%d+%d+%d' % (w,h, xcoodinate, ycoodinate))
    background_label.place(x=0, y=0, relwidth=1, relheight=1)
    
    Te = Text(Fenster, font = ('Arial',12,'bold'), padx = 20, pady = 10)
    Te.place(x = 35, y = 50, width = 270, height = 80)
    Te.insert(END, "This function is not available\n  right now. We are doing our \n        best to make it soon.")
    Fenster.resizable(False, False)
   
#===========================================================================================
#===========================================================================================
#===========================================================================================
#Build main Window
#===========================================================================================
#===========================================================================================
#===========================================================================================

root = Tk()
root.title("Polarographica")     
root.iconbitmap(r"IMAGES\PG.ico")                    
BGRimage = PhotoImage(file = r"IMAGES\BGRMAIN2.png")
w        = BGRimage.width()
h        = BGRimage.height()
root.geometry('%dx%d+0+0' % (w,h))
#colorbgr = Label(root, text= "", bg = '#99ebff')
#colorbgr.place(x = 0, y = 0, width =7000, height = 2000)
background_label = Label(root, image=BGRimage)
background_label.place(x=0, y=0, relwidth=1, relheight=1)
root.resizable(False, False)


menu = Menu(root)                              
root.config(menu=menu)

#=================================================================================
#=================================================================================

filemenu = Menu(menu)                                           
menu.add_cascade(label="File", menu=filemenu)           
filemenu.add_command(label="Open CV/LSV for fitting", command=CV_Data)
filemenu.add_command(label="Open AC-CV/AC-LSV for FFT-Analysis", command=ACCV_Data)
filemenu.add_command(label="Open n-Step CA for fitting", command=CA_Data)
filemenu.add_command(label="Open LSA-CV for fitting", command=LSACV_Data)
filemenu.add_command(label="Open STC-CV for fitting", command=STCCV_Data)
filemenu.add_command(label="Open Impedance/Nyquist for fitting", command=ImpFit_Data)
filemenu.add_command(label="Open Impedance/Nyquist for DRT analysis", command=ImpDRT_Data)


#=================================================================================
#=================================================================================
#=================================================================================
#Classical Evaluation techniques
#=================================================================================
#=================================================================================
#=================================================================================

RSREVphoto        = PhotoImage(file = r"IMAGES\RSrev.png")
RSREVphotoimage   = RSREVphoto.subsample(2, 2) 
RSIRRphoto        = PhotoImage(file = r"IMAGES\RSirr.png")
RSIRRphotoimage   = RSIRRphoto.subsample(2, 2) 
KLLSVphoto        = PhotoImage(file = r"IMAGES\KL_LSV.png")
KLLSVphotoimage   = KLLSVphoto.subsample(2, 2) 
TAFELphoto        = PhotoImage(file = r"IMAGES\TAFEL.png")
TAFELphotoimage   = TAFELphoto.subsample(2, 2) 
COTphoto          = PhotoImage(file = r"IMAGES\COT.png")
COTphotoimage     = COTphoto.subsample(2, 2) 


SemiInfDiff   = Menu(menu)
menu.add_cascade(label="Classical Evaluations", menu=SemiInfDiff)
SemiInfDiff.add_command(label="Randles Sevcik Reversible", command=Get_RSREV_Data,  image=RSREVphotoimage,  compound=tk.LEFT)
SemiInfDiff.add_command(label="Randles Sevcik Irreversible", command=Get_RSIRR_Data,  image=RSIRRphotoimage,  compound=tk.LEFT)
SemiInfDiff.add_command(label="Koutecky-Levich Analysis", command=Get_KL_Data, image=KLLSVphotoimage,  compound=tk.LEFT)
SemiInfDiff.add_command(label="Tafel and Shape-Analysis", command=Get_Tafel_Data, image=TAFELphotoimage,  compound=tk.LEFT)
SemiInfDiff.add_command(label="Cottrell-Analysis", command=Get_Cottrell_Data, image=COTphotoimage,  compound=tk.LEFT)

#=================================================================================
#=================================================================================
#=================================================================================
#Voltamperometric Simulation-techniques
#=================================================================================
#=================================================================================
#=================================================================================

SimVoltamperometric = Menu(menu)
menu.add_cascade(label="Voltamperometric Simulators", menu=SimVoltamperometric)

#----------------------------
#Cyclic Voltammetry simulators
#----------------------------
CVSimulators   = Menu(SimVoltamperometric)
CVphoto        = PhotoImage(file = r"IMAGES\CV.png")
CVphotoimage   = CVphoto.subsample(2, 2) 
CVSimulators.add_command(label="Planar Semi-Infinite Diffusion CV Simulator",               command=Semi_Inf_Planar)
CVSimulators.add_command(label="Planar Finite Reflective Diffusion CV Simulator",           command=Finit_Planar)
CVSimulators.add_command(label="Planar Finite Transmissive Diffusion CV Simulator",         command=Finit_Planar_Trans)
CVSimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion CV Simulator", command=Semi_Inf_Zyl_Ext)
CVSimulators.add_command(label="External Cylindrical Finite Diffusion CV Simulator",        command=Finit_Zyl_Ext)
CVSimulators.add_command(label="Internal Cylindrical Finite Diffusion CV Simulator",        command=Finit_Zyl_Int)
CVSimulators.add_command(label="External Spherical Semi-Infinite Diffusion CV Simulator",   command=Semi_Inf_Sphere_Ext)
CVSimulators.add_command(label="Internal Spherical Finite Diffusion CV Simulator",          command=Finit_Sphere_Int)
CVSimulators.add_command(label="--------------------------------------------------------------------")
CVSimulators.add_command(label="Statistically Weighted Planar Finite Diffusion CV Simulator",               command=Statistical_Finit_Planar)
CVSimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion CV Simulator", command=Statistical_Finit_Zyl_Ext)
CVSimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion CV Simulator", command=Statistical_Finit_Zyl_Int)
CVSimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion CV Simulator",   command=Statistical_Finit_Sphere_Int)
CVSimulators.add_command(label="--------------------------------------------------------------------")
CVSimulators.add_command(label = 'Superimpose different CVs', command=CVSIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'Cyclic Voltammetry', menu = CVSimulators,  image=CVphotoimage,  compound=tk.LEFT)
#----------------------------
#AC-Cyclic Voltammetry Simulators
#----------------------------
ACCVphoto        = PhotoImage(file = r"IMAGES\ACCV.png")
ACCVphotoimage   = ACCVphoto.subsample(2, 2) 
ACCVSimulators   = Menu(SimVoltamperometric)
ACCVSimulators.add_command(label="Planar Semi-Infinite Diffusion AC-CV Simulator", command=AC_Semi_Inf_Planar)
ACCVSimulators.add_command(label="Planar Finite Reflective Diffusion AC-CV Simulator", command=AC_Finit_Planar)
ACCVSimulators.add_command(label="Planar Finite Transmissive Diffusion AC-CV Simulator", command=AC_Finit_Transm_Planar)
ACCVSimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion AC-CV Simulator", command=AC_Semi_Inf_Zyl_Ext)
ACCVSimulators.add_command(label="External Cylindrical Finite Diffusion AC-CV Simulator", command=AC_Finit_Zyl_Ext)
ACCVSimulators.add_command(label="Internal Cylindrical Finite Diffusion AC-CV Simulator", command=AC_Finit_Zyl_Int)
ACCVSimulators.add_command(label="External Spherical Semi-Infinite Diffusion AC-CV Simulator", command=AC_Semi_Inf_Sphere_Ext)
ACCVSimulators.add_command(label="Internal Spherical Finite Diffusion AC-CV Simulator", command=AC_Finit_Sphere_Int)
ACCVSimulators.add_command(label="--------------------------------------------------------------------")
ACCVSimulators.add_command(label="Statistically Weighted Planar Finite Diffusion AC-CV Simulator", command=AC_Statistical_Finit_Planar)
ACCVSimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion AC-CV Simulator", command=AC_Statistical_Finit_Zyl_Ext)
ACCVSimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion AC-CV Simulator", command=AC_Statistical_Finit_Zyl_Int)
ACCVSimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion AC-CV Simulator", command=AC_Statistical_Finit_Sphere_Int)
ACCVSimulators.add_command(label="--------------------------------------------------------------------")
ACCVSimulators.add_command(label = 'Superimpose different AC-CVs', command=ACCV_SIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'AC-Cyclic Voltammetry', menu = ACCVSimulators,  image=ACCVphotoimage,  compound=tk.LEFT)
#----------------------------
#Large Sine Amplitude Cyclic Voltammetry simulators
#----------------------------
LSACVphoto         = PhotoImage(file = r"IMAGES\LSACV.png")
LSACVphotoimage    = LSACVphoto.subsample(2, 2) 
LSA_CVSimulators   = Menu(SimVoltamperometric)
LSA_CVSimulators.add_command(label="Planar Semi-Infinite Diffusion LSA-CV Simulator",               command=LSA_Semi_Inf_Planar)
LSA_CVSimulators.add_command(label="Planar Finite Reflective Diffusion LSA-CV Simulator",           command=LSA_Finit_Planar)
LSA_CVSimulators.add_command(label="Planar Finite Transmissive Diffusion LSA-CV Simulator",         command=LSA_Finit_Planar_Trans)
LSA_CVSimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion LSA-CV Simulator", command=LSA_Semi_Inf_Zyl_Ext)
LSA_CVSimulators.add_command(label="External Cylindrical Finite Diffusion LSA-CV Simulator",        command=LSA_Finit_Zyl_Ext)
LSA_CVSimulators.add_command(label="Internal Cylindrical Finite Diffusion LSA-CV Simulator",        command=LSA_Finit_Zyl_Int)
LSA_CVSimulators.add_command(label="External Spherical Semi-Infinite Diffusion LSA-CV Simulator",   command=LSA_Semi_Inf_Sphere_Ext)
LSA_CVSimulators.add_command(label="Internal Spherical Finite Diffusion LSA-CV Simulator",          command=LSA_Finit_Sphere_Int)
LSA_CVSimulators.add_command(label="--------------------------------------------------------------------")
LSA_CVSimulators.add_command(label="Statistically Weighted Planar Finite Diffusion LSA-CV Simulator",               command=LSA_Statistical_Finit_Planar)
LSA_CVSimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion LSA-CV Simulator", command=LSA_Statistical_Finit_Zyl_Ext)
LSA_CVSimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion LSA-CV Simulator", command=LSA_Statistical_Finit_Zyl_Int)
LSA_CVSimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion LSA-CV Simulator",   command=LSA_Statistical_Finit_Sphere_Int)
LSA_CVSimulators.add_command(label="--------------------------------------------------------------------")
LSA_CVSimulators.add_command(label = 'Superimpose different LSA-CVs', command=LSA_CVSIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'Large Sine Amplitude Voltammetry', menu = LSA_CVSimulators, image=LSACVphotoimage,  compound=tk.LEFT)
#----------------------------
#n-Step Chronoamperometry Simulators
#----------------------------
CASimulators   = Menu(SimVoltamperometric)
CAphoto        = PhotoImage(file = r"IMAGES\CA.png")
CAphotoimage   = CAphoto.subsample(2, 2) 
CASimulators.add_command(label="Planar Semi-Infinite Diffusion n-Step CA Simulator",               command=NST_Semi_Inf_Planar)
CASimulators.add_command(label="Planar Finite Reflective Diffusion n-Step CA Simulator",           command=NST_Finit_Planar)
CASimulators.add_command(label="Planar Finite Transmissive Diffusion n-Step CA Simulator",         command=NST_Finit_Planar_Trans)
CASimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion n-Step CA Simulator", command=NST_Semi_Inf_Zyl_Ext)
CASimulators.add_command(label="External Cylindrical Finite Diffusion n-Step CA Simulator",        command=NST_Finit_Zyl_Ext)
CASimulators.add_command(label="Internal Cylindrical Finite Diffusion n-Step CA Simulator",        command=NST_Semi_Inf_Sphere_Ext)
CASimulators.add_command(label="External Spherical Semi-Infinite Diffusion n-Step CA Simulator",   command=NST_Semi_Inf_Sphere_Ext)
CASimulators.add_command(label="Internal Spherical Finite Diffusion n-Step CA Simulator",          command=NST_Finit_Sphere_Int)
CASimulators.add_command(label="--------------------------------------------------------------------")
CASimulators.add_command(label="Statistically Weighted Planar Finite Diffusion n-Step CA Simulator",                command=NST_Statistical_Finit_Planar)
CASimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion n-Step CA Simulator",  command=NST_Statistical_Finit_Zyl_Ext)
CASimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion n-Step CA Simulator",  command=NST_Statistical_Finit_Zyl_Int)
CASimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion n-Step CA Simulator",    command=NST_Statistical_Finit_Sphere_Int)
CASimulators.add_command(label="--------------------------------------------------------------------")
CASimulators.add_command(label = 'Superimpose different LSA-CVs', command=NST_CASIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'n-Step Chronoamperometry', menu = CASimulators, image=CAphotoimage,  compound=tk.LEFT)
#----------------------------
#Cyclic Staicase Voltammetry simulators
#----------------------------


STC_CVSimulators   = Menu(SimVoltamperometric)
STCVphoto          = PhotoImage(file = r"IMAGES\STCV.png")
STCVphotoimage     = STCVphoto.subsample(2, 2) 
STC_CVSimulators.add_command(label="Planar Semi-Infinite Diffusion Staircase-CV Simulator",               command=STC_Semi_Inf_Planar)
STC_CVSimulators.add_command(label="Planar Finite Reflective Diffusion Staircase-CV Simulator",           command=STC_Finit_Planar)
STC_CVSimulators.add_command(label="Planar Finite Transmissive Diffusion Staircase-CV Simulator",         command=STC_Finit_Planar_Trans)
STC_CVSimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion Staircase-CV Simulator", command=STC_Semi_Inf_Zyl_Ext)
STC_CVSimulators.add_command(label="External Cylindrical Finite Diffusion Staircase-CV Simulator",        command=STC_Finit_Zyl_Ext)
STC_CVSimulators.add_command(label="Internal Cylindrical Finite Diffusion Staircase-CV Simulator",        command=STC_Finit_Zyl_Int)
STC_CVSimulators.add_command(label="External Spherical Semi-Infinite Diffusion Staircase-CV Simulator",   command=STC_Semi_Inf_Sphere_Ext)
STC_CVSimulators.add_command(label="Internal Spherical Finite Diffusion Staircase-CV Simulator",          command=STC_Finit_Sphere_Int)
STC_CVSimulators.add_command(label="--------------------------------------------------------------------")
STC_CVSimulators.add_command(label="Statistically Weighted Planar Finite Diffusion Staircase-CV Simulator",               command=STC_Statistical_Finit_Planar)
STC_CVSimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion Staircase-CV Simulator", command=STC_Statistical_Finit_Zyl_Ext)
STC_CVSimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion Staircase-CV Simulator", command=STC_Statistical_Finit_Zyl_Int)
STC_CVSimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion Staircase-CV Simulator",   command=STC_Statistical_Finit_Sphere_Int)
STC_CVSimulators.add_command(label="--------------------------------------------------------------------")
STC_CVSimulators.add_command(label = 'Superimpose different Staircase-CVs', command=STC_CVSIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'Staircase Cyclic Voltammetry', menu = STC_CVSimulators, image=STCVphotoimage,  compound=tk.LEFT)


#----------------------------
#RIS_Cyclic Voltammetry simulators
#----------------------------
RIS_CVSimulators   = Menu(SimVoltamperometric)
RISCVphoto         = PhotoImage(file = r"IMAGES\RISCV.png")
RISCVphotoimage    = RISCVphoto.subsample(2, 2) 
RIS_CVSimulators.add_command(label="Planar Semi-Infinite Diffusion RIS-CV Simulator",               command=RIS_Semi_Inf_Planar)
RIS_CVSimulators.add_command(label="Planar Finite Reflective Diffusion RIS-CV Simulator",           command=RIS_Finit_Planar)
RIS_CVSimulators.add_command(label="Planar Finite Transmissive Diffusion RIS-CV Simulator",         command=RIS_Finit_Planar_Trans)
RIS_CVSimulators.add_command(label="External Cylindrical Semi-Infinite Diffusion RIS-CV Simulator", command=RIS_Semi_Inf_Zyl_Ext)
RIS_CVSimulators.add_command(label="External Cylindrical Finite Diffusion RIS-CV Simulator",        command=RIS_Finit_Zyl_Ext)
RIS_CVSimulators.add_command(label="Internal Cylindrical Finite Diffusion RIS-CV Simulator",        command=RIS_Finit_Zyl_Int)
RIS_CVSimulators.add_command(label="External Spherical Semi-Infinite Diffusion RIS-CV Simulator",   command=RIS_Semi_Inf_Sphere_Ext)
RIS_CVSimulators.add_command(label="Internal Spherical Finite Diffusion RIS-CV Simulator",          command=RIS_Finit_Sphere_Int)
RIS_CVSimulators.add_command(label="--------------------------------------------------------------------")
RIS_CVSimulators.add_command(label="Statistically Weighted Planar Finite Diffusion RIS-CV Simulator",               command=RIS_Statistical_Finit_Planar)
RIS_CVSimulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion RIS-CV Simulator", command=RIS_Statistical_Finit_Zyl_Ext)
RIS_CVSimulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion RIS-CV Simulator", command=RIS_Statistical_Finit_Zyl_Int)
RIS_CVSimulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion RIS-CV Simulator",   command=RIS_Statistical_Finit_Sphere_Int)
RIS_CVSimulators.add_command(label="--------------------------------------------------------------------")
RIS_CVSimulators.add_command(label = 'Superimpose different RIS-CVs', command=RISCV_SIM_Superimposer)
SimVoltamperometric.add_cascade(label = 'Random-input Signal Voltammetry', menu = RIS_CVSimulators, image=RISCVphotoimage,  compound=tk.LEFT)



#=================================================================================
#=================================================================================
#=================================================================================
#Voltamperometric Evaluation-techniques
#=================================================================================
#=================================================================================
#=================================================================================
EvalVoltamperometric = Menu(menu)
menu.add_cascade(label="Voltamperometric Evaluation", menu=EvalVoltamperometric)
#----------------------------
#Cyclic Voltammetry fitters
#----------------------------
CVFitters         = Menu(menu)
CVFITphoto        = PhotoImage(file = r"IMAGES\CVFIT.png")
CVFITphotoimage   = CVFITphoto.subsample(2, 2) 
CVFitters.add_command(label="Planar Semi-Infinite Diffusion CV Fitter", command=Semi_Inf_Planar_FITTER)
CVFitters.add_command(label="Planar Finite Reflective Diffusion CV Fitter", command=Finit_Planar_FITTER)
CVFitters.add_command(label="Planar Finite Transmissive Diffusion CV Fitter",         command=Finit_Planar_Trans_FITTER)
CVFitters.add_command(label="External Cylindrical Semi-Infinite Diffusion CV Fitter", command=Semi_Inf_Zyl_Ext_FITTER)
CVFitters.add_command(label="External Cylindrical Finite Diffusion CV Fitter", command=Finit_Zyl_Ext_FITTER)
CVFitters.add_command(label="Internal Cylindrical Finite Diffusion CV Fitter", command=Finit_Zyl_Int_FITTER)
CVFitters.add_command(label="External Spherical Semi-Infinite Diffusion CV Fitter", command=Semi_Inf_Sphere_Ext_FITTER)
CVFitters.add_command(label="Internal Spherical Finite Diffusion CV Fitter", command=Finit_Sphere_Int_FITTER)
CVFitters.add_command(label="--------------------------------------------------------------------")
CVFitters.add_command(label="Statistically Weighted Planar Finite Reflective Diffusion CV Fitter", command=Statistical_Finit_Planar_FITTER)
CVFitters.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion CV Fitter", command=Statistical_Finit_Zyl_Ext_FITTER)
CVFitters.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion CV Fitter", command=Statistical_Finit_Zyl_Int_FITTER)
CVFitters.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion CV Fitter", command=Statistical_Finit_Sphere_Int_FITTER)
CVFitters.add_command(label="--------------------------------------------------------------------")
CVFitters.add_command(label = 'Superimpose different CVs for fitting', command=CVFIT_Superimposer)
EvalVoltamperometric.add_cascade(label = 'Cyclic Voltammetry Fitters', menu = CVFitters, image=CVFITphotoimage,  compound=tk.LEFT)
#----------------------------
#ACCV-Fourier-Transformer
#----------------------------
ACCVFFTphoto        = PhotoImage(file = r"IMAGES\ACCVFFT.png")
ACCVFFTphotoimage   = ACCVFFTphoto.subsample(2, 2) 
EvalVoltamperometric.add_command(label = 'FFT-ACCV-Analysis', command = FFTACCV_Experimental, image=ACCVFFTphotoimage,  compound=tk.LEFT)

#----------------------------
#n-Step Chronoamperometry fitters
#----------------------------
CAFitters   = Menu(menu)
CAFITphoto        = PhotoImage(file = r"IMAGES\nSTEPCA.png")
CAFITphotoimage   = CAFITphoto.subsample(2, 2) 
CAFitters.add_command(label="Planar Semi-Infinite Diffusion n-Step CA Fitter",               command=NST_Semi_Inf_Planar_FITTER)
CAFitters.add_command(label="Planar Finite Diffusion n-Step CA Fitter",                      command=NST_Finit_Planar_FITTER)
CAFitters.add_command(label="Planar Finite Transmissive Diffusion n-Step CA Fitter",         command=NST_Finit_Planar_Trans_FITTER)
CAFitters.add_command(label="External Cylindrical Semi-Infinite Diffusion n-Step CA Fitter", command=NST_Semi_Inf_Zyl_Ext_FITTER)
CAFitters.add_command(label="External Cylindrical Finite Diffusion n-Step CA Fitter",        command=NST_Finit_Zyl_Ext_FITTER)
CAFitters.add_command(label="Internal Cylindrical Finite Diffusion n-Step CA Fitter",        command=NST_Finit_Zyl_Int_FITTER)
CAFitters.add_command(label="External Spherical Semi-Infinite Diffusion n-Step CA Fitter",   command=NST_Semi_Inf_Sphere_Ext_FITTER)
CAFitters.add_command(label="Internal Spherical Finite Diffusion n-Step CA Fitter",          command=NST_Finit_Sphere_Int_FITTER)
CAFitters.add_command(label="--------------------------------------------------------------------")
CAFitters.add_command(label="Statistically Weighted Planar Finite Diffusion n-Step CA Fitter",               command=NST_Statistical_Finit_Planar_FITTER)
CAFitters.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion n-Step CA Fitter", command=NST_Statistical_Finit_Zyl_Ext_FITTER)
CAFitters.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion n-Step CA Fitter", command=NST_Statistical_Finit_Zyl_Int_FITTER)
CAFitters.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion n-Step CA Fitter",   command=NST_Statistical_Finit_Sphere_Int_FITTER)
CAFitters.add_command(label="--------------------------------------------------------------------")
CAFitters.add_command(label = 'Superimpose different n-Step CAs for fitting', command=NST_CAFIT_Superimposer)
EvalVoltamperometric.add_cascade(label = 'n-Step Chronoamperometry Fitters', menu = CAFitters,  image=CAFITphotoimage,  compound=tk.LEFT)
#----------------------------
#Large Sine Amplitude Cyclic Voltammetry simulators
#----------------------------
LSA_CVFitters        = Menu(menu)
LSACVFITphoto        = PhotoImage(file = r"IMAGES\LSACVFIT.png")
LSACVFITphotoimage   = LSACVFITphoto.subsample(2, 2) 
LSA_CVFitters.add_command(label="Planar Semi-Infinite Diffusion LSA-CV Fitter",               command=LSA_Semi_Inf_Planar_FITTER)
LSA_CVFitters.add_command(label="Planar Finite Reflective Diffusion LSA-CV Fitter",           command=LSA_Finit_Planar_FITTER)
LSA_CVFitters.add_command(label="Planar Finite Transmissive Diffusion LSA-CV Fitter",         command=LSA_Finit_Planar_Trans_FITTER)
LSA_CVFitters.add_command(label="External Cylindrical Semi-Infinite Diffusion LSA-CV Fitter", command=LSA_Semi_Inf_Zyl_Ext_FITTER)
LSA_CVFitters.add_command(label="External Cylindrical Finite Diffusion LSA-CV Fitter",        command=LSA_Finit_Zyl_Ext_FITTER)
LSA_CVFitters.add_command(label="Internal Cylindrical Finite Diffusion LSA-CV Fitter",        command=LSA_Finit_Zyl_Int_FITTER)
LSA_CVFitters.add_command(label="External Spherical Semi-Infinite Diffusion LSA-CV Fitter",   command=LSA_Semi_Inf_Sphere_Ext_FITTER)
LSA_CVFitters.add_command(label="Internal Spherical Finite Diffusion LSA-CV Fitter",          command=LSA_Finit_Sphere_Int_FITTER)
LSA_CVFitters.add_command(label="--------------------------------------------------------------------")
LSA_CVFitters.add_command(label="Statistically Weighted Planar Finite Diffusion LSA-CV Fitter",               command=LSA_Statistical_Finit_Planar_FITTER)
LSA_CVFitters.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion LSA-CV Fitter", command=LSA_Statistical_Finit_Zyl_Ext_FITTER)
LSA_CVFitters.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion LSA-CV Fitter", command=LSA_Statistical_Finit_Zyl_Int_FITTER)
LSA_CVFitters.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion LSA-CV Fitter",   command=LSA_Statistical_Finit_Sphere_Int_FITTER)
LSA_CVFitters.add_command(label="--------------------------------------------------------------------")
LSA_CVFitters.add_command(label = 'Superimpose different LSA-CVs for fitting', command=LSA_CVFIT_Superimposer)
EvalVoltamperometric.add_cascade(label = 'Large Sine Amplitude Voltammetry Fitters', menu = LSA_CVFitters, image=LSACVFITphotoimage,  compound=tk.LEFT)

#----------------------------
#Staircase Cyclic Voltammetry simulators
#----------------------------
STC_CVFitters       = Menu(menu)
STCVFITphoto        = PhotoImage(file = r"IMAGES\STCVFIT.png")
STCVFITphotoimage   = STCVFITphoto.subsample(2, 2) 
STC_CVFitters.add_command(label="Planar Semi-Infinite Diffusion STC-CV Fitter",               command=STC_Semi_Inf_Planar_FITTER)
STC_CVFitters.add_command(label="Planar Finite Reflective Diffusion STC-CV Fitter",           command=STC_Finit_Planar_FITTER)
STC_CVFitters.add_command(label="Planar Finite Transmissive Diffusion STC-CV Fitter",         command=STC_Finit_Planar_Trans_FITTER)
STC_CVFitters.add_command(label="External Cylindrical Semi-Infinite Diffusion STC-CV Fitter", command=STC_Semi_Inf_Zyl_Ext_FITTER)
STC_CVFitters.add_command(label="External Cylindrical Finite Diffusion STC-CV Fitter",        command=STC_Finit_Zyl_Ext_FITTER)
STC_CVFitters.add_command(label="Internal Cylindrical Finite Diffusion STC-CV Fitter",        command=STC_Finit_Zyl_Int_FITTER)
STC_CVFitters.add_command(label="External Spherical Semi-Infinite Diffusion STC-CV Fitter",   command=STC_Semi_Inf_Sphere_Ext_FITTER)
STC_CVFitters.add_command(label="Internal Spherical Finite Diffusion STC-CV Fitter",          command=STC_Finit_Sphere_Int_FITTER)
STC_CVFitters.add_command(label="--------------------------------------------------------------------")
STC_CVFitters.add_command(label="Statistically Weighted Planar Finite Diffusion STC-CV Fitter",               command=STC_Statistical_Finit_Planar_FITTER)
STC_CVFitters.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion STC-CV Fitter", command=STC_Statistical_Finit_Zyl_Ext_FITTER)
STC_CVFitters.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion STC-CV Fitter", command=STC_Statistical_Finit_Zyl_Int_FITTER)
STC_CVFitters.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion STC-CV Fitter",   command=STC_Statistical_Finit_Sphere_Int_FITTER)
STC_CVFitters.add_command(label="--------------------------------------------------------------------")
STC_CVFitters.add_command(label = 'Superimpose different STC-CVs for fitting', command=STC_CVFIT_Superimposer)
EvalVoltamperometric.add_cascade(label = 'Staircase Cyclic Voltammetry Fitters', menu = STC_CVFitters, image=STCVFITphotoimage,  compound=tk.LEFT)





#=================================================================================
#=================================================================================
#=================================================================================
#Electrrochemical Impedance Spectroscopy
#=================================================================================
#=================================================================================
#=================================================================================


Impedance_Techniques = Menu(menu)
menu.add_cascade(label="Impedance Techniques",  menu=Impedance_Techniques)

PEIS_simulators    = Menu(menu)
EISSIMphoto        = PhotoImage(file = r"IMAGES\EISSIM.png")
EISSIMphotoimage   = EISSIMphoto.subsample(2, 2) 
PEIS_simulators.add_command(label="Planar Semi-Infinite Diffusion Impedance Simulator", command=Semi_Inf_Plan_Imp_Sim)
PEIS_simulators.add_command(label="Planar Finite Transmissive Diffusion Impedance Simulator", command=Fin_Plan_Trans_Imp_Sim)
PEIS_simulators.add_command(label="Planar Finite Reflective Diffusion Impedance Simulator", command=Fin_Plan_Ref_Imp_Sim)
PEIS_simulators.add_command(label="Zylindrical Semi-Infinite Diffusion Impedance Simulator", command=Semi_Inf_Zyl_Imp_Sim)
PEIS_simulators.add_command(label="Spherical Semi-Infinite Diffusion Impedance Simulator", command=Semi_Inf_Sphe_Imp_Sim)
PEIS_simulators.add_command(label="Spherical Internal Finite Diffusion Impedance Simulator", command=Finit_Int_Sphe_Imp_Sim)
Impedance_Techniques.add_cascade(label = 'Potentiostatic EIS Simulators', menu = PEIS_simulators, image=EISSIMphotoimage,  compound=tk.LEFT)

PEIS_fitters   = Menu(menu)
EISFITphoto        = PhotoImage(file = r"IMAGES\EISFIT.png")
EISFITphotoimage   = EISFITphoto.subsample(2, 2) 
PEIS_fitters.add_command(label="Planar Semi-Infinite Diffusion Impedance Fitter", command=Semi_Inf_Plan_Imp_Fit)
PEIS_fitters.add_command(label="Planar Finite Transmissive Diffusion Impedance Fitter", command=Fin_Plan_Trans_Imp_Fit)
PEIS_fitters.add_command(label="Planar Finite Reflective Diffusion Impedance Fitter", command=Fin_Plan_Ref_Imp_Fit)
PEIS_fitters.add_command(label="Zylindrical Semi-Infinite Diffusion Impedance Fitter", command=Semi_Inf_Zyl_Imp_Fit)
PEIS_fitters.add_command(label="Spherical Semi-Infinite Diffusion Impedance Fitter", command=Semi_Inf_Sphe_Imp_Fit)
PEIS_fitters.add_command(label="Spherical Internal Finite Diffusion Impedance Fitter", command=Finit_Int_Sphe_Imp_Fit)
Impedance_Techniques.add_cascade(label = 'Potentiostatic EIS Fitters', menu = PEIS_fitters, image=EISFITphotoimage,  compound=tk.LEFT)

DRT_Transformation  = Menu(menu)
DRTphoto            = PhotoImage(file = r"IMAGES\DRT.png")
DRTphotoimage       = DRTphoto.subsample(2, 2) 
DRT_Transformation.add_command(label="DRT-Tools-NNLS-DRT", command=DRT_Tools_NNLS_DRT)
DRT_Transformation.add_command(label="Native spike-DRT", command=Native_NNLS_DRT)
DRT_Transformation.add_command(label="Gaussian RBF-DRT", command=Gaussian_NNLS_DRT)
DRT_Transformation.add_command(label="Cole-Cole CPE RBF-DRT", command=ColeCole_NNLS_DRT)
DRT_Transformation.add_command(label="Havriliak-Negami DRT", command=HavriliaNegami_NNLS_DRT)
Impedance_Techniques.add_cascade(label = 'DRT-Analysis', menu = DRT_Transformation, image=DRTphotoimage,  compound=tk.LEFT)



#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
#                     BUTTONS ON MAIN WINDOW
#=================================================================================
#=================================================================================
#=================================================================================

button=Button(root,text="Run",bg = '#9BA9C5', command = Get_CT_Data).place(x = 215, y = 90, width = 60, height = 22)
button=Button(root,text="Training",bg = '#9BA9C5', command = NotNow).place(x = 215, y = 117, width = 60, height = 22)
button=Button(root,text="Help",bg = '#9BA9C5', command = NotNow).place(    x = 215, y = 144, width = 60, height = 22)

button=Button(root,text="Run",bg = '#9BA9C5', command = DouglasGunnWithOrWithoutNumba).place(     x = 215, y = 190, width = 60, height = 22)
button=Button(root,text="Training",bg = '#9BA9C5', command = NotNow).place(x = 215, y = 217, width = 60, height = 22)
button=Button(root,text="Help",bg = '#9BA9C5', command = NotNow).place(    x = 215, y = 244, width = 60, height = 22)

button=Button(root,text="Run",bg = '#9BA9C5', command = Load_Computed_CA).place(     x = 215, y = 290, width = 60, height = 22)
button=Button(root,text="Training",bg = '#9BA9C5', command = NotNow).place(x = 215, y = 317, width = 60, height = 22)
button=Button(root,text="Help",bg = '#9BA9C5', command = NotNow).place(    x = 215, y = 344, width = 60, height = 22)

button=Button(root,text="Sim. CV",bg = '#9BA9C5', command = External_3DMT).place(     x = 215, y = 390, width = 60, height = 22)
button=Button(root,text="Fit CV",bg = '#9BA9C5', command = External_3DMT_FITTER).place(x = 215, y = 417, width = 60, height = 22)
button=Button(root,text="Help",bg = '#9BA9C5', command = NotNow).place(    x = 215, y = 444, width = 60, height = 22)


button=Button(root,text="Kramers-Kronig",bg = '#9BA9C5', command = Get_Imp_KK_Data).place(x = 740, y = 100, width = 120, height = 25)
button=Button(root,text="S-Curve Method",bg = '#9BA9C5', command = NotNow).place(x = 740, y = 128, width = 120, height = 25)
button=Button(root,text="Multi-DRT",bg = '#9BA9C5', command = NotNow).place(     x = 740, y = 156, width = 120, height = 25)
button=Button(root,text="Help",bg = '#9BA9C5', command = NotNow).place(          x = 740, y = 184, width = 120, height = 25)


button=Button(root,text="Charge/\nDischarge",bg = '#9BA9C5', command = CDC_FirstLevel).place(x = 467, y = 100, width = 100, height = 50)
button=Button(root,text="Help",bg = '#9BA9C5', command = NotNow).place(              x = 467, y = 156, width = 100, height = 50)


button=Button(root,text="CV-Autofit",bg = '#9BA9C5', command = NotNow).place(x = 467, y = 340, width = 100, height = 50)
button=Button(root,text="Help",bg = '#9BA9C5', command = NotNow).place(      x = 467, y = 396, width = 100, height = 50)





button=Button(root,text="Cyclic Voltammetry",bg = '#9BA9C5', command = PotentiostatScript_CV).place(    x = 740, y = 340,  width = 120, height = 25)
button=Button(root,text="Chronoamperometry ",bg = '#9BA9C5', command = PotentiostatScript_CA).place(    x = 740, y = 368, width = 120, height = 25)
button=Button(root,text="Calibration",       bg = '#9BA9C5', command = PotentiostatScript_Calib).place( x = 740, y = 396, width = 120, height = 25)























mainloop()

