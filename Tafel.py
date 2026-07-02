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


#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================
def BaseGetter():
    global BaseData
    Delim = Delimiter
    
    if Delim == 1:
       
        BasePath = askopenfilename()
        BaseData = np.genfromtxt(BasePath, delimiter='\t')

        
    if Delim == 0:
       
        BasePath = askopenfilename()
        BaseData = np.genfromtxt(BasePath, delimiter='')
        
        
    global BasePotenzial
    global BaseStrom
    
    
    if FirstComesxx  == 1:
        BasePotenzial     =  BaseData[FromRowxx:ToRowxx:Readeveryxx,0::2]  * UmrPot              
        BaseStrom         =  BaseData[FromRowxx:ToRowxx:Readeveryxx,1::2]  * UmRStrom

    if FirstComesxx  == 0:
        BasePotenzial     =  BaseData[FromRowxx:ToRowxx:Readeveryxx,1::2]  * UmrPot
        BaseStrom         =  BaseData[FromRowxx:ToRowxx:Readeveryxx,0::2]  * UmRStrom
        
    if BasePotenzial[0,0] > BasePotenzial[1,0]:
        BasePotenzial     = BasePotenzial[::-1]
        BaseStrom         = BaseStrom[::-1]


# In[ ]:

def Open_Tafel_File():
    root = Toplevel()
    root.title("Your Data")
    global File_Was_Loaded
    File_Was_Loaded = 1
    global ExpData   #needed later to differentiate experimental and simulated data of ACCV
    ExpData      = 1
    Delim        = Delimiter
    global data
    if Delim == 1: 
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='\t')    
    if Delim == 0:   
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='')
    #============================================================================
    #Plotting of loaded file
    #============================================================================
    #============================================================================
    #define entry window
    #============================================================================
    f = Figure(figsize=(5, 5))
    b = f.add_subplot(111)
    #============================================================================
    global Potenzial
    global Strom
    if FirstComesxx  == 1:
        Potenzial     =  data[FromRowxx:ToRowxx:Readeveryxx,0::2] * UmrPot            
        Strom         =  data[FromRowxx:ToRowxx:Readeveryxx,1::2] * UmRStrom
    if FirstComesxx  == 0:
        Potenzial     =  data[FromRowxx:ToRowxx:Readeveryxx,1::2] * UmrPot
        Strom         =  data[FromRowxx:ToRowxx:Readeveryxx,0::2] * UmRStrom
    
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
    TafelWindow()
    
    


# In[ ]:

def Get_Tafel_Data():
    Fenster = Toplevel()                                                         
    Fenster.title("Get-Data for Tafel Analysis")                         
    Fenster.geometry("400x270")                                            
    FromRowxxx_label      = Label(Fenster,text="Read from row")
    FromRowxxx_label.grid(row=0, column=0)
    FromRowxxx_Eingabe    = Entry(Fenster)                                               
    FromRowxxx_Eingabe.grid(row=0, column=1)
    ToRowxxx_label       = Label(Fenster,text="Read to row")
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
    def AcceptTafel():
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
        Delimiter     = var5.get()
    def NextTafel():
        Fenster.destroy()
        Open_Tafel_File()
    Accept = Button(Fenster, text="Accept",command=AcceptTafel)
    Accept.grid(row=8, column=0) 
    Next = Button(Fenster, text="Next",command=NextTafel)
    Next.grid(row=9, column=0)  
    
#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================


def TafelWindow():
    Fenster = Toplevel()                                                         
    Fenster.title("Tafel Curve-Parameter-Getter")                         
    Fenster.geometry("400x400")
    OCP_Label = Label(Fenster,text="OCPs vs E \nRef in V")
    OCP_Label.grid(row=0, column=0)
    RotRates_Label = Label(Fenster,text="Rotation rates \n in rpm")
    RotRates_Label.grid(row=0, column=1)
    Concent_Label = Label(Fenster,text="Concentrations as fract.\n of the highest one")
    Concent_Label.grid(row=0, column=2)
    #IF ONLY ONE OCP
    #______________________________________
    var1 = IntVar()
    Checkbutton(Fenster, text="One OCP", variable=var1).grid(row=80, column=0, sticky=W)
    OneOCP_Eingabe = Entry(Fenster)                                               
    OneOCP_Eingabe.grid(row=81, column=0)
    #IF ONLY ONE Rotation rate
    #______________________________________
    var2 = IntVar()
    Checkbutton(Fenster, text="One Rotation", variable=var2).grid(row=80, column=2, sticky=W)
    OneRot_Eingabe = Entry(Fenster)                                               
    OneRot_Eingabe.grid(row=81, column=2)
    OCPs     = []
    RotRates = []
    Concents = []
    for i in range(NumMeas):  
        en = Entry(Fenster)
        en.grid(row=i+1, column=0)
        OCPs.append(en)
        en1 = Entry(Fenster)
        en1.grid(row=i+1, column=1)
        RotRates.append(en1)    
        en2 = Entry(Fenster)
        en2.grid(row=i+1, column=2)
        Concents.append(en2)
    def GetOCPs():    
        OnlyOneOCP = var1.get()
        OCPots = []    
        if OnlyOneOCP ==0:
            for entry in OCPs:
                OCPots.append(float(entry.get()))      
        if OnlyOneOCP ==1:
            for j in range(NumMeas):
                OCPots.append(float(OneOCP_Eingabe.get()))
        global OCPArray
        OCPArray = np.array(OCPots) 
    def GetRots():    
        Rots= []
        for entry in RotRates:
            Rots.append(float(entry.get()))
        global RotRatesArray
        RotRatesArray = np.array(Rots) 
        global VarParRot
        global VarParCon
        VarParRot = 1
        VarParCon = 0
    def GetConcents():    
        OnlyOneRot = var2.get()  
        if OnlyOneRot ==1:
            global RotRate
            RotRate = (float(OneRot_Eingabe.get()))  
        Concs= []
        for entry in Concents:
            Concs.append(float(entry.get()))
        global ConcentrationsArray
        ConcentrationsArray = np.array(Concs) 
        global VarParRot
        global VarParCon
        VarParRot = 0
        VarParCon = 1
    def Next():
        Fenster.destroy()
        TafelWindowLevel2()
    
    button=Button(Fenster,text="Accept OCPs",command=GetOCPs).grid(row=82,column=0)
    button=Button(Fenster,text="Accept RotRates",command=GetRots).grid(row=82,column=1)
    button=Button(Fenster,text="Accept Concentrations",command=GetConcents).grid(row=82,column=2)
    button=Button(Fenster,text="Next",command=Next).grid(row=83,column=1)







def TafelWindowLevel2():
    Fenster = Toplevel()                                                         
    Fenster.title("Tafel-Parameter Getter Level 2")                         
    Fenster.geometry("400x680")
    n_Label = Label(Fenster, text="n")
    n_Label.grid(row=0, column=0)
    n_Eingabe = Entry(Fenster)                                               
    n_Eingabe.grid(row=0, column=1)
    c_Label = Label(Fenster, text="c")
    c_Label.grid(row=1, column=0)
    c_Eingabe = Entry(Fenster)                                               
    c_Eingabe.grid(row=1, column=1)
    A_Label = Label(Fenster,text="A in cm^2")
    A_Label.grid(row=2, column=0)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.grid(row=2, column=1)
    T_Label = Label(Fenster,text="T in °C")
    T_Label.grid(row=3, column=0)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.grid(row=3, column=1)
    Ru_Label = Label(Fenster,text="Ru in Ohm")
    Ru_Label.grid(row=4, column=0)
    Ru_Eingabe = Entry(Fenster)
    Ru_Eingabe.grid(row=4, column=1)
    #Smoothing
    #________________
    var1 = IntVar()
    Checkbutton(Fenster, text="Smooth data", variable=var1).grid(row=5, column=0, sticky=W)
    SmoothFact_Label = Label(Fenster,text="Smoothing factor")
    SmoothFact_Label.grid(row=6, column=0)
    SmoothFact_Eingabe = Entry(Fenster)
    SmoothFact_Eingabe.grid(row=6, column=1)
    SmoothFrom_Label = Label(Fenster,text="Smooth from")
    SmoothFrom_Label.grid(row=7, column=0)
    SmoothFrom_Eingabe = Entry(Fenster)
    SmoothFrom_Eingabe.grid(row=7, column=1)
    SmoothTo_Label = Label(Fenster,text="Smooth to")
    SmoothTo_Label.grid(row=8, column=0)
    SmoothTo_Eingabe = Entry(Fenster)
    SmoothTo_Eingabe.grid(row=8, column=1)
    #Corrections
    #__________________________
    var2 = IntVar()
    Checkbutton(Fenster, text="Auto Correction", variable=var2).grid(row=9, column=0, sticky=W)
    var3 = IntVar()
    Checkbutton(Fenster, text="Manual Correction", variable=var3).grid(row=10, column=0, sticky=W)
    ManuCorr_Label    = Label(Fenster,text="Cap. Curr in microamp.")
    ManuCorr_Label.grid(row=11, column=0)
    ManuCorr_Eingabe  = Entry(Fenster)
    ManuCorr_Eingabe.grid(row=11, column=1)
    var4 = IntVar()
    Checkbutton(Fenster, text="Base Correction", variable=var4).grid(row=12, column=0, sticky=W)
    ManLimPot_Label = Label(Fenster,text="E of I Lim")
    ManLimPot_Label.grid(row=13, column=0)
    ManLimPot_Eingabe = Entry(Fenster)
    ManLimPot_Eingabe.grid(row=13, column=1)
    Tafelreg_Label = Label(Fenster,text="Tafel Region")
    Tafelreg_Label.grid(row=14, column=0)
    Tafelreg_Eingabe = Entry(Fenster)
    Tafelreg_Eingabe.grid(row=14, column=1)
    TKSTP_Label = Label(Fenster,text="Shifter")
    TKSTP_Label.grid(row=15, column=0)
    TKSTP_Eingabe = Entry(Fenster)
    TKSTP_Eingabe.grid(row=15, column=1)
    TafyLimLow_Label = Label(Fenster,text="Plot y-Lim low")
    TafyLimLow_Label.grid(row=16, column=0)
    TafyLimLow_Eingabe = Entry(Fenster)
    TafyLimLow_Eingabe.grid(row=16, column=1)
    TafyLimUp_Label = Label(Fenster,text="Plot y-Lim up")
    TafyLimUp_Label.grid(row=17, column=0)
    TafyLimUp_Eingabe = Entry(Fenster)
    TafyLimUp_Eingabe.grid(row=17, column=1)
    ShowTafTo_Label = Label(Fenster,text="Plot x-Lim up")
    ShowTafTo_Label.grid(row=18, column=0)
    ShowTafTo_Eingabe = Entry(Fenster)
    ShowTafTo_Eingabe.grid(row=18, column=1)
    REFF_Label = Label(Fenster,text="Refining")
    REFF_Label.grid(row=19, column=0)
    REFF_Eingabe = Entry(Fenster)
    REFF_Eingabe.grid(row=19, column=1)
    EZero_Label = Label(Fenster,text="E_Zero in V")
    EZero_Label.grid(row=24, column=0)
    EZero_Eingabe = Entry(Fenster)
    EZero_Eingabe.grid(row=24, column=1)
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Show Uncorr LSV", variable=var5).grid(row=20, column=0, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Separate LSV Plot", variable=var6).grid(row=22, column=0, sticky=W)
    var7 = IntVar()
    Checkbutton(Fenster, text="Separate Tafel Plot", variable=var7).grid(row=22, column=1, sticky=W)
    var8 = IntVar()
    Checkbutton(Fenster, text="Show Shape Plot", variable=var8).grid(row=20, column=1, sticky=W)
    var9 = IntVar()
    Checkbutton(Fenster, text="Draw Tafel Taker", variable=var9).grid(row=21, column=0, sticky=W)
    var10 = IntVar()
    Checkbutton(Fenster, text="Draw Shape Taker", variable=var10).grid(row=21, column=1, sticky=W)
    var11 = IntVar()
    Checkbutton(Fenster, text="Save as txt", variable=var11).grid(row=25, column=0, sticky=W)
    var12 = IntVar()
    Checkbutton(Fenster, text="Calculate k_zero", variable=var12).grid(row=23, column=0, sticky=W)

    def AcceptParams():
        global n              ;  global c               ;  global A   
        global T              ;  global Ru              ;  global Smooth_Data
        global AutoCorr       ;  global ManLimPot       ;  global Tafelreg
        global TKSTP          ;  global TafyLimLow      ;  global TafyLimUp
        global ShowTafTo      ;  global REFF            ;  global ManuCorr
        global BaseCorr       ;  global ShowUncorrLSV   ;  global SeparateLSV
        global SeparateTafel  ;  global WaveShapeShower ;  global WaveShapeTaker
        global TafelTaker     ;  global AsTxtSaver      ;  global KZeroCalcer
        global EZero
        n                  = float(n_Eingabe.get())
        c                  = float(c_Eingabe.get())
        A                  = float(A_Eingabe.get())
        T                  = float(T_Eingabe.get()) + 273.15
        Ru                 = float(Ru_Eingabe.get())
        Smooth_Data        = var1.get()
        if Smooth_Data     == 1:
            global SmoothFact
            global SmoothFrom
            global SmoothTo
            SmoothFact     = float(SmoothFact_Eingabe.get())
            SmoothFrom     = float(SmoothFrom_Eingabe.get())
            SmoothTo       = float(SmoothTo_Eingabe.get())
        ManuCorr           = var3.get()
        if ManuCorr        == 1:
            global ManuCorrCurr
            ManuCorrCurr   = float(ManuCorr_Eingabe.get())
        AutoCorr         = var2.get()
        BaseCorr         = var4.get()
        ShowUncorrLSV    = var5.get()
        SeparateLSV      = var6.get()
        SeparateTafel    = var7.get()
        WaveShapeShower  = var8.get()
        TafelTaker       = var9.get()
        WaveShapeTaker   = var10.get()
        AsTxtSaver       = var11.get()
        KZeroCalcer      = var12.get()
        ManLimPot        = float(ManLimPot_Eingabe.get())
        Tafelreg         = float(Tafelreg_Eingabe.get())
        TKSTP            = float(TKSTP_Eingabe.get())
        TafyLimLow       = float(TafyLimLow_Eingabe.get())
        TafyLimUp        = float(TafyLimUp_Eingabe.get())
        ShowTafTo        = float(ShowTafTo_Eingabe.get())
        REFF             = float(REFF_Eingabe.get())
        if KZeroCalcer == 1:
            EZero = float(EZero_Eingabe.get()) 
    def Next():    
        if BaseCorr == 0:
            TafelWindowLevel3()
        if BaseCorr == 1:
            BaseGetter()
            TafelWindowLevel3()       
    def Back():
        Fenster.destroy()
        TafelWindow()
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=26,column=1)
    button=Button(Fenster,text="Next",command=Next).grid(row=27,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=28,column=1)
    

def TafelWindowLevel3():
    Fenster = Toplevel()                                                         
    Fenster.title("Tafel-Plot")                         
    Fenster.geometry("1200x600") 
    #DIESE ZWEI FUNKTIONEN MÜSSEN HIER DRIN DEFINIERT WERDEN!!!       
    def GeradenFit(xWert,ySchnitt,Steig):
        return  ySchnitt + Steig*xWert
    def EchtesPotFinden(WoSieSuchenSoll,WasSieSuchenSoll):                         
        idxStart = (np.abs(WoSieSuchenSoll-WasSieSuchenSoll)).argmin()               
        return WoSieSuchenSoll[idxStart] 
     
    global Length
    global ShortLength 
    Length = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:])))*REFF)
    ShortLength = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:]))))   
    global KORRPOTENZIALSUPPERARRAY
    global KORRSTROMSUPPERARRAY
    global UNKORRPOTENZIALSUPPERARRAY
    global UNKORRSTROMSUPPERARRAY 
    KORRPOTENZIALSUPPERARRAY   = np.empty([NumMeas,Length])
    KORRSTROMSUPPERARRAY       = np.empty([NumMeas,Length]) 
    UNKORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,ShortLength])
    UNKORRSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    
    f = Figure(figsize=(12, 6))
    bild1 = f.add_subplot(121)
    if SeparateLSV == 1:
        LSVFenster = Toplevel()                                                         
        LSVFenster.title("LSV-Plot")                         
        LSVFenster.geometry("700x700")    
        LSVPlot = Figure(figsize=(6, 6))   
    if SeparateTafel == 1:
        SepTafelFenster = Toplevel()                                                         
        SepTafelFenster.title("Tafel-Plot")                         
        SepTafelFenster.geometry("700x700")  
        SepTafelPlot = Figure(figsize=(6, 6))   
    if WaveShapeShower == 1:
        WaveShapeFenster = Toplevel()                                                         
        WaveShapeFenster.title("Wave-Shape-Plot")                         
        WaveShapeFenster.geometry("700x700")  
        WaveShapePlot = Figure(figsize=(6, 6))
        
    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #ERSTE FOR SCHLEIFE FÜR TAFEL PLOT DATENEXTRAKTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    
    OCPs = OCPArray

    for i in range(int(NumMeas)):
        StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
        UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
        PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
                   
        if DesiReactxx == 0:
            Stromarrays                = StromarraysROH[::]
            UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::]
            Potenzialarrays            = PotenzialarraysROH[::]   -  Ru*Stromarrays/1000000
        
        if DesiReactxx == 1:
            Stromarrays                = StromarraysROH
            UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
            Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
     
        if BaseCorr == 1:
            BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
            if DesiReactxx == 0:
                BaseStromArrays        = BaseStromArraysROH[::]
            if DesiReactxx == 1:
                BaseStromArrays        = BaseStromArraysROH[::]
            Stromarrays            = Stromarrays - BaseStromArrays 
        LenPotArrays           = len(Potenzialarrays)
        AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=int(LenPotArrays*REFF))#, endpoint=True)
        Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
        AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        if Smooth_Data == 1:  
            from scipy.interpolate import UnivariateSpline
            SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
            SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo) 
            IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
            IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
  
            spl = UnivariateSpline(Potenzialarrays, Stromarrays)
            spl.set_smoothing_factor(SmoothFact)
            AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , int(REFF*len(Potenzialarrays)))
            AuflgelStrarrays = spl(AufgelPotArrays)
    
        #___________________________________________________________
        #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
        #___________________________________________________________
        KorrAufgPotArr         = AufgelPotArrays #Ohmsche Korrektur erfolgte weiter oben

        KorrAufgStrArrays      = AuflgelStrarrays  
        if AutoCorr == 1:
                
            OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
            OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
            KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
        
        if ManuCorr == 1:
            KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr
        
        for j in range(Length):
            KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
            KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
        for j in range(ShortLength):         
            UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
            UNKORRSTROMSUPPERARRAY[i,j]     = StromarraysROH[j]
  
        bild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
        if ShowUncorrLSV == 1:
            bild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='', color='k')
        bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
        bild1.set_ylabel('I [$\mu$A]', fontsize=12)
        for axis in ['top','bottom','left','right']:
            bild1.spines[axis].set_linewidth(2)
        bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        bild1.axvline(OCPs[i],color='k')
        bild1.axhline(0, color='k')
        
        if SeparateLSV == 1:
            LSVbild1 = LSVPlot.add_subplot(111)
            LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
            if ShowUncorrLSV == 1:
                LSVbild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='', color='k')
            LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            LSVbild1.set_ylabel('I [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                LSVbild1.spines[axis].set_linewidth(2)
            LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            LSVbild1.axvline(OCPs[i],color='k')
            LSVbild1.axhline(0, color='k')
        #________________________________________
        #________________________________________
                   #JETZT TAFEL!!!!!
        #________________________________________
        #________________________________________
        
        
    global LogFitStromArray
    global LogKinStromArray
    global EtaArray
    
    
    LogFitStromArray = np.empty([len(OCPs),len(KorrAufgPotArr)])
    LogKinStromArray = np.empty([len(OCPs),len(KorrAufgPotArr)])
    EtaArray         = np.empty([len(OCPs),len(KorrAufgPotArr)])
    
    global TafelbetaArray
    global TafelIEQArray
    TafelbetaArray = np.empty(len(OCPs))
    TafelIEQArray  = np.empty(len(OCPs))

    global FormBetaArray
    global FormIEQArray
    FormBetaArray  = np.empty(len(OCPs))
    FormIEQArray   = np.empty(len(OCPs))

    global PotVonHalbLim
    global PotVonDreivLim
    PotVonHalbLim  = np.empty(len(OCPs))
    PotVonDreivLim = np.empty(len(OCPs))
    
    global EtaStartEchtARRAY
    EtaStartEchtARRAY  = np.empty(len(OCPs))
    
    global EtaEndEchtARRAY
    EtaEndEchtARRAY  = np.empty(len(OCPs))
    
    global TafelStartImLSVARRAY
    TafelStartImLSVARRAY  = np.empty(len(OCPs))
    
    global TafelEndImLSVARRAY
    TafelEndImLSVARRAY  = np.empty(len(OCPs))

    if DesiReactxx   == 0:
        EtaStart       = np.log(0.01)*(R*T)/(n*F)  +  TKSTP
        
    if DesiReactxx   == 1:
        EtaStart       = -np.log(0.01)*(R*T)/(n*F)  +  TKSTP
    
      
    EtaEnd         = EtaStart   +    Tafelreg
    
    if TafelTaker == 1:
        for i in range(int(NumMeas)):
            bild1.axvline(OCPs[i]+EtaEnd, color='k', linestyle='--')
            bild1.axvline(OCPs[i]+EtaStart, color='k', linestyle='--')
            bild1.axvline(ManLimPot, color='k', linestyle='-')
    
    
    if SeparateLSV == 1:
        if TafelTaker == 1:
            for i in range(int(NumMeas)):
                LSVbild1.axvline(OCPs[i]+EtaEnd, color='k', linestyle='--')
                LSVbild1.axvline(OCPs[i]+EtaStart, color='k', linestyle='--')
                LSVbild1.axvline(ManLimPot, color='k', linestyle='-')
    

    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #Zweite FOR SCHLEIFE FÜR TAFEL PLOT Berechnungen von allem Tafeligen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    
    
    
      
    for i in range(int(NumMeas)):
        StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
        UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
        PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
        if DesiReactxx == 0:
            Stromarrays                = StromarraysROH[::]
            UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::]
            Potenzialarrays            = PotenzialarraysROH[::]   -  Ru*Stromarrays/1000000
            
            
        if DesiReactxx == 1:
            Stromarrays                = StromarraysROH
            UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
            Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
        if BaseCorr == 1:
            BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
            if DesiReactxx == 0:
                BaseStromArrays        = BaseStromArraysROH[::]
            if DesiReactxx == 1:
                BaseStromArrays        = BaseStromArraysROH[::]
            Stromarrays            = Stromarrays - BaseStromArrays
        
        LenPotArrays            = len(Potenzialarrays)
        
        AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=int(LenPotArrays*REFF))#, endpoint=True)
        Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
        AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        
        if Smooth_Data == 1:
            
            from scipy.interpolate import UnivariateSpline
            
            SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
            SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo)
            
            IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
            IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
            
            
            spl = UnivariateSpline(Potenzialarrays, Stromarrays)
            spl.set_smoothing_factor(SmoothFact)
            AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , int(REFF*len(Potenzialarrays)))
            AuflgelStrarrays = spl(AufgelPotArrays)
            
            
        #___________________________________________________________
        #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
        #___________________________________________________________
        KorrAufgPotArr         = AufgelPotArrays # Ohmsche Korrektur erfolgte schon weiter oben

        KorrAufgStrArrays      = AuflgelStrarrays  
        if AutoCorr == 1:
                
            OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
            OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
            KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
        
        if ManuCorr == 1:
            KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr
        
        
        Eta = KorrAufgPotArr - OCPs[i]
        
        HydroStrom     = KorrAufgStrArrays
        
        EtaTafCurLim    = ManLimPot-OCPs[i]                                          #Überspannung von Limiting Pot
        EtaTafCurLimEcht= EchtesPotFinden(Eta,EtaTafCurLim)  
        EtaTafCurLiEIdx = np.asscalar(np.where(Eta == EtaTafCurLimEcht) [0])

        LimStrom        = np.absolute(HydroStrom[EtaTafCurLiEIdx])
        
        #print LimStrom #BIS HIER HIN GEHT ALLES
        
        
        if DesiReactxx        == 0:
            LimStrom        = -LimStrom
        if DesiReactxx        == 1:
            LimStrom        =  LimStrom
        
        
        #BETA AUS KURVENFORM
        
        
        LimStromHalb        = LimStrom*0.5
        LimStromDrVi        = LimStrom*0.75
        EchtLimStromHalb    = EchtesPotFinden(HydroStrom,LimStromHalb)
        EchtLimStromDrVi    = EchtesPotFinden(HydroStrom,LimStromDrVi)
        IdxEchtLimStromHalb = np.asscalar(np.where(HydroStrom == EchtLimStromHalb) [0])
        IdxEchtLimStromDrVi = np.asscalar(np.where(HydroStrom == EchtLimStromDrVi) [0])
    
        EtaEHalb            = Eta[IdxEchtLimStromHalb]
        EtaDreiviertel      = Eta[IdxEchtLimStromDrVi]
        KurvenFormbeta      = np.log(3)*R*T/(n*F*(EtaEHalb-EtaDreiviertel))
        
        FormBetaArray[i]    = np.absolute(KurvenFormbeta)
        
        #print KurvenFormbeta
        
    
        
        PotVonHalbLim[i]        = EtaEHalb+OCPs[i]
        PotVonDreivLim[i]       = EtaDreiviertel+OCPs[i]
    
        EtaStartEcht            = EchtesPotFinden(Eta,EtaStart)
        EtaStartEchtARRAY[i]    = EtaStartEcht
        EtaStartIndex           = np.asscalar(np.where(Eta == EtaStartEcht) [0])
        TafelStartImLSVARRAY[i] = EtaStartEcht + OCPs[i]
    
    
        EtaEndEcht              = EchtesPotFinden(Eta,EtaEnd)
        EtaEndEchtARRAY[i]      = EtaEndEcht
        EtaEndIndex             = np.asscalar(np.where(Eta == EtaEndEcht) [0])
        TafelEndImLSVARRAY[i]   = EtaEndEcht + OCPs[i]
     
        LogKinStrom         = np.log(np.absolute(LimStrom*HydroStrom/(LimStrom-HydroStrom+10**(-20))))
           
        if DesiReactxx          == 0:
            TafelFitting, pcov  = curve_fit(GeradenFit, Eta[EtaEndIndex:EtaStartIndex], LogKinStrom[EtaEndIndex:EtaStartIndex])
            Tafelsteigung       = TafelFitting[1]
            Tafelbeta           = -Tafelsteigung*R*T/(n*F)
        
        if DesiReactxx          == 1:
            TafelFitting, pcov  = curve_fit(GeradenFit, Eta[EtaStartIndex:EtaEndIndex], LogKinStrom[EtaStartIndex:EtaEndIndex])
            Tafelsteigung       = TafelFitting[1]
            Tafelbeta           = -Tafelsteigung*R*T/(n*F)
        
        
        TafelbetaArray[i]   = np.absolute(Tafelbeta)
    
        IEQ_Tafel           = np.exp(TafelFitting[0])
        TafelIEQArray[i]    = IEQ_Tafel*0.000001/A
    
        IEQFORM             = np.absolute(LimStrom*np.exp(KurvenFormbeta*n*F*(EtaEHalb)/(R*T)))
        FormIEQArray[i]     = IEQFORM*0.000001/A
        

        LogKinStromArray[i] = LogKinStrom
        LogFitStromArray[i] = GeradenFit(Eta, *TafelFitting)
        EtaArray[i]         = Eta
        


    #_______________________________________________
    #MITTELWERTE DER BETAS BERECHNEN
    #_______________________________________________
    
    global BETA_T
    global BETA_F
    global IEQT
    global IEQF
    global KNULLTARRAY
    global KNULLFARRAY
    global KNULLT
    global KNULLF
    
    if VarParCon ==1:
        global CDecayArray
        CDecayArray = ConcentrationsArray*c
    
    BETA_T = np.absolute(np.mean(TafelbetaArray))
    BETA_F = np.absolute(np.mean(FormBetaArray))

    IEQT   = np.mean(TafelIEQArray)
    IEQTPLOTANNO = IEQT*1000000
    IEQF   = np.mean(FormIEQArray)
    IEQFPLOTANNO = IEQF*1000000
    

    if KZeroCalcer ==1:
        if DesiReactxx == 0:
            if VarParRot ==1:
                KNULLTARRAY = (TafelIEQArray/(n*F*c))*np.exp(BETA_T*n*F*(OCPArray-EZero)/(R*T))
                KNULLT      = np.mean(KNULLTARRAY)
                KNULLFARRAY = (FormIEQArray/(n*F*c))*np.exp(BETA_F*n*F*(OCPArray-EZero)/(R*T))
                KNULLF      = np.mean(KNULLFARRAY)
            if VarParCon ==1:
                KNULLTARRAY = (TafelIEQArray/(n*F*CDecayArray))*np.exp(BETA_T*n*F*(OCPArray-EZero)/(R*T))
                KNULLT      = np.mean(KNULLTARRAY)
                KNULLFARRAY = (FormIEQArray/(n*F*CDecayArray))*np.exp(BETA_F*n*F*(OCPArray-EZero)/(R*T))
                KNULLF      = np.mean(KNULLFARRAY)
            
        
        if DesiReactxx == 1:
            if VarParRot ==1:
                KNULLTARRAY = (TafelIEQArray/(n*F*c))*np.exp(-BETA_T*n*F*(OCPArray-EZero)/(R*T))
                KNULLT      = np.mean(KNULLTARRAY)
                KNULLFARRAY = (FormIEQArray/(n*F*c))*np.exp(-BETA_F*n*F*(OCPArray-EZero)/(R*T))
                KNULLF      = np.mean(KNULLFARRAY)
            if VarParCon ==1:
                KNULLTARRAY = (TafelIEQArray/(n*F*CDecayArray))*np.exp(-BETA_T*n*F*(OCPArray-EZero)/(R*T))
                KNULLT      = np.mean(KNULLTARRAY)
                KNULLFARRAY = (FormIEQArray/(n*F*CDecayArray))*np.exp(-BETA_F*n*F*(OCPArray-EZero)/(R*T))
                KNULLF      = np.mean(KNULLFARRAY)

   

    
    #________________________________________________
    #WAVESHAPESHOWER IM LSV PLOT
    #________________________________________________
    if WaveShapeShower == 1:
        if WaveShapeTaker == 1:
            for i in range(int(NumMeas)):
                bild1.axvline(PotVonHalbLim[i], color='g', linestyle='-')
                bild1.axvline(PotVonDreivLim[i], color='g', linestyle='-')
                
                
    if SeparateLSV == 1:
        if WaveShapeShower == 1:
            if WaveShapeTaker == 1:
                for i in range(int(NumMeas)):
                    LSVbild1.axvline(PotVonHalbLim[i], color='g', linestyle='-')
                    LSVbild1.axvline(PotVonDreivLim[i], color='g', linestyle='-')
        
    
    #__________________________________________________
    #TAFELPLOTS
    #__________________________________________________
    
    Tafelbild = f.add_subplot(122)
    for i in range(NumMeas):
            Tafelbild.plot(EtaArray[i],LogKinStromArray[i], color='k',linestyle='-',marker='')
            Tafelbild.plot(EtaArray[i],LogFitStromArray[i],color='r',linestyle='-',marker='')
    
    Tafelbild.set_xlabel('$\eta$'  '[V]', fontsize=12)
    Tafelbild.set_ylabel('ln(I$_k$)', fontsize=12)
    for axis in ['top','bottom','left','right']:
        Tafelbild.spines[axis].set_linewidth(2)
    Tafelbild.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    if DesiReactxx          == 0:
        Tafelbild.set_xlim([ShowTafTo,0]) 
    if DesiReactxx          == 1:
        Tafelbild.set_xlim([0,ShowTafTo]) 
    Tafelbild.set_ylim([TafyLimLow,TafyLimUp])
    
    if DesiReactxx == 0:
        Tafelbild.annotate(u'\u03B2' '$_{Av}$' '=%5.3f'    % BETA_T, xy=(0.4, 0.93), xycoords='axes fraction',fontsize=12)
        Tafelbild.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % IEQTPLOTANNO, xy=(0.4, 0.85), xycoords='axes fraction',fontsize=12)
        
    if DesiReactxx == 1:
        Tafelbild.annotate(u'\u03B2' '$_{Av}$' '=%5.3f'  % BETA_T, xy=(0.4, 0.15), xycoords='axes fraction',fontsize=12)  
        Tafelbild.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % IEQTPLOTANNO, xy=(0.4, 0.10), xycoords='axes fraction',fontsize=12)
    
    

    
    if SeparateTafel == 1:
        SepTafelbild = SepTafelPlot.add_subplot(111)
        for i in range(NumMeas):
                SepTafelbild.plot(EtaArray[i],LogKinStromArray[i], color='k',linestyle='-',marker='')
                SepTafelbild.plot(EtaArray[i],LogFitStromArray[i],color='r',linestyle='-',marker='')
    
        SepTafelbild.set_xlabel('$\eta$'  '[V]', fontsize=12)
        SepTafelbild.set_ylabel('ln(I$_k$)', fontsize=12)
        for axis in ['top','bottom','left','right']:
            SepTafelbild.spines[axis].set_linewidth(2)
        SepTafelbild.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        if DesiReactxx          == 0:
            SepTafelbild.set_xlim([ShowTafTo,0]) 
        if DesiReactxx          == 1:
            SepTafelbild.set_xlim([0,ShowTafTo]) 
        SepTafelbild.set_ylim([TafyLimLow,TafyLimUp])
        
        if DesiReactxx == 0:
            SepTafelbild.annotate(u'\u03B2' '$_{Av}$' '=%5.3f'    % BETA_T, xy=(0.4, 0.93), xycoords='axes fraction',fontsize=12) 
            SepTafelbild.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % IEQTPLOTANNO, xy=(0.4, 0.85), xycoords='axes fraction',fontsize=12)
        if DesiReactxx == 1:
            SepTafelbild.annotate(u'\u03B2' '$_{Av}$' '=%5.3f'  % BETA_T, xy=(0.4, 0.13), xycoords='axes fraction',fontsize=12)  
            SepTafelbild.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % IEQTPLOTANNO, xy=(0.4, 0.10), xycoords='axes fraction',fontsize=12)
    
    
    
    
    if WaveShapeShower == 1: 
        WaveShapebild = WaveShapePlot.add_subplot(111)
        
        if VarParRot == 1:
            WaveShapebild.plot((1/(RotRatesArray**0.5)),np.absolute(FormBetaArray),linestyle='-',marker='.',color='g')
            WaveShapebild.set_xlabel('1/' r'$\sqrt{\omega} $' ' ' r'[1/$\sqrt{min} $]' , fontsize=12)
            
        if VarParCon == 1:
            WaveShapebild.plot(ConcentrationsArray,np.absolute(FormBetaArray),linestyle='-',marker='.',color='g')
            WaveShapebild.set_xlabel('c' '/C' '$_{max}$'  , fontsize=12)
        
        
        
        for axis in ['top','bottom','left','right']:
            WaveShapebild.spines[axis].set_linewidth(2)
        WaveShapebild.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        WaveShapebild.set_ylim([0,1])
        WaveShapebild.set_ylabel(u'\u03B2' '$_{shape}$' , fontsize=12)
        WaveShapebild.annotate(u'\u03B2' '$_{Av}$' '=%5.3f'    % BETA_F, xy=(0.4, 0.93), xycoords='axes fraction',fontsize=12) 
        WaveShapebild.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % IEQFPLOTANNO, xy=(0.4, 0.85), xycoords='axes fraction',fontsize=12)
    

        
    #ALLE CANVAS EINSTELLUNGEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    canvas = FigureCanvasTkAgg(f, master=Fenster)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, Fenster)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
       
    
    if SeparateLSV == 1:
        canvas = FigureCanvasTkAgg(LSVPlot, master=LSVFenster)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, LSVFenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    if SeparateTafel == 1:
        canvas = FigureCanvasTkAgg(SepTafelPlot, master=SepTafelFenster)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, SepTafelFenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
    if WaveShapeShower == 1:
        canvas = FigureCanvasTkAgg(WaveShapePlot, master=WaveShapeFenster)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, WaveShapeFenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
    
    if AsTxtSaver ==1:
    
        TAFELASTXTsaver()
    


# In[ ]:

def TAFELASTXTsaver():

 
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        

        f.write("Beta-Tafel-Average")
        f.write("\t")
        f.write(str(np.asscalar(BETA_T)))
        f.write("\n")
        f.write("Ieq-Tafel-Average in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(IEQT)))
        f.write("\n")
        if KZeroCalcer ==1:
            f.write("k_zero Tafel in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(KNULLT)))
            f.write("\n")
        f.write("Beta-Shape-Average")
        f.write("\t")
        f.write(str(np.asscalar(BETA_F)))
        f.write("\n")
        f.write("Ieq-Shape-Average A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(IEQF)))
        f.write("\n")
        if KZeroCalcer ==1:
            f.write("k_zero Shape in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(KNULLF)))
            f.write("\n")
        f.write("\n")
  
        f.write("Parameters")
        f.write("\n")
        f.write("A in cm^2")
        f.write("\t")
        f.write(str(A))
        f.write("\n")
        f.write("T in K")
        f.write("\t")
        f.write(str(T))
        f.write("\n")
        if VarParCon ==1:
            f.write("Rotation speed in rpm")
            f.write("\t")
            f.write(str(RotRate))
            f.write("\n")
        

        f.write("Tafel region in V")
        f.write("\t")
        f.write(str(Tafelreg))
        f.write("\n")
        f.write("Pot of Lim-Current in V")
        f.write("\t")
        f.write(str(ManLimPot))
        f.write("\n")
        f.write("Shift of ideal Tafelstart in V")
        f.write("\t")
        f.write(str(TKSTP))
        f.write("\n")
        f.write("Refining Factor")
        f.write("\t")
        f.write(str(REFF))
        f.write("\n")
        f.write("Ohmic Resistance in Ohm")
        f.write("\t")
        f.write(str(Ru))
        f.write("\n")
        if VarParCon ==1:
            f.write("Rotation Rate")
            f.write("\t")
            f.write(str(RotRate))
            f.write("\n")
            
        f.write("Correction Type")
        f.write("\t")
        if AutoCorr == 1:
            f.write("Auto")
            f.write("\n")
        if BaseCorr == 1:
            f.write("Base")
            f.write("\n")
        if ManuCorr == 1:
            f.write("Manual by")
            f.write("\t")
            f.write(str(ManuCorrCurr))
            f.write("\n")
        if Smooth_Data == 1:
            f.write("Data got smoothed")
            f.write("\n")
            f.write("Smoothing factor")
            f.write("\t")
            f.write(str(SmoothFact))
            f.write("\n")
            f.write("Smoothing Startpotential in V")
            f.write("\t")
            f.write(str(SmoothFrom))
            f.write("\n")
            f.write("Smoothing Endpotential in V")
            f.write("\t")
            f.write(str(SmoothTo))
            f.write("\n")
            



            

#______________________________________________________________________________________________________________________
#
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("Individual Tafel-Plot Results")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        if VarParRot ==1:
            f.write("Rotation Rates in rpm")   
            f.write("\t")
        if VarParCon ==1:
            f.write("Concentrations in mol/cm^3")   
            f.write("\t")
        f.write("Individual Tafelbeta")   
        f.write("\t")
        f.write("Individual Tafel-Ieq in A/cm^2")   
        f.write("\t")
        if KZeroCalcer ==1:
            f.write("Individual k_zero Tafel in cm/s")   
            f.write("\t")
        f.write("Individual Shapebeta")   
        f.write("\t")
        f.write("Individual Shape-Ieq in A/cm^2")   
        f.write("\t")
        if KZeroCalcer ==1:
            f.write("Individual k_zero Shape in cm/s")   
            f.write("\t")
        f.write("Potential of 0.5max Current in V")   
        f.write("\t")
        f.write("Potential of 0.75max Current in V")   
        f.write("\t")
        f.write("Fit-Startpot. in Tafelplot in V")   
        f.write("\t")
        f.write("Fit-Endpot. in Tafelplot in V")   
        f.write("\t")
        f.write("Tafel-Startpot in LSV-Plot in V")   
        f.write("\t")
        f.write("Tafel-Endpot in LSV-Plot in V")   
        f.write("\t")
        f.write("\n")
        f.write("\n")

        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            if VarParRot ==1:
                f.write(str(np.asscalar(RotRatesArray[i])))
                f.write("\t")
            if VarParCon ==1:
                f.write(str(np.asscalar(c*ConcentrationsArray[i])))
                f.write("\t")   
            f.write(str(np.asscalar(TafelbetaArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(TafelIEQArray[i])))
            f.write("\t")
            if KZeroCalcer ==1:
                f.write(str(np.asscalar(KNULLTARRAY[i])))
                f.write("\t")
            f.write(str(np.asscalar(FormBetaArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(FormIEQArray[i])))
            f.write("\t")
            if KZeroCalcer ==1:
                f.write(str(np.asscalar(KNULLFARRAY[i])))
                f.write("\t")
            f.write(str(np.asscalar(PotVonHalbLim[i])))
            f.write("\t")
            f.write(str(np.asscalar(PotVonDreivLim[i])))
            f.write("\t")
            f.write(str(np.asscalar(EtaStartEchtARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(EtaEndEchtARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(TafelStartImLSVARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(TafelEndImLSVARRAY[i])))
            f.write("\t")
            f.write("\n")
            

            

        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        
        if VarParRot == 1:
            f.write("Rotation Rates in rpm")
            f.write("\t")
            for jj in range(int(NumMeas)):
                f.write(str(np.asscalar(RotRatesArray[jj])))
                f.write("\t")
                f.write("\t")
            f.write("\n")
            f.write("\t")
        
        if VarParCon == 1:
            f.write("Concentrations as Fract of highest")
            f.write("\t")
            for jj in range(int(NumMeas)):
                f.write(str(np.asscalar(ConcentrationsArray[jj])))
                f.write("\t")
                f.write("\t")
            f.write("\n")
            f.write("\t")
    
    
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(KORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(KORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")
#__________________________________________________________________________________________________-        
        for i in range(3):
            f.write("\n")
            
        f.write("Uncorrected LSVs")
        f.write("\n")
        f.write("\n")
        if VarParRot ==1:
            f.write("Rotation rates in rpm")
            f.write("\t")
        if VarParCon ==1:
            f.write("Concentrations as Fract of highest")
            f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            if VarParRot ==1:
                f.write(str(np.asscalar(RotRatesArray[jj])))
                f.write("\t")
                f.write("\t")
            
            
            if VarParCon ==1:
                f.write(str(np.asscalar(ConcentrationsArray[jj])))
                f.write("\t")
                f.write("\t")
        
        f.write("\n")
        f.write("\t")
     
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        for i in range(3):
            f.write("\n")
    
        for i in range(int(ShortLength)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(UNKORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(UNKORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")

#_________________________________________                
                
        f.write("Tafel Plots Data")
        f.write("\n")
        f.write("\n")
        
        if VarParRot == 1:
            f.write("Rotation Rates in rpm")
            f.write("\t")
            for jj in range(int(NumMeas)):
                f.write(str(np.asscalar(RotRatesArray[jj])))
                f.write("\t")
                f.write("\t")
                f.write("\t")
            f.write("\n")
            f.write("\t")
        
        if VarParCon == 1:
            f.write("Concentrations as Fract of highest")
            f.write("\t")
            for jj in range(int(NumMeas)):
                f.write(str(np.asscalar(ConcentrationsArray[jj])))
                f.write("\t")
                f.write("\t")
                f.write("\t")
            f.write("\n")
            f.write("\t")

        for jj in range(int(NumMeas)):
        
            f.write("Eta in V")
            f.write("\t")
            f.write("LN(I*Ilim/(I-ILim))")
            f.write("\t")
            f.write("Tafelfit LN(I*Ilim/(I-ILim))")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(EtaArray[j:j+1:,i:i+1:]))
                d = str(np.asscalar(LogKinStromArray[j:j+1:,i:i+1:]))
                e = str(np.asscalar(LogFitStromArray[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
                f.write(e)
                f.write("\t")
            f.write("\n")

    root.destroy()



# In[ ]:



