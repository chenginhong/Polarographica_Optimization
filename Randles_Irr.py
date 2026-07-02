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

def Open_RSIRR_File():
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
    
    f = Figure(figsize=(5, 5), dpi=100)
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
    
    RS_IRR_Window()
    
    
    
    


# In[ ]:

def Get_RSIRR_Data():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Get-Irreversible Randles-Sevcik Data File")                         
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
    
    
    
    #First Comes
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



    def AcceptRSIRR():
              
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
    
    
    def NextRSIRR():
        Open_RSIRR_File()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    Accept = Button(Fenster, text="Accept",command=AcceptRSIRR)
    Accept.grid(row=8, column=0) 
    
    Next = Button(Fenster, text="Next",command=NextRSIRR)
    Next.grid(row=9, column=0)  
    
    
    
    

# In[ ]:

def RS_IRR_Window():
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Curve-Getter")                         
    Fenster.geometry("400x400")
    
    OCP_Label = Label(Fenster,text="OCPs vs E \nRef in V")
    OCP_Label.grid(row=0, column=0)
    
    Scanrates_Label = Label(Fenster,text="Scanrates \n in mV/s")
    Scanrates_Label.grid(row=0, column=1)
    
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
    Checkbutton(Fenster, text="One Scanrate", variable=var2).grid(row=80, column=2, sticky=W)
    OneSca_Eingabe = Entry(Fenster)                                               
    OneSca_Eingabe.grid(row=81, column=2)
    

    OCPs     = []
    Scanrates = []
    Concents = []

    for i in range(NumMeas):
        
        en = Entry(Fenster)
        en.grid(row=i+1, column=0)
        OCPs.append(en)

        en1 = Entry(Fenster)
        en1.grid(row=i+1, column=1)
        Scanrates.append(en1)
        
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

    
    
    def GetScanrates():    
        Scanrs= []
        for entry in Scanrates:
            Scanrs.append(float(entry.get()))
        global ScanratesArray
        ScanratesArray = np.array(Scanrs) 
        global VarParSca
        global VarParCon
        VarParSca = 1
        VarParCon = 0
        
        
    def GetConcents():    
        OnlyOneSca = var2.get()
        
        if OnlyOneSca ==1:
            global Scanrate
            Scanrate = (float(OneSca_Eingabe.get()))
        
        Concs= []
        for entry in Concents:
            Concs.append(float(entry.get()))
        global ConcentrationsArray
        ConcentrationsArray = np.array(Concs) 
        global VarParSca
        global VarParCon
        VarParSca = 0
        VarParCon = 1
        
   
        
    def Next():
        RS_IRR_WindowLevel2()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    button=Button(Fenster,text="Accept OCPs",command=GetOCPs).grid(row=82,column=0)
    button=Button(Fenster,text="Accept Scanrates",command=GetScanrates).grid(row=82,column=1)
    button=Button(Fenster,text="Accept Concentrations",command=GetConcents).grid(row=82,column=2)
    button=Button(Fenster,text="Next",command=Next).grid(row=83,column=1)


# In[ ]:

def RS_IRR_WindowLevel2():
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Irreversible Parameter-Getter")                         
    Fenster.geometry("400x400")

    
    n_Label = Label(Fenster, text="n")
    n_Label.grid(row=1, column=0)
    n_Eingabe = Entry(Fenster)                                               
    n_Eingabe.grid(row=1, column=1)

    A_Label = Label(Fenster,text="A in cm^2")
    A_Label.grid(row=2, column=0)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.grid(row=2, column=1)
    
    c_Label = Label(Fenster,text="c in mol/cm^3")
    c_Label.grid(row=3, column=0)
    c_Eingabe = Entry(Fenster)
    c_Eingabe.grid(row=3, column=1)
    
    D_Label = Label(Fenster,text="D in cm^2/s")
    D_Label.grid(row=4, column=0)
    D_Eingabe = Entry(Fenster)
    D_Eingabe.grid(row=4, column=1)
    
    
    T_Label = Label(Fenster,text="T in °C")
    T_Label.grid(row=6, column=0)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.grid(row=6, column=1)
    

    
    def AcceptParams():
        global n
        global A
        global T 
        global D
        global c


        n         = (float(n_Eingabe.get()))
        A         = (float(A_Eingabe.get()))
        T         = (float(T_Eingabe.get())) + 273.15
        D         = (float(D_Eingabe.get()))
        c         = (float(c_Eingabe.get()))
        

    def Next():
        RS_IRR_WindowLevel3()
        
        def quit():
            Fenster.destroy()
        quit()   
        
    def Back():
        RS_IRR_Window()
        
        def quit():
            Fenster.destroy()
        quit()   
    
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=8,column=0)
    button=Button(Fenster,text="Next",command=Next).grid(row=9,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=9,column=0)
    


# In[ ]:

def RS_IRR_WindowLevel3():
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Irreversible Parameter-Getter-2")                         
    Fenster.geometry("400x520")
    
    
    REFF_Label = Label(Fenster,text="Refining factor")
    REFF_Label.grid(row=2, column=0)
    REFF_Eingabe = Entry(Fenster)
    REFF_Eingabe.grid(row=2, column=1)
    
    Ru_Label = Label(Fenster,text="Ru in Ohm")
    Ru_Label.grid(row=3, column=0)
    Ru_Eingabe = Entry(Fenster)
    Ru_Eingabe.grid(row=3, column=1)
    
    
    #RefToCotti
    #_________________________

    if VarParSca  == 1:
        var1 = IntVar()
        Checkbutton(Fenster, text="Ref to Cottrell", variable=var1).grid(row=4, column=0, sticky=W)

    
        CottiSlope_Label = Label(Fenster,text="Cott Slope microamp/s^0.5")
        CottiSlope_Label.grid(row=5, column=0)
        CottiSlope_Eingabe = Entry(Fenster)
        CottiSlope_Eingabe.grid(row=5, column=1)
    
      
    #Smoothing
    #________________
    var4 = IntVar()
    Checkbutton(Fenster, text="Smooth data", variable=var4).grid(row=6, column=0, sticky=W)
    
    SmoothFact_Label = Label(Fenster,text="Smoothing factor")
    SmoothFact_Label.grid(row=7, column=0)
    SmoothFact_Eingabe = Entry(Fenster)
    SmoothFact_Eingabe.grid(row=7, column=1)
    
    SmoothFrom_Label = Label(Fenster,text="Smooth from")
    SmoothFrom_Label.grid(row=8, column=0)
    SmoothFrom_Eingabe = Entry(Fenster)
    SmoothFrom_Eingabe.grid(row=8, column=1)
    
    SmoothTo_Label = Label(Fenster,text="Smooth to")
    SmoothTo_Label.grid(row=9, column=0)
    SmoothTo_Eingabe = Entry(Fenster)
    SmoothTo_Eingabe.grid(row=9, column=1)
    

     
    
    #Corrections
    #__________________________
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Auto Correction", variable=var5).grid(row=10, column=0, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Separate RS-Plot", variable=var6).grid(row=13, column=1, sticky=W)
    
    var7 = IntVar()
    Checkbutton(Fenster, text="Manual Correction", variable=var7).grid(row=11, column=0, sticky=W)
    ManuCorr_Label    = Label(Fenster,text="DS-cap.(+ox or -red)\nin (microamp.*s)/mV")
    ManuCorr_Label.grid(row=12, column=0)
    ManuCorr_Eingabe  = Entry(Fenster)
    ManuCorr_Eingabe.grid(row=12, column=1)
    
    var8 = IntVar()
    Checkbutton(Fenster, text="Base Correction", variable=var8).grid(row=10, column=1, sticky=W)
    var9 = IntVar()
    Checkbutton(Fenster, text="Separate LSV-Plot", variable=var9).grid(row=13, column=0, sticky=W)
    var10 = IntVar()
    Checkbutton(Fenster, text="Show Uncorr LSV", variable=var10).grid(row=15, column=0, sticky=W)
    var11 = IntVar()
    Checkbutton(Fenster, text="Matsuda-Ayabe-1", variable=var11).grid(row=14, column=0, sticky=W)
    var12 = IntVar()
    Checkbutton(Fenster, text="Matsuda-Ayabe-2", variable=var12).grid(row=14, column=1, sticky=W)
    var13 = IntVar()
    Checkbutton(Fenster, text="Save Results as Txt", variable=var13).grid(row=18, column=0, sticky=W)
    var14 = IntVar()
    Checkbutton(Fenster, text="Calculate k_zero", variable=var14).grid(row=16, column=0, sticky=W)
    
    EZero_Label = Label(Fenster,text="E_0 vs ref. in V (for calc k_0)")
    EZero_Label.grid(row=17, column=0)
    EZero_Eingabe = Entry(Fenster)
    EZero_Eingabe.grid(row=17, column=1)
    
    
    
    

    def AcceptParams():
               
        global REFF
        global Ru 
        if VarParSca  == 1:
            global RefToCotti    
            global CottiSlope    
            global Cotti_D     
            global Cotti_c 
        global Smooth_Data  
        global SmoothFact     
        global SmoothFrom    
        global SmoothTo     
        global AutoCorr 
        global OnlyRSPlot
        global SeparateLSV
        global SeparateMA1
        global SeparateMA2
        global BaseCorr
        global ManuCorr
        global ManuCorrCurr
        global ShowUncorrLSV
        global AsTxtSaver
        global KZeroCalcer
        global EZero
        
        
        
        REFF        = (float(REFF_Eingabe.get()))
        Ru          = (float(Ru_Eingabe.get()))


        Smooth_Data   = var4.get()
        AutoCorr      = var5.get()
        OnlyRSPlot    = var6.get()
        ManuCorr      = var7.get()
        BaseCorr      = var8.get()
        SeparateLSV   = var9.get()
        ShowUncorrLSV = var10.get()
        SeparateMA1   = var11.get()
        SeparateMA2   = var12.get()
        AsTxtSaver    = var13.get()
        KZeroCalcer   = var14.get()
        
        if KZeroCalcer ==1:
            EZero = (float(EZero_Eingabe.get()))

        
        if ManuCorr == 1:
            ManuCorrCurr = (float(ManuCorr_Eingabe.get()))
        
        
        if  Smooth_Data == 1:
            SmoothFact  = (float(SmoothFact_Eingabe.get()))
            SmoothFrom  = (float(SmoothFrom_Eingabe.get()))
            SmoothTo    = (float(SmoothTo_Eingabe.get()))
        
        
        if VarParSca  == 1:    
            RefToCotti = var1.get()
            
            if RefToCotti  == 1:
                CottiSlope  = (float(CottiSlope_Eingabe.get()))

                
    def Next():
        
        if BaseCorr == 0:
            RS_IRR_WindowLevel4()
        if BaseCorr == 1:
            BaseGetter()
            RS_IRR_WindowLevel4()
        

        
        
    def Back():
        RS_IRR_WindowLevel2()
        
        def quit():
            Fenster.destroy()
        quit()   
    

    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=19,column=1)
    button=Button(Fenster,text="Next",command=Next).grid(row=20 ,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=21,column=1)
    


# In[ ]:

def RS_IRR_WindowLevel4():
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Irreversible Level4")                         
    
    Fenster.geometry("700x700")
    
    
    
    #DIESE FUNKTIONEN MÜSSEN HIER DRIN DEFINIERT WERDEN!!!
        
    def GeradenFit(xWert,ySchnitt,Steig):
        return  ySchnitt + Steig*xWert

    def EchtesPotFinden(WoSieSuchenSoll,WasSieSuchenSoll):                         
        idxStart = (np.abs(WoSieSuchenSoll-WasSieSuchenSoll)).argmin()               
        return WoSieSuchenSoll[idxStart]     
    
    def find_nearest(AuflStromBisPeak,HalbMaxAuflStrom):                         
        idx = (np.abs(AuflStromBisPeak-HalbMaxAuflStrom)).argmin()               
        return AuflStromBisPeak[idx]

    
    global Length
    global ShortLength
    
    Length = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:])))*REFF)
    ShortLength = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:]))))
    
    
    global KORRPOTENZIALSUPPERARRAY
    global KORRSTROMSUPPERARRAY
    global UNKORRPOTENZIALSUPPERARRAY
    global UNKORRSTROMSUPPERARRAY
    global I_Peak
    global FITTED_I_Peak
    global WurzelScanratesarray
    global LnScanratesarray 
    global PeakPotentiale
    global FIT_PeakPotentiale
    global HalbPeakPotenziale
    
    
    
    
    KORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,Length])
    KORRSTROMSUPPERARRAY     = np.empty([NumMeas,Length])
    
    UNKORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,ShortLength])
    UNKORRSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    
    
        
    if VarParSca == 1:
        f = Figure(figsize=(7, 7), dpi=80)
    
    if VarParCon == 1:
        f = Figure(figsize=(7, 7), dpi=80)
        
    
    if SeparateLSV == 1:
        SepLSVFenster = Toplevel()                                                         
        SepLSVFenster.title("LSV-Plot")                         
        SepLSVFenster.geometry("700x700")
        
        SepLSVPlot = Figure(figsize=(6, 6), dpi=100)
    
    
    if SeparateMA1 == 1:
        MA1Fenster = Toplevel()                                                         
        MA1Fenster.title("Matsuda-Ayabe-Plot 1")                         
        MA1Fenster.geometry("700x700")
        
        MA1Plot = Figure(figsize=(6, 6), dpi=100)
    
    if SeparateMA2 == 1:
        MA2Fenster = Toplevel()                                                         
        MA2Fenster.title("Matsuda-Ayabe-Plot 2")                         
        MA2Fenster.geometry("700x700")
        
        MA2Plot = Figure(figsize=(6, 6), dpi=100)
        
    

    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #FOR SCHLEIFE WENN SCANRATES GEÄNDERT WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    
    if VarParSca == 1:
    
        global SHAPE_BETAARRAY
        global RS_IEQARRAY
        global RS_KNULLARRAY
        global LNPLOT_IEQARRAY
        global LNPLOT_KNULLARRAY
        global SHAPE_IEQARRAY
        global SHAPE_KNULLARRAY

   
    
        PeakPotentiale     = np.empty(int(NumMeas))
        HalbPeakPotenziale = np.empty(int(NumMeas))
        SHAPE_BETAARRAY    = np.empty(int(NumMeas))
        LNPLOT_IEQARRAY    = np.empty(int(NumMeas))
        RS_IEQARRAY        = np.empty(int(NumMeas))
        
        
    
        
        OCPs = OCPArray
        
        #________________________________________________________________      
        
               
        I_Peak                     = np.empty(int(NumMeas))
        Potenzialmaximumsarray     = np.empty(int(NumMeas))
        Potenzialhalbmaximumsarray = np.empty(int(NumMeas))
        

        for i in range(int(NumMeas)):
            
            StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
            UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
            if DesiReactxx == 0:
                Stromarrays                = StromarraysROH[::-1]
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::-1]
                Potenzialarrays            = PotenzialarraysROH[::-1]   -  Ru*Stromarrays/1000000
            
            
            if DesiReactxx == 1:
                Stromarrays                = StromarraysROH
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
                Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
            if BaseCorr == 1:
                BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
                if DesiReactxx == 0:
                    BaseStromArrays        = BaseStromArraysROH[::-1]
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
            KorrAufgPotArr         = AufgelPotArrays #Korrektur war schon weiter oben
            
 
            KorrAufgStrArrays      = AuflgelStrarrays  
            if AutoCorr == 1:
                
                OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
                OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
                KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
             
            if ManuCorr == 1:
                KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr*ScanratesArray[i] #MauCorrCurr ist hier eine Kapazität von oben
                #Es ist nur noch als Curr benannt, weil ich zu faul zum Umdefinieren bin                                                                    
            
            for j in range(Length):
                
                KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
                KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
            for j in range(ShortLength):
                
                UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
                UNKORRSTROMSUPPERARRAY[i,j]     = Stromarrays[j]
            
            
            
            
            
            
            if DesiReactxx == 0:
                PeakMax    = np.min(KorrAufgStrArrays) 
            
            if DesiReactxx == 1:
                PeakMax    = np.max(KorrAufgStrArrays) 
            
            
    
            I_Peak[i]       = PeakMax
            
            #AB HIER NEUES
        
            IndexPeak                    = np.asscalar(np.where(KorrAufgStrArrays == PeakMax) [0]) 
            Peakpotenzial                = KorrAufgPotArr[IndexPeak]
            Potenzialmaximumsarray[i]    = Peakpotenzial
            PeakPotentiale[i]            = Peakpotenzial  #globales zum Rausnehmen
            HalbMaxAuflStrom             = PeakMax*0.5
            AuflStromBisPeak             = KorrAufgStrArrays[0:IndexPeak]
            HalbMaxAuflStromBisPeak      = (find_nearest(AuflStromBisPeak, HalbMaxAuflStrom))
            IndexHalbMaxAuflStromBisPeak = np.asscalar(np.where(AuflStromBisPeak ==    HalbMaxAuflStromBisPeak) [0])
            Halbpeakpotenzial            = (np.asscalar(KorrAufgPotArr[IndexHalbMaxAuflStromBisPeak]))
            Potenzialhalbmaximumsarray[i]= Halbpeakpotenzial
            HalbPeakPotenziale[i]        = Halbpeakpotenzial #globales zum Rausnehmen    
            
            #JETZT MÜSSTE ALLES WICHTIGE DEFINIERT SEIN
        
        
        
            bild1 = f.add_subplot(221)
            if ShowUncorrLSV ==1:
                bild1.plot(UnkorrPotenzialarrays,Stromarrays*0.001, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays*0.001, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [mA]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            bild1.axvline(Potenzialhalbmaximumsarray[i],color='r')
            bild1.axvline(Potenzialmaximumsarray[i],color='b')    
            
            
            if SeparateLSV == 1:
                LSVbild1 = SepLSVPlot.add_subplot(111)
                if ShowUncorrLSV ==1:
                    LSVbild1.plot(UnkorrPotenzialarrays,Stromarrays*0.001, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays*0.001, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [mA]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
                LSVbild1.axvline(Potenzialhalbmaximumsarray[i],color='r')
                LSVbild1.axvline(Potenzialmaximumsarray[i],color='b')    
        

        

        WurzelScanratesarray = (ScanratesArray*0.001)**0.5
        LnScanratesarray     = np.log(ScanratesArray*0.001)
        
        #FITFUNKTIONEN NOCH OHNE SCANRATELIMITIERUNG
        
        RSIRRFitting, pcov = curve_fit(GeradenFit, WurzelScanratesarray,I_Peak)  
        
        FITTED_I_Peak = GeradenFit(WurzelScanratesarray, *RSIRRFitting)
        
        #RSIRR-FIT Beta und IEQ
        #____________________________________________________________________________________________________________
        
        beta_RS = (RSIRRFitting[1]*0.000001*(R*T)**0.5 / (0.4967*n*F*A*c*(n*F*D)**0.5))**2
        global RS_BETA
        RS_BETA = beta_RS
        
        #Wenn Cottiref
        if RefToCotti == 1:
            beta_n_RS = np.absolute((RSIRRFitting[1]*(R*T)**0.5)/(np.pi**0.5 * 0.4967 * CottiSlope * F**0.5))**2      
            global nBETA_RS
            nBETA_RS = beta_n_RS  #globales Rausholen
            
        #IEQs Berechnen
        global RS_IEQ
        global RS_KNULL
        if DesiReactxx == 0:
            for i in range(int(NumMeas)):
                RS_IEQARRAY[i] = np.absolute(2.1815*n*F*c*(beta_RS*n*F*D*0.001*ScanratesArray[i]/(R*T))**0.5 * np.exp(beta_RS*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            RS_IEQ      = np.mean(RS_IEQARRAY)
            RS_IEQma    = RS_IEQ*1000000
            
            if KZeroCalcer ==1:
                RS_KNULLARRAY = (RS_IEQARRAY/(n*F*c))*np.exp(RS_BETA*n*F*(OCPArray-EZero)/(R*T))
                RS_KNULL      = np.mean(RS_KNULLARRAY)
            
        if DesiReactxx == 1:
            for i in range(int(NumMeas)):
                RS_IEQARRAY[i] = np.absolute(2.1815*n*F*c*(beta_RS*n*F*D*0.001*ScanratesArray[i]/(R*T))**0.5 * np.exp(-beta_RS*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            RS_IEQ      = np.mean(RS_IEQARRAY)
            RS_IEQma    = RS_IEQ*1000000
            
            if KZeroCalcer ==1:
                RS_KNULLARRAY = (RS_IEQARRAY/(n*F*c))*np.exp(-RS_BETA*n*F*(OCPArray-EZero)/(R*T))
                RS_KNULL      = np.mean(RS_KNULLARRAY)
   
        #___________________________________________________________________________________________________________
        
        
        #AYABE1-FIT Beta und IEQ
        #____________________________________________________________________________________________________________
        AyabeFit,pcov     = curve_fit(GeradenFit, LnScanratesarray,Potenzialmaximumsarray)
        betaVonAY = -(R*T/(2*n*F*AyabeFit[1])) 
        betaAy1 = np.absolute(betaVonAY)
        global LNPLOT_BETA
        LNPLOT_BETA = betaAy1 #globales Rausholen
        FIT_PeakPotentiale = GeradenFit(LnScanratesarray, *AyabeFit*1000)
        global LNPLOT_IEQ
        global LNPLOT_KNULL
        if DesiReactxx == 0:
            for i in range(int(NumMeas)):
                LNPLOT_IEQARRAY[i]= np.absolute(2.1815*n*F*c*(betaAy1*n*F*D*0.001*ScanratesArray[i]/(R*T))**0.5 * np.exp(betaAy1*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            LNPLOT_IEQ    = np.mean(LNPLOT_IEQARRAY)
            LNPLOT_IEQma  = LNPLOT_IEQ*1000000
            
            if KZeroCalcer ==1:
                LNPLOT_KNULLARRAY = (LNPLOT_IEQARRAY/(n*F*c))*np.exp(LNPLOT_BETA*n*F*(OCPArray-EZero)/(R*T))
                LNPLOT_KNULL    = np.mean(LNPLOT_KNULLARRAY)
        
        if DesiReactxx == 1:
            for i in range(int(NumMeas)):
                LNPLOT_IEQARRAY[i] = np.absolute(2.1815*n*F*c*(betaAy1*n*F*D*0.001*ScanratesArray[i]/(R*T))**0.5 * np.exp(-betaAy1*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            LNPLOT_IEQ    = np.mean(LNPLOT_IEQARRAY)
            LNPLOT_IEQma  = LNPLOT_IEQ*1000000
            
            if KZeroCalcer ==1:
                LNPLOT_KNULLARRAY = (LNPLOT_IEQARRAY/(n*F*c))*np.exp(-LNPLOT_BETA*n*F*(OCPArray-EZero)/(R*T))
                LNPLOT_KNULL    = np.mean(LNPLOT_KNULLARRAY)
        #_____________________________________________________________________________________________________________
        
        
        
        
        
        #AYABE2-Beta und IEQ
        #_____________________________________________________________________________________________________________
        Delta_EP_EPHalb = np.absolute(1.85*R*T/(F*n*(Potenzialhalbmaximumsarray - Potenzialmaximumsarray)))
        SHAPE_BETAARRAY = Delta_EP_EPHalb
        beta = np.absolute(np.mean(Delta_EP_EPHalb))
        betaAy2 = beta
        global SHAPE_BETA
        SHAPE_BETA = betaAy2 #globales Rausholen
        global SHAPE_IEQ
        global SHAPE_KNULL
        if DesiReactxx == 0:
            SHAPE_IEQARRAY = 2.1815*n*F*c*(Delta_EP_EPHalb*n*F*D*0.001*ScanratesArray/(R*T))**0.5 * np.exp(Delta_EP_EPHalb*n*F*(Potenzialmaximumsarray-OCPArray)/(R*T))
            SHAPE_IEQ      = np.mean(SHAPE_IEQARRAY)
            SHAPE_IEQma    = SHAPE_IEQ*1000000
 
            if KZeroCalcer ==1:
                SHAPE_KNULLARRAY = (SHAPE_IEQARRAY/(n*F*c))*np.exp(SHAPE_BETA*n*F*(OCPArray-EZero)/(R*T))
                SHAPE_KNULL    = np.mean(SHAPE_KNULLARRAY)
        
        if DesiReactxx == 1:
            SHAPE_IEQARRAY = 2.1815*n*F*c*(Delta_EP_EPHalb*n*F*D*0.001*ScanratesArray/(R*T))**0.5 * np.exp(-Delta_EP_EPHalb*n*F*(Potenzialmaximumsarray-OCPArray)/(R*T))
            SHAPE_IEQ      = np.mean(SHAPE_IEQARRAY)
            SHAPE_IEQma    = SHAPE_IEQ*1000000
            

            if KZeroCalcer ==1:
                SHAPE_KNULLARRAY = (SHAPE_IEQARRAY/(n*F*c))*np.exp(-SHAPE_BETA*n*F*(OCPArray-EZero)/(R*T))
                SHAPE_KNULL    = np.mean(SHAPE_KNULLARRAY)
            
        #_____________________________________________________________________________________________________________
        
        
        
        
        
        
        
        bild3 = f.add_subplot(223)
        bild3.plot(LnScanratesarray,1000*Potenzialmaximumsarray,linestyle='',marker='.')
        bild3.plot(LnScanratesarray,FIT_PeakPotentiale,linestyle='-',marker='')
        bild3.set_xlabel(r'$\ln{\nu} $'  , fontsize=12)
        bild3.set_ylabel('E$_p$' '-E$_{p/2}$' '[mV]' , fontsize=12)
        for axis in ['top','bottom','left','right']:
            bild3.spines[axis].set_linewidth(2)
        bild3.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        if DesiReactxx == 0:
            bild3.annotate(u'\u03B2=%5.3f'    % betaAy1, xy=(0.4, 0.93), xycoords='axes fraction')
            bild3.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % LNPLOT_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
            
        if DesiReactxx == 1:
            bild3.annotate(u'\u03B2=%5.3f'  % betaAy1, xy=(0.4, 0.13), xycoords='axes fraction')
            bild3.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % LNPLOT_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
        
        
        if SeparateMA1 ==1:
            MA1bild3 = MA1Plot.add_subplot(111)
            MA1bild3.plot(LnScanratesarray,1000*Potenzialmaximumsarray,linestyle='',marker='.')
            MA1bild3.plot(LnScanratesarray,FIT_PeakPotentiale,linestyle='-',marker='')
            MA1bild3.set_xlabel(r'$\ln{\nu} $'  , fontsize=12)
            MA1bild3.set_ylabel('E$_p$' '-E$_{p/2}$' '[mV]' , fontsize=12)
            for axis in ['top','bottom','left','right']:
                MA1bild3.spines[axis].set_linewidth(2)
            MA1bild3.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
            if DesiReactxx == 0:
                MA1bild3.annotate(u'\u03B2=%5.3f'    % betaAy1, xy=(0.4, 0.93), xycoords='axes fraction')
                MA1bild3.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % LNPLOT_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
            if DesiReactxx == 1:
                MA1bild3.annotate(u'\u03B2=%5.3f'  % betaAy1, xy=(0.4, 0.13), xycoords='axes fraction')
                MA1bild3.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'   % LNPLOT_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
            
            canvas = FigureCanvasTkAgg(MA1Plot, master=MA1Fenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, MA1Fenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        
        
        bild4  = f.add_subplot(224)
        bild4.plot(WurzelScanratesarray,Delta_EP_EPHalb,linestyle='-',marker='.')
        bild4.set_xlabel(r'$\sqrt{\nu} $' ' ' r'[$\sqrt{V} $]' , fontsize=12)
        bild4.set_ylabel(u'\u03B2' , fontsize=12)
        for axis in ['top','bottom','left','right']:
            bild4.spines[axis].set_linewidth(2)
        bild4.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        bild4.set_ylim([0,1])
        
        if DesiReactxx == 0:
            bild4.annotate(u'\u03B2=%5.3f'    % betaAy2, xy=(0.4, 0.93), xycoords='axes fraction') 
            bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
        if DesiReactxx == 1:
            bild4.annotate(u'\u03B2=%5.3f'  % betaAy2, xy=(0.4, 0.13), xycoords='axes fraction')  
            bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
        
        if SeparateMA2 ==1:
            MA2bild4  = MA2Plot.add_subplot(111)
            MA2bild4.plot(WurzelScanratesarray,Delta_EP_EPHalb,linestyle='-',marker='.')
            MA2bild4.set_xlabel(r'$\sqrt{\nu} $' ' ' r'[$\sqrt{V} $]' , fontsize=12)
            MA2bild4.set_ylabel(u'\u03B2' , fontsize=12)
            for axis in ['top','bottom','left','right']:
                MA2bild4.spines[axis].set_linewidth(2)
            MA2bild4.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            MA2bild4.set_ylim([0,1])
            
            if DesiReactxx == 0:
                MA2bild4.annotate(u'\u03B2=%5.3f'    % betaAy2, xy=(0.4, 0.93), xycoords='axes fraction')
                MA2bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
            if DesiReactxx == 1:
                MA2bild4.annotate(u'\u03B2=%5.3f'  % betaAy2, xy=(0.4, 0.13), xycoords='axes fraction')  
                MA2bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
            
            canvas = FigureCanvasTkAgg(MA2Plot, master=MA2Fenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, MA2Fenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        
        
        bild2 = f.add_subplot(222)
        
        
        
        
        def RSIRRSCPLOTTER():
        
            bild2.plot(WurzelScanratesarray,0.001*GeradenFit(WurzelScanratesarray, *RSIRRFitting),color='r',linestyle='-',marker='')
            bild2.plot(WurzelScanratesarray,0.001*I_Peak,linestyle='',marker='.')
            bild2.set_xlabel(r'$\sqrt{\nu} $' ' ' r'[$\sqrt{V} $]' , fontsize=12)
            bild2.set_ylabel('I$_p$ [mA]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
               
            
            #___________________________________________________________________________________________________________________
            
            if RefToCotti == 0:
            
                if DesiReactxx == 0:
                    bild2.annotate(u'\u03B2=%8.3f' % beta_RS, xy=(0.4, 0.93), xycoords='axes fraction')
                    bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
                
                if DesiReactxx == 1:
                    bild2.annotate(u'\u03B2=%8.3f' % beta_RS, xy=(0.4, 0.13), xycoords='axes fraction')
                    bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
            
            if RefToCotti == 1:
               
                BeschreibefensterCotti = Toplevel()                                                         
                BeschreibefensterCotti.title("Information")                         
                BeschreibefensterCotti.geometry("200x200")    
                Te = Text(BeschreibefensterCotti, height=6, width=30)
                Te.pack()
                Te.insert(END, "Randles Slope got\ndivided by Cottrell Slope")
                    

                if DesiReactxx == 0:
                    bild2.annotate('n'  u'\u03B2=%5.3f'    % beta_n_RS, xy=(0.4, 0.93), xycoords='axes fraction') 
                    bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
                if DesiReactxx == 1:
                    bild2.annotate('n'  u'\u03B2=%5.3f'  % beta_n_RS, xy=(0.4, 0.13), xycoords='axes fraction') 
                    bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
        

           
        
        RSIRRSCPLOTTER()

        
        
        canvas = FigureCanvasTkAgg(f, master=Fenster)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, Fenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        


        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #ONLY RSIRR PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyRSPlot == 1:
            
            RSIRRFenster = Toplevel()                                                         
            RSIRRFenster.title("Irreversible Randles-Sevcik Plot")                         
            RSIRRFenster.geometry("700x700")
        
            RSIRRPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = RSIRRPlot.add_subplot(111)
            
            RSIRRSCPLOTTER()
                
            canvas = FigureCanvasTkAgg(RSIRRPlot, master=RSIRRFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, RSIRRFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        #Das ist das Canvas für den Separaten LSV Plot, was vllt in der for schleife nicht geht???
        
        if SeparateLSV == 1:
            canvas = FigureCanvasTkAgg(SepLSVPlot, master=SepLSVFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, SepLSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
        if AsTxtSaver ==1:
            RSIRRSCASTXTsaver()
    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #FOR SCHLEIFE WENN CONCENTRATIONS GEÄNDERT WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________ 
    

    
    if VarParCon == 1:
    
        
        global SHAPE_BETAARRAYCONC
        global RS_IEQARRAYCONC
        global RS_KNULLARRAYCONC
        global SHAPE_IEQARRAYCONC
        global SHAPE_KNULLARRAYCONC

   
    
        PeakPotentiale     = np.empty(int(NumMeas))
        HalbPeakPotenziale = np.empty(int(NumMeas))
        SHAPE_BETAARRAYCONC= np.empty(int(NumMeas))
        RS_IEQARRAYCONC    = np.empty(int(NumMeas))
        
        
        OCPs = OCPArray
        
        #________________________________________________________________      
        
               
        I_Peak                     = np.empty(int(NumMeas))
        Potenzialmaximumsarray     = np.empty(int(NumMeas))
        Potenzialhalbmaximumsarray = np.empty(int(NumMeas))
        
        

        for i in range(int(NumMeas)):
            
            StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
            UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
            if DesiReactxx == 0:
                Stromarrays                = StromarraysROH[::-1]
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::-1]
                Potenzialarrays            = PotenzialarraysROH[::-1]   -  Ru*Stromarrays/1000000
            
            
            if DesiReactxx == 1:
                Stromarrays                = StromarraysROH
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
                Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
            if BaseCorr == 1:
                BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
                if DesiReactxx == 0:
                    BaseStromArrays        = BaseStromArraysROH[::-1]
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
            KorrAufgPotArr         = AufgelPotArrays #Korrektur war schon weiter oben
            
 
            KorrAufgStrArrays      = AuflgelStrarrays  
            if AutoCorr == 1:
                
                OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
                OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
                KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
             
            if ManuCorr == 1:
                KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr*Scanrate#MauCorrCurr ist hier eine Kapazität von oben
                #Es ist nur noch als Curr benannt, weil ich zu faul zum Umdefinieren bin 
                
            
            
            for j in range(Length):
                
                KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
                KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
            for j in range(ShortLength):
                
                UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
                UNKORRSTROMSUPPERARRAY[i,j]     = Stromarrays[j]
            
 
            
            
            if DesiReactxx == 0:
                PeakMax    = np.min(KorrAufgStrArrays) 
            
            if DesiReactxx == 1:
                PeakMax    = np.max(KorrAufgStrArrays) 
            
    
            I_Peak[i]       = PeakMax
            
            
            #AB HIER NEUES
        
            IndexPeak                    = np.asscalar(np.where(KorrAufgStrArrays == PeakMax) [0]) 
            Peakpotenzial                = KorrAufgPotArr[IndexPeak]
            Potenzialmaximumsarray[i]    = Peakpotenzial
            PeakPotentiale[i]            = Peakpotenzial  #globales zum Rausnehmen
            HalbMaxAuflStrom             = PeakMax*0.5
            AuflStromBisPeak             = KorrAufgStrArrays[0:IndexPeak]
            HalbMaxAuflStromBisPeak      = (find_nearest(AuflStromBisPeak, HalbMaxAuflStrom))
            IndexHalbMaxAuflStromBisPeak = np.asscalar(np.where(AuflStromBisPeak ==    HalbMaxAuflStromBisPeak) [0])
            Halbpeakpotenzial            = (np.asscalar(KorrAufgPotArr[IndexHalbMaxAuflStromBisPeak]))
            Potenzialhalbmaximumsarray[i]= Halbpeakpotenzial
            HalbPeakPotenziale[i]        = Halbpeakpotenzial #globales zum Rausnehmen    
        
            #JETZT MÜSSTE ALLES WICHTIGE DEFINIERT SEIN
        
        
        
            bild1 = f.add_subplot(221)
            if ShowUncorrLSV ==1:
                bild1.plot(UnkorrPotenzialarrays,Stromarrays*0.001, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays*0.001, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [mA]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            
            if SeparateLSV == 1:
                LSVbild1 = SepLSVPlot.add_subplot(111)
                if ShowUncorrLSV ==1:
                    LSVbild1.plot(UnkorrPotenzialarrays,Stromarrays*0.001, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays*0.001, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [mA]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
                

        global CDecayArray

        CDecayArray = ConcentrationsArray*c
       

        #FITFUNKTIONEN NOCH OHNE SCANRATELIMITIERUNG
        
        RSIRRFitting, pcov = curve_fit(GeradenFit, CDecayArray,I_Peak)  
         
        
        FITTED_I_Peak = GeradenFit(CDecayArray, *RSIRRFitting)
        
        #RSIRR-FIT Beta und IEQ
        #____________________________________________________________________________________________________________
        
        beta_RS = (RSIRRFitting[1]*0.000001*(R*T)**0.5 / (0.4967*n*F*A*(n*F*D*Scanrate*0.001)**0.5))**2
        global RS_BETACONC
        RS_BETACONC = beta_RS
        
        #HIER GEHT KEIN REF TO COTTI
            
        #IEQs Berechnen
        global RS_IEQCONC
        global RS_KNULLCONC
        if DesiReactxx == 0:
            for i in range(int(NumMeas)):
                RS_IEQARRAYCONC[i] = np.absolute(2.1815*n*F*CDecayArray[i]*(beta_RS*n*F*D*0.001*Scanrate/(R*T))**0.5 * np.exp(beta_RS*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            RS_IEQCONC      = np.mean(RS_IEQARRAYCONC)
            RS_IEQma    = RS_IEQCONC*1000000
            
            if KZeroCalcer ==1:
                RS_KNULLARRAYCONC = (RS_IEQARRAYCONC/(n*F*CDecayArray[i]))*np.exp(RS_BETACONC*n*F*(OCPArray-EZero)/(R*T))
                RS_KNULLCONC      = np.mean(RS_KNULLARRAYCONC)
            
        if DesiReactxx == 1:
            for i in range(int(NumMeas)):
                RS_IEQARRAYCONC[i] = np.absolute(2.1815*n*F*CDecayArray[i]*(beta_RS*n*F*D*0.001*Scanrate/(R*T))**0.5 * np.exp(-beta_RS*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            RS_IEQCONC      = np.mean(RS_IEQARRAYCONC)
            RS_IEQma    = RS_IEQCONC*1000000
            
            if KZeroCalcer ==1:
                RS_KNULLARRAYCONC = (RS_IEQARRAYCONC/(n*F*CDecayArray[i]))*np.exp(-RS_BETACONC*n*F*(OCPArray-EZero)/(R*T))
                RS_KNULLCONC      = np.mean(RS_KNULLARRAYCONC)
   
        #___________________________________________________________________________________________________________
        
        
        #AYABE2-Beta und IEQ
        #_____________________________________________________________________________________________________________
        Delta_EP_EPHalb = np.absolute(1.85*R*T/(F*n*(Potenzialhalbmaximumsarray - Potenzialmaximumsarray)))
        SHAPE_BETAARRAYCONC = Delta_EP_EPHalb
        beta = np.absolute(np.mean(Delta_EP_EPHalb))
        betaAy2 = beta
        global SHAPE_BETACONC
        SHAPE_BETACONC = betaAy2 #globales Rausholen
        global SHAPE_IEQCONC
        global SHAPE_KNULLCONC
        if DesiReactxx == 0:
            SHAPE_IEQARRAYCONC = 2.1815*n*F*CDecayArray*(Delta_EP_EPHalb*n*F*D*0.001*Scanrate/(R*T))**0.5 * np.exp(Delta_EP_EPHalb*n*F*(Potenzialmaximumsarray-OCPArray)/(R*T))
            SHAPE_IEQCONC      = np.mean(SHAPE_IEQARRAYCONC)
            SHAPE_IEQma    = SHAPE_IEQCONC*1000000
 
            if KZeroCalcer ==1:
                SHAPE_KNULLARRAYCONC = (SHAPE_IEQARRAYCONC/(n*F*CDecayArray))*np.exp(SHAPE_BETACONC*n*F*(OCPArray-EZero)/(R*T))
                SHAPE_KNULLCONC    = np.mean(SHAPE_KNULLARRAYCONC)
        
        if DesiReactxx == 1:
            SHAPE_IEQARRAYCONC = 2.1815*n*F*CDecayArray*(Delta_EP_EPHalb*n*F*D*0.001*Scanrate/(R*T))**0.5 * np.exp(-Delta_EP_EPHalb*n*F*(Potenzialmaximumsarray-OCPArray)/(R*T))
            SHAPE_IEQCONC      = np.mean(SHAPE_IEQARRAYCONC)
            SHAPE_IEQma    = SHAPE_IEQCONC*1000000
            

            if KZeroCalcer ==1:
                SHAPE_KNULLARRAYCONC = (SHAPE_IEQARRAYCONC/(n*F*CDecayArray))*np.exp(-SHAPE_BETACONC*n*F*(OCPArray-EZero)/(R*T))
                SHAPE_KNULLCONC    = np.mean(SHAPE_KNULLARRAYCONC)
            
        #MATSUDA AYABE 1 --> GEHT HIER NICHT
        #_____________________________________________________________________________________________________________
        
        
        bild3 = f.add_subplot(223)
        bild3.annotate('Matsuda-Ayabe1 NOT possible here' , xy=(0.1, 0.93), xycoords='axes fraction')
        for axis in ['top','bottom','left','right']:
            bild3.spines[axis].set_linewidth(2)
        bild3.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
        
        if SeparateMA1 ==1:
            MA1bild3 = MA1Plot.add_subplot(111)
            MA1bild3.annotate('Matsuda-Ayabe1 NOT possible here', xy=(0.1, 0.93), xycoords='axes fraction')
            for axis in ['top','bottom','left','right']:
                MA1bild3.spines[axis].set_linewidth(2)
            MA1bild3.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
            canvas = FigureCanvasTkAgg(MA1Plot, master=MA1Fenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, MA1Fenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        #_______________________________________________________________________________________________________________    
        
        
                           
        #MATSUDA AYABE 2
        #_________________________________________________________________________________________________________________
                           
        
        bild4  = f.add_subplot(224)
        bild4.plot(ConcentrationsArray,Delta_EP_EPHalb,linestyle='-',marker='.')
        bild4.set_xlabel('c/''c$_{max}$' , fontsize=12)
        bild4.set_ylabel(u'\u03B2' , fontsize=12)
        for axis in ['top','bottom','left','right']:
            bild4.spines[axis].set_linewidth(2)
        bild4.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        bild4.set_ylim([0,1])
        
        if DesiReactxx == 0:
            bild4.annotate(u'\u03B2=%5.3f'    % betaAy2, xy=(0.4, 0.93), xycoords='axes fraction') 
            bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
        if DesiReactxx == 1:
            bild4.annotate(u'\u03B2=%5.3f'  % betaAy2, xy=(0.4, 0.13), xycoords='axes fraction')  
            bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
        
        if SeparateMA2 ==1:
            MA2bild4  = MA2Plot.add_subplot(111)
            MA2bild4.plot(ConcentrationsArray,Delta_EP_EPHalb,linestyle='-',marker='.')
            MA2bild4.set_xlabel('c/''c$_{max}$' , fontsize=12)
            MA2bild4.set_ylabel(u'\u03B2' , fontsize=12)
            for axis in ['top','bottom','left','right']:
                MA2bild4.spines[axis].set_linewidth(2)
            MA2bild4.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            MA2bild4.set_ylim([0,1])
            
            if DesiReactxx == 0:
                MA2bild4.annotate(u'\u03B2=%5.3f'    % betaAy2, xy=(0.4, 0.93), xycoords='axes fraction')
                MA2bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
            if DesiReactxx == 1:
                MA2bild4.annotate(u'\u03B2=%5.3f'  % betaAy2, xy=(0.4, 0.13), xycoords='axes fraction')  
                MA2bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
            
            canvas = FigureCanvasTkAgg(MA2Plot, master=MA2Fenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, MA2Fenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        

        
        
        bild2 = f.add_subplot(222)
        
        def RSIRRCONCPLOTTER():
        
            bild2.plot(ConcentrationsArray,FITTED_I_Peak,color='r',linestyle='-',marker='')
            bild2.plot(ConcentrationsArray,I_Peak,linestyle='',marker='.')
            bild2.set_xlabel('c/''c$_{max}$' , fontsize=12)
            bild2.set_ylabel('I$_p$ [mA]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
            #________________________________________________________
            #________________________________________________________
            #ANNOTOATIONS FOR PARAMETERS AND CALCULATIONS
            #________________________________________________________
            #________________________________________________________
            #__________________________________________________________________________________________________________________
         
            beta_RS = (np.absolute(RSIRRFitting[1]*0.000001*(R*T)**0.5 / (0.4967*n*F*A*(Scanrate*0.001)**0.5 *(n*F*D)**0.5)))**2
            

            if DesiReactxx == 0:
                bild2.annotate(u'\u03B2=%8.3f' % beta_RS, xy=(0.4, 0.93), xycoords='axes fraction')
                bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
                
            if DesiReactxx == 1:
                bild2.annotate(u'\u03B2=%8.3f' % beta_RS, xy=(0.4, 0.13), xycoords='axes fraction')
                bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
            
            
           
        RSIRRCONCPLOTTER()

       
        
        canvas = FigureCanvasTkAgg(f, master=Fenster)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, Fenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        


        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #ONLY RSIRR PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyRSPlot == 1:
            
            RSIRRFenster = Toplevel()                                                         
            RSIRRFenster.title("Irreversible Randles-Sevcik Plot")                         
            RSIRRFenster.geometry("700x700")
        
            RSIRRPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = RSIRRPlot.add_subplot(111)
            
            RSIRRCONCPLOTTER()
                
            canvas = FigureCanvasTkAgg(RSIRRPlot, master=RSIRRFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, RSIRRFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        #Das ist das Canvas für den Separaten LSV Plot, was vllt in der for schleife nicht geht???
        
        if SeparateLSV == 1:
            canvas = FigureCanvasTkAgg(SepLSVPlot, master=SepLSVFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, SepLSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
            
        if AsTxtSaver ==1:
            RSIRRCONCASTXTsaver()
    
    


# In[ ]:

def RSIRRSCASTXTsaver():
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")

        f.write("Beta-RS")
        f.write("\t")
        f.write(str(np.asscalar(RS_BETA)))
        f.write("\n")
        if RefToCotti  ==1:
            f.write("nBeta-RS from Cottrellref.")
            f.write("\t")
            f.write(str(np.asscalar(nBETA_RS)))
            f.write("\n")
        
        f.write("Beta-LN-Plot")
        f.write("\t")
        f.write(str(np.asscalar(LNPLOT_BETA)))
        f.write("\n")
        f.write("Average Beta-Shape-Plot")
        f.write("\t")
        f.write(str(np.asscalar(SHAPE_BETA)))
        f.write("\n")
        f.write("Average IEQ-RS in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(RS_IEQ)))
        f.write("\n")
        f.write("Average IEQ-LN-Plot in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(LNPLOT_IEQ)))
        f.write("\n")
        f.write("Average IEQ-Shape-Plot in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(SHAPE_IEQ)))
        f.write("\n")
        if KZeroCalcer ==1:
            f.write("Average KZero-RS in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(RS_KNULL)))
            f.write("\n")
            f.write("Average KZero-LN-Plot in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(LNPLOT_KNULL)))
            f.write("\n")
            f.write("Average KZero-Shape-Plot in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_KNULL)))
            f.write("\n")
        f.write("\n")
       

        
        
        f.write("Parameters")
        f.write("\n")
        f.write("D in cm^2/s")
        f.write("\t")
        f.write(str(D))
        f.write("\n")
        f.write("A in cm^2")
        f.write("\t")
        f.write(str(A))
        f.write("\n")
        f.write("T in K")
        f.write("\t")
        f.write(str(T))
        f.write("\n")
        f.write("c in mol/cm^3")
        f.write("\t")
        f.write(str(c))
        f.write("\n")
        f.write("Refining Factor")
        f.write("\t")
        f.write(str(REFF))
        f.write("\n")
        f.write("Ohmic Resistance in Ohm")
        f.write("\t")
        f.write(str(Ru))
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
            f.write("\t")
            f.write("(microams*s)/mV")
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
        if RefToCotti  ==1:
            f.write("D ref. to Cottrell")
            f.write("\n")
                
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("RS-IRR-Data")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        f.write("Scanrates in mV/s")   
        f.write("\t")
        f.write("Root of Scanrates in (V/s)^0.5")   
        f.write("\t")
        f.write("I-Peak in microampere")   
        f.write("\t")
        f.write("Fit I-Peak in microampere")   
        f.write("\t")
        f.write("Nat. Logarithm of Scanrates")   
        f.write("\t")
        f.write("E-Peak in V")   
        f.write("\t")
        f.write("Fit E-Peak in V")   
        f.write("\t")
        f.write("E-HalfPeak in V")   
        f.write("\t")
        f.write("beta-Shape-Plot")   
        f.write("\t")
        f.write("Ieq-Randles-Plot")   
        f.write("\t")
        f.write("Ieq-LN-Plot")   
        f.write("\t")
        f.write("Ieq-Shape-Plot")   
        f.write("\t")
        if KZeroCalcer ==1:
            f.write("K_Zero-Randles-Plot")   
            f.write("\t")
            f.write("K_Zero-LN-Plot")   
            f.write("\t")
            f.write("K_Zero-Shape-Plot")   
            f.write("\t")
        f.write("\n")
        
        
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(ScanratesArray[i])))
            f.write("\t")
            f.write(str(np.asscalar((WurzelScanratesarray[i]))))
            f.write("\t")
            f.write(str(np.asscalar(I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(FITTED_I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(LnScanratesarray[i])))
            f.write("\t")
            f.write(str(np.asscalar(PeakPotentiale[i])))
            f.write("\t")
            f.write(str(np.asscalar(0.001*FIT_PeakPotentiale[i])))
            f.write("\t")
            f.write(str(np.asscalar(HalbPeakPotenziale[i])))
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_BETAARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(RS_IEQARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(LNPLOT_IEQARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_IEQARRAY[i])))
            f.write("\t")
            if KZeroCalcer ==1:
                f.write(str(np.asscalar(RS_KNULLARRAY[i])))
                f.write("\t")
                f.write(str(np.asscalar(LNPLOT_KNULLARRAY[i])))
                f.write("\t")
                f.write(str(np.asscalar(SHAPE_KNULLARRAY[i])))
                f.write("\t")
            
            f.write("\n")
            

        
        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Scanrate in mV/s")
        f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            f.write(str(np.asscalar(ScanratesArray[jj])))
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
        f.write("Scanrate in mV/s")
        f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            f.write(str(np.asscalar(ScanratesArray[jj])))
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

    root.destroy()



# In[ ]:

def RSIRRCONCASTXTsaver():
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")

        f.write("Beta-RS")
        f.write("\t")
        f.write(str(np.asscalar(RS_BETACONC)))
        f.write("\n")
        f.write("Average Beta-Shape-Plot")
        f.write("\t")
        f.write(str(np.asscalar(SHAPE_BETACONC)))
        f.write("\n")
        f.write("Average IEQ-RS in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(RS_IEQCONC)))
        f.write("\n")
        f.write("Average IEQ-Shape-Plot in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(SHAPE_IEQCONC)))
        f.write("\n")
        if KZeroCalcer ==1:
            f.write("Average KZero-RS in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(RS_KNULLCONC)))
            f.write("\n")
            f.write("Average KZero-Shape-Plot in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_KNULLCONC)))
            f.write("\n")
        f.write("\n")
       

        
        
        f.write("Parameters")
        f.write("\n")
        f.write("D in cm^2/s")
        f.write("\t")
        f.write(str(D))
        f.write("\n")
        f.write("A in cm^2")
        f.write("\t")
        f.write(str(A))
        f.write("\n")
        f.write("T in K")
        f.write("\t")
        f.write(str(T))
        f.write("\n")
        f.write("Highest c in mol/cm^3")
        f.write("\t")
        f.write(str(c))
        f.write("\n")
        f.write("Refining Factor")
        f.write("\t")
        f.write(str(REFF))
        f.write("\n")
        f.write("Ohmic Resistance in Ohm")
        f.write("\t")
        f.write(str(Ru))
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
            f.write("\t")
            f.write("(microams*s)/mV")
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
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("RS-IRR-Data")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        f.write("Concentrations in mol/cm^3")   
        f.write("\t")
        f.write("I-Peak in microampere")   
        f.write("\t")
        f.write("Fit I-Peak in microampere")   
        f.write("\t")
        f.write("E-Peak in V")    
        f.write("\t")
        f.write("E-HalfPeak in V")   
        f.write("\t")
        f.write("beta-Shape-Plot")   
        f.write("\t")
        f.write("Ieq-Randles-Plot")   
        f.write("\t")
        f.write("Ieq-Shape-Plot")   
        f.write("\t")
        if KZeroCalcer ==1:
            f.write("K_Zero-Randles-Plot")   
            f.write("\t")
            f.write("K_Zero-Shape-Plot")   
            f.write("\t")
        f.write("\n")
        
        
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(CDecayArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(FITTED_I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(PeakPotentiale[i])))
            f.write("\t")
            f.write(str(np.asscalar(HalbPeakPotenziale[i])))
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_BETAARRAYCONC[i])))
            f.write("\t")
            f.write(str(np.asscalar(RS_IEQARRAYCONC[i])))
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_IEQARRAYCONC[i])))
            f.write("\t")
            if KZeroCalcer ==1:
                f.write(str(np.asscalar(RS_KNULLARRAYCONC[i])))
                f.write("\t")
                f.write(str(np.asscalar(SHAPE_KNULLARRAYCONC[i])))
                f.write("\t")
            
            f.write("\n")
            

        
        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Concentrations as fract. of highest")
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
        f.write("Scanrate in mV/s")
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

    root.destroy()


# In[ ]:



