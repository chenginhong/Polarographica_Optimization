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


# In[ ]:

def Open_RSREV_File():
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
    
    
    RandlesRevWindow()
    


# In[ ]:

def Get_RSREV_Data():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Get-Randles-Sevcik File")                         
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



    def AcceptRSREV():
              
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
    
    
    def NextRSREV():
        Open_RSREV_File()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    Accept = Button(Fenster, text="Accept",command=AcceptRSREV)
    Accept.grid(row=8, column=0) 
    
    Next = Button(Fenster, text="Next",command=NextRSREV)
    Next.grid(row=9, column=0)  
    


# In[ ]:

def RandlesRevWindow():
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
        RS_REV_NextLevel2()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    button=Button(Fenster,text="Accept OCPs",command=GetOCPs).grid(row=82,column=0)
    button=Button(Fenster,text="Accept Scanrates",command=GetScanrates).grid(row=82,column=1)
    button=Button(Fenster,text="Accept Concentrations",command=GetConcents).grid(row=82,column=2)
    button=Button(Fenster,text="Next",command=Next).grid(row=83,column=1)

    


# In[ ]:

def RS_REV_NextLevel2():
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Reversible Parameter-Getter")                         
    Fenster.geometry("400x400")
    
    
    Getter_Label = Label(Fenster, text="Get From RS Rev")
    Getter_Label.grid(row=0, column=0)
    var1 = IntVar()
    Checkbutton(Fenster, text="get n", variable=var1).grid(row=0, column=1, sticky=W)
    var2 = IntVar()
    Checkbutton(Fenster, text="get c", variable=var2).grid(row=0,column=2, sticky=W)
    var3 = IntVar()
    Checkbutton(Fenster, text="get D", variable=var3).grid(row=0,column=3, sticky=W)
    
        
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
        global get_n
        global get_c
        global get_D
        
    
        get_n = var1.get()
        get_c = var2.get()
        get_D = var3.get()
        
        if get_n ==1:
            if get_c ==0:
                if get_D ==0:
              
                    #n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    D         = (float(D_Eingabe.get()))
                    c         = (float(c_Eingabe.get()))
        
        if get_n ==0:
            if get_c ==1:
                if get_D ==0:
              
                    n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    D         = (float(D_Eingabe.get()))
                    #c         = (float(c_Eingabe.get()))
        
        if get_n ==0:
            if get_c ==0:
                if get_D ==1:
              
                    n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    #D         = (float(D_Eingabe.get()))
                    c         = (float(c_Eingabe.get()))
    
    
    
    def Next():
        RS_REV_NextLevel3()
        
        def quit():
            Fenster.destroy()
        quit()   
        
    def Back():
        RandlesRevWindow()
        
        def quit():
            Fenster.destroy()
        quit()   
    
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=8,column=0)
    button=Button(Fenster,text="Next",command=Next).grid(row=9,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=9,column=0)
    


# In[ ]:

def RS_REV_NextLevel3():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Reversible Parameter-Getter-2")                         
    Fenster.geometry("400x400")
    
    
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
    
    
    if get_n == 1:
    
    
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
    ManuCorr_Label    = Label(Fenster,text="Cap. Curr in microamp.")
    ManuCorr_Label.grid(row=12, column=0)
    ManuCorr_Eingabe  = Entry(Fenster)
    ManuCorr_Eingabe.grid(row=12, column=1)
    
    var8 = IntVar()
    Checkbutton(Fenster, text="Base Correction", variable=var8).grid(row=10, column=1, sticky=W)
    var9 = IntVar()
    Checkbutton(Fenster, text="Separate LSV-Plot", variable=var9).grid(row=13, column=0, sticky=W)
    var10 = IntVar()
    Checkbutton(Fenster, text="Show Uncorr LSV", variable=var10).grid(row=14, column=0, sticky=W)
    var11 = IntVar()
    Checkbutton(Fenster, text="Save Results as Txt", variable=var11).grid(row=14, column=1, sticky=W)
    

    def AcceptParams():
               
        global REFF
        global Ru 
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
        global Instead_n_D
        global Instead_n_c
        global SeparateLSV
        global BaseCorr
        global ManuCorr
        global ManuCorrCurr
        global ShowUncorrLSV
        global AsTxtSaver
        
        
        REFF            = (float(REFF_Eingabe.get()))
        Ru              = (float(Ru_Eingabe.get()))


        Smooth_Data     = var4.get()
        AutoCorr        = var5.get()
        OnlyRSPlot      = var6.get()
        ManuCorr        = var7.get()
        BaseCorr        = var8.get()
        SeparateLSV     = var9.get()
        ShowUncorrLSV   = var10.get()
        AsTxtSaver      = var11.get()
        
        
        
        if ManuCorr == 1:
            ManuCorrCurr   =(float(ManuCorr_Eingabe.get()))
        
        
        if  Smooth_Data == 1:
            SmoothFact  = (float(SmoothFact_Eingabe.get()))
            SmoothFrom  = (float(SmoothFrom_Eingabe.get()))
            SmoothTo    = (float(SmoothTo_Eingabe.get()))
        
        
        if get_n == 1:
            
            RefToCotti = var1.get()
            
            if RefToCotti  == 1:
                CottiSlope  = (float(CottiSlope_Eingabe.get()))

                
    def Next():
        
        if BaseCorr == 0:
            RS_REV_NextLevel4()
        if BaseCorr == 1:
            BaseGetter()
            RS_REV_NextLevel4()
             
        
    def Back():
        RS_REV_NextLevel2()
        
        def quit():
            Fenster.destroy()
        quit()   
    

    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=15,column=1)
    button=Button(Fenster,text="Next",command=Next).grid(row=16 ,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=17,column=1)
    


# In[ ]:

def RS_REV_NextLevel4():
    Fenster = Toplevel()                                                         
    Fenster.title("Reversible Randles-Sevcik Analysis")                         
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
    global I_Peak
    global FITTED_I_Peak
    global WurzelScanratesarray
    
    
    KORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,Length])
    KORRSTROMSUPPERARRAY     = np.empty([NumMeas,Length])
    
    UNKORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,ShortLength])
    UNKORRSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
 
    
    f = Figure(figsize=(12, 6), dpi=100)
    
    if SeparateLSV == 1:
        LSVFenster = Toplevel()                                                         
        LSVFenster.title("LSV-Plot")                         
        LSVFenster.geometry("700x700")
        
        LSVPlot = Figure(figsize=(6, 6), dpi=100)
    
    #f.subplots_adjust(hspace=0.4)
    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #FOR SCHLEIFE WENN ROT RATES GEÄNDERT WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    
    if VarParSca == 1:
    
        OCPs = OCPArray
        
        #________________________________________________________________      
        
               
        I_Peak    = np.empty(int(NumMeas))
        
        
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
            KorrAufgPotArr         = AufgelPotArrays #Ohmsche Korrektur erfolgte schon weiter oben
            
 
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
                UNKORRSTROMSUPPERARRAY[i,j]     = Stromarrays[j]
            
            
            if DesiReactxx == 0:
                PeakMax    = np.min(KorrAufgStrArrays) 
            
            if DesiReactxx == 1:
                PeakMax    = np.max(KorrAufgStrArrays) 
            
            

    
            I_Peak[i] = PeakMax
        
            
        
            bild1 = f.add_subplot(121)
            if ShowUncorrLSV ==1:
                bild1.plot(UnkorrPotenzialarrays,Stromarrays, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            
            if SeparateLSV == 1:
                LSVbild1 = LSVPlot.add_subplot(111)
                if ShowUncorrLSV ==1:
                    LSVbild1.plot(UnkorrPotenzialarrays,Stromarrays, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [$\mu$A]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
        
        
        

        WurzelScanratesarray = (ScanratesArray*0.001)**0.5

        RSREVFitting, pcov = curve_fit(GeradenFit, WurzelScanratesarray,I_Peak)   
        
        FITTED_I_Peak = GeradenFit(WurzelScanratesarray, *RSREVFitting)
        
        bild2 = f.add_subplot(122)
        
        def RSREVSCPLOTTER():
        
            bild2.plot(WurzelScanratesarray,GeradenFit(WurzelScanratesarray, *RSREVFitting),color='r',linestyle='-',marker='')
            bild2.plot(WurzelScanratesarray,I_Peak,linestyle='',marker='.')
            bild2.set_xlabel(r'$\sqrt{\nu} $' ' ' r'[$\sqrt{V} $]' , fontsize=12)
            bild2.set_ylabel('I$_p$ [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
            #________________________________________________________
            #________________________________________________________
            #ANNOTOATIONS FOR PARAMETERS AND CALCULATION OF n D or c
            #________________________________________________________
            #________________________________________________________
        
        
            #IF SCHLEIFE WENN n BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
            if get_n ==1:
         
                global n_RS
                n_RS = (np.absolute(RSREVFitting[1]*0.000001*(R*T)**0.5 / (0.4463*F*A*c*(F*D)**0.5)))**0.666666667

                if RefToCotti == 1:
               
                    BeschreibefensterCotti = Toplevel()                                                         
                    BeschreibefensterCotti.title("Information")                         
                    BeschreibefensterCotti.geometry("200x200")    
                    Te = Text(BeschreibefensterCotti, height=6, width=30)
                    Te.pack()
                    Te.insert(END, "Randles Slope got\ndivided by Cottrell Slope")
                    
                    
                
                    n_RS = np.absolute((RSREVFitting[1]*(R*T)**0.5)/(np.pi**0.5 * 0.4463 * CottiSlope * F**0.5))**2
      
                
                if DesiReactxx == 0:
                    bild2.annotate('n     =%8.3f' % n_RS, xy=(0.5, 0.93), xycoords='axes fraction')     
                if DesiReactxx == 1:
                    bild2.annotate('n     =%8.3f' % n_RS, xy=(0.5, 0.13), xycoords='axes fraction')  
        
        
            #IF SCHLEIFE WENN D BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_D ==1:
            
                global D_RS
                D_RS = 1000000*(0.000001*RSREVFitting[1]*(R*T)**0.5 / (0.4463*n*F*A*c*(n*F)**0.5))**2
            
                
                if DesiReactxx == 0:
                    bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_RS, xy=(0.5, 0.93), xycoords='axes fraction') 
                if DesiReactxx == 1:
                    bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_RS, xy=(0.5, 0.13), xycoords='axes fraction') 
                

            #IF SCHLEIFE WENN c BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_c ==1:
            
                global c_RS
                c_RS = np.absolute(1000000*(RSREVFitting[1]*0.000001*(R*T)**0.5) / (0.4463*n*F*A*(n*F*D)**0.5))
            
                if DesiReactxx == 0:
                    bild2.annotate('c =%5.3f*10$^{-6} mol/cm^{3}$' % c_RS, xy=(0.5, 0.93), xycoords='axes fraction') 
                if DesiReactxx == 1:
                    bild2.annotate('c =%5.3f*10$^{-6} mol/cm^{3}$' % c_RS, xy=(0.5, 0.13), xycoords='axes fraction') 
                

        
        RSREVSCPLOTTER()

        
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
        #ONLY RS PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyRSPlot == 1:
            
            RSFenster = Toplevel()                                                         
            RSFenster.title("Reversible Randles-Sevcik Plot")                         
            RSFenster.geometry("700x700")
        
            RSREVPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = RSREVPlot.add_subplot(111)
            
            RSREVSCPLOTTER()
                
            canvas = FigureCanvasTkAgg(RSREVPlot, master=RSFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, RSFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        if SeparateLSV == 1:
            canvas = FigureCanvasTkAgg(LSVPlot, master=LSVFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, LSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        if AsTxtSaver ==1:
            RSREVSCASTXTsaver()
    
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
    
        OCPs = OCPArray
        
        #________________________________________________________________      
        
               
        I_Peak    = np.empty(int(NumMeas))
        

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
            KorrAufgPotArr         = AufgelPotArrays #Ohmsche Korrektur erfolgte schon weiter oben
            
 
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
                UNKORRSTROMSUPPERARRAY[i,j]     = Stromarrays[j]
            
            
            if DesiReactxx == 0:
                PeakMax    = np.min(KorrAufgStrArrays) 
            
            if DesiReactxx == 1:
                PeakMax    = np.max(KorrAufgStrArrays) 
            

            I_Peak[i] = PeakMax
        

            bild1 = f.add_subplot(121)
            if ShowUncorrLSV ==1:
                bild1.plot(UnkorrPotenzialarrays,Stromarrays, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            
        
            if SeparateLSV == 1:
                LSVbild1 = LSVPlot.add_subplot(111)
                if ShowUncorrLSV ==1:
                    LSVbild1.plot(UnkorrPotenzialarrays,Stromarrays, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [$\mu$A]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
                
                
                
                
                

        Concentrationsarray      = ConcentrationsArray * c
        

        RSREVFitting, pcov = curve_fit(GeradenFit, Concentrationsarray,I_Peak)   
        
        cMax = c*1000000
        
        FITTED_I_Peak = GeradenFit(Concentrationsarray, *RSREVFitting)
        
        bild2 = f.add_subplot(122)
        
        def RSREVCONCPLOTTER():
        
            bild2.plot(ConcentrationsArray,GeradenFit(Concentrationsarray, *RSREVFitting),color='r',linestyle='-',marker='')
            bild2.plot(ConcentrationsArray,I_Peak,linestyle='',marker='.')
            bild2.set_xlabel('fract.(c$_{max}$)', fontsize=12)
            bild2.set_ylabel('I$_p$ [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
            if DesiReactxx == 0:
                bild2.annotate( 'c$_{max}$ =%5.3f*10$^{-6} mol/cm^{3}$' % cMax , xy=(0.5, 0.83), xycoords='axes fraction')
            if DesiReactxx == 1:
                bild2.annotate( 'c$_{max}$ =%5.3f*10$^{-6} mol/cm^{3}$' % cMax , xy=(0.5, 0.23), xycoords='axes fraction')          
            
        
        
            #________________________________________________________
            #________________________________________________________
            #ANNOTOATIONS FOR PARAMETERS AND CALCULATION OF n D or c
            #________________________________________________________
            #________________________________________________________
        
            
            #IF SCHLEIFE WENN n BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
            if get_n ==1:
         
                global n_RS
                n_RS = (np.absolute(RSREVFitting[1]*0.000001*(R*T)**0.5 / (0.4463*F*A*((Scanrate*0.001)**0.5)*(F*D)**0.5)))**0.666666667

                if RefToCotti == 1:
               
                    BeschreibefensterCotti = Toplevel()                                                         
                    BeschreibefensterCotti.title("Information")                         
                    BeschreibefensterCotti.geometry("200x200")    
                    Te = Text(BeschreibefensterCotti, height=6, width=30)
                    Te.pack()
                    Te.insert(END, "Randles Slope got\ndivided by Cottrell Slope")
                    
                    
                
                    n_RS = (np.absolute(RSREVFitting[1]*c*(R*T)**0.5/(0.4463*np.pi**0.5 *(Scanrate*0.001)**0.5 *CottiSlope *F**0.5)))**2
      
                if DesiReactxx == 0:
                    bild2.annotate('n     =%8.3f' % n_RS, xy=(0.5, 0.93), xycoords='axes fraction')     
                if DesiReactxx == 1:
                    bild2.annotate('n     =%8.3f' % n_RS, xy=(0.5, 0.13), xycoords='axes fraction')           
        
        
            #IF SCHLEIFE WENN D BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_D ==1:
            
                global D_RS
                D_RS = 1000000*(0.000001*RSREVFitting[1]*(R*T)**0.5 / (0.4463*n*F*A*(n*F*Scanrate*0.001)**0.5))**2
            
                if DesiReactxx == 0:
                    bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_RS, xy=(0.5, 0.93), xycoords='axes fraction') 
                if DesiReactxx == 1:
                    bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_RS, xy=(0.5, 0.13), xycoords='axes fraction') 
        
            #IF SCHLEIFE WENN c BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_c ==1:

                if DesiReactxx == 0:
                    bild2.annotate('GETTING c NOT POSSIBLE', xy=(0.5, 0.93), xycoords='axes fraction') 
                if DesiReactxx == 1:
                    bild2.annotate('GETTING c NOT POSSIBLE', xy=(0.5, 0.13), xycoords='axes fraction')  
        
        
        
        
        
        RSREVCONCPLOTTER()

        
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
        #ONLY RSREV PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyRSPlot == 1:
            
            RSFenster = Toplevel()                                                         
            RSFenster.title("Reversible Randles-Sevcik Plot")                         
            RSFenster.geometry("700x700")
        
            RSREVPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = RSREVPlot.add_subplot(111)
            
            RSREVCONCPLOTTER()
                
            canvas = FigureCanvasTkAgg(RSREVPlot, master=RSFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, RSFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        if SeparateLSV == 1:
            canvas = FigureCanvasTkAgg(LSVPlot, master=LSVFenster)
            canvas.draw()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2Tk(canvas, LSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        
        if AsTxtSaver ==1:
            RSREVCONCASTXTsaver()
      


# In[1]:

def RSREVSCASTXTsaver():
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        if get_n == 1:
            f.write("n-RS-Rev")
            f.write("\t")
            f.write(str(np.asscalar(n_RS)))
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
                
            
  

            
        if get_c == 1:
            f.write("c-RS-REV")
            f.write("\t")
            f.write(str(0.000001*np.asscalar(c_RS)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
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
                
                
                
        if get_D == 1:
            f.write("D-RS-REV")  
            f.write("\t")
            f.write(str(0.000001*np.asscalar(D_RS)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
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
        
        f.write("RS-REV-Plot Data")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        f.write("Scanrates in mV/s")   
        f.write("\t")
        f.write("Root of Scanrates in (mV/s)^0.5")   
        f.write("\t")
        f.write("I-Peak in microampere")   
        f.write("\t")
        f.write("Fit I-Peak in microampere")   
        f.write("\t")
        f.write("\n")
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(ScanratesArray[i])))
            f.write("\t")
            f.write(str(np.asscalar((ScanratesArray[i])**0.5)))
            f.write("\t")
            f.write(str(np.asscalar(I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(FITTED_I_Peak[i])))
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

def RSREVCONCASTXTsaver():
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        if get_n == 1:
            f.write("n-RS-Rev")
            f.write("\t")
            f.write(str(np.asscalar(n_RS)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("Scanrate in mV/s")
            f.write("\t")
            f.write(str(Scanrate))
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
            f.write("c_max in mol/cm^3")
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
                
            
  

            
        if get_c == 1:
            f.write("c-RS-REV")
            f.write("\t")
            f.write("getting c not possible here")
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("Scanrate in mV/s")
            f.write("\t")
            f.write(str(Scanrate))
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
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
                
                
                
        if get_D == 1:
            f.write("D-RS-REV")  
            f.write("\t")
            f.write(str(0.000001*np.asscalar(D_RS)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("Scanrate in mV/s")
            f.write("\t")
            f.write(str(Scanrate))
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
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
        
        f.write("RS-REV-Plot Data")
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
        f.write("\n")
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(c*ConcentrationsArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(FITTED_I_Peak[i])))
            f.write("\t")
            f.write("\n")
            

        
        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Concentr. as Fract. of highest")
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
        f.write("Concentr. as Fract. of highest")
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

