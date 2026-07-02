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

def Open_Cottrell_File():
    root = Toplevel()
    root.title("Your Data")
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
    global Time
    
    
            
    
    if FirstComesxx  == 1:
        Potenzial     =  data[FromRowxx:ToRowxx:Readeveryxx,0::2] * UmrPot            
        Strom         =  data[FromRowxx:ToRowxx:Readeveryxx,1::2] * UmRStrom
    if FirstComesxx  == 0:
        Potenzial     =  data[FromRowxx:ToRowxx:Readeveryxx,1::2] * UmrPot
        Strom         =  data[FromRowxx:ToRowxx:Readeveryxx,0::2] * UmRStrom
            
    #==================================================================================
    
    global Weite
    Weite = Potenzial.shape[0]
    global NumMeas
    NumMeas = Potenzial.shape[1]
    
    for i in range(Weite):
        Stromarrays             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))                          
        Potenzialarrays         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:])) 
        LenPotArrays            = len(Potenzialarrays)
    
        
        b.plot(Potenzial, 0.001*Strom, linestyle='-',marker='',color='k')
        b.set_xlabel('t / s', fontsize=12)
        b.set_ylabel('I / mA', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b.spines[axis].set_linewidth(2)
            b.spines[axis].set_color('k')
        b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
   
    #=========================================================  
    
         
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    CottrellWindow()




def Get_Cottrell_Data():
    Fenster = Toplevel()                                                         
    Fenster.title("Get-Data for Cottrell Analysis")                         
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

    UmrPot_Label = Label(Fenster,text="Time Factor to be seconds")
    UmrPot_Label.grid(row=4, column=0)
    UmrPot_Eingabe = Entry(Fenster)
    UmrPot_Eingabe.grid(row=4, column=1)
    
    #Desired Reaction
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
    Checkbutton(Fenster, text="Time", variable=var3).grid(row=6,column=1, sticky=W)
    var4 = IntVar()
    Checkbutton(Fenster, text="Current", variable=var4).grid(row=6,column=2, sticky=W)
    
    Delimiterxxx_Label  = Label(Fenster,text="Delimiter")
    Delimiterxxx_Label.grid(row=7, column=0)
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Tab", variable=var5).grid(row=7,column=1, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Space", variable=var6).grid(row=7,column=2, sticky=W)

    def AcceptCottrell():
              
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
    
    
    def NextCottrell():
        Open_Cottrell_File()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    Accept = Button(Fenster, text="Accept",command=AcceptCottrell)
    Accept.grid(row=8, column=0) 
    
    Next = Button(Fenster, text="Next",command=NextCottrell)
    Next.grid(row=9, column=0)  


# In[ ]:

def CottrellWindow():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Cottrell-Analysis")                         
    #Fenster.geometry("150x200")
    

    var3 = IntVar()
    Checkbutton(Fenster, text="Baseline Correction", variable=var3).grid(row=1,column=0, sticky=W)
    
    c_Label = Label(Fenster,text="c in mol/cm^3")
    c_Label.grid(row=2, column=0)
    
    Concents = []

    for i in range(NumMeas):
        
        en = Entry(Fenster)
        en.grid(row=i+3, column=0)
        Concents.append(en)    
        
    
    def AcceptParams():

        global BaseCorr
        
        Concs= []
        for entry in Concents:
            Concs.append(float(entry.get()))
        global ConcentrationsArray
        ConcentrationsArray = np.array(Concs) 
        
        BaseCorr = var3.get()
        
    def Next():
        
        def quit():
            Fenster.destroy()
        quit() 
        
        if BaseCorr == 0:
            CottrellWindowLevel2()
        if BaseCorr == 1:
            BaseGetter()
            CottrellWindowLevel2()
        
        
    
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=30,column=0)
    button=Button(Fenster,text="Next",command=Next).grid(row=31,column=0)


# In[ ]:

def CottrellWindowLevel2():
    Fenster = Toplevel()                                                         
    Fenster.title("Cottrell-Analysis")                         
    #Fenster.geometry("1200x600")
    
    
    f = Figure(figsize=(12, 6), dpi=100)

    global ShortLength
    ShortLength = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:]))))
    
    global ZEITSUPERARRAY
    global INVWURZZEITSUPPERARRAY
    global UNKORRSTROMSUPPERARRAY
    global STROMSUPPERARRAY
    
    ZEITSUPERARRAY             = np.empty([NumMeas,ShortLength])
    INVWURZZEITSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    UNKORRSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    STROMSUPPERARRAY           = np.empty([NumMeas,ShortLength])
    
    
    for i in range(int(NumMeas)):
        StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
        ZeitArraysROH              = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
        InvWurzZeit                = 1/((ZeitArraysROH**0.5))
        
        if DesiReactxx == 0:
            Stromarrays            = -StromarraysROH
            
        if DesiReactxx == 1:
            Stromarrays            = StromarraysROH
   
            
        if BaseCorr == 1:
            
            BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
            if DesiReactxx == 0:
                BaseStromArrays    = -BaseStromArraysROH[::]
            if DesiReactxx == 1:
                BaseStromArrays    = BaseStromArraysROH[::]
            Stromarrays            = Stromarrays - BaseStromArrays
                                              

        for j in range(ShortLength):
                
            ZEITSUPERARRAY[i,j]          = ZeitArraysROH[j] 
            INVWURZZEITSUPPERARRAY[i,j]  = InvWurzZeit[j]
            UNKORRSTROMSUPPERARRAY[i,j]  = StromarraysROH[j]
            STROMSUPPERARRAY[i,j]        = Stromarrays[j]

        
        
    #Chronoplot
    #______________________________________________________
    Chronoplot = f.add_subplot(121)
    for i in range(NumMeas):
        Chronoplot.plot(ZEITSUPERARRAY[i],UNKORRSTROMSUPPERARRAY[i], color='k',linestyle='-',marker='')
        Chronoplot.plot(ZEITSUPERARRAY[i],STROMSUPPERARRAY[i],color='r',linestyle='-',marker='')
    Chronoplot.set_xlabel('t in [s]', fontsize=12)
    Chronoplot.set_ylabel('I [$\mu$A]', fontsize=12)
    for axis in ['top','bottom','left','right']:
        Chronoplot.spines[axis].set_linewidth(2)
    Chronoplot.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
    #Cottiplot für Fit
    #_______________________________________________________
        
    Cottiplot = f.add_subplot(122)
    for i in range(NumMeas):
        Cottiplot.plot(INVWURZZEITSUPPERARRAY[i],STROMSUPPERARRAY[i],color='r',linestyle='',marker='.')
    Cottiplot.set_xlabel(r'$\sqrt{t^{-1}}$' '[s $^{-0.5}$]''t in [s]', fontsize=12)
    Cottiplot.set_ylabel('I [$\mu$A]', fontsize=12)
    for axis in ['top','bottom','left','right']:
        Cottiplot.spines[axis].set_linewidth(2)
    Cottiplot.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
        
        
    
    canvas = FigureCanvasTkAgg(f, master=Fenster)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, Fenster)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    

    CottrellWindowLevel3()
    


# In[ ]:

def CottrellWindowLevel3():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Cottrell-Analysis")                         
    #Fenster.geometry("150x200")
    
    FitFrom_Label = Label(Fenster, text="Fit From")
    FitFrom_Label.grid(row=1, column=0)
    FitFrom_Eingabe = Entry(Fenster)                                               
    FitFrom_Eingabe.grid(row=2, column=0)

    FitTo_Label = Label(Fenster,text="Fit To")
    FitTo_Label.grid(row=3, column=0)
    FitTo_Eingabe = Entry(Fenster)
    FitTo_Eingabe.grid(row=4, column=0)
    
    var1 = IntVar()
    Checkbutton(Fenster, text="Save as txt", variable=var1).grid(row=5,column=0, sticky=W)
    
    def Next():
        
        global FitFrom
        global FitTo
        global AsTxtSaver
        
        FitFrom    = float(FitFrom_Eingabe.get())
        FitTo      = float(FitTo_Eingabe.get())
        AsTxtSaver = int(var1.get())
        
        def quit():
            Fenster.destroy()
        quit() 
        
        CottrellWindowLevel4()
        
        
        
    button=Button(Fenster,text="Next",command=Next).grid(row=6,column=0)
    


# In[ ]:

def CottrellWindowLevel4():  
    
    Fenster = Toplevel()                                                         
    Fenster.title("Cottrell-Plot")                         
    #Fenster.geometry("700x700")
    
    #DIESE ZWEI FUNKTIONEN MÜSSEN HIER DRIN DEFINIERT WERDEN!!!
        
    def GeradenFit(xWert,ySchnitt,Steig):
        return  ySchnitt + Steig*xWert

    def EchtesPotFinden(WoSieSuchenSoll,WasSieSuchenSoll):                         
        idxStart = (np.abs(WoSieSuchenSoll-WasSieSuchenSoll)).argmin()               
        return WoSieSuchenSoll[idxStart] 
    
    

    global FITSTROMSUPPERARRAY
    FITSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    global CottiSlopesArray
    CottiSlopesArray        = np.empty(int(NumMeas))
    global CottiCutsArray
    CottiCutsArray          = np.empty(int(NumMeas))
    
    
    f = Figure(figsize=(6, 6), dpi=100)

          
    
    for i in range(int(NumMeas)):
    
        FitFromPotenzial       = EchtesPotFinden(INVWURZZEITSUPPERARRAY[i,::],FitFrom)
        FitToPotenzial         = EchtesPotFinden(INVWURZZEITSUPPERARRAY[i,::],FitTo)
        
        #print FitFromPotenzial
        #print FitToPotenzial
        
        IdxFitFrom = np.asscalar(np.where(INVWURZZEITSUPPERARRAY[i,::] == FitFromPotenzial) [0])
        IdxFitTo   = np.asscalar(np.where(INVWURZZEITSUPPERARRAY[i,::] == FitToPotenzial) [0])
        
        if FitFromPotenzial<FitToPotenzial:
            CottiFitting, pcov  = curve_fit(GeradenFit, INVWURZZEITSUPPERARRAY[i,IdxFitTo:IdxFitFrom], STROMSUPPERARRAY[i,IdxFitTo:IdxFitFrom])
        
        if FitFromPotenzial>FitToPotenzial:
            CottiFitting, pcov  = curve_fit(GeradenFit, INVWURZZEITSUPPERARRAY[i,IdxFitFrom:IdxFitTo], STROMSUPPERARRAY[i,IdxFitFrom:IdxFitTo])
        
        
        CottiSlopesArray[i] = CottiFitting[1]
        CottiCutsArray[i]   = CottiFitting[0]
        
    #Fitstromarrays berechnen
    #________________________________________________________________
    for i in range(len(CottiSlopesArray)):
 
        FITSTROMSUPPERARRAY[i] = INVWURZZEITSUPPERARRAY[i,::]*CottiSlopesArray[i] + CottiCutsArray[i]       
            

    #PLOTTEN
    #_________________________________________________________________
    #Cottiplot für Fit
    #_________________________________________________________________
        
    Cottiplot = f.add_subplot(111)
    for i in range(NumMeas):
        Cottiplot.plot(INVWURZZEITSUPPERARRAY[i],STROMSUPPERARRAY[i],color='r',linestyle='',marker='.')
        Cottiplot.plot(INVWURZZEITSUPPERARRAY[i],FITSTROMSUPPERARRAY[i],color='r',linestyle='-',marker='')
    
        Cottiplot.annotate('m=%5.3f'    % CottiSlopesArray[i], xy=(0.15, 0.93 - 0.05*i), xycoords='axes fraction')
    
    
    Cottiplot.set_xlabel(r'$\sqrt{t^{-1}}$' '[s $^{-0.5}$]''t in [s]', fontsize=12)
    Cottiplot.set_ylabel('I [$\mu$A]', fontsize=12)
    for axis in ['top','bottom','left','right']:
        Cottiplot.spines[axis].set_linewidth(2)
    Cottiplot.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
        
        
    
    canvas = FigureCanvasTkAgg(f, master=Fenster)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, Fenster)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)    
    
    
    if AsTxtSaver ==1:
        CottiAsTxtSaver()
       


# In[ ]:

def CottiAsTxtSaver():
    
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        f.write("c in mol/cm^3")
        f.write("\t")
        f.write("Cottrellslope in (µA cm^3)/(s^0.5 mol)")
        f.write("\n")
        f.write("\n")
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(ConcentrationsArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(CottiSlopesArray[i])))
            f.write("\t")
            f.write("\n")
        
        f.write("\n")
        
        f.write("Correction type")
        f.write("\t")
        if BaseCorr == 0:
            f.write("none")
        if BaseCorr == 1:
            f.write("Base")
        
        f.write("\n")
        f.write("Fit from")
        f.write("\t")
        f.write(str((FitFrom)))
        f.write("\t")
        f.write("\n")
        f.write("Fit to")
        f.write("\t")
        f.write(str((FitTo)))
        f.write("\t")
        
        f.write("\n")
        f.write("\n")
        
        f.write("Data")
        f.write("\n")
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(ConcentrationsArray[i])))
            f.write("\t")
            f.write("\t")
            f.write("\t")
            f.write("\t")
        f.write("\n")
        
        for i in range(int(NumMeas)):
            f.write("t in s")
            f.write("\t")
            f.write(" 1/t^0.5 in 1/s^0.5")
            f.write("\t")
            f.write("I in microamps")
            f.write("\t")
            f.write("Fit-I in microamps")
            f.write("\t")
        f.write("\n")
        
        for j in range(int(ShortLength)):
            for i in range(int(NumMeas)):
                f.write(str(np.asscalar(ZEITSUPERARRAY[i,j])))
                f.write("\t")
                f.write(str(np.asscalar(INVWURZZEITSUPPERARRAY[i,j])))
                f.write("\t")
                f.write(str(np.asscalar(STROMSUPPERARRAY[i,j])))
                f.write("\t")
                f.write(str(np.asscalar(FITSTROMSUPPERARRAY[i,j])))
                f.write("\t")
            
            f.write("\n")
            

    root.destroy()


# In[ ]:



