

#==========================================================================================================
#==========================================================================================================
#Import all the required modules
#==========================================================================================================
#==========================================================================================================
import        numpy as np
import        matplotlib.pyplot as plt
from          tqdm import tqdm
import        importlib
import        tkinter               as tk
import        time
from          tkinter               import *
from          tkinter               import ttk
from          tkinter.scrolledtext  import ScrolledText




def CDCwrite(Text): 
    CDCtext_area.insert(INSERT, "%s" %Text )
 

def CDCMainConsole(NCycles, I, ECutLow, ECutUp, 
                         alphapos, kzeropos, alphaneg, kzeroneg, FlowRate, 
                         Poro, lenz, leny, lenx, Incz, CurrFindIter, 
                         Res, Rincr, credneg, coxneg, credpos, coxpos, ActRetainneg, ActRetainpos, T,
                         Ezeropos, Ezeroneg, Area, TankVolneg, TankVolpos):
    global CDCMainwin
    global CDCtext_area
    CDCMainwin = tk.Toplevel() 
    CDCMainwin.geometry('%dx%d+0+0' % (600,475))
    CDCMainwin.title("Charge/Discharge Console") 
    def FINALGO():
        CYCLELOPPERFUNC_CDC(NCycles, I, ECutLow, ECutUp, 
                         alphapos, kzeropos, alphaneg, kzeroneg, FlowRate, 
                         Poro, lenz, leny, lenx, Incz, CurrFindIter, 
                         Res, Rincr, credneg, coxneg, credpos, coxpos, ActRetainneg, ActRetainpos, T,
                         Ezeropos, Ezeroneg, Area, TankVolneg, TankVolpos)
    # Title Label 
    ttk.Label(CDCMainwin,  text = "Charge/Dischare Console", font = ("Times New Roman", 15),  background = 'green',  foreground = "white").place(x = 75, y = 10, width = 300, height = 25)
    CDCtext_area = tk.scrolledtext.ScrolledText(CDCMainwin,  wrap = tk.WORD,  font = ("Times New Roman", 9)) 
    CDCtext_area.place(x = 75, y = 100,  width = 450,  height = 300) 
    # Placing cursor in the text area 
    CDCtext_area.focus() 
    button=Button(CDCMainwin,text="Start Computations",bg = '#9BA9C5', command = FINALGO).place(x = 75, y = 40, width = 300, height = 40)
    #button=Button(CDCMainwin,text="Start Computations",bg = '#9BA9C5', command = CDClooper).place(x = 75, y = 40, width = 300, height = 40)

      
    


def CDC_FirstLevel():
    Fenster = tk.Toplevel()                                                         
    Fenster.title("Rudimentary Charge/Discharge")                         
    Fenster.geometry("500x325")
    colorbgr = Label(Fenster, text= "", bg = '#FF8566')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000)
    Te = Text(Fenster, font = ('Arial',12,'bold'), padx = 20, pady = 10)
    Te.place(x = 50, y = 50, width = 400, height = 175)
    Te.insert(END, "CAUTION! This is a rudimentary version of \nthe charge/discharge function which is \nnot fully checked. Despite there is\nan entry for different values of n (transferred\nelectrons) this field is not active and always \nn=1. Imporved versions will follow in next \nversions of Polarographica.")
    def NoMatterGo():
        Fenster.destroy()
        ChargeDischargeWithOrWithoutNumba()
    button=Button(Fenster,text="No matter, start!", bg = '#9BA9C5', command = NoMatterGo).place(x = 150, y = 250, width = 200, height = 50)
    Fenster.resizable(False, False)

 


def ChargeDischargeWithOrWithoutNumba():
    NumbaExister = importlib.util.find_spec('numba')
    if NumbaExister == None:
        VerySlow_CDC()
    if NumbaExister is not None:
        InputParamsChargeDischarge()
        



def VerySlow_CDC():
    Fenster = tk.Toplevel()                                                         
    Fenster.title("Painfully slow warner")                         
    WaitingDuck_image      = tk.PhotoImage(file = r"IMAGES\Waitingduck.png")       #Image as background in toplevel
    background_label       = tk.Label(Fenster, image=WaitingDuck_image)            #Image as background in toplevel
    background_label.image = WaitingDuck_image                                     #Image as background in toplevel
    w                      = WaitingDuck_image.width()                             #Image as background in toplevel
    h                      = WaitingDuck_image.height()                            #Image as background in toplevel
    xcoodinate             = (Fenster.winfo_screenwidth()  -  w)/2                #Pop up window in center of screen
    ycoodinate             = (Fenster.winfo_screenheight()  - h)/2                #Pop up window in center of screen
    Fenster.geometry('%dx%d+%d+%d' % (w,h, xcoodinate, ycoodinate))
    background_label.place(x=0, y=0, relwidth=1, relheight=1)
    Te = Text(Fenster, font = ('Arial',12,'bold'), padx = 20, pady = 10)
    Te.place(x = 320, y = 38, width = 250, height = 320)
    Te.insert(END, "WARNING! You do NOT \nhave installed the module \nNumba. Your computation\nwill take ages! \nWe highly recommend to \ninstall numba. If you are \nusing the Anaconda \npackage manager, you \ncan go for the command \n\nconda install numba \n\nMore info at\n\nhttps://numba.pydata.org/ ")
    def RunWithoutNumba():
        Fenster.destroy()
        InputParamsChargeDischarge()
    button=Button(Fenster,text="Continue without\nnumba\n(not recommended)", bg = '#9BA9C5', command = RunWithoutNumba).place(x = 320, y = 385, width = 250, height = 75)
    Fenster.resizable(False, False)




def InputParamsChargeDischarge():
    Fenster_CDC_in = tk.Toplevel()                                                         
    Fenster_CDC_in.title("Input for Charge-Discharge")   
    Fenster_CDC_in.geometry("900x500")
    colorbgr = Label(Fenster_CDC_in, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000)
    Fenster_CDC_in.resizable(False, False)
    
    #=========================================================================================================================
    #     NEGATIVE SIDE INPUTS
    #=========================================================================================================================
    neg_Label = Label(Fenster_CDC_in, text="Inputs for negative side", bg = '#80ffbf', font = ('Arial',16,'bold'))
    neg_Label.place(x = 25, y = 20)
    
    nneg_Label = Label(Fenster_CDC_in, text="n... keep it 1 ", bg = '#80ffbf')
    nneg_Label.place(x = 25, y = 65)
    nneg_Eingabe = Entry(Fenster_CDC_in)  
    nneg_Eingabe.insert(END, 1)
    #nneg_Eingabe.pack()                                             
    nneg_Eingabe.place(x = 150, y = 65)
    
    alphaneg_Label = Label(Fenster_CDC_in, text= u"\u03B1*",bg = '#80ffbf')
    alphaneg_Label.place(x = 25, y = 90)
    alphaneg_Eingabe = Entry(Fenster_CDC_in)                                               
    alphaneg_Eingabe.insert(END, 0.5)
    #alphaneg_Eingabe.pack()                                             
    alphaneg_Eingabe.place(x = 150, y = 90)
    
    kzeroneg_Label = Label(Fenster_CDC_in,text="k°[cm/s]*",bg = '#80ffbf')
    kzeroneg_Label.place(x = 25, y = 115)
    kzeroneg_Eingabe = Entry(Fenster_CDC_in)
    kzeroneg_Eingabe.insert(END, 1e-8)
    #kzeroneg_Eingabe.pack()    
    kzeroneg_Eingabe.place(x = 150, y = 115)
    
    credneg_Label = Label(Fenster_CDC_in,text="c_red [mol/L]*",bg = '#80ffbf')
    credneg_Label.place(x = 25, y = 140)
    credneg_Eingabe = Entry(Fenster_CDC_in)
    credneg_Eingabe.insert(END, 0.99)
    #credneg_Eingabe.pack() 
    credneg_Eingabe.place(x = 150, y = 140) 
    
    coxneg_Label = Label(Fenster_CDC_in,text="c_ox [mol/L]*",bg = '#80ffbf')
    coxneg_Label.place(x = 25, y = 165)
    coxneg_Eingabe = Entry(Fenster_CDC_in)
    coxneg_Eingabe.insert(END, 0.01)
    #coxneg_Eingabe.pack() 
    coxneg_Eingabe.place(x = 150, y = 165)
    
    TankVolneg_Label = Label(Fenster_CDC_in,text="Tank Volume [L]*",bg = '#80ffbf')
    TankVolneg_Label.place(x = 25, y = 190)
    TankVolneg_Eingabe = Entry(Fenster_CDC_in)
    TankVolneg_Eingabe.insert(END, 0.1)
    #TankVolneg_Eingabe.pack() 
    TankVolneg_Eingabe.place(x = 150, y = 190)
    
    Ezeroneg_Label = Label(Fenster_CDC_in,text="E°(vs. SHE)[V]*",bg = '#80ffbf')
    Ezeroneg_Label.place(x = 25, y = 215)
    Ezeroneg_Eingabe = Entry(Fenster_CDC_in)
    Ezeroneg_Eingabe.insert(END, -0.255)
    #Ezeroneg_Eingabe.pack() 
    Ezeroneg_Eingabe.place(x = 150, y = 215)
    
    ActRetainneg_Label = Label(Fenster_CDC_in,text= u"\u0394 k° [cm/(s^2)]*",bg = '#80ffbf')
    ActRetainneg_Label.place(x = 25, y = 240)
    ActRetainneg_Eingabe = Entry(Fenster_CDC_in)
    ActRetainneg_Eingabe.insert(END, 0)
    #ActRetainneg_Eingabe.pack() 
    ActRetainneg_Eingabe.place(x = 150, y = 240) 
    
    
    #=========================================================================================================================
    #     POSITIVE SIDE INPUTS
    #=========================================================================================================================
    
    pos_Label = Label(Fenster_CDC_in, text="Inputs for positive side", bg = '#80ffbf', font = ('Arial',16,'bold'))
    pos_Label.place(x = 315, y = 20)
    
    npos_Label = Label(Fenster_CDC_in, text="n... keep it 1", bg = '#80ffbf')
    npos_Label.place(x = 315, y = 65)
    npos_Eingabe = Entry(Fenster_CDC_in)          
    npos_Eingabe.insert(END, 1)
    #npos_Eingabe.pack()
    npos_Eingabe.place(x = 440, y = 65)
    
    alphapos_Label = Label(Fenster_CDC_in, text= u"\u03B1*",bg = '#80ffbf')
    alphapos_Label.place(x = 315, y = 90)
    alphapos_Eingabe = Entry(Fenster_CDC_in)  
    alphapos_Eingabe.insert(END, 0.5)
    #alphapos_Eingabe.pack()                                             
    alphapos_Eingabe.place(x = 440, y = 90)
    
    kzeropos_Label = Label(Fenster_CDC_in,text="k°[cm/s]*",bg = '#80ffbf')
    kzeropos_Label.place(x = 315, y = 115)
    kzeropos_Eingabe = Entry(Fenster_CDC_in)
    kzeropos_Eingabe.insert(END, 1e-8)
    #kzeropos_Eingabe.pack() 
    kzeropos_Eingabe.place(x = 440, y = 115)
    
    credpos_Label = Label(Fenster_CDC_in,text="c_red in mol/L*",bg = '#80ffbf')
    credpos_Label.place(x = 315, y = 140)
    credpos_Eingabe = Entry(Fenster_CDC_in)
    credpos_Eingabe.insert(END, 0.01)
    #credpos_Eingabe.pack() 
    credpos_Eingabe.place(x = 440, y = 140) 
    
    coxpos_Label = Label(Fenster_CDC_in,text="c_ox in mol/L*",bg = '#80ffbf')
    coxpos_Label.place(x = 315, y = 165)
    coxpos_Eingabe = Entry(Fenster_CDC_in)
    coxpos_Eingabe.insert(END, 0.99)
    #coxpos_Eingabe.pack() 
    coxpos_Eingabe.place(x = 440, y = 165) 
    
    TankVolpos_Label = Label(Fenster_CDC_in,text="Tank Volume [L]*",bg = '#80ffbf')
    TankVolpos_Label.place(x = 315, y = 190)
    TankVolpos_Eingabe = Entry(Fenster_CDC_in)
    TankVolpos_Eingabe.insert(END, 0.1)
    #TankVolpos_Eingabe.pack()
    TankVolpos_Eingabe.place(x = 440, y = 190)
    
    Ezeropos_Label = Label(Fenster_CDC_in,text="E°(vs. SHE)[V]*",bg = '#80ffbf')
    Ezeropos_Label.place(x = 315, y = 215)
    Ezeropos_Eingabe = Entry(Fenster_CDC_in)
    Ezeropos_Eingabe.insert(END, 0.955)
    #Ezeropos_Eingabe.pack()
    Ezeropos_Eingabe.place(x = 440, y = 215)
    
    ActRetainpos_Label = Label(Fenster_CDC_in,text= u"\u0394 k° [cm/(s^2)]*",bg = '#80ffbf')
    ActRetainpos_Label.place(x = 315, y = 240)
    ActRetainpos_Eingabe = Entry(Fenster_CDC_in)
    ActRetainpos_Eingabe.insert(END, 0)
    #ActRetainpos_Eingabe.pack()
    ActRetainpos_Eingabe.place(x = 440, y = 240) 
 
    #=========================================================================================================================
    #     GENERAL INPUTS
    #=========================================================================================================================
    General_Label = Label(Fenster_CDC_in, text="General Inputs", bg = '#80ffbf', font = ('Arial',16,'bold'))
    General_Label.place(x = 600, y = 20)
    
    T_Label = Label(Fenster_CDC_in,text="T [°C]*",bg = '#80ffbf')
    T_Label.place(x = 600, y = 65)
    T_Eingabe = Entry(Fenster_CDC_in)
    T_Eingabe.insert(END, 25)
    #T_Eingabe.pack()
    T_Eingabe.place(x = 725, y = 65)
    
    I_Label = Label(Fenster_CDC_in,text="Drawn Current [A]*",bg = '#80ffbf')
    I_Label.place(x = 600, y = 90)
    I_Eingabe = Entry(Fenster_CDC_in)
    I_Eingabe.insert(END, 1)
    #I_Eingabe.pack()
    I_Eingabe.place(x = 725, y = 90)
    
    R_Label = Label(Fenster_CDC_in,text= u"Cell Resistance [\u03A9]*",bg = '#80ffbf')
    R_Label.place(x = 600, y = 115)
    R_Eingabe = Entry(Fenster_CDC_in)
    R_Eingabe.insert(END, 0)
    #R_Eingabe.pack()
    R_Eingabe.place(x = 725, y = 115)
    
    Rincr_Label = Label(Fenster_CDC_in,text= "R increase [\u03A9/s]*",bg = '#80ffbf')
    Rincr_Label.place(x = 600, y = 140)
    Rincr_Eingabe = Entry(Fenster_CDC_in)
    Rincr_Eingabe.insert(END, 0)
    #Rincr_Eingabe.pack()
    Rincr_Eingabe.place(x = 725, y = 140)
    
    ECutUp_Label = Label(Fenster_CDC_in,text= "Max cell voltage [V]*",bg = '#80ffbf')
    ECutUp_Label.place(x = 600, y = 165)
    ECutUp_Eingabe = Entry(Fenster_CDC_in)
    ECutUp_Eingabe.insert(END, 1.7)
    #ECutUp_Eingabe.pack()
    ECutUp_Eingabe.place(x = 725, y = 165)
    
    ECutLow_Label = Label(Fenster_CDC_in,text= "Min cell voltage [V]*",bg = '#80ffbf')
    ECutLow_Label.place(x = 600, y = 190)
    ECutLow_Eingabe = Entry(Fenster_CDC_in)
    ECutLow_Eingabe.insert(END, 0.8)
    #ECutLow_Eingabe.pack()
    ECutLow_Eingabe.place(x = 725, y = 190)
    
    FlowRate_Label = Label(Fenster_CDC_in,text= "Flowrate [mL/min]*",bg = '#80ffbf')
    FlowRate_Label.place(x = 600, y = 215)
    FlowRate_Eingabe = Entry(Fenster_CDC_in)
    FlowRate_Eingabe.insert(END, 50)
    #FlowRate_Eingabe.pack()
    FlowRate_Eingabe.place(x = 725, y = 215)
    
    NCycles_Label = Label(Fenster_CDC_in,text= "n-Cycles*",bg = '#80ffbf')
    NCycles_Label.place(x = 600, y = 240)
    NCycles_Eingabe = Entry(Fenster_CDC_in)
    NCycles_Eingabe.insert(END, 4)
    #NCycles_Eingabe.pack()
    NCycles_Eingabe.place(x = 725, y = 240)
    
    
    Border1_Label = Label(Fenster_CDC_in,text="-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    Border1_Label.place(x = 0, y = 265)
    
    #=========================================================================================================================
    #     ELECTRODE INPUTS
    #=========================================================================================================================
    
    
    Electrode_Label = Label(Fenster_CDC_in, text="Electrode Parameters", bg = '#80ffbf', font = ('Arial',16,'bold'))
    Electrode_Label.place(x = 25, y = 300)
    
    A_Label = Label(Fenster_CDC_in,text="Common A [cm^2]*",bg = '#80ffbf')
    A_Label.place(x = 25, y = 340)
    A_Eingabe = Entry(Fenster_CDC_in)
    A_Eingabe.insert(END, 120)
    #A_Eingabe.pack()
    A_Eingabe.place(x = 150, y = 340)
    
    lenx_Label = Label(Fenster_CDC_in,text="Len in x [cm]*",bg = '#80ffbf')
    lenx_Label.place(x = 25, y = 365)
    lenx_Eingabe = Entry(Fenster_CDC_in)
    lenx_Eingabe.insert(END, 0.4)
    #lenx_Eingabe.pack()
    lenx_Eingabe.place(x = 150, y = 365)
    
    leny_Label = Label(Fenster_CDC_in,text="Len in y [cm]*",bg = '#80ffbf')
    leny_Label.place(x = 25, y = 390)
    leny_Eingabe = Entry(Fenster_CDC_in)
    leny_Eingabe.insert(END, 3)
    #leny_Eingabe.pack()
    leny_Eingabe.place(x = 150, y = 390)
    
    lenz_Label = Label(Fenster_CDC_in,text="Len in z [cm]*",bg = '#80ffbf')
    lenz_Label.place(x = 25, y = 415)
    lenz_Eingabe = Entry(Fenster_CDC_in)
    lenz_Eingabe.insert(END, 5)
    #lenz_Eingabe.pack()
    lenz_Eingabe.place(x = 150, y = 415)
    
    Incz_Label = Label(Fenster_CDC_in,text="Increments in z*",bg = '#80ffbf')
    Incz_Label.place(x = 25, y = 440)
    Incz_Eingabe = Entry(Fenster_CDC_in)
    Incz_Eingabe.insert(END, 100)
    #Incz_Eingabe.pack()
    Incz_Eingabe.place(x = 150, y = 440)
    
    Poro_Label = Label(Fenster_CDC_in,text="Open porosity [%]*",bg = '#80ffbf')
    Poro_Label.place(x = 25, y = 465)
    Poro_Eingabe = Entry(Fenster_CDC_in)
    Poro_Eingabe.insert(END, 0.95)
    #Poro_Eingabe.pack()
    Poro_Eingabe.place(x = 150, y = 465)
    
    
    #=========================================================================================================================
    #     OTHER INPUTS
    #=========================================================================================================================
    
    
    Other_Label = Label(Fenster_CDC_in, text="Miscellaneous Inputs", bg = '#80ffbf', font = ('Arial',16,'bold'))
    Other_Label.place(x = 315, y = 300)
    
    CurrIncMod = IntVar()
    Checkbutton(Fenster_CDC_in, text="Change I-find Interations", variable=CurrIncMod ,bg = '#80ffbf').place(x = 315, y = 340)
    CurrFindIter_Eingabe = Entry(Fenster_CDC_in)
    CurrFindIter_Eingabe.place(x = 492, y = 340, width = 70, height = 20)
    
    
    
    
    
    def Go_ChargeDischarge():
        nneg           = (float(nneg_Eingabe.get()))
        alphaneg       = (float(alphaneg_Eingabe.get()))
        kzeroneg       = (float(kzeroneg_Eingabe.get()))
        credneg        = 0.001*(float(credneg_Eingabe.get()))
        coxneg         = 0.001*(float(coxneg_Eingabe.get()))
        TankVolneg     = 1000*(float(TankVolneg_Eingabe.get()))
        Ezeroneg       = (float(Ezeroneg_Eingabe.get()))
        ActRetainneg   = (float(ActRetainneg_Eingabe.get()))
        npos           = (float(npos_Eingabe.get()))
        alphapos       = (float(alphapos_Eingabe.get()))
        kzeropos       = (float(kzeropos_Eingabe.get()))
        credpos        = 0.001*(float(credpos_Eingabe.get()))
        coxpos         = 0.001*(float(coxpos_Eingabe.get()))
        TankVolpos     = 1000*(float(TankVolpos_Eingabe.get()))
        Ezeropos       = (float(Ezeropos_Eingabe.get()))
        ActRetainpos   = (float(ActRetainpos_Eingabe.get()))
        T              = (float(T_Eingabe.get())) + 273.15
        I              = (float(I_Eingabe.get()))
        Res            = (float(R_Eingabe.get()))
        Rincr          = (float(Rincr_Eingabe.get()))
        ECutUp         = (float(ECutUp_Eingabe.get()))
        ECutLow        = (float(ECutLow_Eingabe.get()))
        FlowRate       = (float(FlowRate_Eingabe.get()))
        NCycles        = (int(NCycles_Eingabe.get()))
        Area           = (float(A_Eingabe.get()))
        lenx           = (float(lenx_Eingabe.get()))
        leny           = (float(leny_Eingabe.get()))
        lenz           = (float(lenz_Eingabe.get()))
        lenz           = (float(lenz_Eingabe.get()))
        Incz           = (int(Incz_Eingabe.get()))
        Poro           = (float(Poro_Eingabe.get()))
        CurrIncModder  = CurrIncMod.get()
        CurrFindIter   = 50
        
        if CurrIncModder == 1:
            CurrFindIter = (int(CurrFindIter_Eingabe.get()))
            

        CDCMainConsole(NCycles, I, ECutLow, ECutUp, 
                         alphapos, kzeropos, alphaneg, kzeroneg, FlowRate, 
                         Poro, lenz, leny, lenx, Incz, CurrFindIter, 
                         Res, Rincr, credneg, coxneg, credpos, coxpos, ActRetainneg, ActRetainpos, T,
                         Ezeropos, Ezeroneg, Area, TankVolneg, TankVolpos)
        
        
        
    #=================================================================================================================
    
    button=Button(Fenster_CDC_in,text="Next", bg = '#9BA9C5', command = Go_ChargeDischarge).place(x = 600, y = 340, width = 250, height = 75)
    
    






def CYCLELOPPERFUNC_CDC(N_Cycless, Drawn_Currentt, LowerPotLimm, UpperPotlimm, 
                         alphaPosSitee, kzeroPosSitee, alphaNegSitee, kzeroNegSitee, FlowRate_mLPerMinn, 
                         OpenPorosityy, cm_in_zz, cm_in_yy, cm_in_xx, Incr_per_zz, CurrFindIterationss, 
                         Initial_Ress, RES_INCRR, Cred_negg, Cox_negg, Cred_poss, Cox_poss, Decayernegg, Decayerposs,
                         Temperaturee, Ezeroposs, Ezeronegg, Areaa, TankVolnegg, TankVolposs):

    
    NumbaExister = importlib.util.find_spec('numba')
    if NumbaExister is not None:
        from numba import jit 
        #========================================================================================================
        #========================================================================================================
        #IMPORTANT FUNCTIONS HIDDEN IN DROPDOWN
        
        @jit(nopython=True) 
        def Current_Butler_Volmer(cs_red, cs_ox, alpha, kzero, A, Eappl, Ezero):
            Eeq = Ezero + (R*T/(n*F))*np.log(cs_ox/cs_red)
            Ieq = n*F*A*kzero*cs_ox*np.exp(-(1-alpha)*n*F*(Eeq-Ezero)/(R*T))
            return Ieq*( np.exp(alpha*n*F*(Eappl - Eeq)/(R*T))   -  np.exp(-(1-alpha)*n*F*(Eappl - Eeq)/(R*T)) )
        
        @jit(nopython=True) 
        def Galvanopotfind(cs_red_array, cs_ox_array, alpha, kzero, A, Eappl, Ezero, Target_Current, MaxIter, dt, dV, printwhichloop):
            I_Array          = np.zeros(len(cs_red_array))
            cred_Array       = np.zeros(len(cs_red_array))
            cox_Array        = np.zeros(len(cs_red_array))
            #print("OOOOOOOOXXXXXXX", cox_Array[0], cox_Array[-1])
            #print("RRREEEEEEEEEDDD", cred_Array[0], cred_Array[-1])
            dN               = 0
            for i in range(len(I_Array)):
                I_Array[i]     = Current_Butler_Volmer(cs_red_array[i], cs_ox_array[i], alpha, kzero, A, Eappl, Ezero)/(float(len(I_Array)))
                dn             = I_Array[i]*dt/float(n)   # change in the amount of concentration per vloumeelement, when passing the electrode
                cred_Array[i]  = cs_red_array[i]-dn/dV
                cox_Array[i]   = cs_ox_array[i] +dn/dV
                dN            += dn                # Compute total dN at time instance by summing dn of all volume elements
            Current            =(dN/dt)*n*F           # Compute total curent at given time instance
            #print("CURRRRRR", Current)
            #Current            = ((cred_Array[0]-cred_Array[-1])*ActVol/Passing_time)*F     # bad approximate...
            NUMMM = 0
            Skip_Pre_Loop = "no"
            while np.abs(Current - Target_Current) > 0.0000000000001:
                NUMMM +=1
                if np.abs(Current - Target_Current) > 0.1 and Skip_Pre_Loop == "no":
                    if printwhichloop == "yes":
                        print("Pre-loop entered")
                    if Current-Target_Current > 0.1:
                        Eappl          = Eappl - 0.01
                    if Current-Target_Current < 0.1:
                        Eappl          = Eappl + 0.01
                    if NUMMM == MaxIter:
                        Skip_Pre_Loop = "yes"
                        NUMMM    = 0
                        MaxInter = 20*MaxIter
                if np.abs(Current - Target_Current) > 0.01 and np.abs(Current - Target_Current) <= 0.1 or Skip_Pre_Loop == "yes":
                    if printwhichloop == "yes":
                        print("Fine-loop entered")
                    if Current - Target_Current > 0:
                        Eappl          = Eappl - 0.00001
                    if Current - Target_Current < 0:
                        Eappl          = Eappl + 0.00001  
                if np.abs(Current - Target_Current) > 0.001 and np.abs(Current - Target_Current) <= 0.01:
                    if printwhichloop == "yes":
                        print("Fine-loop-2 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.00001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.00001  
                if np.abs(Current - Target_Current) > 0.0001 and np.abs(Current - Target_Current) <= 0.001:
                    if printwhichloop == "yes":
                        print("Fine-loop-3 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.000001  
                if np.abs(Current - Target_Current) > 0.00001 and np.abs(Current - Target_Current) <= 0.0001:
                    if printwhichloop == "yes":
                        print("Fine-loop-4 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.0000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.0000001  
                if np.abs(Current - Target_Current) > 0.000001 and np.abs(Current - Target_Current) <= 0.00001:
                    if printwhichloop == "yes":
                        print("Fine-loop-5 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.00000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.00000001  
                if np.abs(Current - Target_Current) > 0.0000001 and np.abs(Current - Target_Current) <= 0.000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-6 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.000000001  
                if np.abs(Current - Target_Current) > 0.00000001 and np.abs(Current - Target_Current) <= 0.0000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-7 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.0000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.0000000001  
                if np.abs(Current - Target_Current) > 0.000000001 and np.abs(Current - Target_Current) <= 0.00000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-8 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.00000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.00000000001 
                if np.abs(Current - Target_Current) > 0.0000000001 and np.abs(Current - Target_Current) <= 0.000000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-9 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.000000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.000000000001 
                if np.abs(Current - Target_Current) > 0.00000000001 and np.abs(Current - Target_Current) <= 0.0000000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-10 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.0000000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.0000000000001 
                if np.abs(Current - Target_Current) > 0.000000000001 and np.abs(Current - Target_Current) <= 0.00000000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-11 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.00000000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.00000000000001 
                if np.abs(Current - Target_Current) > 0.0000000000001 and np.abs(Current - Target_Current) <= 0.000000000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-12 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.000000000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.000000000000001 
             
        
                I_Array          = np.zeros(len(cs_red_array))
                dN               = 0
                for i in range(len(I_Array)):
                    I_Array[i]     = Current_Butler_Volmer(cs_red_array[i], cs_ox_array[i], alpha, kzero, A, Eappl, Ezero)/(float(len(I_Array)))
                    dn             = I_Array[i]*dt/float(n)   # change in the amount of V2+ and V3+ per vloumeelement, when passing the electrode
                    #print("DDDDDDDDDDNNNNNNNNNN", dn)
                    cred_Array[i]  = cs_red_array[i]-dn/dV
                    cox_Array[i]   = cs_ox_array[i] +dn/dV
                    dN            += dn                # Compute total dN at time instance by summing dn of all volume elements
                Current            =(dN/dt)*n*F           # Compute total curent at given time instance 
                #Current            = ((cred_Array[0]-cred_Array[-1])*ActVol/Passing_time)*F
                if printwhichloop == "yes":
                    print("Iteration number:",NUMMM,  "Current Approximate:", Current, "A")
                if NUMMM == MaxIter:
                    break
            
            return Eappl, cred_Array, cox_Array, Current
        
        
        
        
        
    if NumbaExister == None:
        
        #========================================================================================================
        #========================================================================================================
        #IMPORTANT FUNCTIONS HIDDEN IN DROPDOWN
         
        def Current_Butler_Volmer(cs_red, cs_ox, alpha, kzero, A, Eappl, Ezero):
            Eeq = Ezero + (R*T/(n*F))*np.log(cs_ox/cs_red)
            #Ieq = n*F*A*kzero*cs_red*np.exp(alpha*n*F*(Eeq-Ezero)/(R*T))
            Ieq = n*F*A*kzero*cs_ox*np.exp(-(1-alpha)*n*F*(Eeq-Ezero)/(R*T))
            return Ieq*( np.exp(alpha*n*F*(Eappl - Eeq)/(R*T))   -  np.exp(-(1-alpha)*n*F*(Eappl - Eeq)/(R*T)) )
         
        def Galvanopotfind(cs_red_array, cs_ox_array, alpha, kzero, A, Eappl, Ezero, Target_Current, MaxIter, dt, dV, printwhichloop):
            I_Array          = np.zeros(len(cs_red_array))
            cred_Array       = np.zeros(len(cs_red_array))
            cox_Array        = np.zeros(len(cs_red_array))
            #print("OOOOOOOOXXXXXXX", cox_Array[0], cox_Array[-1])
            #print("RRREEEEEEEEEDDD", cred_Array[0], cred_Array[-1])
            dN               = 0
            for i in range(len(I_Array)):
                I_Array[i]     = Current_Butler_Volmer(cs_red_array[i], cs_ox_array[i], alpha, kzero, A, Eappl, Ezero)/(float(len(I_Array)))
                dn             = I_Array[i]*dt/float(n)   # change in the amount of concentration per vloumeelement, when passing the electrode
                cred_Array[i]  = cs_red_array[i]-dn/dV
                cox_Array[i]   = cs_ox_array[i] +dn/dV
                dN            += dn                # Compute total dN at time instance by summing dn of all volume elements
            Current            =(dN/dt)*n*F           # Compute total curent at given time instance
            #print("CURRRRRR", Current)
            #Current            = ((cred_Array[0]-cred_Array[-1])*ActVol/Passing_time)*F     # bad approximate...
            NUMMM = 0
            Skip_Pre_Loop = "no"
            while np.abs(Current - Target_Current) > 0.0000000000001:
                NUMMM +=1
                if np.abs(Current - Target_Current) > 0.1 and Skip_Pre_Loop == "no":
                    if printwhichloop == "yes":
                        print("Pre-loop entered")
                    if Current-Target_Current > 0.1:
                        Eappl          = Eappl - 0.01
                    if Current-Target_Current < 0.1:
                        Eappl          = Eappl + 0.01
                    if NUMMM == MaxIter:
                        Skip_Pre_Loop = "yes"
                        NUMMM    = 0
                        MaxInter = 20*MaxIter
                if np.abs(Current - Target_Current) > 0.01 and np.abs(Current - Target_Current) <= 0.1 or Skip_Pre_Loop == "yes":
                    if printwhichloop == "yes":
                        print("Fine-loop entered")
                    if Current - Target_Current > 0:
                        Eappl          = Eappl - 0.00001
                    if Current - Target_Current < 0:
                        Eappl          = Eappl + 0.00001  
                if np.abs(Current - Target_Current) > 0.001 and np.abs(Current - Target_Current) <= 0.01:
                    if printwhichloop == "yes":
                        print("Fine-loop-2 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.00001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.00001  
                if np.abs(Current - Target_Current) > 0.0001 and np.abs(Current - Target_Current) <= 0.001:
                    if printwhichloop == "yes":
                        print("Fine-loop-3 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.000001  
                if np.abs(Current - Target_Current) > 0.00001 and np.abs(Current - Target_Current) <= 0.0001:
                    if printwhichloop == "yes":
                        print("Fine-loop-4 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.0000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.0000001  
                if np.abs(Current - Target_Current) > 0.000001 and np.abs(Current - Target_Current) <= 0.00001:
                    if printwhichloop == "yes":
                        print("Fine-loop-5 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.00000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.00000001  
                if np.abs(Current - Target_Current) > 0.0000001 and np.abs(Current - Target_Current) <= 0.000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-6 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.000000001  
                if np.abs(Current - Target_Current) > 0.00000001 and np.abs(Current - Target_Current) <= 0.0000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-7 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.0000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.0000000001  
                if np.abs(Current - Target_Current) > 0.000000001 and np.abs(Current - Target_Current) <= 0.00000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-8 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.00000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.00000000001 
                if np.abs(Current - Target_Current) > 0.0000000001 and np.abs(Current - Target_Current) <= 0.000000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-9 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.000000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.000000000001 
                if np.abs(Current - Target_Current) > 0.00000000001 and np.abs(Current - Target_Current) <= 0.0000000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-10 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.0000000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.0000000000001 
                if np.abs(Current - Target_Current) > 0.000000000001 and np.abs(Current - Target_Current) <= 0.00000000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-11 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.00000000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.00000000000001 
                if np.abs(Current - Target_Current) > 0.0000000000001 and np.abs(Current - Target_Current) <= 0.000000000001:
                    if printwhichloop == "yes":
                        print("Fine-loop-12 entered")
                    if Current-Target_Current > 0:
                        Eappl          = Eappl - 0.000000000000001
                    if Current-Target_Current < 0:
                        Eappl          = Eappl + 0.000000000000001 
             
        
                I_Array          = np.zeros(len(cs_red_array))
                dN               = 0
                for i in range(len(I_Array)):
                    I_Array[i]     = Current_Butler_Volmer(cs_red_array[i], cs_ox_array[i], alpha, kzero, A, Eappl, Ezero)/(float(len(I_Array)))
                    dn             = I_Array[i]*dt/float(n)   # change in the amount of V2+ and V3+ per vloumeelement, when passing the electrode
                    #print("DDDDDDDDDDNNNNNNNNNN", dn)
                    cred_Array[i]  = cs_red_array[i]-dn/dV
                    cox_Array[i]   = cs_ox_array[i] +dn/dV
                    dN            += dn                # Compute total dN at time instance by summing dn of all volume elements
                Current            =(dN/dt)*n*F           # Compute total curent at given time instance 
                #Current            = ((cred_Array[0]-cred_Array[-1])*ActVol/Passing_time)*F
                if printwhichloop == "yes":
                    print("Iteration number:",NUMMM,  "Current Approximate:", Current, "A")
                if NUMMM == MaxIter:
                    break
            
            return Eappl, cred_Array, cox_Array, Current
        
        
        
    #==========================================================================================================
    #==========================================================================================================
    #Define the Initiol Conditions
    #==========================================================================================================
    n =  1.0   ;   R   =  8.314   ;  F   =  96485.0   ;   T   = Temperaturee
    #========================================================================================================== 
    
    #_____negative size________;____________positive site_____
    cV2       = Cred_negg      ;       cV4       = Cred_poss      
    cV3       = Cox_negg       ;       cV5       = Cox_poss
    E0V2V3    = Ezeronegg      ;       E0V4V5    = Ezeroposs
    #=========================================================
    Vneg      =  TankVolnegg   ;       Vpos      =  TankVolposs             
    nV2       =  cV2*Vneg      ;       nV4       =  cV4*Vpos
    nV3       =  cV3*Vneg      ;       nV5       =  cV5*Vpos    
    #=========================================================
    alpha_neg = alphaNegSitee  ;       alpha_pos  = alphaPosSitee  
    A_neg     = Areaa          ;       A_pos      = Areaa           
    #=========================================================
    #=========================================================
    E_ges_neg = E0V2V3 + (R*T/(n*F))*np.log(cV3/cV2)
    E_ges_pos = E0V4V5 + (R*T/(n*F))*np.log(cV5/cV4)
    E_tot     = E_ges_pos - E_ges_neg
    

    
    #=============================================================================================================
    #=============================================================================================================
    #=============================================================================================================
    
    
    def ChargeDischarge(N_Cycles, Drawn_Current, LowerPotLim, UpperPotlim, 
                         alphaPosSite, kzeroPosSite, alphaNegSite, kzeroNegSite, FlowRate_mLPerMin, 
                         OpenPorosity, cm_in_z, cm_in_y, cm_in_x, Incr_per_z, CurrFindIterations, 
                         Initial_Res, RES_INCR, Cred_neg, Cox_neg, Cred_pos, Cox_pos, Decayerpos, Decayerneg):
        
        z_in_Electrode   = np.linspace(0,cm_in_z,Incr_per_z)      # lenght of the electrode cut in 10000 slices
        ylen             = cm_in_y                                # width of the electrode in y-direction
        xlen             = cm_in_x                                # Thickness of the electrode
        OpPor            = OpenPorosity                           # Open Porosite: percent of the electrode which is electrolyte (if filled)
        ActVol           = max(z_in_Electrode)*ylen*xlen*OpPor    # Volume of the electrode filled with electrolyte   
        dV               = ActVol/float(len(z_in_Electrode))      # Effective differential Volume element filled with electrolyte
        FlowRate         = FlowRate_mLPerMin                      # cm^3/min
        ProjVol          = OpPor*ylen*xlen                        # Projected Volume area
        Velocity         = FlowRate/(60*ProjVol)                  # flowrate in cm/s
        dt               = (z_in_Electrode[1]-z_in_Electrode[0])/Velocity
        DRAWNCURR        = Drawn_Current    
        N_Steps_in_loop  = CurrFindIterations    
        MaxTime          = nV2*F/DRAWNCURR
        UPLIM            = UpperPotlim
        LOWLIM           = LowerPotLim
        cV2              = Cred_neg
        cV3              = Cox_neg
        cV4              = Cred_pos
        cV5              = Cox_pos
        cV2_init_array   = cV2*np.ones(len(z_in_Electrode))
        cV3_init_array   = cV3*np.ones(len(z_in_Electrode))
        cV4_init_array   = cV4*np.ones(len(z_in_Electrode))
        cV5_init_array   = cV5*np.ones(len(z_in_Electrode))
        
        
        #Initialize Output Arrays
        Exp_Durat_final     = np.array([0])
        Pots_neg_final      = np.array([0])
        Pots_pos_final      = np.array([0])
        THECURRENT_n_final  = np.array([0])
        THECURRENT_p_final  = np.array([0])
        cV2_of_t_final      = np.array([0])
        cV3_of_t_final      = np.array([0])
        cV4_of_t_final      = np.array([0])
        cV5_of_t_final      = np.array([0])
        CellVoltage_final   = np.array([0])
        NernstPot_Pos_final = np.array([0])
        NernstPot_Neg_final = np.array([0])
        NernstPot_Tot_final = np.array([0])
        
        RESISTANCE          = Initial_Res
        kzero_neg           = kzeroNegSite
        kzero_pos           = kzeroPosSite
        
        for cyy in range(N_Cycles):
            cy = cyy + 1
           
            
            CDCwrite("==================================================\n")
            CDCwrite("Cycle  ")
            CDCwrite(cy)
            CDCwrite("  of  ")
            CDCwrite(N_Cycles)
            CDCwrite("\n==================================================\n")
        
            #====================================================================================
            ###########                FIRST DISCHARGE CYCLE            #########################
            #====================================================================================
            #Initialize calculations
            #====================================================================================

            if cy%2 != 0 :    
                
                
                CDCwrite("!!!!!!!!!  DISCHARGE !!!!!!!!!\n")
                CDCwrite("==================================================\n")
                
                AAA = Galvanopotfind(cV2_init_array, cV3_init_array, alpha_neg, kzero_neg, A_neg, Eappl = -0.2, Ezero = E0V2V3, 
                                     Target_Current =DRAWNCURR, MaxIter = 5000, dt = dt, dV = dV, printwhichloop = "no")
          
    
                BBB = Galvanopotfind(cV4_init_array, cV5_init_array, alpha_pos, kzero_pos, A_pos, Eappl =  0.9, Ezero = E0V4V5, 
                                     Target_Current =-DRAWNCURR, MaxIter = 5000, dt = dt, dV = dV, printwhichloop = "no")
        
                #====================================================================================
                #Return and print first results
                #====================================================================================
                print_it = "yes"
                if print_it == "yes":
                    CDCwrite("Effective flow velocity                                                       = ")
                    CDCwrite(round(Velocity, 3))
                    CDCwrite(" cm/s\n")
                    CDCwrite("Contact time per z-element                                                = ")
                    CDCwrite(round(dt, 4))
                    CDCwrite(" s \n")
                    CDCwrite("Initial cell voltage under current                                          = ")
                    CDCwrite( round(BBB[0]-AAA[0],3))
                    CDCwrite(" V\n")
                    CDCwrite("Total Nernstian Potential Difference at start (charged)      = ")
                    CDCwrite(round(E_tot,3))
                    CDCwrite(" V\n")
                    CDCwrite("Negative Nernstian Potential Difference at start (charged) = ")
                    CDCwrite(round(E_ges_neg,3))
                    CDCwrite(" V\n")
                    CDCwrite("Positive Nernstian Potential Difference at start (charged)  = ")
                    CDCwrite(round(E_ges_pos,3))
                    CDCwrite(" V\n")
                    CDCwrite("Maximum time for full conversion                                      = ")
                    CDCwrite(round(MaxTime,2))
                    CDCwrite(" s\n")
                    CDCwrite("SOC neg      = ")
                    CDCwrite(round((cV2/(cV2+cV3))*100,3))
                    CDCwrite("\n")
                    CDCwrite("SOC pos     = ")
                    CDCwrite(round((cV5/(cV4+cV5))*100,3))
                    CDCwrite("\n")
                    CDCtext_area.see(tk.END)
                    CDCMainwin.update()
    
                #====================================================================================
                #Initialize empty arrays to store the solution(s)
                #====================================================================================
                initialie_calcs = "yes"
                if initialie_calcs == "yes":
                    Exp_Durat     = np.arange(dt, MaxTime, dt )  + np.max(Exp_Durat_final)
                    Pots_neg      = np.zeros(len(Exp_Durat))
                    Pots_pos      = np.zeros(len(Exp_Durat))
                    THECURRENT_n  = np.zeros(len(Exp_Durat))
                    THECURRENT_p  = np.zeros(len(Exp_Durat))
                    cV2_of_t      = np.zeros(len(Exp_Durat))
                    cV3_of_t      = np.zeros(len(Exp_Durat))
                    cV4_of_t      = np.zeros(len(Exp_Durat))
                    cV5_of_t      = np.zeros(len(Exp_Durat))
                    CellVoltage   = np.zeros(len(Exp_Durat))
                    NernstPot_Pos = np.zeros(len(Exp_Durat))
                    NernstPot_Neg = np.zeros(len(Exp_Durat))
                    NernstPot_Tot = np.zeros(len(Exp_Durat))
                    
                    
                    CDCwrite("Timepoints = ")
                    CDCwrite(len(Exp_Durat))
                    CDCwrite("\n")
                    CDCtext_area.see(tk.END)
                    CDCMainwin.update()
                #====================================================================================
                #ESSENTIAL For-Loop for the calculations + slicing the resuts
                #====================================================================================
                
                #for k in tqdm(range(len(Exp_Durat))):
                lennn0  = 0
                lennn1  = int(0.05*len(Exp_Durat))
                lennn2  = int(0.10*len(Exp_Durat))
                lennn3  = int(0.15*len(Exp_Durat))
                lennn4  = int(0.20*len(Exp_Durat))
                lennn5  = int(0.25*len(Exp_Durat))
                lennn6  = int(0.30*len(Exp_Durat))
                lennn7  = int(0.35*len(Exp_Durat))
                lennn8  = int(0.40*len(Exp_Durat))
                lennn9  = int(0.45*len(Exp_Durat))
                lennn10 = int(0.50*len(Exp_Durat))
                lennn11 = int(0.55*len(Exp_Durat))
                lennn12 = int(0.60*len(Exp_Durat))
                lennn13 = int(0.65*len(Exp_Durat))
                lennn14 = int(0.70*len(Exp_Durat))
                lennn15 = int(0.75*len(Exp_Durat))
                lennn16 = int(0.80*len(Exp_Durat))
                lennn17 = int(0.85*len(Exp_Durat))
                lennn18 = int(0.90*len(Exp_Durat))
                lennn19 = int(0.95*len(Exp_Durat))
                lennn20 = int(1.00*len(Exp_Durat))
                
                
                CDCwrite("Computing, please wait:")
                CDCtext_area.see(tk.END)
                CDCMainwin.update()
                
                for k in range(len(Exp_Durat)):
                    if k == 0:
                        Result_n     = Galvanopotfind(cV2_init_array, cV3_init_array, alpha_neg, kzero_neg, A_neg, Eappl = AAA[0], Ezero = E0V2V3, 
                                                      Target_Current = DRAWNCURR, MaxIter = 10000, dt = dt, dV = dV, printwhichloop = "no")
                        Result_p     = Galvanopotfind(cV4_init_array, cV5_init_array, alpha_pos, kzero_pos, A_pos, Eappl = BBB[0], Ezero = E0V4V5, 
                                                      Target_Current =-DRAWNCURR, MaxIter = 10000,  dt = dt, dV = dV, printwhichloop = "no")
                    
                    if k > 0:
                        Result_n     = Galvanopotfind(cV2_init_array, cV3_init_array, alpha_neg, kzero_neg, A_neg, Eappl = Result_n[0], Ezero = E0V2V3, 
                                                      Target_Current = DRAWNCURR, MaxIter = N_Steps_in_loop,  dt = dt, dV = dV, printwhichloop = "no")  
                        Result_p     = Galvanopotfind(cV4_init_array, cV5_init_array, alpha_pos, kzero_pos, A_pos, Eappl = Result_p[0], Ezero = E0V4V5, 
                                                      Target_Current =-DRAWNCURR, MaxIter = N_Steps_in_loop,  dt = dt, dV = dV, printwhichloop = "no")
                     
                    
                    if k == lennn0 or k == lennn1 or k == lennn2 or k == lennn3 or k == lennn4 or k == lennn5 or k == lennn6 or k == lennn7 or k == lennn8 or k == lennn9 or k == lennn10 or k == lennn11 or k == lennn12 or k == lennn13 or k == lennn14 or k == lennn15 or k == lennn16 or k == lennn17 or k == lennn18 or k == lennn19 or k == lennn20 :
                        CDCwrite("#")
                        CDCtext_area.see(tk.END)
                        CDCMainwin.update()
                                        
                    
                    cV2                  = (cV2*(100-ActVol-dV) + Result_n[1][-1]*dV)/(100-ActVol) 
                    cV3                  = (cV3*(100-ActVol-dV) + Result_n[2][-1]*dV)/(100-ActVol) 
                    cV2_of_t[k]          = cV2
                    cV3_of_t[k]          = cV3
                    cV2_init_array[1::]  = Result_n[1][0:-1:]
                    cV3_init_array[1::]  = Result_n[2][0:-1:]
                    cV2_init_array[0]    = cV2_of_t[k] 
                    cV3_init_array[0]    = cV3_of_t[k] 
                    RESISTANCE          += RES_INCR*dt
                    kzero_neg           += Decayerneg*dt
                    kzero_pos           += Decayerpos*dt
                    THECURRENT_n[k]      = Result_n[3] 
                    Pots_neg[k]          = Result_n[0]  
                    cV4                  = (cV4*(100-ActVol-dV) + Result_p[1][-1]*dV)/(100-ActVol) 
                    cV5                  = (cV5*(100-ActVol-dV) + Result_p[2][-1]*dV)/(100-ActVol) 
                    cV4_of_t[k]          = cV4
                    cV5_of_t[k]          = cV5
                    cV4_init_array[1::]  = Result_p[1][0:-1:]
                    cV5_init_array[1::]  = Result_p[2][0:-1:]
                    cV4_init_array[0]    = cV4_of_t[k] 
                    cV5_init_array[0]    = cV5_of_t[k] 
                    THECURRENT_p[k]      = Result_p[3]
                    Pots_pos[k]          = Result_p[0] 
                    NernstPot_Neg[k]     = E0V2V3 + (R*T/(n*F))*np.log(cV3/cV2)
                    NernstPot_Pos[k]     = E0V4V5 + (R*T/(n*F))*np.log(cV5/cV4)
                    NernstPot_Tot[k]     = NernstPot_Pos[k] - NernstPot_Neg[k]
                    CellVoltage[k]       = (Pots_pos[k] - Pots_neg[k]) - RESISTANCE*THECURRENT_n[k]   # has to be - because a positive currents forces a more negative potential during discharge
                    STOP                 = k
                    if np.abs(CellVoltage[k]) < LOWLIM or np.abs(CellVoltage[k]) > UPLIM:
                        STOP = k
                        CDCwrite("\n")
                        CDCwrite("Stop criterion achieved by limiting voltage. Loop left at index: ")
                        CDCwrite(k)
                        CDCwrite("\n")
                        CDCtext_area.see(tk.END)
                        CDCMainwin.update()
                        break
                    
                    if np.isnan(THECURRENT_n[k]) == True  or np.isnan(THECURRENT_p[k]) == True or np.isnan(cV2)== True or np.isnan(cV3)== True or np.isnan(cV4)== True  or np.isnan(cV5)== True or np.isnan(Pots_neg[k]) == True or np.isnan(Pots_pos[k]) == True or np.any(np.isnan(cV2_init_array)) == True or np.any(np.isnan(cV3_init_array)) == True or np.any(np.isnan(cV4_init_array)) == True or np.any(np.isnan(cV5_init_array)) == True:
                        STOP =          k-5
                        
                        time.sleep(2)
                        cV2_init_array = cV2_of_t[STOP]*np.ones(len(z_in_Electrode))
                        cV3_init_array = cV3_of_t[STOP]*np.ones(len(z_in_Electrode))
                        cV4_init_array = cV4_of_t[STOP]*np.ones(len(z_in_Electrode))
                        cV5_init_array = cV5_of_t[STOP]*np.ones(len(z_in_Electrode))
                        
                        CDCwrite("\n")
                        CDCwrite("Stop criterion because current cannot be held at index: ")
                        CDCwrite(k)
                        CDCwrite("\n")
                        CDCtext_area.see(tk.END)
                        CDCMainwin.update()
                        break
               
                Exp_Durat           = Exp_Durat[0:STOP]
                Pots_neg            = Pots_neg[0:STOP] 
                Pots_pos            = Pots_pos[0:STOP]
                THECURRENT_n        = THECURRENT_n[0:STOP]
                THECURRENT_p        = THECURRENT_p[0:STOP]
                cV2_of_t            = cV2_of_t[0:STOP]
                cV3_of_t            = cV3_of_t[0:STOP]
                cV4_of_t            = cV4_of_t[0:STOP]
                cV5_of_t            = cV5_of_t[0:STOP]
                CellVoltage         = CellVoltage[0:STOP]
                NernstPot_Pos       = NernstPot_Pos[0:STOP]
                NernstPot_Neg       = NernstPot_Neg[0:STOP]
                NernstPot_Tot       = NernstPot_Tot[0:STOP]
        
                Exp_Durat_final     = np.concatenate([Exp_Durat_final, Exp_Durat])
                Pots_neg_final      = np.concatenate([Pots_neg_final, Pots_neg])
                Pots_pos_final      = np.concatenate([Pots_pos_final, Pots_pos])
                THECURRENT_n_final  = np.concatenate([THECURRENT_n_final, THECURRENT_n])
                THECURRENT_p_final  = np.concatenate([THECURRENT_p_final, THECURRENT_p])
                cV2_of_t_final      = np.concatenate([cV2_of_t_final, cV2_of_t])
                cV3_of_t_final      = np.concatenate([cV3_of_t_final, cV3_of_t])
                cV4_of_t_final      = np.concatenate([cV4_of_t_final, cV4_of_t])
                cV5_of_t_final      = np.concatenate([cV5_of_t_final, cV5_of_t])
                CellVoltage_final   = np.concatenate([CellVoltage_final, CellVoltage])
                NernstPot_Pos_final = np.concatenate([NernstPot_Pos_final, NernstPot_Pos])
                NernstPot_Neg_final = np.concatenate([NernstPot_Neg_final, NernstPot_Neg])
                NernstPot_Tot_final = np.concatenate([NernstPot_Tot_final, NernstPot_Tot])
    
    
            if cy%2 == 0:
                
                CDCwrite("!!!!!!!!!  CHARGE !!!!!!!!!\n")
                CDCwrite("==================================================\n")
                AAA2 = Galvanopotfind(cV2_init_array, cV3_init_array, alpha_neg, kzero_neg, A_neg, Eappl = -0.2, Ezero = E0V2V3, 
                         Target_Current =-DRAWNCURR, MaxIter = 5000,  dt = dt, dV = dV,  printwhichloop = "no")
                BBB2 = Galvanopotfind(cV4_init_array, cV5_init_array, alpha_pos, kzero_pos, A_pos, Eappl =  0.9, Ezero = E0V4V5, 
                                     Target_Current =DRAWNCURR, MaxIter = 5000,  dt = dt, dV = dV,  printwhichloop = "no")
            
                E_ges_neg2 = E0V2V3 + (R*T/(n*F))*np.log(cV3/cV2)
                E_ges_pos2 = E0V4V5 + (R*T/(n*F))*np.log(cV5/cV4)
                E_tot2     = E_ges_pos2 - E_ges_neg2
                
                #====================================================================================
                #Return and print first results
                #====================================================================================
                print_it = "yes"
                if print_it == "yes":
                    CDCwrite("Effective flow velocity                                                       = ")
                    CDCwrite(round(Velocity, 3))
                    CDCwrite(" cm/s\n")
                    CDCwrite("Contact time per z-element                                                = ")
                    CDCwrite(round(dt, 4))
                    CDCwrite(" s \n")
                    CDCwrite("Initial cell voltage under current                                          = ")
                    CDCwrite( round(BBB2[0]-AAA2[0],3))
                    CDCwrite(" V\n")
                    CDCwrite("Total Nernstian Potential Difference at start (charged)      = ")
                    CDCwrite(round(E_tot2,3))
                    CDCwrite(" V\n")
                    CDCwrite("Negative Nernstian Potential Difference at start (charged) = ")
                    CDCwrite(round(E_ges_neg2,3))
                    CDCwrite(" V\n")
                    CDCwrite("Positive Nernstian Potential Difference at start (charged)  = ")
                    CDCwrite(round(E_ges_pos2,3))
                    CDCwrite(" V\n")
                    CDCwrite("Maximum time for full conversion                                      = ")
                    CDCwrite(round(MaxTime,2))
                    CDCwrite(" s\n")
                    CDCwrite("SOC neg      = ")
                    CDCwrite(round((cV2/(cV2+cV3))*100,3))
                    CDCwrite("\n")
                    CDCwrite("SOC pos     = ")
                    CDCwrite(round((cV5/(cV4+cV5))*100,3))
                    CDCwrite("\n")
                    CDCtext_area.see(tk.END)
                    CDCMainwin.update()
    
                #====================================================================================
                #Initialize empty arrays to store the solution(s)
                #====================================================================================
                initialze_calcs = "yes"
                if initialze_calcs == "yes":
                    Exp_Durat2     = np.arange(dt, MaxTime, dt ) + np.max(Exp_Durat_final)
                    Pots_neg2      = np.zeros(len(Exp_Durat2))
                    Pots_pos2      = np.zeros(len(Exp_Durat2))
                    THECURRENT_n2  = np.zeros(len(Exp_Durat2))
                    THECURRENT_p2  = np.zeros(len(Exp_Durat2))
                    cV2_of_t2      = np.zeros(len(Exp_Durat2))
                    cV3_of_t2      = np.zeros(len(Exp_Durat2))
                    cV4_of_t2      = np.zeros(len(Exp_Durat2))
                    cV5_of_t2      = np.zeros(len(Exp_Durat2))
                    CellVoltage2   = np.zeros(len(Exp_Durat2))
                    NernstPot_Pos2 = np.zeros(len(Exp_Durat2))
                    NernstPot_Neg2 = np.zeros(len(Exp_Durat2))
                    NernstPot_Tot2 = np.zeros(len(Exp_Durat2))
                    CDCwrite("Timepoints = ")
                    CDCwrite(len(Exp_Durat2))
                    CDCwrite("\n")
                    CDCtext_area.see(tk.END)
                    CDCMainwin.update()
                    
                #====================================================================================
                #ESSENTIAL For-Loop for the calculations + slicing the resuts
                #====================================================================================
               
                
               
                
                lennn0  = 0
                lennn1  = int(0.05*len(Exp_Durat2))
                lennn2  = int(0.10*len(Exp_Durat2))
                lennn3  = int(0.15*len(Exp_Durat2))
                lennn4  = int(0.20*len(Exp_Durat2))
                lennn5  = int(0.25*len(Exp_Durat2))
                lennn6  = int(0.30*len(Exp_Durat2))
                lennn7  = int(0.35*len(Exp_Durat2))
                lennn8  = int(0.40*len(Exp_Durat2))
                lennn9  = int(0.45*len(Exp_Durat2))
                lennn10 = int(0.50*len(Exp_Durat2))
                lennn11 = int(0.55*len(Exp_Durat2))
                lennn12 = int(0.60*len(Exp_Durat2))
                lennn13 = int(0.65*len(Exp_Durat2))
                lennn14 = int(0.70*len(Exp_Durat2))
                lennn15 = int(0.75*len(Exp_Durat2))
                lennn16 = int(0.80*len(Exp_Durat2))
                lennn17 = int(0.85*len(Exp_Durat2))
                lennn18 = int(0.90*len(Exp_Durat2))
                lennn19 = int(0.95*len(Exp_Durat2))
                lennn20 = int(1.00*len(Exp_Durat2))
                #for k in tqdm(range(len(Exp_Durat2))):
                CDCwrite("Computing, please wait:")
                CDCtext_area.see(tk.END)
                CDCMainwin.update()
                
                for k in range(len(Exp_Durat2)):
                    
                    if k == 0:
                        Result_n2     = Galvanopotfind(cV2_init_array, cV3_init_array, alpha_neg, kzero_neg, A_neg, Eappl = AAA2[0], Ezero = E0V2V3, 
                                                       Target_Current =-DRAWNCURR, MaxIter = 10000,  dt = dt, dV = dV, printwhichloop = "no")
                        Result_p2     = Galvanopotfind(cV4_init_array, cV5_init_array, alpha_pos, kzero_pos, A_pos, Eappl = BBB2[0], Ezero = E0V4V5, 
                                                       Target_Current =DRAWNCURR, MaxIter = 10000,  dt = dt, dV = dV, printwhichloop = "no")
                    
                    if k > 0:
                        Result_n2     = Galvanopotfind(cV2_init_array, cV3_init_array, alpha_neg, kzero_neg, A_neg, Eappl = Result_n2[0], Ezero = E0V2V3, 
                                                      Target_Current =-DRAWNCURR, MaxIter = N_Steps_in_loop,  dt = dt, dV = dV, printwhichloop = "no")  
                        Result_p2     = Galvanopotfind(cV4_init_array, cV5_init_array, alpha_pos, kzero_pos, A_pos, Eappl = Result_p2[0], Ezero = E0V4V5, 
                                                      Target_Current =DRAWNCURR, MaxIter = N_Steps_in_loop,  dt = dt, dV = dV, printwhichloop = "no")  
                                                  
                    if k == lennn0 or k == lennn1 or k == lennn2 or k == lennn3 or k == lennn4 or k == lennn5 or k == lennn6 or k == lennn7 or k == lennn8 or k == lennn9 or k == lennn10 or k == lennn11 or k == lennn12 or k == lennn13 or k == lennn14 or k == lennn15 or k == lennn16 or k == lennn17 or k == lennn18 or k == lennn19 or k == lennn20 :
                        CDCwrite("#")
                        CDCtext_area.see(tk.END)
                        CDCMainwin.update()
                    
                    cV2                  = (cV2*(100-ActVol-dV) + Result_n2[1][-1]*dV)/(100-ActVol) 
                    cV3                  = (cV3*(100-ActVol-dV) + Result_n2[2][-1]*dV)/(100-ActVol) 
                    cV2_of_t2[k]         = cV2
                    cV3_of_t2[k]         = cV3
                    cV2_init_array[1::]  = Result_n2[1][0:-1:]
                    cV3_init_array[1::]  = Result_n2[2][0:-1:]
                    cV2_init_array[0]    = cV2_of_t2[k] 
                    cV3_init_array[0]    = cV3_of_t2[k] 
                    RESISTANCE          += RES_INCR*dt
                    kzero_neg           += Decayerneg*dt
                    kzero_pos           += Decayerpos*dt
                    THECURRENT_n2[k]     = Result_n2[3] 
                    Pots_neg2[k]         = Result_n2[0]
                    cV4                  = (cV4*(100-ActVol-dV) + Result_p2[1][-1]*dV)/(100-ActVol) 
                    cV5                  = (cV5*(100-ActVol-dV) + Result_p2[2][-1]*dV)/(100-ActVol) 
                    cV4_of_t2[k]         = cV4
                    cV5_of_t2[k]         = cV5
                    cV4_init_array[1::]  = Result_p2[1][0:-1:]
                    cV5_init_array[1::]  = Result_p2[2][0:-1:]
                    cV4_init_array[0]    = cV4_of_t2[k] 
                    cV5_init_array[0]    = cV5_of_t2[k] 
                    THECURRENT_p2[k]     = Result_p2[3]
                    Pots_pos2[k]         = Result_p2[0]
                    NernstPot_Neg2[k]    = E0V2V3 + (R*T/(n*F))*np.log(cV3/cV2)
                    NernstPot_Pos2[k]    = E0V4V5 + (R*T/(n*F))*np.log(cV5/cV4)
                    NernstPot_Tot2[k]    = NernstPot_Pos2[k] - NernstPot_Neg2[k]
                    CellVoltage2[k]      = (Pots_pos2[k] - Pots_neg2[k]) - RESISTANCE*THECURRENT_n2[k]   #also - RES..., since the current flips the sign
                    STOP                 = k
                    if np.abs(CellVoltage2[k]) < LOWLIM or np.abs(CellVoltage2[k]) > UPLIM:
                        STOP = k-1
                        CDCwrite("\n")
                        CDCwrite("Stop criterion achieved by limiting voltage. Loop left at index: ")
                        CDCwrite(k)
                        CDCwrite("\n")
                        CDCtext_area.see(tk.END)
                        CDCMainwin.update()
                        break
                    
                    if np.isnan(THECURRENT_n2[k]) == True  or np.isnan(THECURRENT_p2[k]) == True or np.isnan(cV2)== True or np.isnan(cV3)== True or np.isnan(cV4)== True  or np.isnan(cV5)== True or np.isnan(Pots_neg2[k]) == True or np.isnan(Pots_pos2[k]) == True or np.any(np.isnan(cV2_init_array)) == True or np.any(np.isnan(cV3_init_array)) == True or np.any(np.isnan(cV4_init_array)) == True or np.any(np.isnan(cV5_init_array)) == True:
                        
                        STOP =           k-5

                        time.sleep(2)
                        cV2_init_array = cV2_of_t2[STOP]*np.ones(len(z_in_Electrode))
                        cV3_init_array = cV3_of_t2[STOP]*np.ones(len(z_in_Electrode))
                        cV4_init_array = cV4_of_t2[STOP]*np.ones(len(z_in_Electrode))
                        cV5_init_array = cV5_of_t2[STOP]*np.ones(len(z_in_Electrode))

                        CDCwrite("\n")
                        CDCwrite("Stop criterion because current cannot be held at index: ")
                        CDCwrite(k)
                        CDCwrite("\n")
                        CDCtext_area.see(tk.END)
                        CDCMainwin.update()
                        break
                
                Exp_Durat2     = Exp_Durat2[0:STOP]
                Pots_neg2      = Pots_neg2[0:STOP] 
                Pots_pos2      = Pots_pos2[0:STOP]
                THECURRENT_n2  = THECURRENT_n2[0:STOP]
                THECURRENT_p2  = THECURRENT_p2[0:STOP]
                cV2_of_t2      = cV2_of_t2[0:STOP]
                cV3_of_t2      = cV3_of_t2[0:STOP]
                cV4_of_t2      = cV4_of_t2[0:STOP]
                cV5_of_t2      = cV5_of_t2[0:STOP]
                CellVoltage2   = CellVoltage2[0:STOP]
                NernstPot_Pos2 = NernstPot_Pos2[0:STOP]
                NernstPot_Neg2 = NernstPot_Neg2[0:STOP]
                NernstPot_Tot2 = NernstPot_Tot2[0:STOP]
                
                Exp_Durat_final     = np.concatenate([Exp_Durat_final, Exp_Durat2])
                Pots_neg_final      = np.concatenate([Pots_neg_final, Pots_neg2])
                Pots_pos_final      = np.concatenate([Pots_pos_final, Pots_pos2])
                THECURRENT_n_final  = np.concatenate([THECURRENT_n_final, THECURRENT_n2])
                THECURRENT_p_final  = np.concatenate([THECURRENT_p_final, THECURRENT_p2])
                cV2_of_t_final      = np.concatenate([cV2_of_t_final, cV2_of_t2])
                cV3_of_t_final      = np.concatenate([cV3_of_t_final, cV3_of_t2])
                cV4_of_t_final      = np.concatenate([cV4_of_t_final, cV4_of_t2])
                cV5_of_t_final      = np.concatenate([cV5_of_t_final, cV5_of_t2])
                CellVoltage_final   = np.concatenate([CellVoltage_final, CellVoltage2])
                NernstPot_Pos_final = np.concatenate([NernstPot_Pos_final, NernstPot_Pos2])
                NernstPot_Neg_final = np.concatenate([NernstPot_Neg_final, NernstPot_Neg2])
                NernstPot_Tot_final = np.concatenate([NernstPot_Tot_final, NernstPot_Tot2])
    
    
        
        #====================================================================================
        #Plotting of the Results
        #====================================================================================
        plot_it = "yes"
        if plot_it == "yes":
            plt.figure(figsize=(6,6), dpi = 100)
            
            plt.subplot(221)
            plt.plot(Exp_Durat_final[1::], Pots_neg_final[1::], color = 'blue', label = "neg site pots BV only")
            plt.plot(Exp_Durat_final[1::], Pots_pos_final[1::], color = 'red', label = "pos site pots BV only")
            plt.plot(Exp_Durat_final[1::], NernstPot_Neg_final[1::], color = 'black', label = "neg site pots Nernst", linestyle = ':', linewidth = 1)
            plt.plot(Exp_Durat_final[1::], NernstPot_Pos_final[1::], color = 'black', label = "pos site pots Nernst", linestyle = ':', linewidth = 1)
            plt.axhline(E0V2V3, color = 'black', linewidth = 0.5)
            plt.axhline(E0V4V5, color = 'black', linewidth = 0.5)
            plt.ylim(-0.55,1.29)
            plt.legend(frameon = False)
            
            plt.subplot(222)
            plt.plot(Exp_Durat_final[1::], CellVoltage_final[1::], color = 'black', label = "Cell voltage BV + Ohm")
            plt.plot(Exp_Durat_final[1::], NernstPot_Tot_final[1::], color = 'black', label = "Cell voltage Nernst", linestyle = ':', linewidth = 1)     
            plt.ylim(0.52,1.8)
            plt.legend(frameon = False)
            
            plt.subplot(223)
            #plt.ylim(0.999,1.001)
            plt.plot(Exp_Durat_final[1::], THECURRENT_n_final[1::], color = 'blue', label = "neg current")
            plt.plot(Exp_Durat_final[1::], -THECURRENT_p_final[1::], color = 'red', label = "-pos current", linestyle = ':', linewidth = 2)
            plt.ylim(-1.6,1.1)
            plt.legend(frameon = False)
            
            plt.subplot(224)
            plt.plot(Exp_Durat_final[1::], 1000*cV2_of_t_final[1::], color ='magenta', label = "c_V2+")
            plt.plot(Exp_Durat_final[1::], 1000*cV3_of_t_final[1::], color ='green', label = "c_V3+")
            plt.plot(Exp_Durat_final[1::], 1000*cV4_of_t_final[1::], color ='blue', label = "c_V4+", linestyle = ':', linewidth = 1)
            plt.plot(Exp_Durat_final[1::], 1000*cV5_of_t_final[1::], color ='orange', label = "c_V5+", linestyle = ':', linewidth = 1)
            plt.axhline(0, color = 'black', linewidth = 0.5)
            plt.axhline(1, color = 'black', linewidth = 0.5)
            plt.ylim(-0.1,1.4)
            plt.legend(ncol = 2, frameon = False)
            
            
            plt.tight_layout()
            plt.show()
    
    
    #=============================================================================================================
    #=============================================================================================================
    #=============================================================================================================
    
    
    Result = ChargeDischarge(N_Cycles = N_Cycless, Drawn_Current = Drawn_Currentt, LowerPotLim = LowerPotLimm, 
                             UpperPotlim = UpperPotlimm, alphaPosSite =  alphaPosSitee, kzeroPosSite = kzeroPosSitee, 
                             alphaNegSite = alphaNegSitee, kzeroNegSite = kzeroNegSitee, 
                             FlowRate_mLPerMin = FlowRate_mLPerMinn, OpenPorosity = OpenPorosityy, 
                             cm_in_z = cm_in_zz, cm_in_y = cm_in_yy, cm_in_x = cm_in_xx, Incr_per_z = Incr_per_zz, 
                             Initial_Res = Initial_Ress, RES_INCR = RES_INCRR, CurrFindIterations = CurrFindIterationss, 
                             Cred_neg = Cred_negg, Cox_neg = Cox_negg, Cred_pos = Cred_poss, Cox_pos = Cox_poss, 
                             Decayerneg = Decayernegg, Decayerpos = Decayerposs )
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
