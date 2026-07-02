# coding: utf-8
'''

#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================
#Author:       Tim Tichter
#Years:        2019-2022
#Function:     POLAROGRAPHICAS DRT Function - here, the GUI is set up.
#=======================================================================================================================
#=======================================================================================================================
#importing all required modules from Python
#=======================================================================================================================
#=======================================================================================================================
'''

from tkinter                           import *
from tkinter.filedialog                import askopenfilename
from tkinter.filedialog                import asksaveasfilename
from tkinter                           import messagebox                    # Python 3
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases          import key_press_handler
from matplotlib.figure                 import Figure
from matplotlib.backend_bases          import key_press_handler
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
########################################################################################################################
########################################################################################################################
#  Re-Load all functions required for the actual computationn of the DRTs
########################################################################################################################
########################################################################################################################
from DRT_Native_Gauss_ColeCole_HavNeg  import *
from DRT_Tools_translation             import *
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
#  Define all warnings if inputs are ill-conditioned or do not match
#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================

def DRT_Tools_NNLS_DRT():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    if File_Was_Loaded ==1:
        DRT_Tools_NNLS_DRT_Entry_Window()

def No_Loaded_File_Warner():
    Fenster = Toplevel()                                                         
    Fenster.title("No loaded file found")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "Caution! You loaded no file that\nshould be analyzed. Go to\nopen file first and select\na file that should be analyzed.")


def No_KK_check_possible():
    messagebox.showerror(title="Invalid Inputs", 
                         message="Caution! You want to perform a Kramers-Kronig \ncheck. However, you did not select to compute the \nDRT from both - real and imaginary part. \nThis is, however, mandatory for a KK-check. Your \nrequest for KK-check will be ignored until valid \ninputs are provided.")
     
def No_AntiKK_check_possible():
    messagebox.showerror(title="Invalid Inputs", 
                         message="Caution! You want to do an Anti Kramers-Kronig \ncheck. However, you did not select to compute the \nDRT from both - real and imaginary part. \nThis is, however, mandatory for an anti KK-check. \nYour request for anti KK-check will be ignored \nuntil valid inputs are provided.")

def No_CombiKK_check_possible():
     messagebox.showerror(title="Invalid Inputs", 
                         message="Caution! You want to do an combi Kramers-Kronig \ncheck. However, you did not select to compute the \nDRT from the combined Re-Im fit. \nThis is, however, mandatory for an combi KK-check. \nYour request for combi KK-check will be ignored until \nvalid inputs are provided.")

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Define a function for coloring labels of y-axis individually. This function is copied from
#  https://stackoverflow.com/questions/33159134/matplotlib-y-axis-label-with-multiple-colors
#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================

def multicolor_ylabel(ax,list_of_strings,list_of_colors,axis='x',anchorpad=0,**kw):
    """this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',rotation=90,**kw)) 
                     for text,color in zip(list_of_strings[::-1],list_of_colors) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(1.10, 0.1), 
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_ybox)

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Define the function for opening an EIS file for DRT analysis
#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================

def Open_ImpDRT_File():
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
    global Time
    #==================================================================================
    Frequency         =  data[::Readeveryxx,0:1:1] 
    Z_REAL         =  data[::Readeveryxx,1:2:1] * UmrPot            
    Z_IMAG             =  data[::Readeveryxx,2:3:1] * UmRZ_IMAG
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
        Z_IMAGarrays             = np.squeeze(np.transpose(Z_IMAG[0::,i:i+1:]))                          
        Z_REALarrays         = np.squeeze(np.transpose(Z_REAL[0::,i:i+1:])) 
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

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Define the entry-window how to open an EIS file for DRT analysis
#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================
    
def Get_ImpDRT_Data():
    Fenster = Toplevel()  
    Fenster.geometry("300x225")   
    colorbgr = Label(Fenster, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000) 
    Fenster.resizable(False, False)
    Fenster.title("Get-Data")  
                                                          
    FromRowxxx_label      = Label(Fenster,text="Skip Header", bg = '#80ffbf')
    FromRowxxx_label.grid(row=0, column=0)
    FromRowxxx_Eingabe    = Entry(Fenster)    
    FromRowxxx_Eingabe.insert(END, 0)                                           
    FromRowxxx_Eingabe.grid(row=0, column=1)

    ToRowxxx_label       = Label(Fenster,text="Skip Footer", bg = '#80ffbf')
    ToRowxxx_label.grid(row=1, column=0)
    ToRowxxx_Eingabe     = Entry(Fenster)
    ToRowxxx_Eingabe.insert(END, 0)
    ToRowxxx_Eingabe.grid(row=1, column=1)
    
    Readeveryxxx_label   = Label(Fenster,text="Read every", bg = '#80ffbf')
    Readeveryxxx_label.grid(row=2, column=0)
    Readeveryxxx_Eingabe = Entry(Fenster)
    Readeveryxxx_Eingabe.insert(END, 1)
    Readeveryxxx_Eingabe.grid(row=2, column=1)
    
    UmrPot_Label = Label(Fenster,text="Z_real Factor to be Ohm", bg = '#80ffbf')
    UmrPot_Label.grid(row=3, column=0)
    UmrPot_Eingabe = Entry(Fenster)
    UmrPot_Eingabe.insert(END, 1)
    UmrPot_Eingabe.grid(row=3, column=1)
        
    UmRZ_IMAG_Label = Label(Fenster,text="Z_imag Factor to be Ohm", bg = '#80ffbf')
    UmRZ_IMAG_Label.grid(row=4, column=0)
    UmRZ_IMAG_Eingabe = Entry(Fenster)
    UmRZ_IMAG_Eingabe.insert(END, 1)
    UmRZ_IMAG_Eingabe.grid(row=4, column=1)
    
    HasToBe_Label  = Label(Fenster,text="Data Order", bg = '#80ffbf')
    HasToBe_Label.grid(row=6, column=0)
    HasToBe_Label1  = Label(Fenster,text="Freq___Z.real___Z.Imag", bg = '#80ffbf')
    HasToBe_Label1.grid(row=6, column=1)
        
    Delimiterxxx_Label  = Label(Fenster,text="Delimiter", bg = '#80ffbf')
    Delimiterxxx_Label.grid(row=7, column=0)
    
    var0 = IntVar()
    var0.set(1)
    
    def ChangeDelimiter():
        var0.get()
    
    Radiobutton(Fenster, text='Tab',    padx = 20, variable=var0, value=1,    command = ChangeDelimiter, bg = '#80ffbf').place(x = 118, y = 125, height = 25)
    Radiobutton(Fenster, text='Space',  padx = 20, variable=var0, value=0,  command = ChangeDelimiter, bg = '#80ffbf').place(x = 185, y = 125, height = 25)
    

    def Accept_ImpDRT():
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
        Delimiter         = var0.get()
    def Next_ImpDRT():
        Open_ImpDRT_File()
        def quit():
            Fenster.destroy()
        quit()
    Accept = Button(Fenster, text="Accept",command=Accept_ImpDRT)
    Accept.place( x = 25, y = 160, width = 115, height = 50)  
    Next = Button(Fenster, text="Next",command=Next_ImpDRT)
    Next.place( x = 150, y = 160, width = 115, height = 50)

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Define the inputs for EIS analysis via DRT-tools translation
#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================  
    
def DRT_Tools_NNLS_DRT_Entry_Window():
    Fenster = Toplevel()  
    Fenster.geometry("600x500")
    Fenster.title("DRT_Tools_NNLS_DRT")
    
    Te = Text(Fenster, height=50, width=70, bg = '#ffffff')
    Te.place(x = 230, y = 75, width = 320, height = 350)
    Te.insert(END, "IMPORTANT NOTE:\nThis function is a pythonic \nversion/translation of the great \nopen-source software DRTtools \n(originally written in Matlab). The\nmain difference is, that we use the \nLawson-Hanson (NNLS) algorithm for \noptimization instead of quadprog.\n \nWHEN USING THIS FUNCTION, CHECK OUT\nTHE FOLLOWING ORIGINAL WORK BY!\n\nT. H. Wan, M. Saccoccio, \nC. Chen and F. Ciucci\n\nDOI: 10.1016/j.electacta.2015.09.097\n\nas well as the following web-pages\n\nhttps://ciucci.org/project/drt/\nhttps://github.com/ciuccislab/DRTtools")
    
    Resol_Label = Label(Fenster, text="Increase Data res. by fact.*")
    Resol_Label.grid(row=0, column=0)
    Resol_Eingabe = Entry(Fenster)      
    Resol_Eingabe.insert(END, 10)                                        
    Resol_Eingabe.grid(row=0, column=1)
        
    reg_par_Label = Label(Fenster, text="Tikhonov regul. parameter*")
    reg_par_Label.grid(row=1, column=0)
    reg_par_Eingabe = Entry(Fenster)   
    reg_par_Eingabe.insert(END, 0.01)                                             
    reg_par_Eingabe.grid(row=1, column=1)
    
    epsilon_Label = Label(Fenster,text="epsilon*")
    epsilon_Label.grid(row=2, column=0)
    epsilon_Eingabe = Entry(Fenster)
    epsilon_Eingabe.insert(END, 2.5) 
    epsilon_Eingabe.grid(row=2, column=1)
    
    var1 = IntVar()
    Checkbutton(Fenster, text="Use imaginary part", variable=var1).grid(row=3, column=0, sticky=W)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Use real part", variable=var2).grid(row=4, column=0, sticky=W)
    
    var3 = IntVar(value = 1)
    Checkbutton(Fenster, text="Combined re_im", variable=var3).grid(row=5, column=0, sticky=W)
    
    Sep_Label = Label(Fenster,text="Choose ONE type of radial basis function*")
    Sep_Label.grid(row=6, column=0)
    
    var4 = IntVar(value = 1)
    Checkbutton(Fenster, text="gaussian", variable=var4).grid(row=7, column=0, sticky=W)
    
    var5 = IntVar()
    Checkbutton(Fenster, text="C2_Matern", variable=var5).grid(row=8, column=0, sticky=W)
    
    var6 = IntVar()
    Checkbutton(Fenster, text="C4_Matern", variable=var6).grid(row=9, column=0, sticky=W)
    
    var7 = IntVar()
    Checkbutton(Fenster, text="C6_Matern", variable=var7).grid(row=10, column=0, sticky=W)
    
    var8 = IntVar()
    Checkbutton(Fenster, text="Inv._quadratic", variable=var8).grid(row=11, column=0, sticky=W)
    
    Sep2_Label = Label(Fenster,text="Choose ONE type of derivative order*")
    Sep2_Label.grid(row=12, column=0)
    
    var9 = IntVar(value = 1)
    Checkbutton(Fenster, text="1st_order", variable=var9).grid(row=13, column=0, sticky=W)
    
    var10 = IntVar()
    Checkbutton(Fenster, text="2nd_order", variable=var10).grid(row=14, column=0, sticky=W)
    
    var11 = IntVar()
    Checkbutton(Fenster, text="Save_Data", variable=var11).grid(row=15, column=0, sticky=W)
      
    def AcceptParams():
        global Resol
        global Tikh_Par
        global eps
        global Imag_User
        global Real_User
        global Comb_User
        global deru
        global RBF_Type
        global As_txt_saver
    
        Resol             = int(float(Resol_Eingabe.get()))
        Tikh_Par          = (float(reg_par_Eingabe.get()))
        eps               = (float(epsilon_Eingabe.get()))
        Imag_User         = var1.get()
        Real_User         = var2.get()
        Comb_User         = var3.get()
        gaussrbf          = var4.get()
        c2rbf             = var5.get()
        c4rbf             = var6.get()
        c6rbf             = var7.get()
        invqurbf          = var8.get()
        deru1             = var9.get()
        deru2             = var10.get()
        #rbf_type
        if gaussrbf == 1:
            RBF_Type = 'gaussian'
        if c2rbf == 1:
            RBF_Type = 'C2_matern'
        if c4rbf == 1:
            RBF_Type = 'C4_matern'
        if c6rbf == 1:
            RBF_Type = 'C6_matern'
        if invqurbf  == 1:
            RBF_Type = 'inverse_quadratic'  
        #der_used
        if deru1  == 1:
            deru = '1st-order'
        if deru2  == 1:
            deru = '2nd-order'
        
        As_txt_saver      = var11.get()
        
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=16,column=0)
    button=Button(Fenster,text="Next",command=DRT_Tools_NNLS_DRT_Transformation).grid(row=17,column=0)

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Run the EIS analysis via DRT-tools translation
#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================  

def DRT_Tools_NNLS_DRT_Transformation():
   
    root = Toplevel()
    root.title("DRT_Tools_DRT-Function")
    global Z_imag_Meas
    global Z_real_Meas
    global freq
    Z_imag_Meas     = np.squeeze(Z_IMAG)        #Nach dem Laden heißen die Arrays nur so wie früher
    Z_real_Meas     = np.squeeze(Z_REAL)    #deshalb hier die Umbenennung
    freq            = np.squeeze(Frequency_Array)                
    Z               = Z_real_Meas + 1j*Z_imag_Meas
    freq            = freq[np.imag(Z) < 0]
    Z               = Z[np.imag(Z) < 0]
    ##############################################################################
    #initialise parameters and bounds
    epsilon       = eps
    rbf_type      = RBF_Type
    der_used      = deru
    damp          = Tikh_Par
    taumax        = ceil(np.max(np.log10(1./freq)))+1
    taumin        = floor(np.min(np.log10(1./freq)))-1
    freq_out      = np.logspace(-taumin,-taumax,Resol*len(freq))
    ##############################################################################
    A_im          = assemble_A_im(freq,epsilon,rbf_type)
    M_im          = assemble_M_im(freq,epsilon,rbf_type,der_used)
    H_im, f_im    = quad_format(A_im,-Z.imag,M_im,damp)
    A_re          = assemble_A_re(freq,epsilon,rbf_type)
    M_re          = assemble_M_re(freq,epsilon,rbf_type,der_used)
    H_re, f_re    = quad_format(A_re,(Z.real),M_re,damp)
    H_co, f_co    = quad_format_combined(A_re, A_im,(Z.real), -Z.imag, M_re, M_im,damp)
    x_imag_2      = nnls(H_im,np.abs(f_im))[0]
    x_real_2      = nnls(H_re,np.abs(f_re))[0]
    x_comb_2      = nnls(H_co,np.abs(f_co))[0] 
    global gamma_imag_fine
    global gamma_real_fine
    global gamma_comb_fine
    global TAUARRAY
    TAUARRAY         = 1./freq_out 

    gamma_imag_fine  = map_array_to_gamma(freq_out,freq,x_imag_2[2:],epsilon,rbf_type)
    gamma_real_fine  = map_array_to_gamma(freq_out,freq,x_real_2[2:],epsilon,rbf_type)
    gamma_comb_fine  = map_array_to_gamma(freq_out,freq,x_comb_2[2:],epsilon,rbf_type)
    
    Abbildung = Figure(figsize=(7, 5), dpi=100)
    b = Abbildung.add_subplot(111)
    if Imag_User == 1:
        b.plot(np.log10(TAUARRAY),gamma_imag_fine,color='b', linestyle ='-',label = 'Im_DRT')
    if Real_User == 1:
        b.plot(np.log10(TAUARRAY),gamma_real_fine,color='r', linestyle ='-', label = 'Re_DRT')
    if Comb_User == 1:
        b.plot(np.log10(TAUARRAY),gamma_comb_fine,color='k', linestyle ='-', label = 'Comb_DRT')
    b.legend(frameon = False, fontsize = 15)    
    b.set_xlabel('log$_{10}\,(\\tau\,/\,\mathrm{s})$', fontsize=15)
    b.set_ylabel('$\Delta Z_{\mathrm{Real}}\,\\tau\,\gamma(\\tau)\,/\,\Omega$', fontsize=15)
    for axis in ['top','bottom','left','right']:
        b.spines[axis].set_linewidth(2)
        b.spines[axis].set_color('k')
    b.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
    secax = b.secondary_xaxis('top', functions=(log_to_ln, ln_to_log))
    secax.set_xlabel('ln$\,(\\tau\,/\,\mathrm{s})$', fontsize = 15)
    secax.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    if As_txt_saver == 1:
        DRT_Tools_NNLS_DRT_as_Txt_Saver()

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Save the result from EIS analysis via DRT-tools translation
#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================     
        
def DRT_Tools_NNLS_DRT_as_Txt_Saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))
    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n") 
        fi.write("DRT calc. from")
        fi.write("\n") 
        fi.write("Imaginary-Part:")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("yes")
        if Imag_User != 1:
            fi.write("no")
        fi.write("\n") 
        fi.write("Real-Part:")
        fi.write("\t")
        if Real_User == 1:
            fi.write("yes")
        if Real_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Combined Real_Im:")
        fi.write("\t")
        if Comb_User == 1:
            fi.write("yes")
        if Comb_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Resolution increase")
        fi.write("\t")
        fi.write(str(Resol))
        fi.write("\n")
        fi.write("Reg_Par")
        fi.write("\t")
        fi.write(str(Tikh_Par))
        fi.write("\n")
        fi.write("epsilon")
        fi.write("\t")
        fi.write(str(eps))
        fi.write("\n")
        fi.write("RBF-Type:")
        fi.write("\t")
        if RBF_Type == 'gaussian':
            fi.write("gaussian")
        if RBF_Type == 'C2_matern':
            fi.write("C2_matern")
        if RBF_Type == 'C4_matern':
            fi.write("C4_matern")
        if RBF_Type == 'C6_matern':
            fi.write("C6_matern")
        if RBF_Type == 'inverse_quadratic':
            fi.write("inverse_quadratic")
        fi.write("\n")
        fi.write("Derivative-Order:")
        fi.write("\t")
        if deru == '1st-order':
            fi.write("1st-order")
        if deru == '2nd-order':
            fi.write("2nd-order")
        fi.write("\n")
        fi.write("\n")
        fi.write("log10(tau / s)")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("gamma(tau)_Imag")
            fi.write("\t")
        if Real_User == 1:
            fi.write("gamma(tau)_Real")
            fi.write("\t")
        if Comb_User == 1:
            fi.write("gamma(tau)_Combined")
        fi.write("\n")
        for i in range(len(TAUARRAY)):
            fi.write(str(np.asscalar(np.log10(TAUARRAY[i]))))
            fi.write("\t")
            if Imag_User == 1:
                fi.write(str(np.asscalar(gamma_imag_fine[i])))
                fi.write("\t")
            if Real_User == 1:
                fi.write(str(np.asscalar(gamma_real_fine[i])))
                fi.write("\t")
            if Comb_User == 1:
                fi.write(str(np.asscalar(gamma_comb_fine[i])))
                fi.write("\t")
            fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency [Hz]")
        fi.write("\t")
        fi.write("Z_real_Meas [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_Meas [Ohm]")
        fi.write("\t")
        fi.write("\n")
        for i in range(len(Z_real_Meas)):
            fi.write(str(np.asscalar(Frequency_Array[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_real_Meas[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_imag_Meas[i])))
            fi.write("\t")
            fi.write("\n")
    root.destroy()

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
#     Here, everything for the GUI of the Native, the Gaussian and the ColeCole DRT is defined
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################


#============================================================================================================
#============================================================================================================
#============================================================================================================   
#  Here, everything for the GUI of the native DRT is defined
#============================================================================================================ 
#============================================================================================================
#============================================================================================================
    
def Native_NNLS_DRT():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    if File_Was_Loaded ==1:
        Native_DRT_Entry_Window()

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Define the inputs for EIS analysis via Native DRT
#=======================================================================================================================    
#=======================================================================================================================
#======================================================================================================================= 

def Native_DRT_Entry_Window():
    Fenster = Toplevel()  
    Fenster.geometry("550x300")
    colorbgr = Label(Fenster, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000) 
    Fenster.resizable(False, False)
    Fenster.title("Native NNLS DRT")
     
    Resol_Label = Label(Fenster, text="Increase Data res. by fact.*", bg = '#80ffbf')
    Resol_Label.place(x = 25, y = 10)
    Resol_Eingabe = Entry(Fenster)      
    Resol_Eingabe.insert(END, 10)                                        
    Resol_Eingabe.place(x = 185, y = 10, width = 80, height = 19)
    
    extend_low_Label = Label(Fenster,text="Extend low*", bg = '#80ffbf')
    extend_low_Label.place(x = 25, y = 35)
    extend_low_Eingabe = Entry(Fenster)
    extend_low_Eingabe.insert(END, 0) 
    extend_low_Eingabe.place(x = 185, y = 35, width = 80, height = 19)
        
    extend_high_Label = Label(Fenster,text="Extend high*", bg = '#80ffbf')
    extend_high_Label.place(x = 25, y = 60)
    extend_high_Eingabe = Entry(Fenster)
    extend_high_Eingabe.insert(END, 0) 
    extend_high_Eingabe.place(x = 185, y = 60, width = 80, height = 19)
    
    reg_par_Label = Label(Fenster, text="Tikhonov parameter*", bg = '#80ffbf')
    reg_par_Label.place(x = 300, y = 10)
    reg_par_Eingabe = Entry(Fenster)   
    reg_par_Eingabe.insert(END, 0.01)                                             
    reg_par_Eingabe.place(x = 440, y = 10, width = 80, height = 19)
    
    SepLine0 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine0.place(x = 0, y = 85, height = 10)
    
    Lab1 = Label(Fenster,text="DRT from:", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab1.place(x = 25, y = 100 )

    var1 = IntVar()
    Checkbutton(Fenster, text="Imaginary part", variable=var1, bg = '#80ffbf').place(x = 25, y = 180)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Real part", variable=var2, bg = '#80ffbf').place(x = 25, y = 155)
    
    var3 = IntVar(value = 1)
    Checkbutton(Fenster, text="Combined fit", variable=var3, bg = '#80ffbf').place(x = 25, y = 130)
    
    var0 = StringVar()
    var0.set('linear')  # initializing the choice, i.e. Python
    
    def ChangePenalty():
        var0.get()

    Lab2 = Label(Fenster,text="Penalty-Type", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab2.place(x = 182, y = 100 )
    
    Radiobutton(Fenster, text='linear',    padx = 20, variable=var0, value='linear',    command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 130, height = 25)
    Radiobutton(Fenster, text='1st-order', padx = 20, variable=var0, value='1st-order', command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 155, height = 25)
    Radiobutton(Fenster, text='2nd-order', padx = 20, variable=var0, value='2nd-order', command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 180, height = 25)
    
    var5 = StringVar()
    var5.set(None)  # initializing the choice, i.e. Python
    
    def GetBackCalckType():
        var5.get()

    Lab2 = Label(Fenster,text="Calc. back to Z", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab2.place(x = 330, y = 100 )
    
    Radiobutton(Fenster, text='Kramers-Kronig like', padx = 20, variable=var5, value='KramersKronig',        command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 130, height = 25)
    Radiobutton(Fenster, text='Anti Kramers-Kronig', padx = 20, variable=var5, value='AntiKramersKronig',    command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 155, height = 25)
    Radiobutton(Fenster, text='Combi Kramers-Kronig', padx = 20, variable=var5, value='CombinedBack',         command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 180, height = 25)
    
    SepLine1 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine1.place(x = 0, y = 205, height = 10)
                
    var6 = IntVar()
    Checkbutton(Fenster, text="Remove Real-offset", variable=var6, bg = '#80ffbf').place(x = 25, y = 230)
    
    var4 = IntVar()
    Checkbutton(Fenster, text="Save Data", variable=var4, bg = '#80ffbf').place(x = 25, y = 255)         
                
    def AcceptParams():
        global Resol
        global Tikh_Par
        global ext_rangelow
        global ext_rangehigh
        global Imag_User
        global Real_User
        global Comb_User
        global As_txt_saver
        global PenaltyType
        global RemoveRealOffset
        global BackCalculation
        global KK_enable
        global AntiKK_enable
        global CombiKK_enable
        Resol             = int(float(Resol_Eingabe.get()))
        Tikh_Par          = (float(reg_par_Eingabe.get()))
        ext_rangelow      = (float(extend_low_Eingabe.get()))
        ext_rangehigh     = (float(extend_high_Eingabe.get()))
        Imag_User         = var1.get()
        Real_User         = var2.get()
        Comb_User         = var3.get()
        As_txt_saver      = var4.get()
        PenaltyType       = var0.get()
        RemoveRealOffset  = var6.get()
        BackCalculation   = var5.get()
        KK_enable         = 0
        AntiKK_enable     = 0
        CombiKK_enable    = 0
        if BackCalculation == 'KramersKronig' and Imag_User == 0 or BackCalculation == 'KramersKronig' and Real_User == 0 :
            No_KK_check_possible()
        if BackCalculation == 'KramersKronig' and Imag_User == 1 and Real_User == 1:
            KK_enable         = 1
        if BackCalculation == 'AntiKramersKronig' and Imag_User == 0 or BackCalculation == 'AntiKramersKronig' and Real_User == 0:
            No_AntiKK_check_possible()
        if BackCalculation == 'AntiKramersKronig' and Imag_User == 1 and Real_User == 1:
            AntiKK_enable     = 1
        if BackCalculation == 'CombinedBack' and Comb_User == 0:
            No_CombiKK_check_possible()
        if BackCalculation == 'CombinedBack' and Comb_User == 1:
            CombiKK_enable    = 1
                
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).place( x = 185, y = 230, width = 130, height = 50)
    button=Button(Fenster,text="Next",command= Native_DRT_Transformation).place(x = 335, y = 230, width = 130, height = 50)

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Run the EIS analysis via Native DRT
#=======================================================================================================================    
#=======================================================================================================================
#======================================================================================================================= 

def Native_DRT_Transformation():
    root = Toplevel()
    root.title("Result native DRT-Transformation")
    global TAUARRAY
    global gamma_imag_spike
    global gamma_real_spike
    global gamma_comb_spike
    global Z_imag_Meas
    global Z_real_Meas
    global freq
    global Z_Back
    global Z_meas
    Z_imag_Meas     = np.squeeze(Z_IMAG)               #Nach dem Laden heißen die Arrays nur so wie früher
    Z_real_Meas     = np.squeeze(Z_REAL)           #deshalb hier die Umbenennung
    freq            = np.squeeze(Frequency_Array) 
    Z               = Z_real_Meas + 1j*Z_imag_Meas
    freq            = freq[np.imag(Z) < 0]
    Z               = Z[np.imag(Z) < 0]
    Z_meas          = Z.copy()
    Real_Offset = np.min(Z.real).copy()
    if RemoveRealOffset == 1:
        Z = Z - np.min(Z.real)
    if Imag_User == 1:
        IM_Output        = DRT_Native(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "imag", Reg_Type = PenaltyType)
        TAUARRAY         = 10**IM_Output[0]
        gamma_imag_spike = IM_Output[1]
    if Real_User == 1:
        RE_Output        = DRT_Native(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "real", Reg_Type = PenaltyType)
        TAUARRAY         = 10**RE_Output[0]
        gamma_real_spike = RE_Output[1]
    if Comb_User == 1:
        CO_Output        = DRT_Native(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "comb", Reg_Type = PenaltyType)
        TAUARRAY         = 10**CO_Output[0]
        gamma_comb_spike = CO_Output[1]   
    ##########################################################################################
    # Back-calculation from DRTs to impedances in any permutation, KK, anti KK, combi
    ########################################################################################## 
    if KK_enable == 1:
        Z_Back           = KramersKronig_Multiplication(MATRIX = IM_Output[4], Real_DRT = RE_Output[1], Imag_DRT = IM_Output[1])
        Z_Back = Z_Back + Real_Offset     
    if AntiKK_enable == 1:
        Z_Back           = AntiKramersKronig_Multiplication(MATRIX = IM_Output[4], Real_DRT = RE_Output[1], Imag_DRT = IM_Output[1])
        if RemoveRealOffset == 1:
            Z_Back = Z_Back + Real_Offset
    if CombiKK_enable == 1:
        Z_Back           = CombiKramersKronig_Multiplication(MATRIX = CO_Output[4], Comb_DRT = CO_Output[1])
        if RemoveRealOffset == 1:
            Z_Back = Z_Back + Real_Offset
    ##########################################################################################
    #  Plot the DRT
    ##########################################################################################
    if BackCalculation == "None" or BackCalculation == 'KramersKronig' and KK_enable == 0 or BackCalculation == 'AntiKramersKronig' and AntiKK_enable == 0 or BackCalculation == 'CombinedBack' and CombiKK_enable == 0:
        Abbildung = Figure(figsize=(7, 5), dpi=100)
        b = Abbildung.add_subplot(111)
        if Imag_User == 1:
            b.plot(IM_Output[0],IM_Output[1],color='blue', linestyle ='-' ,label = 'Imag-DRT')
        if Real_User == 1:
            b.plot(RE_Output[0],RE_Output[1],color='red', linestyle ='-'  ,label = 'Real-DRT')
        if Comb_User == 1:
            b.plot(CO_Output[0],CO_Output[1],color='black', linestyle ='-',label = 'Combi-DRT')   
        b.legend(frameon = False, fontsize = 15)    
        b.set_xlabel('log$_{10}\,(\\tau\,/\,\mathrm{s})$', fontsize=15)
        b.set_ylabel('$\Delta Z_{\mathrm{Real}}\,\\tau\,\gamma(\\tau)\,/\,\Omega$', fontsize=15)
        for axis in ['top','bottom','left','right']:
            b.spines[axis].set_linewidth(2)
            b.spines[axis].set_color('k')
        b.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        secax = b.secondary_xaxis('top', functions=(log_to_ln, ln_to_log))
        secax.set_xlabel('ln$\,(\\tau\,/\,\mathrm{s})$', fontsize = 15)
        secax.tick_params(tickdir='in', width = 2, length=6, labelsize=15)  
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    ##########################################################################################
    #  Plot the back-calculation result (in case it was activated)
    ##########################################################################################
    if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
        #============================================================================
        #Plotting of DRT and back-computation
        #============================================================================ 
        f  = Figure(figsize=(11, 6), dpi=100)        
        b1 = f.add_subplot(121)
        b2 = f.add_subplot(122) 
        #=========================================================
        #final plotting
        #=========================================================
        if Imag_User == 1:
            b1.plot(IM_Output[0],IM_Output[1],color='blue', linestyle ='-' ,label = 'Imag-DRT')
        if Real_User == 1:
            b1.plot(RE_Output[0],RE_Output[1],color='red', linestyle ='-'  ,label = 'Real-DRT')
        if Comb_User == 1:
            b1.plot(CO_Output[0],CO_Output[1],color='black', linestyle ='-',label = 'Combi-DRT')
        b1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=2, fontsize = 15)
        b1.set_xlabel('log$_{10}\,(\\tau\,/\,\mathrm{s})$', fontsize=15)
        b1.set_ylabel('$\Delta Z_{\mathrm{Real}}\,\\tau\,\gamma(\\tau)\,/\,\Omega$', fontsize=15)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        secax = b1.secondary_xaxis('top', functions=(log_to_ln, ln_to_log))
        secax.set_xlabel('ln$\,(\\tau\,/\,\mathrm{s})$', fontsize = 15)
        secax.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        if KK_enable == 1:
            LABEL_y_a = "Im from Re-DRT"
            LABEL_y_b = "Re from Im-DRT"
        if AntiKK_enable == 1:
            LABEL_y_a = "Im from Im-DRT"
            LABEL_y_b = "Re from Re-DRT"
        if CombiKK_enable == 1:
            LABEL_y_a = "Im from Comb-DRT"
            LABEL_y_b = "Re from Comb-DRT"      
        b2.plot(np.log10(freq), -Z_Back.imag, linestyle='-',marker='', linewidth = 0.5, color='blue', label = LABEL_y_a)
        b2.plot(np.log10(freq),  Z_Back.real, linestyle='-',marker='', linewidth = 0.5, color='red',  label = LABEL_y_b)
        b2.plot(np.log10(freq), -Z_meas.imag, linestyle='',marker='.', markersize = 3,  color='blue', label = 'Im-meas')
        b2.plot(np.log10(freq),  Z_meas.real, linestyle='',marker='.', markersize = 3,  color='red',  label = 'Re-meas')
        ax2 = b2.twinx()
        ax2.plot(np.log10(freq), 100*np.abs((Z_Back.imag-Z_meas.imag)/Z_meas.imag), color = 'cyan',   linewidth = 0.5, label = '$\Delta$Im')
        ax2.plot(np.log10(freq), 100*np.abs((Z_Back.real-Z_meas.real)/Z_meas.real), color = 'orange', linewidth = 0.5, label = '$\Delta$Re')
        ax2.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        b2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=2, fontsize = 15)
        b2.set_xlabel('log$_{10}\,(\\nu\,/\,\mathrm{s}^{-1})$', fontsize=15)
        b2.set_ylabel('$-\mathfrak{IM}(Z(j\omega)/\Omega)$ and $\mathfrak{RE}(Z(j\omega)/\Omega)$', fontsize = 15)
        multicolor_ylabel(ax2,('$\Delta Z_{\mathrm{Re}}(j\omega))$','and','$\Delta Z_{\mathrm{Im}}(j\omega))$', 'in %'),('black', 'cyan','black','orange'), axis='y',size=15, weight='bold')
        for axis in ['top','bottom','left','right']:
            b2.spines[axis].set_linewidth(2)
            b2.spines[axis].set_color('k')
        b2.tick_params(tickdir='in', width = 2, length=6, labelsize=15)       
        #========================================================= 
        canvas = FigureCanvasTkAgg(f, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        f.tight_layout()
    if As_txt_saver == 1:
        Native_DRT_as_Txt_Saver()

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Save the result from the EIS analysis via Native DRT
#=======================================================================================================================    
#=======================================================================================================================
#======================================================================================================================= 
        
def Native_DRT_as_Txt_Saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))
    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")   
        fi.write("Native (spike)-DRT calc. from")
        fi.write("\n") 
        fi.write("Imaginary-Part:")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("yes")
        if Imag_User != 1:
            fi.write("no")
        fi.write("\n") 
        fi.write("Real-Part:")
        fi.write("\t")
        if Real_User == 1:
            fi.write("yes")
        if Real_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Combined Real_Im:")
        fi.write("\t")
        if Comb_User == 1:
            fi.write("yes")
        if Comb_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Real-Offset-removal")
        fi.write("\t")
        if RemoveRealOffset == 0:
            fi.write("no")
        if RemoveRealOffset == 1:
            fi.write("yes")
        fi.write("\n")
        fi.write("Resolution increase")
        fi.write("\t")
        fi.write(str(Resol))
        fi.write("\n")
        fi.write("Lower Range Extension")
        fi.write("\t")
        fi.write(str(ext_rangelow))
        fi.write("\n")
        fi.write("Upper Range Extension")
        fi.write("\t")
        fi.write(str(ext_rangehigh))
        fi.write("\n")
        fi.write("Reg_Par")
        fi.write("\t")
        fi.write(str(Tikh_Par))
        fi.write("\n")
        fi.write("Regularization type")
        fi.write("\t")
        if PenaltyType == "linear":
            fi.write("linear regularization")
            fi.write("\n")
        if PenaltyType == "1st-order":
            fi.write("1st-derivative regularization")
            fi.write("\n")
        if PenaltyType == "2nd-order":
            fi.write("2nd-derivative regularization")
            fi.write("\n")       
        fi.write("Back-calculation")
        fi.write("\t")
        if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
            if KK_enable == 1:
                fi.write("Kramers-Kronig")
                fi.write("\n")
            if AntiKK_enable == 1:
                fi.write("Anti-Kramers-Kronig")
                fi.write("\n")
            if CombiKK_enable == 1:
                fi.write("Combi-Kramers-Kronig")
                fi.write("\n")    
        if KK_enable == 0 and AntiKK_enable == 0 and CombiKK_enable == 0:
            fi.write("None")
            fi.write("\n") 
        fi.write("\n")
        fi.write("log10(tau / s)")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("gamma(tau)_Imag")
            fi.write("\t")
        if Real_User == 1:
            fi.write("gamma(tau)_Real")
            fi.write("\t")
        if Comb_User == 1:
            fi.write("gamma(tau)_Combined")
        fi.write("\n")       
        for i in range(len(TAUARRAY)):
            fi.write(str(np.asscalar(np.log10(TAUARRAY[i]))))
            fi.write("\t")
            if Imag_User == 1:
                fi.write(str(np.asscalar(gamma_imag_spike[i])))
                fi.write("\t")
            if Real_User == 1:
                fi.write(str(np.asscalar(gamma_real_spike[i])))
                fi.write("\t")
            if Comb_User == 1:
                fi.write(str(np.asscalar(gamma_comb_spike[i])))
                fi.write("\t")
            fi.write("\n")  
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency [Hz]")
        fi.write("\t")
        fi.write("Z_real_used [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_used [Ohm]")
        fi.write("\t")
        if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
            if KK_enable == 1:
                fi.write("Z_real_from_ImDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_ReDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")
            if AntiKK_enable == 1:
                fi.write("Z_real_from_ReDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_ImDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")
            if CombiKK_enable == 1:
                fi.write("Z_real_from_CombDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_CombDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")       
        fi.write("\n")        
        for i in range(len(freq)):
            fi.write(str(np.asscalar(freq[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_meas[i].real)))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_meas[i].imag)))
            fi.write("\t")
            if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
                fi.write(str(np.asscalar(Z_Back[i].real)))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_Back[i].imag)))
                fi.write("\t")
                fi.write(str(100*np.asscalar(np.abs((Z_Back[i].real-Z_meas[i].real)/Z_meas[i].real))))
                fi.write("\t")
                fi.write(str(100*np.asscalar(np.abs((Z_Back[i].imag-Z_meas[i].imag)/Z_meas[i].imag))))
                fi.write("\t")
            fi.write("\n") 
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency_original [Hz]")
        fi.write("\t")
        fi.write("Z_real_original [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_originnal [Ohm]")
        fi.write("\n")
        for i in range(len(Frequency_Array)):
            fi.write(str(np.asscalar(Frequency_Array[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_real_Meas[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_imag_Meas[i])))
            fi.write("\n")
    root.destroy()
    
#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Define the inputs for EIS analysis via Gaussian DRT
#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================      
    
def Gaussian_NNLS_DRT():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    if File_Was_Loaded ==1:
        Gaussian_DRT_Entry_Window()

def Gaussian_DRT_Entry_Window():
    Fenster = Toplevel()  
    Fenster.geometry("550x300")
    colorbgr = Label(Fenster, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000) 
    Fenster.resizable(False, False)
    Fenster.title("Gaussian NNLS DRT")
     
    Resol_Label = Label(Fenster, text="Increase Data res. by fact.*", bg = '#80ffbf')
    Resol_Label.place(x = 25, y = 10)
    Resol_Eingabe = Entry(Fenster)      
    Resol_Eingabe.insert(END, 10)                                        
    Resol_Eingabe.place(x = 185, y = 10, width = 80, height = 19)
    
    extend_low_Label = Label(Fenster,text="Extend low*", bg = '#80ffbf')
    extend_low_Label.place(x = 25, y = 35)
    extend_low_Eingabe = Entry(Fenster)
    extend_low_Eingabe.insert(END, 0) 
    extend_low_Eingabe.place(x = 185, y = 35, width = 80, height = 19)
        
    extend_high_Label = Label(Fenster,text="Extend high*", bg = '#80ffbf')
    extend_high_Label.place(x = 25, y = 60)
    extend_high_Eingabe = Entry(Fenster)
    extend_high_Eingabe.insert(END, 0) 
    extend_high_Eingabe.place(x = 185, y = 60, width = 80, height = 19)
    
    reg_par_Label = Label(Fenster, text="Tikhonov parameter*", bg = '#80ffbf')
    reg_par_Label.place(x = 300, y = 10)
    reg_par_Eingabe = Entry(Fenster)   
    reg_par_Eingabe.insert(END, 0.01)                                             
    reg_par_Eingabe.place(x = 440, y = 10, width = 80, height = 19)
    
    FWHM_Label = Label(Fenster, text="FWHM of RBF*", bg = '#80ffbf')
    FWHM_Label.place(x = 300, y = 35)
    FWHM_Eingabe = Entry(Fenster)   
    FWHM_Eingabe.insert(END, 0.25)                                             
    FWHM_Eingabe.place(x = 440, y = 35, width = 80, height = 19)
    
    SepLine0 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine0.place(x = 0, y = 85, height = 10)
    
    Lab1 = Label(Fenster,text="DRT from:", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab1.place(x = 25, y = 100 )

    var1 = IntVar()
    Checkbutton(Fenster, text="Imaginary part", variable=var1, bg = '#80ffbf').place(x = 25, y = 180)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Real part", variable=var2, bg = '#80ffbf').place(x = 25, y = 155)
    
    var3 = IntVar(value = 1)
    Checkbutton(Fenster, text="Combined fit", variable=var3, bg = '#80ffbf').place(x = 25, y = 130)
    
    var0 = StringVar()
    var0.set('linear')  # initializing the choice, i.e. Python
    
    def ChangePenalty():
        var0.get()

    Lab2 = Label(Fenster,text="Penalty-Type", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab2.place(x = 182, y = 100 )
    
    Radiobutton(Fenster, text='linear',    padx = 20, variable=var0, value='linear',    command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 130, height = 25)
    Radiobutton(Fenster, text='1st-order', padx = 20, variable=var0, value='1st-order', command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 155, height = 25)
    Radiobutton(Fenster, text='2nd-order', padx = 20, variable=var0, value='2nd-order', command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 180, height = 25)
    
    var5 = StringVar()
    var5.set(None)  # initializing the choice, i.e. Python
    
    def GetBackCalckType():
        var5.get()

    Lab2 = Label(Fenster,text="Calc. back to Z", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab2.place(x = 330, y = 100 )
    
    Radiobutton(Fenster, text='Kramers-Kronig like', padx = 20, variable=var5, value='KramersKronig',        command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 130, height = 25)
    Radiobutton(Fenster, text='Anti Kramers-Kronig', padx = 20, variable=var5, value='AntiKramersKronig',    command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 155, height = 25)
    Radiobutton(Fenster, text='Combi Kramers-Kronig', padx = 20, variable=var5, value='CombinedBack',         command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 180, height = 25)
    
    SepLine1 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine1.place(x = 0, y = 205, height = 10)
    
    var6 = IntVar()
    Checkbutton(Fenster, text="Remove Real-offset", variable=var6, bg = '#80ffbf').place(x = 25, y = 230)
    
    var4 = IntVar()
    Checkbutton(Fenster, text="Save Data", variable=var4, bg = '#80ffbf').place(x = 25, y = 255)    
    
    def AcceptParams():
        global Resol
        global Tikh_Par
        global ext_rangelow
        global ext_rangehigh
        global Imag_User
        global Real_User
        global Comb_User
        global As_txt_saver
        global Gaussain_FWHM
        global PenaltyType
        global RemoveRealOffset
        global BackCalculation
        global KK_enable
        global AntiKK_enable
        global CombiKK_enable
        Resol             = int(float(Resol_Eingabe.get()))
        Tikh_Par          = (float(reg_par_Eingabe.get()))
        ext_rangelow      = (float(extend_low_Eingabe.get()))
        ext_rangehigh     = (float(extend_high_Eingabe.get()))
        Imag_User         = var1.get()
        Real_User         = var2.get()
        Comb_User         = var3.get()
        As_txt_saver      = var4.get()   
        Gaussain_FWHM     = (float(FWHM_Eingabe.get()))
        PenaltyType       = var0.get()
        RemoveRealOffset  = var6.get()
        BackCalculation   = var5.get()
        KK_enable         = 0
        AntiKK_enable     = 0
        CombiKK_enable    = 0
        if BackCalculation == 'KramersKronig' and Imag_User == 0 or BackCalculation == 'KramersKronig' and Real_User == 0 :
            No_KK_check_possible()
        if BackCalculation == 'KramersKronig' and Imag_User == 1 and Real_User == 1:
            KK_enable         = 1
        if BackCalculation == 'AntiKramersKronig' and Imag_User == 0 or BackCalculation == 'AntiKramersKronig' and Real_User == 0:
            No_AntiKK_check_possible()
        if BackCalculation == 'AntiKramersKronig' and Imag_User == 1 and Real_User == 1:
            AntiKK_enable     = 1
        if BackCalculation == 'CombinedBack' and Comb_User == 0:
            No_CombiKK_check_possible()
        if BackCalculation == 'CombinedBack' and Comb_User == 1:
            CombiKK_enable    = 1
         
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).place( x = 185, y = 230, width = 130, height = 50)
    button=Button(Fenster,text="Next",command=Gaussian_DRT_Transformation).place(x = 335, y = 230, width = 130, height = 50)

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Run the EIS analysis via Gaussian DRT
#=======================================================================================================================    
#=======================================================================================================================
#======================================================================================================================= 

def Gaussian_DRT_Transformation():
    root = Toplevel()
    root.title("Result Gaussian DRT-Transformation")
    global TAUARRAY
    global gamma_imag_gauss
    global gamma_real_gauss
    global gamma_comb_gauss
    global Z_imag_Meas
    global Z_real_Meas
    global freq
    global Z_Back   
    global Z_meas
    Z_imag_Meas     = np.squeeze(Z_IMAG)               #Nach dem Laden heißen die Arrays nur so wie früher
    Z_real_Meas     = np.squeeze(Z_REAL)           #deshalb hier die Umbenennung
    freq            = np.squeeze(Frequency_Array)
    Z               = Z_real_Meas + 1j*Z_imag_Meas
    freq            = freq[np.imag(Z) < 0]
    Z               = Z[np.imag(Z) < 0]
    Z_meas          = Z.copy()
    Real_Offset = np.min(Z.real).copy()
    if RemoveRealOffset == 1:
        Z = Z - np.min(Z.real)  
    if Imag_User == 1:
        IM_Output        = DRT_Gaussian_RBF(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "imag", FWHM = Gaussain_FWHM, Reg_Type = PenaltyType)
        TAUARRAY         = 10**IM_Output[0]
        gamma_imag_gauss = IM_Output[2]
    if Real_User == 1:
        RE_Output        = DRT_Gaussian_RBF(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "real", FWHM = Gaussain_FWHM, Reg_Type = PenaltyType)
        TAUARRAY         = 10**RE_Output[0]
        gamma_real_gauss = RE_Output[2]
    if Comb_User == 1:
        CO_Output        = DRT_Gaussian_RBF(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "comb", FWHM = Gaussain_FWHM, Reg_Type = PenaltyType)
        TAUARRAY         = 10**CO_Output[0]
        gamma_comb_gauss = CO_Output[2]     
    ##########################################################################################
    # Back-calculation from DRTs to impedances in any permutation, KK, anti KK, combi
    ########################################################################################## 
    if KK_enable == 1:
        Z_Back           = KramersKronig_Multiplication(MATRIX = IM_Output[5], Real_DRT = RE_Output[2], Imag_DRT = IM_Output[2])
        Z_Back = Z_Back + Real_Offset     
    if AntiKK_enable == 1:
        Z_Back           = AntiKramersKronig_Multiplication(MATRIX = IM_Output[5], Real_DRT = RE_Output[2], Imag_DRT = IM_Output[2])
        if RemoveRealOffset == 1:
            Z_Back = Z_Back + Real_Offset
    if CombiKK_enable == 1:
        Z_Back           = CombiKramersKronig_Multiplication(MATRIX = CO_Output[5], Comb_DRT = CO_Output[2])
        if RemoveRealOffset == 1:
            Z_Back = Z_Back + Real_Offset
    ##########################################################################################
    #  Plot the DRT
    ##########################################################################################
    if BackCalculation == "None" or BackCalculation == 'KramersKronig' and KK_enable == 0 or BackCalculation == 'AntiKramersKronig' and AntiKK_enable == 0 or BackCalculation == 'CombinedBack' and CombiKK_enable == 0:
        Abbildung = Figure(figsize=(7, 5), dpi=100)
        b = Abbildung.add_subplot(111)
        if Imag_User == 1:
            b.plot(IM_Output[0],IM_Output[2],color='blue', linestyle ='-' ,label = 'Imag-DRT')
        if Real_User == 1:
            b.plot(RE_Output[0],RE_Output[2],color='red', linestyle ='-'  ,label = 'Real-DRT')
        if Comb_User == 1:
            b.plot(CO_Output[0],CO_Output[2],color='black', linestyle ='-',label = 'Combi-DRT')     
        b.legend(frameon = False, fontsize = 15)    
        b.set_xlabel('log$_{10}\,(\\tau\,/\,\mathrm{s})$', fontsize=15)
        b.set_ylabel('$\Delta Z_{\mathrm{Real}}\,\\tau\,\gamma(\\tau)\,/\,\Omega$', fontsize=15)
        for axis in ['top','bottom','left','right']:
            b.spines[axis].set_linewidth(2)
            b.spines[axis].set_color('k')
        b.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        secax = b.secondary_xaxis('top', functions=(log_to_ln, ln_to_log))
        secax.set_xlabel('ln$\,(\\tau\,/\,\mathrm{s})$', fontsize = 15)
        secax.tick_params(tickdir='in', width = 2, length=6, labelsize=15)   
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1) 
    ##########################################################################################
    #  Plot the back-calculation result (in case it was activated)
    ##########################################################################################
    if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
        #============================================================================
        #Plotting of DRT and back-computation
        #============================================================================ 
        f  = Figure(figsize=(11, 6), dpi=100)        
        b1 = f.add_subplot(121)
        b2 = f.add_subplot(122) 
        #=========================================================
        #final plotting
        #=========================================================
        if Imag_User == 1:
            b1.plot(IM_Output[0],IM_Output[2],color='blue', linestyle ='-' ,label = 'Imag-DRT')
        if Real_User == 1:
            b1.plot(RE_Output[0],RE_Output[2],color='red', linestyle ='-'  ,label = 'Real-DRT')
        if Comb_User == 1:
            b1.plot(CO_Output[0],CO_Output[2],color='black', linestyle ='-',label = 'Combi-DRT')
        b1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=2, fontsize = 15)   
        b1.set_xlabel('log$_{10}\,(\\tau\,/\,\mathrm{s})$', fontsize=15)
        b1.set_ylabel('$\Delta Z_{\mathrm{Real}}\,\\tau\,\gamma(\\tau)\,/\,\Omega$', fontsize=15)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        secax = b1.secondary_xaxis('top', functions=(log_to_ln, ln_to_log))
        secax.set_xlabel('ln$\,(\\tau\,/\,\mathrm{s})$', fontsize = 15)
        secax.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        if KK_enable == 1:
            LABEL_y_a = "Im from Re-DRT"
            LABEL_y_b = "Re from Im-DRT"
        if AntiKK_enable == 1:
            LABEL_y_a = "Im from Im-DRT"
            LABEL_y_b = "Re from Re-DRT"
        if CombiKK_enable == 1:
            LABEL_y_a = "Im from Comb-DRT"
            LABEL_y_b = "Re from Comb-DRT"  
        b2.plot(np.log10(freq), -Z_Back.imag, linestyle='-',marker='', linewidth = 0.5, color='blue', label = LABEL_y_a)
        b2.plot(np.log10(freq),  Z_Back.real, linestyle='-',marker='', linewidth = 0.5, color='red',  label = LABEL_y_b)
        b2.plot(np.log10(freq), -Z_meas.imag, linestyle='',marker='.', markersize = 3,  color='blue', label = 'Im-meas')
        b2.plot(np.log10(freq),  Z_meas.real, linestyle='',marker='.', markersize = 3,  color='red',  label = 'Re-meas')
        ax2 = b2.twinx()
        ax2.plot(np.log10(freq), 100*np.abs((Z_Back.imag-Z_meas.imag)/Z_meas.imag), color = 'cyan',   linewidth = 0.5, label = '$\Delta$Im')
        ax2.plot(np.log10(freq), 100*np.abs((Z_Back.real-Z_meas.real)/Z_meas.real), color = 'orange', linewidth = 0.5, label = '$\Delta$Re')
        ax2.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        b2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=2, fontsize = 15)
        b2.set_xlabel('log$_{10}\,(\\nu\,/\,\mathrm{s}^{-1})$', fontsize=15)
        b2.set_ylabel('$-\mathfrak{IM}(Z(j\omega)/\Omega)$ and $\mathfrak{RE}(Z(j\omega)/\Omega)$', fontsize = 15)
        multicolor_ylabel(ax2,('$\Delta Z_{\mathrm{Re}}(j\omega))$','and','$\Delta Z_{\mathrm{Im}}(j\omega))$', 'in %'),('black', 'cyan','black','orange'), axis='y',size=15, weight='bold')
        for axis in ['top','bottom','left','right']:
            b2.spines[axis].set_linewidth(2)
            b2.spines[axis].set_color('k')
        b2.tick_params(tickdir='in', width = 2, length=6, labelsize=15)       
        #========================================================= 
        canvas = FigureCanvasTkAgg(f, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        f.tight_layout()
    
    if As_txt_saver == 1:
        Gaussian_DRT_as_Txt_Saver()

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Save the result from the EIS analysis via Gaussian DRT
#=======================================================================================================================    
#=======================================================================================================================
#======================================================================================================================= 

def Gaussian_DRT_as_Txt_Saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))
    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")    
        fi.write("Gaussian-DRT calc. from")
        fi.write("\n") 
        fi.write("Imaginary-Part:")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("yes")
        if Imag_User != 1:
            fi.write("no")
        fi.write("\n") 
        fi.write("Real-Part:")
        fi.write("\t")
        if Real_User == 1:
            fi.write("yes")
        if Real_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Combined Real_Im:")
        fi.write("\t")
        if Comb_User == 1:
            fi.write("yes")
        if Comb_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Real-Offset-removal")
        fi.write("\t")
        if RemoveRealOffset == 0:
            fi.write("no")
        if RemoveRealOffset == 1:
            fi.write("yes")
        fi.write("\n")
        fi.write("Resolution increase")
        fi.write("\t")
        fi.write(str(Resol))
        fi.write("\n")
        fi.write("Lower Range Extension")
        fi.write("\t")
        fi.write(str(ext_rangelow))
        fi.write("\n")
        fi.write("Upper Range Extension")
        fi.write("\t")
        fi.write(str(ext_rangehigh))
        fi.write("\n")
        fi.write("Reg_Par")
        fi.write("\t")
        fi.write(str(Tikh_Par))
        fi.write("\n")
        fi.write("FWHM")
        fi.write("\t")
        fi.write(str(Gaussain_FWHM))
        fi.write("\n")
        fi.write("Regularization type")
        fi.write("\t")
        if PenaltyType == "linear":
            fi.write("linear regularization")
            fi.write("\n")
        if PenaltyType == "1st-order":
            fi.write("1st-derivative regularization")
            fi.write("\n")
        if PenaltyType == "2nd-order":
            fi.write("2nd-derivative regularization")
            fi.write("\n")
        fi.write("Back-calculation")
        fi.write("\t")
        if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
            if KK_enable == 1:
                fi.write("Kramers-Kronig")
                fi.write("\n")
            if AntiKK_enable == 1:
                fi.write("Anti-Kramers-Kronig")
                fi.write("\n")
            if CombiKK_enable == 1:
                fi.write("Combi-Kramers-Kronig")
                fi.write("\n")    
        if KK_enable == 0 and AntiKK_enable == 0 and CombiKK_enable == 0:
            fi.write("None")
            fi.write("\n") 
        fi.write("\n")
        fi.write("log10(tau / s)")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("gamma(tau)_Imag")
            fi.write("\t")
        if Real_User == 1:
            fi.write("gamma(tau)_Real")
            fi.write("\t")
        if Comb_User == 1:
            fi.write("gamma(tau)_Combined")
        fi.write("\n") 
        for i in range(len(TAUARRAY)):
            fi.write(str(np.asscalar(np.log10(TAUARRAY[i]))))
            fi.write("\t")
            if Imag_User == 1:
                fi.write(str(np.asscalar(gamma_imag_gauss[i])))
                fi.write("\t")
            if Real_User == 1:
                fi.write(str(np.asscalar(gamma_real_gauss[i])))
                fi.write("\t")
            if Comb_User == 1:
                fi.write(str(np.asscalar(gamma_comb_gauss[i])))
                fi.write("\t")
            fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency [Hz]")
        fi.write("\t")
        fi.write("Z_real_used [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_used [Ohm]")
        fi.write("\t")
        if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
            if KK_enable == 1:
                fi.write("Z_real_from_ImDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_ReDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")
            if AntiKK_enable == 1:
                fi.write("Z_real_from_ReDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_ImDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")
            if CombiKK_enable == 1:
                fi.write("Z_real_from_CombDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_CombDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")       
        fi.write("\n")        
        for i in range(len(freq)):
            fi.write(str(np.asscalar(freq[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_meas[i].real)))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_meas[i].imag)))
            fi.write("\t")
            if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
                fi.write(str(np.asscalar(Z_Back[i].real)))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_Back[i].imag)))
                fi.write("\t")
                fi.write(str(100*np.asscalar(np.abs((Z_Back[i].real-Z_meas[i].real)/Z_meas[i].real))))
                fi.write("\t")
                fi.write(str(100*np.asscalar(np.abs((Z_Back[i].imag-Z_meas[i].imag)/Z_meas[i].imag))))
                fi.write("\t")
            fi.write("\n") 
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency_original [Hz]")
        fi.write("\t")
        fi.write("Z_real_original [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_originnal [Ohm]")
        fi.write("\n")
        for i in range(len(Frequency_Array)):
            fi.write(str(np.asscalar(Frequency_Array[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_real_Meas[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_imag_Meas[i])))
            fi.write("\n")
    root.destroy()
    
    
#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Define the inputs for EIS analysis via ColeCole DRT
#=======================================================================================================================    
#=======================================================================================================================
#======================================================================================================================= 
    
    
def ColeCole_NNLS_DRT():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    if File_Was_Loaded ==1:
        ColeCole_DRT_Entry_Window()



def ColeCole_DRT_Entry_Window():
    Fenster = Toplevel()  
    Fenster.geometry("550x300")
    colorbgr = Label(Fenster, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000) 
    Fenster.resizable(False, False)
    Fenster.title("Cole-Cole NNLS DRT")
     
    Resol_Label = Label(Fenster, text="Increase Data res. by fact.*", bg = '#80ffbf')
    Resol_Label.place(x = 25, y = 10)
    Resol_Eingabe = Entry(Fenster)      
    Resol_Eingabe.insert(END, 10)                                        
    Resol_Eingabe.place(x = 185, y = 10, width = 80, height = 19)
    
    extend_low_Label = Label(Fenster,text="Extend low*", bg = '#80ffbf')
    extend_low_Label.place(x = 25, y = 35)
    extend_low_Eingabe = Entry(Fenster)
    extend_low_Eingabe.insert(END, 0) 
    extend_low_Eingabe.place(x = 185, y = 35, width = 80, height = 19)
        
    extend_high_Label = Label(Fenster,text="Extend high*", bg = '#80ffbf')
    extend_high_Label.place(x = 25, y = 60)
    extend_high_Eingabe = Entry(Fenster)
    extend_high_Eingabe.insert(END, 0) 
    extend_high_Eingabe.place(x = 185, y = 60, width = 80, height = 19)
    
    reg_par_Label = Label(Fenster, text="Tikhonov parameter*", bg = '#80ffbf')
    reg_par_Label.place(x = 300, y = 10)
    reg_par_Eingabe = Entry(Fenster)   
    reg_par_Eingabe.insert(END, 0.01)                                             
    reg_par_Eingabe.place(x = 440, y = 10, width = 80, height = 19)
    
    ColeColeGamma_Label = Label(Fenster, text="Cole-Cole Parameter*", bg = '#80ffbf')
    ColeColeGamma_Label.place(x = 300, y = 35)
    ColeColeGamma_Eingabe = Entry(Fenster)   
    ColeColeGamma_Eingabe.insert(END, 0.95)                                             
    ColeColeGamma_Eingabe.place(x = 440, y = 35, width = 80, height = 19)
    
    SepLine0 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine0.place(x = 0, y = 85, height = 10)
    
    
    Lab1 = Label(Fenster,text="DRT from:", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab1.place(x = 25, y = 100 )

    var1 = IntVar()
    Checkbutton(Fenster, text="Imaginary part", variable=var1, bg = '#80ffbf').place(x = 25, y = 180)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Real part", variable=var2, bg = '#80ffbf').place(x = 25, y = 155)
    
    var3 = IntVar(value = 1)
    Checkbutton(Fenster, text="Combined fit", variable=var3, bg = '#80ffbf').place(x = 25, y = 130)
    
   
    var0 = StringVar()
    var0.set('linear')  # initializing the choice, i.e. Python
    
    def ChangePenalty():
        var0.get()

    Lab2 = Label(Fenster,text="Penalty-Type", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab2.place(x = 182, y = 100 )
    
    
    Radiobutton(Fenster, text='linear',    padx = 20, variable=var0, value='linear',    command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 130, height = 25)
    Radiobutton(Fenster, text='1st-order', padx = 20, variable=var0, value='1st-order', command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 155, height = 25)
    Radiobutton(Fenster, text='2nd-order', padx = 20, variable=var0, value='2nd-order', command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 180, height = 25)
    
    
    var5 = StringVar()
    var5.set(None)  # initializing the choice, i.e. Python
    
    def GetBackCalckType():
        var5.get()

    Lab2 = Label(Fenster,text="Calc. back to Z", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab2.place(x = 330, y = 100 )
    
    Radiobutton(Fenster, text='Kramers-Kronig like', padx = 20, variable=var5, value='KramersKronig',        command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 130, height = 25)
    Radiobutton(Fenster, text='Anti Kramers-Kronig', padx = 20, variable=var5, value='AntiKramersKronig',    command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 155, height = 25)
    Radiobutton(Fenster, text='Combi Kramers-Kronig', padx = 20, variable=var5, value='CombinedBack',         command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 180, height = 25)
    
    SepLine1 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine1.place(x = 0, y = 205, height = 10)
    
    var6 = IntVar()
    Checkbutton(Fenster, text="Remove Real-offset", variable=var6, bg = '#80ffbf').place(x = 25, y = 230)
    
    var4 = IntVar()
    Checkbutton(Fenster, text="Save Data", variable=var4, bg = '#80ffbf').place(x = 25, y = 255)    
    
    def AcceptParams():
        global Resol
        global Tikh_Par
        global ext_rangelow
        global ext_rangehigh
        global Imag_User
        global Real_User
        global Comb_User
        global As_txt_saver
        global ColeColeGamma
        global PenaltyType
        global RemoveRealOffset
        global BackCalculation
        global KK_enable
        global AntiKK_enable
        global CombiKK_enable
        Resol             = int(float(Resol_Eingabe.get()))
        Tikh_Par          = (float(reg_par_Eingabe.get()))
        ext_rangelow      = (float(extend_low_Eingabe.get()))
        ext_rangehigh     = (float(extend_high_Eingabe.get()))
        Imag_User         = var1.get()
        Real_User         = var2.get()
        Comb_User         = var3.get()
        As_txt_saver      = var4.get()   
        ColeColeGamma     = (float(ColeColeGamma_Eingabe.get()))
        PenaltyType       = var0.get()
        RemoveRealOffset  = var6.get()
        BackCalculation   = var5.get()
        KK_enable         = 0
        AntiKK_enable     = 0
        CombiKK_enable    = 0
        if BackCalculation == 'KramersKronig' and Imag_User == 0 or BackCalculation == 'KramersKronig' and Real_User == 0 :
            No_KK_check_possible()
        if BackCalculation == 'KramersKronig' and Imag_User == 1 and Real_User == 1:
            KK_enable         = 1
        if BackCalculation == 'AntiKramersKronig' and Imag_User == 0 or BackCalculation == 'AntiKramersKronig' and Real_User == 0:
            No_AntiKK_check_possible()
        if BackCalculation == 'AntiKramersKronig' and Imag_User == 1 and Real_User == 1:
            AntiKK_enable     = 1
        if BackCalculation == 'CombinedBack' and Comb_User == 0:
            No_CombiKK_check_possible()
        if BackCalculation == 'CombinedBack' and Comb_User == 1:
            CombiKK_enable    = 1
              
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).place( x = 185, y = 230, width = 130, height = 50)
    button=Button(Fenster,text="Next",command=ColeCole_DRT_Transformation).place(x = 335, y = 230, width = 130, height = 50)
    
#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Run the EIS analysis via ColeCole DRT
#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================     


def ColeCole_DRT_Transformation():
    root = Toplevel()
    root.title("Result Cole-Cole DRT-Transformation")
    global TAUARRAY
    global gamma_imag_ColeCole
    global gamma_real_ColeCole
    global gamma_comb_ColeCole
    global Z_imag_Meas
    global Z_real_Meas
    global freq
    global Z_Back 
    global Z_meas
    Z_imag_Meas     = np.squeeze(Z_IMAG)               #Nach dem Laden heißen die Arrays nur so wie früher
    Z_real_Meas     = np.squeeze(Z_REAL)           #deshalb hier die Umbenennung
    freq            = np.squeeze(Frequency_Array) 
    Z               = Z_real_Meas + 1j*Z_imag_Meas
    freq            = freq[np.imag(Z) < 0]
    Z               = Z[np.imag(Z) < 0]
    Z_meas          = Z.copy()
    Real_Offset = np.min(Z.real).copy()
    if RemoveRealOffset == 1:
        Z = Z - np.min(Z.real)  
    if Imag_User == 1:
        IM_Output           = DRT_ColeCole_RBF(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "imag", GAMMA_Cole = ColeColeGamma, Reg_Type = PenaltyType)
        TAUARRAY            = 10**IM_Output[0]
        gamma_imag_ColeCole = IM_Output[2]
    if Real_User == 1:
        RE_Output           = DRT_ColeCole_RBF(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "real", GAMMA_Cole = ColeColeGamma, Reg_Type = PenaltyType)
        TAUARRAY            = 10**RE_Output[0]
        gamma_real_ColeCole = RE_Output[2]
    if Comb_User == 1:
        CO_Output           = DRT_ColeCole_RBF(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "comb", GAMMA_Cole = ColeColeGamma, Reg_Type = PenaltyType)
        TAUARRAY            = 10**CO_Output[0]
        gamma_comb_ColeCole = CO_Output[2]    
    ##########################################################################################
    # Back-calculation from DRTs to impedances in any permutation, KK, anti KK, combi
    ########################################################################################## 
    if KK_enable == 1:
        Z_Back           = KramersKronig_Multiplication(MATRIX = IM_Output[5], Real_DRT = RE_Output[2], Imag_DRT = IM_Output[2])
        Z_Back = Z_Back + Real_Offset     
    if AntiKK_enable == 1:
        Z_Back           = AntiKramersKronig_Multiplication(MATRIX = IM_Output[5], Real_DRT = RE_Output[2], Imag_DRT = IM_Output[2])
        if RemoveRealOffset == 1:
            Z_Back = Z_Back + Real_Offset
    if CombiKK_enable == 1:
        Z_Back           = CombiKramersKronig_Multiplication(MATRIX = CO_Output[5], Comb_DRT = CO_Output[2])
        if RemoveRealOffset == 1:
            Z_Back = Z_Back + Real_Offset
    ##########################################################################################
    #  Plot the DRT
    ##########################################################################################
    if BackCalculation == "None" or BackCalculation == 'KramersKronig' and KK_enable == 0 or BackCalculation == 'AntiKramersKronig' and AntiKK_enable == 0 or BackCalculation == 'CombinedBack' and CombiKK_enable == 0:
        Abbildung = Figure(figsize=(7, 5), dpi=100)
        b = Abbildung.add_subplot(111)
        if Imag_User == 1:
            b.plot(IM_Output[0],IM_Output[2],color='blue', linestyle ='-' ,label = 'Imag-DRT')
        if Real_User == 1:
            b.plot(RE_Output[0],RE_Output[2],color='red', linestyle ='-'  ,label = 'Real-DRT')
        if Comb_User == 1:
            b.plot(CO_Output[0],CO_Output[2],color='black', linestyle ='-',label = 'Combi-DRT')    
        b.legend(frameon = False, fontsize = 15)    
        b.set_xlabel('log$_{10}\,(\\tau\,/\,\mathrm{s})$', fontsize=15)
        b.set_ylabel('$\Delta Z_{\mathrm{Real}}\,\\tau\,\gamma(\\tau)\,/\,\Omega$', fontsize=15)
        for axis in ['top','bottom','left','right']:
            b.spines[axis].set_linewidth(2)
            b.spines[axis].set_color('k')
        b.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        secax = b.secondary_xaxis('top', functions=(log_to_ln, ln_to_log))
        secax.set_xlabel('ln$\,(\\tau\,/\,\mathrm{s})$', fontsize = 15)
        secax.tick_params(tickdir='in', width = 2, length=6, labelsize=15)  
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    ##########################################################################################
    #  Plot the back-calculation result (in case it was activated)
    ##########################################################################################
    if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
        #============================================================================
        #Plotting of DRT and back-computation
        #============================================================================ 
        f  = Figure(figsize=(11, 6), dpi=100)        
        b1 = f.add_subplot(121)
        b2 = f.add_subplot(122) 
        #=========================================================
        #final plotting
        #=========================================================
        if Imag_User == 1:
            b1.plot(IM_Output[0],IM_Output[2],color='blue', linestyle ='-' ,label = 'Imag-DRT')
        if Real_User == 1:
            b1.plot(RE_Output[0],RE_Output[2],color='red', linestyle ='-'  ,label = 'Real-DRT')
        if Comb_User == 1:
            b1.plot(CO_Output[0],CO_Output[2],color='black', linestyle ='-',label = 'Combi-DRT')
        b1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=2, fontsize = 15)   
        b1.set_xlabel('log$_{10}\,(\\tau\,/\,\mathrm{s})$', fontsize=15)
        b1.set_ylabel('$\Delta Z_{\mathrm{Real}}\,\\tau\,\gamma(\\tau)\,/\,\Omega$', fontsize=15)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        secax = b1.secondary_xaxis('top', functions=(log_to_ln, ln_to_log))
        secax.set_xlabel('ln$\,(\\tau\,/\,\mathrm{s})$', fontsize = 15)
        secax.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        if KK_enable == 1:
            LABEL_y_a = "Im from Re-DRT"
            LABEL_y_b = "Re from Im-DRT"
        if AntiKK_enable == 1:
            LABEL_y_a = "Im from Im-DRT"
            LABEL_y_b = "Re from Re-DRT"
        if CombiKK_enable == 1:
            LABEL_y_a = "Im from Comb-DRT"
            LABEL_y_b = "Re from Comb-DRT"  
        b2.plot(np.log10(freq), -Z_Back.imag, linestyle='-',marker='', linewidth = 0.5, color='blue', label = LABEL_y_a)
        b2.plot(np.log10(freq),  Z_Back.real, linestyle='-',marker='', linewidth = 0.5, color='red',  label = LABEL_y_b)
        b2.plot(np.log10(freq), -Z_meas.imag, linestyle='',marker='.', markersize = 3,  color='blue', label = 'Im-meas')
        b2.plot(np.log10(freq),  Z_meas.real, linestyle='',marker='.', markersize = 3,  color='red',  label = 'Re-meas')
        ax2 = b2.twinx()
        ax2.plot(np.log10(freq), 100*np.abs((Z_Back.imag-Z_meas.imag)/Z_meas.imag), color = 'cyan',   linewidth = 0.5, label = '$\Delta$Im')
        ax2.plot(np.log10(freq), 100*np.abs((Z_Back.real-Z_meas.real)/Z_meas.real), color = 'orange', linewidth = 0.5, label = '$\Delta$Re')
        ax2.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        b2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=2, fontsize = 15)
        b2.set_xlabel('log$_{10}\,(\\nu\,/\,\mathrm{s}^{-1})$', fontsize=15)
        b2.set_ylabel('$-\mathfrak{IM}(Z(j\omega)/\Omega)$ and $\mathfrak{RE}(Z(j\omega)/\Omega)$', fontsize = 15)
        multicolor_ylabel(ax2,('$\Delta Z_{\mathrm{Re}}(j\omega))$','and','$\Delta Z_{\mathrm{Im}}(j\omega))$', 'in %'),('black', 'cyan','black','orange'), axis='y',size=15, weight='bold')
        for axis in ['top','bottom','left','right']:
            b2.spines[axis].set_linewidth(2)
            b2.spines[axis].set_color('k')
        b2.tick_params(tickdir='in', width = 2, length=6, labelsize=15)       
        #========================================================= 
        canvas = FigureCanvasTkAgg(f, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        f.tight_layout()
    if As_txt_saver == 1:
        ColeCole_DRT_as_Txt_Saver()

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Save the result from the EIS analysis via ColeCole DRT
#=======================================================================================================================    
#=======================================================================================================================
#======================================================================================================================= 

def ColeCole_DRT_as_Txt_Saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))
    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")   
        fi.write("Cole-Cole-DRT calc. from")
        fi.write("\n") 
        fi.write("Imaginary-Part:")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("yes")
        if Imag_User != 1:
            fi.write("no")
        fi.write("\n") 
        fi.write("Real-Part:")
        fi.write("\t")
        if Real_User == 1:
            fi.write("yes")
        if Real_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Combined Real_Im:")
        fi.write("\t")
        if Comb_User == 1:
            fi.write("yes")
        if Comb_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Real-Offset-removal")
        fi.write("\t")
        if RemoveRealOffset == 0:
            fi.write("no")
        if RemoveRealOffset == 1:
            fi.write("yes")
        fi.write("\n")
        fi.write("Resolution increase")
        fi.write("\t")
        fi.write(str(Resol))
        fi.write("\n")
        fi.write("Lower Range Extension")
        fi.write("\t")
        fi.write(str(ext_rangelow))
        fi.write("\n")
        fi.write("Upper Range Extension")
        fi.write("\t")
        fi.write(str(ext_rangehigh))
        fi.write("\n")
        fi.write("Reg_Par")
        fi.write("\t")
        fi.write(str(Tikh_Par))
        fi.write("\n")
        fi.write("Cole-Cole-Gamma")
        fi.write("\t")
        fi.write(str(ColeColeGamma))
        fi.write("\n")
        fi.write("Regularization type")
        fi.write("\t")
        if PenaltyType == "linear":
            fi.write("linear regularization")
            fi.write("\n")
        if PenaltyType == "1st-order":
            fi.write("1st-derivative regularization")
            fi.write("\n")
        if PenaltyType == "2nd-order":
            fi.write("2nd-derivative regularization")
            fi.write("\n")   
        fi.write("Back-calculation")
        fi.write("\t")
        if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
            if KK_enable == 1:
                fi.write("Kramers-Kronig")
                fi.write("\n")
            if AntiKK_enable == 1:
                fi.write("Anti-Kramers-Kronig")
                fi.write("\n")
            if CombiKK_enable == 1:
                fi.write("Combi-Kramers-Kronig")
                fi.write("\n")    
        if KK_enable == 0 and AntiKK_enable == 0 and CombiKK_enable == 0:
            fi.write("None")
            fi.write("\n") 
        fi.write("\n")
        fi.write("log10(tau / s)")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("gamma(tau)_Imag")
            fi.write("\t")
        if Real_User == 1:
            fi.write("gamma(tau)_Real")
            fi.write("\t")
        if Comb_User == 1:
            fi.write("gamma(tau)_Combined")
        fi.write("\n")      
        for i in range(len(TAUARRAY)):
            fi.write(str(np.asscalar(np.log10(TAUARRAY[i]))))
            fi.write("\t")
            if Imag_User == 1:
                fi.write(str(np.asscalar(gamma_imag_ColeCole[i])))
                fi.write("\t")
            if Real_User == 1:
                fi.write(str(np.asscalar(gamma_real_ColeCole[i])))
                fi.write("\t")
            if Comb_User == 1:
                fi.write(str(np.asscalar(gamma_comb_ColeCole[i])))
                fi.write("\t")
            fi.write("\n") 
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency [Hz]")
        fi.write("\t")
        fi.write("Z_real_used [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_used [Ohm]")
        fi.write("\t")
        if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
            if KK_enable == 1:
                fi.write("Z_real_from_ImDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_ReDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")
            if AntiKK_enable == 1:
                fi.write("Z_real_from_ReDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_ImDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")
            if CombiKK_enable == 1:
                fi.write("Z_real_from_CombDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_CombDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")       
        fi.write("\n")        
        for i in range(len(freq)):
            fi.write(str(np.asscalar(freq[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_meas[i].real)))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_meas[i].imag)))
            fi.write("\t")
            if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
                fi.write(str(np.asscalar(Z_Back[i].real)))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_Back[i].imag)))
                fi.write("\t")
                fi.write(str(100*np.asscalar(np.abs((Z_Back[i].real-Z_meas[i].real)/Z_meas[i].real))))
                fi.write("\t")
                fi.write(str(100*np.asscalar(np.abs((Z_Back[i].imag-Z_meas[i].imag)/Z_meas[i].imag))))
                fi.write("\t")
            fi.write("\n") 
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency_original [Hz]")
        fi.write("\t")
        fi.write("Z_real_original [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_originnal [Ohm]")
        fi.write("\n")
        for i in range(len(Frequency_Array)):
            fi.write(str(np.asscalar(Frequency_Array[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_real_Meas[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_imag_Meas[i])))
            fi.write("\n")    
    root.destroy()   
    
    
    
    
#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Define the inputs for EIS analysis via Havriliak-Negami DRT
#=======================================================================================================================    
#=======================================================================================================================
#======================================================================================================================= 
    
    
def HavriliaNegami_NNLS_DRT():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    if File_Was_Loaded ==1:
        HavriliaNegami_DRT_Entry_Window()



def HavriliaNegami_DRT_Entry_Window():
    Fenster = Toplevel()  
    Fenster.geometry("550x300")
    colorbgr = Label(Fenster, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000) 
    Fenster.resizable(False, False)
    Fenster.title("Havriliak-Negami NNLS DRT")
     
    Resol_Label = Label(Fenster, text="Increase Data res. by fact.*", bg = '#80ffbf')
    Resol_Label.place(x = 25, y = 10)
    Resol_Eingabe = Entry(Fenster)      
    Resol_Eingabe.insert(END, 10)                                        
    Resol_Eingabe.place(x = 185, y = 10, width = 80, height = 19)
    
    extend_low_Label = Label(Fenster,text="Extend low*", bg = '#80ffbf')
    extend_low_Label.place(x = 25, y = 35)
    extend_low_Eingabe = Entry(Fenster)
    extend_low_Eingabe.insert(END, 0) 
    extend_low_Eingabe.place(x = 185, y = 35, width = 80, height = 19)
        
    extend_high_Label = Label(Fenster,text="Extend high*", bg = '#80ffbf')
    extend_high_Label.place(x = 25, y = 60)
    extend_high_Eingabe = Entry(Fenster)
    extend_high_Eingabe.insert(END, 0) 
    extend_high_Eingabe.place(x = 185, y = 60, width = 80, height = 19)
    
    reg_par_Label = Label(Fenster, text="Tikhonov parameter*", bg = '#80ffbf')
    reg_par_Label.place(x = 300, y = 10)
    reg_par_Eingabe = Entry(Fenster)   
    reg_par_Eingabe.insert(END, 0.01)                                             
    reg_par_Eingabe.place(x = 440, y = 10, width = 80, height = 19)
    
    HavNegAlpha_Label = Label(Fenster, text="Havriliak-Negami alpha*", bg = '#80ffbf')
    HavNegAlpha_Label.place(x = 300, y = 35)
    HavNegAlpha_Eingabe = Entry(Fenster)   
    HavNegAlpha_Eingabe.insert(END, 0.95)                                             
    HavNegAlpha_Eingabe.place(x = 440, y = 35, width = 80, height = 19)
    
    HavNegBeta_Label = Label(Fenster, text="Havriliak-Negami beta*", bg = '#80ffbf')
    HavNegBeta_Label.place(x = 300, y = 60)
    HavNegBeta_Eingabe = Entry(Fenster)   
    HavNegBeta_Eingabe.insert(END, 0.95)                                             
    HavNegBeta_Eingabe.place(x = 440, y = 60, width = 80, height = 19)
    
    SepLine0 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine0.place(x = 0, y = 85, height = 10)
    
    
    Lab1 = Label(Fenster,text="DRT from:", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab1.place(x = 25, y = 100 )

    var1 = IntVar()
    Checkbutton(Fenster, text="Imaginary part", variable=var1, bg = '#80ffbf').place(x = 25, y = 180)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Real part", variable=var2, bg = '#80ffbf').place(x = 25, y = 155)
    
    var3 = IntVar(value = 1)
    Checkbutton(Fenster, text="Combined fit", variable=var3, bg = '#80ffbf').place(x = 25, y = 130)
    
   
    var0 = StringVar()
    var0.set('linear')  # initializing the choice, i.e. Python
    
    def ChangePenalty():
        var0.get()

    Lab2 = Label(Fenster,text="Penalty-Type", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab2.place(x = 182, y = 100 )
    
    
    Radiobutton(Fenster, text='linear',    padx = 20, variable=var0, value='linear',    command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 130, height = 25)
    Radiobutton(Fenster, text='1st-order', padx = 20, variable=var0, value='1st-order', command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 155, height = 25)
    Radiobutton(Fenster, text='2nd-order', padx = 20, variable=var0, value='2nd-order', command = ChangePenalty, bg = '#80ffbf').place(x = 165, y = 180, height = 25)
    
    
    var5 = StringVar()
    var5.set(None)  # initializing the choice, i.e. Python
    
    def GetBackCalckType():
        var5.get()

    Lab2 = Label(Fenster,text="Calc. back to Z", bg = '#80ffbf', font = ('Arial',12,'bold'))
    Lab2.place(x = 330, y = 100 )
    
    Radiobutton(Fenster, text='Kramers-Kronig like',  padx = 20, variable=var5, value='KramersKronig',       command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 130, height = 25)
    Radiobutton(Fenster, text='Anti Kramers-Kronig',  padx = 20, variable=var5, value='AntiKramersKronig',   command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 155, height = 25)
    Radiobutton(Fenster, text='Combi Kramers-Kronig', padx = 20, variable=var5, value='CombinedBack',        command = GetBackCalckType, bg = '#80ffbf').place(x = 313, y = 180, height = 25)
    
    SepLine1 = Label(Fenster,text="-------------------------------------------------------------------------------------------------------------------------------------------------------------------------",bg = '#80ffbf')
    SepLine1.place(x = 0, y = 205, height = 10)
    
    var6 = IntVar()
    Checkbutton(Fenster, text="Remove Real-offset", variable=var6, bg = '#80ffbf').place(x = 25, y = 230)
    
    var4 = IntVar()
    Checkbutton(Fenster, text="Save Data", variable=var4, bg = '#80ffbf').place(x = 25, y = 255)    
    
    def AcceptParams():
        global Resol
        global Tikh_Par
        global ext_rangelow
        global ext_rangehigh
        global Imag_User
        global Real_User
        global Comb_User
        global As_txt_saver
        global HavNegAlpha
        global HavNegBeta
        global PenaltyType
        global RemoveRealOffset
        global BackCalculation
        global KK_enable
        global AntiKK_enable
        global CombiKK_enable
        Resol             = int(float(Resol_Eingabe.get()))
        Tikh_Par          = (float(reg_par_Eingabe.get()))
        ext_rangelow      = (float(extend_low_Eingabe.get()))
        ext_rangehigh     = (float(extend_high_Eingabe.get()))
        Imag_User         = var1.get()
        Real_User         = var2.get()
        Comb_User         = var3.get()
        As_txt_saver      = var4.get()   
        HavNegAlpha       = (float(HavNegAlpha_Eingabe.get()))
        HavNegBeta        = (float(HavNegBeta_Eingabe.get()))
        PenaltyType       = var0.get()
        RemoveRealOffset  = var6.get()
        BackCalculation   = var5.get()
        KK_enable         = 0
        AntiKK_enable     = 0
        CombiKK_enable    = 0
        if BackCalculation == 'KramersKronig' and Imag_User == 0 or BackCalculation == 'KramersKronig' and Real_User == 0 :
            No_KK_check_possible()
        if BackCalculation == 'KramersKronig' and Imag_User == 1 and Real_User == 1:
            KK_enable         = 1
        if BackCalculation == 'AntiKramersKronig' and Imag_User == 0 or BackCalculation == 'AntiKramersKronig' and Real_User == 0:
            No_AntiKK_check_possible()
        if BackCalculation == 'AntiKramersKronig' and Imag_User == 1 and Real_User == 1:
            AntiKK_enable     = 1
        if BackCalculation == 'CombinedBack' and Comb_User == 0:
            No_CombiKK_check_possible()
        if BackCalculation == 'CombinedBack' and Comb_User == 1:
            CombiKK_enable    = 1
              
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).place( x = 185, y = 230, width = 130, height = 50)
    button=Button(Fenster,text="Next",command=HavriliakNegami_DRT_Transformation).place(x = 335, y = 230, width = 130, height = 50)
    
#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Run the EIS analysis via ColeCole DRT
#=======================================================================================================================    
#=======================================================================================================================
#=======================================================================================================================     


def HavriliakNegami_DRT_Transformation():
    root = Toplevel()
    root.title("Result Havriliak-Negami DRT-Transformation")
    global TAUARRAY
    global gamma_imag_HavriliakNegami
    global gamma_real_HavriliakNegami
    global gamma_comb_HavriliakNegami
    global Z_imag_Meas
    global Z_real_Meas
    global freq
    global Z_Back 
    global Z_meas
    Z_imag_Meas     = np.squeeze(Z_IMAG)               #Nach dem Laden heißen die Arrays nur so wie früher
    Z_real_Meas     = np.squeeze(Z_REAL)           #deshalb hier die Umbenennung
    freq            = np.squeeze(Frequency_Array) 
    Z               = Z_real_Meas + 1j*Z_imag_Meas
    freq            = freq[np.imag(Z) < 0]
    Z               = Z[np.imag(Z) < 0]
    Z_meas          = Z.copy()
    Real_Offset = np.min(Z.real).copy()
    if RemoveRealOffset == 1:
        Z = Z - np.min(Z.real)  
    if Imag_User == 1:
        IM_Output                  = DRT_HavriliakNegami_RBF(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "imag", ALPHA_HAVNEG = HavNegAlpha, BETA_HAVNEG = HavNegBeta, Reg_Type = PenaltyType)
        TAUARRAY                   = 10**IM_Output[0]
        gamma_imag_HavriliakNegami = IM_Output[2]
    if Real_User == 1:
        RE_Output                  = DRT_HavriliakNegami_RBF(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "real", ALPHA_HAVNEG = HavNegAlpha, BETA_HAVNEG = HavNegBeta, Reg_Type = PenaltyType)
        TAUARRAY                   = 10**RE_Output[0]
        gamma_real_HavriliakNegami = RE_Output[2]
    if Comb_User == 1:
        CO_Output                  = DRT_HavriliakNegami_RBF(freq_raw = freq, Z_raw = Z, damp = Tikh_Par, low_ext = ext_rangelow, high_ext = ext_rangehigh, resol = Resol, DRT_Type = "comb", ALPHA_HAVNEG = HavNegAlpha, BETA_HAVNEG = HavNegBeta, Reg_Type = PenaltyType)
        TAUARRAY                   = 10**CO_Output[0]
        gamma_comb_HavriliakNegami = CO_Output[2]    
    ##########################################################################################
    # Back-calculation from DRTs to impedances in any permutation, KK, anti KK, combi
    ########################################################################################## 
    if KK_enable == 1:
        Z_Back           = KramersKronig_Multiplication(MATRIX = IM_Output[5], Real_DRT = RE_Output[2], Imag_DRT = IM_Output[2])
        Z_Back = Z_Back + Real_Offset     
    if AntiKK_enable == 1:
        Z_Back           = AntiKramersKronig_Multiplication(MATRIX = IM_Output[5], Real_DRT = RE_Output[2], Imag_DRT = IM_Output[2])
        if RemoveRealOffset == 1:
            Z_Back = Z_Back + Real_Offset
    if CombiKK_enable == 1:
        Z_Back           = CombiKramersKronig_Multiplication(MATRIX = CO_Output[5], Comb_DRT = CO_Output[2])
        if RemoveRealOffset == 1:
            Z_Back = Z_Back + Real_Offset
    ##########################################################################################
    #  Plot the DRT
    ##########################################################################################
    if BackCalculation == "None" or BackCalculation == 'KramersKronig' and KK_enable == 0 or BackCalculation == 'AntiKramersKronig' and AntiKK_enable == 0 or BackCalculation == 'CombinedBack' and CombiKK_enable == 0:
        Abbildung = Figure(figsize=(7, 5), dpi=100)
        b = Abbildung.add_subplot(111)
        if Imag_User == 1:
            b.plot(IM_Output[0],IM_Output[2],color='blue', linestyle ='-' ,label = 'Imag-DRT')
        if Real_User == 1:
            b.plot(RE_Output[0],RE_Output[2],color='red', linestyle ='-'  ,label = 'Real-DRT')
        if Comb_User == 1:
            b.plot(CO_Output[0],CO_Output[2],color='black', linestyle ='-',label = 'Combi-DRT')    
        b.legend(frameon = False, fontsize = 15)    
        b.set_xlabel('log$_{10}\,(\\tau\,/\,\mathrm{s})$', fontsize=15)
        b.set_ylabel('$\Delta Z_{\mathrm{Real}}\,\\tau\,\gamma(\\tau)\,/\,\Omega$', fontsize=15)
        for axis in ['top','bottom','left','right']:
            b.spines[axis].set_linewidth(2)
            b.spines[axis].set_color('k')
        b.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        secax = b.secondary_xaxis('top', functions=(log_to_ln, ln_to_log))
        secax.set_xlabel('ln$\,(\\tau\,/\,\mathrm{s})$', fontsize = 15)
        secax.tick_params(tickdir='in', width = 2, length=6, labelsize=15)  
        canvas = FigureCanvasTkAgg(Abbildung, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    ##########################################################################################
    #  Plot the back-calculation result (in case it was activated)
    ##########################################################################################
    if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
        #============================================================================
        #Plotting of DRT and back-computation
        #============================================================================ 
        f  = Figure(figsize=(11, 6), dpi=100)        
        b1 = f.add_subplot(121)
        b2 = f.add_subplot(122) 
        #=========================================================
        #final plotting
        #=========================================================
        if Imag_User == 1:
            b1.plot(IM_Output[0],IM_Output[2],color='blue', linestyle ='-' ,label = 'Imag-DRT')
        if Real_User == 1:
            b1.plot(RE_Output[0],RE_Output[2],color='red', linestyle ='-'  ,label = 'Real-DRT')
        if Comb_User == 1:
            b1.plot(CO_Output[0],CO_Output[2],color='black', linestyle ='-',label = 'Combi-DRT')
        b1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=2, fontsize = 15)   
        b1.set_xlabel('log$_{10}\,(\\tau\,/\,\mathrm{s})$', fontsize=15)
        b1.set_ylabel('$\Delta Z_{\mathrm{Real}}\,\\tau\,\gamma(\\tau)\,/\,\Omega$', fontsize=15)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        secax = b1.secondary_xaxis('top', functions=(log_to_ln, ln_to_log))
        secax.set_xlabel('ln$\,(\\tau\,/\,\mathrm{s})$', fontsize = 15)
        secax.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        if KK_enable == 1:
            LABEL_y_a = "Im from Re-DRT"
            LABEL_y_b = "Re from Im-DRT"
        if AntiKK_enable == 1:
            LABEL_y_a = "Im from Im-DRT"
            LABEL_y_b = "Re from Re-DRT"
        if CombiKK_enable == 1:
            LABEL_y_a = "Im from Comb-DRT"
            LABEL_y_b = "Re from Comb-DRT"  
        b2.plot(np.log10(freq), -Z_Back.imag, linestyle='-',marker='', linewidth = 0.5, color='blue', label = LABEL_y_a)
        b2.plot(np.log10(freq),  Z_Back.real, linestyle='-',marker='', linewidth = 0.5, color='red',  label = LABEL_y_b)
        b2.plot(np.log10(freq), -Z_meas.imag, linestyle='',marker='.', markersize = 3,  color='blue', label = 'Im-meas')
        b2.plot(np.log10(freq),  Z_meas.real, linestyle='',marker='.', markersize = 3,  color='red',  label = 'Re-meas')
        ax2 = b2.twinx()
        ax2.plot(np.log10(freq), 100*np.abs((Z_Back.imag-Z_meas.imag)/Z_meas.imag), color = 'cyan',   linewidth = 0.5, label = '$\Delta$Im')
        ax2.plot(np.log10(freq), 100*np.abs((Z_Back.real-Z_meas.real)/Z_meas.real), color = 'orange', linewidth = 0.5, label = '$\Delta$Re')
        ax2.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
        b2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=2, fontsize = 15)
        b2.set_xlabel('log$_{10}\,(\\nu\,/\,\mathrm{s}^{-1})$', fontsize=15)
        b2.set_ylabel('$-\mathfrak{IM}(Z(j\omega)/\Omega)$ and $\mathfrak{RE}(Z(j\omega)/\Omega)$', fontsize = 15)
        multicolor_ylabel(ax2,('$\Delta Z_{\mathrm{Re}}(j\omega))$','and','$\Delta Z_{\mathrm{Im}}(j\omega))$', 'in %'),('black', 'cyan','black','orange'), axis='y',size=15, weight='bold')
        for axis in ['top','bottom','left','right']:
            b2.spines[axis].set_linewidth(2)
            b2.spines[axis].set_color('k')
        b2.tick_params(tickdir='in', width = 2, length=6, labelsize=15)       
        #========================================================= 
        canvas = FigureCanvasTkAgg(f, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        f.tight_layout()
    if As_txt_saver == 1:
        HavriliakNegami_DRT_as_Txt_Saver()

#=======================================================================================================================  
#======================================================================================================================= 
#=======================================================================================================================  
#  Save the result from the EIS analysis via ColeCole DRT
#=======================================================================================================================    
#=======================================================================================================================
#======================================================================================================================= 

def HavriliakNegami_DRT_as_Txt_Saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))
    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")   
        fi.write("Havriliak-Negami-DRT calc. from")
        fi.write("\n") 
        fi.write("Imaginary-Part:")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("yes")
        if Imag_User != 1:
            fi.write("no")
        fi.write("\n") 
        fi.write("Real-Part:")
        fi.write("\t")
        if Real_User == 1:
            fi.write("yes")
        if Real_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Combined Real_Im:")
        fi.write("\t")
        if Comb_User == 1:
            fi.write("yes")
        if Comb_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Real-Offset-removal")
        fi.write("\t")
        if RemoveRealOffset == 0:
            fi.write("no")
        if RemoveRealOffset == 1:
            fi.write("yes")
        fi.write("\n")
        fi.write("Resolution increase")
        fi.write("\t")
        fi.write(str(Resol))
        fi.write("\n")
        fi.write("Lower Range Extension")
        fi.write("\t")
        fi.write(str(ext_rangelow))
        fi.write("\n")
        fi.write("Upper Range Extension")
        fi.write("\t")
        fi.write(str(ext_rangehigh))
        fi.write("\n")
        fi.write("Reg_Par")
        fi.write("\t")
        fi.write(str(Tikh_Par))
        fi.write("\n")
        fi.write("Havriliak-Negami-alpha")
        fi.write("\t")
        fi.write(str(HavNegAlpha))
        fi.write("\n")
        fi.write("Havriliak-Negami-beta")
        fi.write("\t")
        fi.write(str(HavNegBeta))
        fi.write("\n")
        fi.write("Regularization type")
        fi.write("\t")
        if PenaltyType == "linear":
            fi.write("linear regularization")
            fi.write("\n")
        if PenaltyType == "1st-order":
            fi.write("1st-derivative regularization")
            fi.write("\n")
        if PenaltyType == "2nd-order":
            fi.write("2nd-derivative regularization")
            fi.write("\n")   
        fi.write("Back-calculation")
        fi.write("\t")
        if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
            if KK_enable == 1:
                fi.write("Kramers-Kronig")
                fi.write("\n")
            if AntiKK_enable == 1:
                fi.write("Anti-Kramers-Kronig")
                fi.write("\n")
            if CombiKK_enable == 1:
                fi.write("Combi-Kramers-Kronig")
                fi.write("\n")    
        if KK_enable == 0 and AntiKK_enable == 0 and CombiKK_enable == 0:
            fi.write("None")
            fi.write("\n") 
        fi.write("\n")
        fi.write("log10(tau / s)")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("gamma(tau)_Imag")
            fi.write("\t")
        if Real_User == 1:
            fi.write("gamma(tau)_Real")
            fi.write("\t")
        if Comb_User == 1:
            fi.write("gamma(tau)_Combined")
        fi.write("\n")      
        for i in range(len(TAUARRAY)):
            fi.write(str(np.asscalar(np.log10(TAUARRAY[i]))))
            fi.write("\t")
            if Imag_User == 1:
                fi.write(str(np.asscalar(gamma_imag_HavriliakNegami[i])))
                fi.write("\t")
            if Real_User == 1:
                fi.write(str(np.asscalar(gamma_real_HavriliakNegami[i])))
                fi.write("\t")
            if Comb_User == 1:
                fi.write(str(np.asscalar(gamma_comb_HavriliakNegami[i])))
                fi.write("\t")
            fi.write("\n") 
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency [Hz]")
        fi.write("\t")
        fi.write("Z_real_used [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_used [Ohm]")
        fi.write("\t")
        if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
            if KK_enable == 1:
                fi.write("Z_real_from_ImDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_ReDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")
            if AntiKK_enable == 1:
                fi.write("Z_real_from_ReDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_ImDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")
            if CombiKK_enable == 1:
                fi.write("Z_real_from_CombDRT [Ohm]")
                fi.write("\t")
                fi.write("Z_imag_from_CombDRT [Ohm]")
                fi.write("\t")
                fi.write("Delta_Z_real_[%]")
                fi.write("\t")
                fi.write("Delta_Z_imag_[%]")
                fi.write("\t")       
        fi.write("\n")        
        for i in range(len(freq)):
            fi.write(str(np.asscalar(freq[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_meas[i].real)))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_meas[i].imag)))
            fi.write("\t")
            if KK_enable == 1 or AntiKK_enable == 1 or CombiKK_enable == 1:
                fi.write(str(np.asscalar(Z_Back[i].real)))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_Back[i].imag)))
                fi.write("\t")
                fi.write(str(100*np.asscalar(np.abs((Z_Back[i].real-Z_meas[i].real)/Z_meas[i].real))))
                fi.write("\t")
                fi.write(str(100*np.asscalar(np.abs((Z_Back[i].imag-Z_meas[i].imag)/Z_meas[i].imag))))
                fi.write("\t")
            fi.write("\n") 
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency_original [Hz]")
        fi.write("\t")
        fi.write("Z_real_original [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_originnal [Ohm]")
        fi.write("\n")
        for i in range(len(Frequency_Array)):
            fi.write(str(np.asscalar(Frequency_Array[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_real_Meas[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_imag_Meas[i])))
            fi.write("\n")    
    root.destroy()   
    
    
    
    
    
    

    
    
    
    
    
    
