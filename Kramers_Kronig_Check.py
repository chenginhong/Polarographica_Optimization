# -*- coding: utf-8 -*-
"""
#Author:       Tim Tichter
#Years:        2022
#Function:     POLAROGRAPHICAS KK-Check function (Based on DRT)
"""

import numpy                           as np
from scipy.optimize                    import nnls
from math                              import ceil,floor
from tkinter                           import *
from tkinter.filedialog                import askopenfilename
from tkinter.filedialog                import asksaveasfilename
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases          import key_press_handler
from matplotlib.figure                 import Figure
from matplotlib.backend_bases          import key_press_handler
from DRT_Native_Gauss_ColeCole_HavNeg  import build_RC, build_RegMat_Template, DRT_Native, KramersKronig_Multiplication
from DRT_Interface                     import multicolor_ylabel



def Get_Imp_KK_Data():
    Fenster = Toplevel()  
    Fenster.geometry("300x250")   
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
    
    Re_to_Ohm_Label = Label(Fenster,text="Z_real Factor to be Ohm", bg = '#80ffbf')
    Re_to_Ohm_Label.grid(row=3, column=0)
    Re_to_Ohm_Eingabe = Entry(Fenster)
    Re_to_Ohm_Eingabe.insert(END, 1)
    Re_to_Ohm_Eingabe.grid(row=3, column=1)
        
    Im_to_Ohm_Label = Label(Fenster,text="Z_imag Factor to be Ohm", bg = '#80ffbf')
    Im_to_Ohm_Label.grid(row=4, column=0)
    Im_to_Ohm_Eingabe = Entry(Fenster)
    Im_to_Ohm_Eingabe.insert(END, 1)
    Im_to_Ohm_Eingabe.grid(row=4, column=1)
    
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
    
    var1 = IntVar()
    Checkbutton(Fenster, text="Save KK-result", variable=var1, bg = '#80ffbf').place(x = 137, y = 150, height = 25)  
                
    def Accept_ImpDRT():
        global FromRowxx
        global ToRowxx
        global Readeveryxx
        global Re_to_Ohm
        global Im_to_Ohm
        global Delimiter
        global Save_KK_Result
        FromRowxx     = (int(FromRowxxx_Eingabe.get()))
        ToRowxx       = (int(ToRowxxx_Eingabe.get()))
        Readeveryxx   = (int(Readeveryxxx_Eingabe.get()))
        Re_to_Ohm     = (float(Re_to_Ohm_Eingabe.get()))
        Im_to_Ohm     = (float(Im_to_Ohm_Eingabe.get()))
        Delimiter     = var0.get()
        Save_KK_Result= var1.get()
    def Next_ImpDRT():
        Do_KK_Check()  
    Accept = Button(Fenster, text="Accept",command=Accept_ImpDRT)
    Accept.place( x = 25, y = 185, width = 115, height = 50)  
    Next = Button(Fenster, text="Next",command=Next_ImpDRT)
    Next.place( x = 150, y = 185, width = 115, height = 50)
    
    
    
    
def Do_KK_Check():
    root = Toplevel()
    root.title("Result from Kramers-Kronig check")
    if Delimiter == 1:
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='\t', skip_header = FromRowxx, skip_footer = ToRowxx)
    if Delimiter == 0:     
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='',   skip_header = FromRowxx, skip_footer = ToRowxx)
    #==================================================================================
    # Define the loaded data a global variables and put it into correct order
    #==================================================================================
    global Frequency
    global Z_Real
    global Z_Imag
    global Freq
    global Z_Real_used
    global Z_Imag_used
    global Z_Back_from_KK
    global Z_real_deviation
    global Z_imag_deviation
    Frequency         =  data[::Readeveryxx,0:1:1] 
    Z_Real            =  data[::Readeveryxx,1:2:1] * Re_to_Ohm            
    Z_Imag            =  data[::Readeveryxx,2:3:1] * Im_to_Ohm
    if Frequency[1] > Frequency[0]:            
        Frequency     = Frequency[::-1]
        Z_Real        = Z_Real[::-1]             
        Z_Imag        = Z_Imag[::-1]
    Freq              = Frequency[Z_Imag < 0]
    Z_Real_used       = Z_Real[Z_Imag < 0]
    Z_Imag_used       = Z_Imag[Z_Imag < 0]
    Real_Offset       = np.min(Z_Real_used)
    Z                 = Z_Real_used - Real_Offset + 1j*Z_Imag_used
    #==================================================================================
    # Compute the native DRT as an intermediate step of KK-correlation
    #==================================================================================
    IM_DRT_Native     = DRT_Native(freq_raw = Freq, Z_raw = Z, damp = 1e-4, low_ext = 5, high_ext = 3, resol = 10, DRT_Type = "imag", Reg_Type = 'linear')
    RE_DRT_Native     = DRT_Native(freq_raw = Freq, Z_raw = Z, damp = 1e-4, low_ext = 5, high_ext = 3, resol = 10, DRT_Type = "real", Reg_Type = 'linear')
    #==================================================================================
    # Perform Kramers_Kronig multiplication
    #==================================================================================
    Z_Back_from_KK    = KramersKronig_Multiplication(MATRIX = IM_DRT_Native[4], Real_DRT = RE_DRT_Native[1], Imag_DRT = IM_DRT_Native[1]) + Real_Offset
    Z_real_deviation  = 100*np.abs((Z_Back_from_KK.real-Z_Real_used)/Z_Real_used)
    Z_imag_deviation  = 100*np.abs((Z_Back_from_KK.imag-Z_Imag_used)/Z_Imag_used)
    #==================================================================================
    #Plotting of loaded file and KK result
    #==================================================================================
    f  = Figure(figsize=(11, 6), dpi=100)        
    b1 = f.add_subplot(121)
    b2 = f.add_subplot(122)  
    b1.plot(Z_Back_from_KK.real, -Z_Back_from_KK.imag, linestyle='-', linewidth = 0.5, marker='', color='red', label = 'KK-result')   
    b1.plot(Z_Real, -Z_Imag, linestyle='', marker='o', markersize = 5, color='lightgrey', label = 'native data')   
    b1.plot(Z_Real_used, -Z_Imag_used, linestyle='', marker='.', markersize = 3, color='black', label = 'used data')  
    b1.set_xlabel('$\mathfrak{RE}(Z(j\omega)/\Omega)$', fontsize = 15)
    b1.set_ylabel('$-\mathfrak{IM}(Z(j\omega)/\Omega)$', fontsize = 15)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=15)
    b1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=2, fontsize = 15)
    b2.plot(np.log10(Freq), -Z_Back_from_KK.imag, linestyle='-',marker='', linewidth = 0.5, color='blue', label = "Im from Re-DRT")
    b2.plot(np.log10(Freq),  Z_Back_from_KK.real, linestyle='-',marker='', linewidth = 0.5, color='red',  label = "Re from Im-DRT")
    b2.plot(np.log10(Frequency), -Z_Imag, linestyle='',marker='.', markersize = 3,  color='blue', label = 'Im-meas')
    b2.plot(np.log10(Frequency),  Z_Real, linestyle='',marker='.', markersize = 3,  color='red',  label = 'Re-meas')
    ax2 = b2.twinx()
    ax2.plot(np.log10(Freq), Z_imag_deviation, color = 'cyan',   linewidth = 0.5, label = '$\Delta$Im')
    ax2.plot(np.log10(Freq), Z_real_deviation, color = 'orange', linewidth = 0.5, label = '$\Delta$Re')
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
    
    if Save_KK_Result == 1:
        KK_Result_As_txt_saver()
    


    
################################################################################################################################
#    Optionally save the data from the KK-check
################################################################################################################################
    
def KK_Result_As_txt_saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Data of KK-check",filetypes = (("save as","*.txt"),("all files","*.*")))
    with open(root.filename,"w") as fi:
        fi.write( "Result from KK-check."
                 +"\nKK-check is essentially performed by using the native DRT on the EIS data provided by the user."
                 +"\nAfter reconstruction, the DRT is back-multiplied according to the Kramers-Kronig relations. This means that the"
                 +"\nreal-DRT is used to compute the imaginary part and the imaginary DRT is used to compute the real part "
                 +"\nof the experimentally measured EIS data. For computing the DRT, any inductive behaviour of the experimental"
                 +"\nEIS data is removed first. Also, the real offset is removed. The DRT-reconstruction follows the Native-DRT"
                 +"\nfunction which can be accessed in the DRT computation module. The lower range of the time-constant scale"
                 +"\nis extended by five decades, the upper range by two decades. The resolution is increased by factor ten."
                 +"\nThe Tikhonov regularization parameter is kept at 0.0001. Though this is not necessarily a good measure for"
                 +"\ncomputing the DRT itself (since it will overfit), it is a good choice for checking the Kramers-Kronig behaviour.")
        fi.write("\n")
        fi.write("\n")
        fi.write("Used Frequency [Hz]")
        fi.write("\t")
        fi.write("Z_real_used [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_used [Ohm]")
        fi.write("\t")
        fi.write("Z_real_back [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_back [Ohm]")
        fi.write("\t")
        fi.write("Z_real_deviation in [%]")
        fi.write("\t")
        fi.write("Z_imag_deviation in [%]")
        fi.write("\n")
        for i in range(len(Freq)):
            fi.write(str(np.asscalar(Freq[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_Real_used[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_Imag_used[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_Back_from_KK[i].real)))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_Back_from_KK[i].imag)))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_real_deviation[i])))
            fi.write("\t")  
            fi.write(str(np.asscalar(Z_imag_deviation[i])))
            fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Original EIS data")
        fi.write("\n")
        fi.write("Frequency_orig [Hz]")
        fi.write("\t")
        fi.write("Z_real_orig [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_orig [Ohm]")
        fi.write("\n")
        for i in range(len(Freq)):
            fi.write(str(np.asscalar(Frequency[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_Real[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_Imag[i])))
            fi.write("\n")
        fi.close()
    root.destroy()
        
    
    
    
