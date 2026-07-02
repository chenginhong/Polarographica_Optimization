# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 19:58:46 2021

@author: gisel
"""

#==========================================================================================================================
#Load all the required modules
#==========================================================================================================================
import matplotlib.pyplot as plt
import time
import numpy as np
import time
import os
import importlib
import tkinter                           as tk
from   tkinter                           import *
from   tkinter.filedialog                import askopenfilename
from   tkinter.filedialog                import asksaveasfilename



#=======================================================================================
#Check if numba for acceleration/multiprocessing is installed
#=======================================================================================
NumbaExister = importlib.util.find_spec('numba')
def DouglasGunnWithOrWithoutNumba():
    if NumbaExister == None:
        VerySlow_DourglasGunn()
    if NumbaExister is not None:
        Grid_Definer()
    
def VerySlow_DourglasGunn():
    Fenster = tk.Toplevel()                                                         
    Fenster.title("Painfully slow warner")                         
    WaitingDuck_image      = tk.PhotoImage(file = r"IMAGES\Waitingduck.png")       #Image as background in toplevel
    background_label       = tk.Label(Fenster, image=WaitingDuck_image)            #Image as background in toplevel
    background_label.image = WaitingDuck_image                                     #Image as background in toplevel
    w                      = WaitingDuck_image.width()                             #Image as background in toplevel
    h                      = WaitingDuck_image.height()                            #Image as background in toplevel
    xcoodinate             = (Fenster.winfo_screenwidth()  -  w)/2                 #Pop up window in center of screen
    ycoodinate             = (Fenster.winfo_screenheight()  - h)/2                 #Pop up window in center of screen
    Fenster.geometry('%dx%d+%d+%d' % (w,h, xcoodinate, ycoodinate))
    background_label.place(x=0, y=0, relwidth=1, relheight=1)
    Te = Text(Fenster, font = ('Arial',12,'bold'), padx = 20, pady = 10)
    Te.place(x = 320, y = 38, width = 250, height = 320)
    Te.insert(END, "WARNING! You do NOT \nhave installed the module \nNumba. Your computation\nwill take ages! \nWe highly recommend to \ninstall numba. If you are \nusing the Anaconda \npackage manager, you \ncan go for the command \n\nconda install numba \n\nMore info at\n\nhttps://numba.pydata.org/ ")
    def Run_DouglasGunn_WithoutNumba():
        Fenster.destroy()
        Grid_Definer()
    button=Button(Fenster,text="Continue without\nnumba\n(not recommended)", bg = '#9BA9C5', command = Run_DouglasGunn_WithoutNumba).place(x = 320, y = 385, width = 250, height = 75)
    Fenster.resizable(False, False)
    
    
    
#=======================================================================================
#This is the function for setting the input of the numerical computations
#======================================================================================= 
def Grid_Definer():
    Fenster = Toplevel()                                                         
    Fenster.title("Grid Definer")                         
    Fenster.geometry("300x250")    
    colorbgr = Label(Fenster, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000)                                        
    SpatialGrid_label      = Label(Fenster, text= u'Voxel size in \u03BCm', bg = '#80ffbf')
    SpatialGrid_label.place(x = 20, y = 20)
    SpatialGrid_Eingabe    = Entry(Fenster)  
    SpatialGrid_Eingabe.insert(END, 5)                                             
    SpatialGrid_Eingabe.place(x = 150, y = 20)
    DiffKoeff_label       = Label(Fenster,text= r'D in 10^-6 cm^2/s', bg = '#80ffbf')
    DiffKoeff_label.place(x =20, y = 45)
    DiffKoeff_Eingabe     = Entry(Fenster)
    DiffKoeff_Eingabe.insert(END, 1) 
    DiffKoeff_Eingabe.place(x = 150, y = 45)
    dt_min_label       = Label(Fenster,text= r'dt(min) / s', bg = '#80ffbf')
    dt_min_label.place(x =20, y = 70)
    dt_min_Eingabe     = Entry(Fenster)
    dt_min_Eingabe.insert(END, 0.001) 
    dt_min_Eingabe.place(x = 150, y = 70)
    dt_max_label       = Label(Fenster,text= r'dt(max) / s', bg = '#80ffbf')
    dt_max_label.place(x =20, y = 95)
    dt_max_Eingabe     = Entry(Fenster)
    dt_max_Eingabe.insert(END, 2) 
    dt_max_Eingabe.place(x = 150, y = 95)
    dt_Exp_label       = Label(Fenster,text= r'dt expansion', bg = '#80ffbf')
    dt_Exp_label.place(x =20, y = 120)
    dt_Exp_Eingabe     = Entry(Fenster)
    dt_Exp_Eingabe.insert(END, 0.05) 
    dt_Exp_Eingabe.place(x = 150, y = 120)
    dt_Num_label       = Label(Fenster,text= r'num of dt', bg = '#80ffbf')
    dt_Num_label.place(x =20, y = 145)
    dt_Num_Eingabe     = Entry(Fenster)
    dt_Num_Eingabe.insert(END, 410) 
    dt_Num_Eingabe.place(x = 150, y = 145)
    
    def NextGrid():
        global DiffKoeff
        DiffKoeff     = (float(DiffKoeff_Eingabe.get()))
        global SpatialGrid
        SpatialGrid       = (float(SpatialGrid_Eingabe.get()))
        global dtmin
        dtmin             = (float(dt_min_Eingabe.get()))
        global dtmax
        dtmax             = (float(dt_max_Eingabe.get()))
        global dtexp
        dtexp             = (float(dt_Exp_Eingabe.get()))
        global dtnum
        dtnum             = (int(dt_Num_Eingabe.get()))
        Fenster.destroy()
        GO_DOUGLAS_GUNN()
    
    Next = Button(Fenster, text="Load your\nnumpy file",command=NextGrid)
    Next.place(x = 150, y = 180, width = 120, height = 50)



#=======================================================================================
#define the Douglas Gunn related functions accordingly - here without numba
#=======================================================================================
if NumbaExister == None:
    #This is a TDMA solver aka the Thomas algorithm. The classical sp solver does not work here which is why I
    #wrote this solver.
    def Thomas(a,b,c,L):
        cs = np.zeros(len(c))
        Ls = np.zeros(len(L))
        x  = np.zeros(len(L))
        cs[0] = c[0]/b[0]
        Ls[0] = L[0]/b[0]
        for i in range(len(L)-1):
            if i>0:
                cs[i]= c[i]/(b[i] - a[i-1]*cs[i-1])
                Ls[i]= (L[i] - a[i-1]*Ls[i-1])/(b[i] - a[i-1]*cs[i-1])
        Ls[-1] = (L[-1] - a[-1]*Ls[-2])/(b[-1] - a[-1]*cs[-1])
        x[-1] = Ls[-1]
        for i in range(len(L)):
            if i > 0:
                x[-(i+1)] = Ls[-(i+1)] - cs[-(i)]*x[-(i)]
        return x
    #This function is the matrixbuilder for the Douglas-Gunn algorithm. It is basically a 1-D Crank Nicolson type of matrix
    #which involves any internal boundary conditions from a Starting-Array. It returns three vectors which will be used by the Thomas
    #function to solve for the next time-(sub)-iteration
    def Matbuilder(Array_Vector, ll):
        lenx       = len(Array_Vector)
        Middle     = (2+2*ll)*np.ones(lenx)
        Middle[0]  = (2+ll)
        Middle[-1] = (2+ll)
        OffDiag    = -ll*np.ones(lenx-1)
        Middle     = Middle*Array_Vector + (Array_Vector-1)**2
        up1        = OffDiag*Array_Vector[0:-1:]
        low1       = OffDiag*Array_Vector[1::]
        return low1,  Middle, up1
    def Deriv_3D(stackk, Startt):
        sta                     = np.zeros((len(stackk[::,0,0])+2, len(stackk[0,::,0])+2, len(stackk[0,0,::])+2))
        sta[1:-1:,1:-1:,1:-1:] += -6*stackk
        sta[0:-2:,1:-1:,1:-1:] += stackk[::,::,::]
        sta[2::,1:-1:,1:-1:]   += stackk[::,::,::]
        sta[1:-1:,0:-2:,1:-1:] += stackk[::,::,::]
        sta[1:-1:,2::,1:-1:]   += stackk[::,::,::]
        sta[1:-1:,1:-1:,0:-2:] += stackk[::,::,::]
        sta[1:-1:,1:-1:,2::]   += stackk[::,::,::]
        sta[1:-1:,1:-1:,1:-1:]  = sta[1:-1:,1:-1:,1:-1:]*(Startt-1)**2
        Result                 = np.sum(sta[1:-1:,1:-1:,1:-1:])
        return Result
 
#=======================================================================================
#define the Douglas Gunn related functions accordingly - here WITH Numba
#=======================================================================================
if NumbaExister is not None:
    from numba import jit
    #This is a TDMA solver aka the Thomas algorithm. The classical sp solver does not work here.
    @jit(nopython=True)     #Function decorator for multiproceccing and multitheading
    def Thomas(a,b,c,L):
        cs = np.zeros(len(c))
        Ls = np.zeros(len(L))
        x  = np.zeros(len(L))
        cs[0] = c[0]/b[0]
        Ls[0] = L[0]/b[0]
        for i in range(len(L)-1):
            if i>0:
                cs[i]= c[i]/(b[i] - a[i-1]*cs[i-1])
                Ls[i]= (L[i] - a[i-1]*Ls[i-1])/(b[i] - a[i-1]*cs[i-1])
        Ls[-1] = (L[-1] - a[-1]*Ls[-2])/(b[-1] - a[-1]*cs[-1])
        x[-1] = Ls[-1]
        for i in range(len(L)):
            if i > 0:
                x[-(i+1)] = Ls[-(i+1)] - cs[-(i)]*x[-(i)]
        return x
    #This function is the matrixbuilder for the Douglas-Gunn algorithm. It is basically a 1-D Crank Nicolson type of matrix
    #which involves any internal boundary conditions from a Starting-Array. It returns three vectors which will be used by the Thomas
    #function to solve for the next time-(sub)-iteration
    @jit(nopython=True)     #Function decorator for multiproceccing and multitheading
    def Matbuilder(Array_Vector, ll):
        lenx       = len(Array_Vector)
        Middle     = (2+2*ll)*np.ones(lenx)
        Middle[0]  = (2+ll)
        Middle[-1] = (2+ll)
        OffDiag    = -ll*np.ones(lenx-1)
        Middle     = Middle*Array_Vector + (Array_Vector-1)**2
        up1        = OffDiag*Array_Vector[0:-1:]
        low1       = OffDiag*Array_Vector[1::]
        return low1,  Middle, up1
    
    @jit(nopython=True)      #Function decorator for multiproceccing and multitheading
    def Deriv_3D(stackk, Startt):
        sta                     = np.zeros((len(stackk[::,0,0])+2, len(stackk[0,::,0])+2, len(stackk[0,0,::])+2))
        sta[1:-1:,1:-1:,1:-1:] += -6*stackk
        sta[0:-2:,1:-1:,1:-1:] += stackk[::,::,::]
        sta[2::,1:-1:,1:-1:]   += stackk[::,::,::]
        sta[1:-1:,0:-2:,1:-1:] += stackk[::,::,::]
        sta[1:-1:,2::,1:-1:]   += stackk[::,::,::]
        sta[1:-1:,1:-1:,0:-2:] += stackk[::,::,::]
        sta[1:-1:,1:-1:,2::]   += stackk[::,::,::]
        sta[1:-1:,1:-1:,1:-1:]  = sta[1:-1:,1:-1:,1:-1:]*(Startt-1)**2
        Result                 = np.sum(sta[1:-1:,1:-1:,1:-1:])
        return Result  


    
#This function is basically the heart of the computations of the 3D diffusion stuff. It is the Douglas gunn modification of the
#3D Crank Nicolson algorithm. It is particularly useful, since it splits each and every time step into three substeps and
#therefore involves tridiagonal matrices only. Since it removes the recursive behaviour of the 3D Crank-Nicolson technique
#by splitting each time iteration into 3N^2 steps, which are almost independent of each other, it suggests to utilize
#Multiprocessing for massive speed up in the computations.
def DouglasGunn(stack, START, ll):
    Interm_Array    = np.ones_like(stack)
    stack_one_Third = np.ones_like(stack)
    stack_two_Third = np.ones_like(stack)
    stack_complete  = np.ones_like(stack)
    #===========================================================================================================================
    #Timestep 1/3
    #===========================================================================================================================
    Interm_Array[1:-1:,1:-1:,1:-1:] = (ll*(stack[0:-2:,1:-1:,1:-1:] - 10*stack[1:-1:,1:-1:,1:-1:] + stack[2::,1:-1:,1:-1:]
                              + 2*(stack[1:-1:,0:-2:,1:-1:] + stack[1:-1:,2::,1:-1:]
                              + stack[1:-1:,1:-1:,0:-2:]    + stack[1:-1:,1:-1:,2::])) + 2*stack[1:-1:,1:-1:,1:-1:])*START
    #=====================================================================================================
    #SOLVE THE 1/3 STEP
    #=====================================================================================================
    for j in range(len(START[0,::,0])):
        for k in range(len(START[0,0,::])):
            aa,bb,cc                          = Matbuilder(START[::,j,k], ll)
            stack_one_Third[1:-1:,j+1,k+1]    = Thomas(aa,bb,cc, Interm_Array[1:-1:,j+1,k+1])
    stack_one_Third[ 0,1:-1:,1:-1:]           = stack_one_Third[ 1,1:-1:,1:-1:]
    stack_one_Third[-1,1:-1:,1:-1:]           = stack_one_Third[-2,1:-1:,1:-1:]
    stack_one_Third[1:-1:, 0,1:-1:]           = stack_one_Third[1:-1:, 1,1:-1:]
    stack_one_Third[1:-1:,-1,1:-1:]           = stack_one_Third[1:-1:,-2,1:-1:]
    stack_one_Third[1:-1:,1:-1:, 0]           = stack_one_Third[1:-1:,1:-1:, 1]
    stack_one_Third[1:-1:,1:-1:,-1]           = stack_one_Third[1:-1:,1:-1:,-2]
    #===========================================================================================================================
    #Timestep 2/3
    #===========================================================================================================================
    Interm_Array[1:-1:,1:-1:,1:-1:] = (ll*(stack_one_Third[0:-2:,1:-1:,1:-1:] - 2*stack_one_Third[1:-1:,1:-1:,1:-1:]
                            + stack_one_Third[2::,1:-1:,1:-1:]
                            +   stack[0:-2:,1:-1:,1:-1:] - 8*stack[1:-1:,1:-1:,1:-1:] + stack[2::,1:-1:,1:-1:]
                            +   stack[1:-1:,0:-2:,1:-1:] + stack[1:-1:,2::,1:-1:]
                            + 2*stack[1:-1:,1:-1:,0:-2:] + 2*stack[1:-1:,1:-1:,2::]) + 2*stack[1:-1:,1:-1:,1:-1:])*START
    #=====================================================================================================
    #SOLVE THE 2/3 STEP
    #=====================================================================================================
    for i in range(len(START[::,0,0])):
        for k in range(len(START[0,0,::])):
            aa,bb,cc                          = Matbuilder(START[i,::,k], ll)
            stack_two_Third[i+1,1:-1:,k+1]    = Thomas(aa,bb,cc, Interm_Array[i+1,1:-1:,k+1])
    stack_two_Third[ 0,1:-1:,1:-1:]           = stack_two_Third[ 1,1:-1:,1:-1:]
    stack_two_Third[-1,1:-1:,1:-1:]           = stack_two_Third[-2,1:-1:,1:-1:]
    stack_two_Third[1:-1:, 0,1:-1:]           = stack_two_Third[1:-1:, 1,1:-1:]
    stack_two_Third[1:-1:,-1,1:-1:]           = stack_two_Third[1:-1:,-2,1:-1:]
    stack_two_Third[1:-1:,1:-1:, 0]           = stack_two_Third[1:-1:,1:-1:, 1]
    stack_two_Third[1:-1:,1:-1:,-1]           = stack_two_Third[1:-1:,1:-1:,-2]
    #===========================================================================================================================
    #Timestep 3/3
    #===========================================================================================================================
    Interm_Array[1:-1:,1:-1:,1:-1:] = (ll*(stack_one_Third[0:-2:,1:-1:,1:-1:] - 2*stack_one_Third[1:-1:,1:-1:,1:-1:]
                            +   stack_one_Third[2::,1:-1:,1:-1:]
                            +   stack_two_Third[1:-1:,0:-2:,1:-1:] - 2*stack_two_Third[1:-1:,1:-1:,1:-1:]
                            +   stack_two_Third[1:-1:,2::,1:-1:]
                            +   stack[0:-2:,1:-1:,1:-1:] - 6*stack[1:-1:,1:-1:,1:-1:] + stack[2::,1:-1:,1:-1:]
                            +   stack[1:-1:,0:-2:,1:-1:] + stack[1:-1:,2::,1:-1:]
                            +   stack[1:-1:,1:-1:,0:-2:] + stack[1:-1:,1:-1:,2::]) + 2*stack[1:-1:,1:-1:,1:-1:])*START
    for i in range(len(START[::,0,0])):
        for j in range(len(START[0,::,0])):
            aa,bb,cc                         = Matbuilder(START[i,j,::], ll)
            stack_complete[i+1,j+1,1:-1:]    = Thomas(aa,bb,cc, Interm_Array[i+1,j+1,1:-1:])
    stack_complete[ 0,1:-1:,1:-1:]           = stack_complete[ 1,1:-1:,1:-1:]
    stack_complete[-1,1:-1:,1:-1:]           = stack_complete[-2,1:-1:,1:-1:]
    stack_complete[1:-1:, 0,1:-1:]           = stack_complete[1:-1:, 1,1:-1:]
    stack_complete[1:-1:,-1,1:-1:]           = stack_complete[1:-1:,-2,1:-1:]
    stack_complete[1:-1:,1:-1:, 0]           = stack_complete[1:-1:,1:-1:, 1]
    stack_complete[1:-1:,1:-1:,-1]           = stack_complete[1:-1:,1:-1:,-2]

    return stack_complete





def GO_DOUGLAS_GUNN():
    #=============================================================================================================================
    #=============================================================================================================================
    #=============================================================================================================================
    #=============================================================================================================================
    print("Compiling Douglas-Gunn related functions")
    dimx, dimy, dimz  = 20,20,20
    
    stack             = np.ones((dimx,dimy,dimz))
    stack[::,::,0]    = np.zeros((dimx,dimy))
    START             = stack
    #print(START)
    ActSites          = Deriv_3D(START, START)
    print ("Initializing functions... wait a moment.")
    dx                = SpatialGrid/10000.0           # calculate from mu m to cm
    D                 = DiffKoeff/1000000.0           # since the input is *1000000
    timepoints        = np.arange(0,10,1)
    dt_zero           = dtmin
    dt_max            = dtmax
    dt_max_arr        = dt_max*np.ones(len(timepoints))
    t_expand          = dtexp
    dt_nonlim         = dt_zero*np.exp(t_expand*timepoints)
    dt_bounded        = dt_nonlim*dt_max_arr/(dt_nonlim + dt_max_arr)
    t                 = np.zeros(len(dt_bounded))
    tt                = 0
    for i in range(len(dt_bounded)):
        tt     += dt_bounded[i]
        t[i]    = tt
    #print ("maximum time is:",max(t))
    #===================================================================
    #Go for Douglas Gunn
    #===================================================================
    ll_Array     = D*dt_bounded/dx**2  #0.5*np.ones(len(dt_bounded))#
    Flux_Array   = np.zeros(len(ll_Array))
    st           = np.zeros((len(stack[::,0,0])+2, len(stack[0,::,0])+2, len(stack[0,0,::])+2))
    st[1:-1:,1:-1:,1:-1:]        = stack
    st[ 0,1:-1:,1:-1:]           = stack[ 0,::,::]
    st[-1,1:-1:,1:-1:]           = stack[-1,::,::]
    st[1:-1:, 0,1:-1:]           = stack[::, 0,::]
    st[1:-1:,-1,1:-1:]           = stack[::,-1,::]
    st[1:-1:,1:-1:, 0]           = stack[::,::, 0]
    st[1:-1:,1:-1:,-1]           = stack[::,::,-1]
    
    time_start = time.perf_counter()
    sttt = DouglasGunn(st, START, ll_Array[0])
    print("Initializing functions took", time.perf_counter() - time_start, "seconds")
    

    #=========================================================================================
    #=========================================================================================
    #=========================================================================================
    #RE-DEFINE to use the compiled functions
    #=========================================================================================
    #=========================================================================================
    #=========================================================================================

    Path                 = askopenfilename()
    stack                = np.load(Path)
    START                = stack.copy()
    
    #=====================================================================================================================
    #Define grid sizes
    #=====================================================================================================================
    
    dx                   = SpatialGrid/10000.0           # calculate from mu m to cm
    D                    = DiffKoeff/1000000.0           # since the input is *1000000
    timepoints           = np.arange(0,dtnum,1)
    dt_zero              = dtmin
    dt_max               = dtmax
    dt_max_arr           = dt_max*np.ones(len(timepoints))
    t_expand             = dtexp
    dt_nonlim            = dt_zero*np.exp(t_expand*timepoints)
    dt_bounded           = dt_nonlim*dt_max_arr/(dt_nonlim + dt_max_arr)
    t                    = np.zeros(len(dt_bounded))
    tt                   = 0
    
    for i in range(len(dt_bounded)):
        tt     += dt_bounded[i]
        t[i]    = tt
    #====================================================================================================================
    #Go for Douglas Gunn
    #====================================================================================================================
    ll_Array              = D*dt_bounded/dx**2
    Flux_Array            = np.zeros(len(ll_Array))
    st                    = np.zeros((len(stack[::,0,0])+2, len(stack[0,::,0])+2, len(stack[0,0,::])+2))
    st[1:-1:,1:-1:,1:-1:] = stack
    st[ 0,1:-1:,1:-1:]    = stack[ 0,::,::]
    st[-1,1:-1:,1:-1:]    = stack[-1,::,::]
    st[1:-1:, 0,1:-1:]    = stack[::, 0,::]
    st[1:-1:,-1,1:-1:]    = stack[::,-1,::]
    st[1:-1:,1:-1:, 0]    = stack[::,::, 0]
    st[1:-1:,1:-1:,-1]    = stack[::,::,-1]
    #====================================================================================================================
    #Douglas-Gunn-Main-Loop
    #====================================================================================================================
    print ("maximum time is:",max(t))
    staaart = time.perf_counter()
    
    
    
    for i in range(len(timepoints)):
        start = time.perf_counter()
        if i == 0:
            ActSites  = Deriv_3D(START, START)
            print ("Number of active sites is:", ActSites)
        if i == 2:
            StartCounting = time.perf_counter()
        if i == 12:
            EndCounting = time.perf_counter()
            print("10 steps 3Dderiv Jit took", EndCounting-StartCounting, "s")
        Flux_Array[i] = D*Deriv_3D(st[1:-1:, 1:-1:, 1:-1:], START)/(ActSites*dx)
        st         = DouglasGunn(st, START, ll_Array[i])
        print ("computation time of current full-loop timestep",i, " of ", len(timepoints), ":", "\t", (time.perf_counter()-start), "s")
        
        if i == 0:
            OUTPATH  = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))
        file = open(OUTPATH, "w")
        for j in range(len(t)):
            file.write(str(t[j]))
            file.write("\t")
            file.write(str(Flux_Array[j]/D**0.5))
            file.write("\n")
        file.close()
    
    print ("Total time was", (time.perf_counter()-staaart), "s")
    
    
    #===================================================================================================================
    #Results
    #===================================================================================================================
    print ("maximum time is:",max(t))
    print ("Number of active sites is:", ActSites)
    plt.figure(figsize = (5,5), dpi = 100)
    plt.plot(np.log10(t),np.log10(Flux_Array), color = 'k', linestyle = ':', linewidth = 2.5, label = 'Your Simulation')
    plt.plot(np.log10(t),np.log10(D**0.5/(t*np.pi)**0.5), color = 'r', label = 'planar Semiinf.')
    plt.xlabel(r'$\mathrm{log}_{10}(t / \mathrm{s})$', fontsize = 15)
    plt.ylabel(r'$\mathrm{log}_{10}(J_{\mathrm{norm}}(t))$', fontsize = 15)
    plt.legend(frameon = False, fontsize = 15)
    plt.ylim(-7,-1)
    plt.xlim(-3,3)
    plt.show()
    
    



