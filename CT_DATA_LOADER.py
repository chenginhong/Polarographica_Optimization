# -*- coding: utf-8 -*-
"""
Created on Tue May 11 10:03:13 2021

@author: gisel
"""
import matplotlib.pyplot as plt
import time
import numpy as np
from skimage import io
import os
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D
from math import ceil
from tkinter                           import *
from tkinter.filedialog                import askopenfilename
from tkinter.filedialog                import asksaveasfilename


#=========================================================================================================================
#=========================================================================================================================
#=========================================================================================================================
#Go for the computations
#=========================================================================================================================
#=========================================================================================================================
#=========================================================================================================================




def Get_CT_Data():
    Fenster = Toplevel()                                                         
    Fenster.title("Get CT Data")                         
    Fenster.geometry("300x200")  
    colorbgr = Label(Fenster, text= "", bg = '#80ffbf')
    colorbgr.place(x = 0, y = 0, width =7000, height = 2000)
    DefPath_Label = Label(Fenster, text="Specify path to image folder", bg = '#80ffbf', font = ('Arial',14,'bold'))                                          
    DefPath_Label.place(x = 20, y = 20)
    INPUTCT_label      = Label(Fenster,text="Path:", bg = '#80ffbf')
    INPUTCT_label.place(x = 20, y = 70)
    INPUTCT_Eingabe    = Entry(Fenster)                                               
    INPUTCT_Eingabe.place(x = 120, y = 70)
    var1 = IntVar()
    Checkbutton(Fenster, text="Plot Stack", variable=var1, bg = '#80ffbf').place(x = 20, y = 120)
    def NextCT():
        Pathstr    = (str(INPUTCT_Eingabe.get()))
        PlotStack  = var1.get()
        Fenster.destroy()
        CT_Data_Opener(Pathname = Pathstr, PLOTSTACK = PlotStack)
    Next = Button(Fenster, text="Next",command=NextCT)
    Next.place(x = 120, y = 110, width = 120, height = 50)







def CT_Data_Opener(Pathname, PLOTSTACK):
    Path = Pathname
    #Path = r'C:\Users\gisel\Desktop\20210401_H5_FINERESOL_EXTENDED1902x1902x300+200\20210402_Downpixeler'
    dir = Path
    #----------------------------------------------------------------------------
    # import dataset
    num_files =[]
    for img_files in os.listdir(dir):
        if img_files.endswith(".tif"):
            num_files.append(img_files)
    first_image = io.imread(dir + '/' + num_files[0])
    xsize = first_image.shape[0]
    ysize = first_image.shape[1]
    files = len(num_files)
    print("\n")
    print("x_size = ", xsize)
    print("y_size = ", ysize)
    print("z_size = ", files)
    
    dimx    = xsize
    dimy    = ysize
    dimz    = files
    stack   = np.zeros((dimx, dimy, dimz), dtype = np.float32)
    COUNTER = 0
    print('reading training dataset...')
    for n in tqdm(range(dimz)):
        if io.imread(dir + '/' + num_files[COUNTER]).ndim == 2:
            stack[:,:,n] = io.imread(dir + '/' + num_files[COUNTER])[::,::]
        elif io.imread(dir + '/' + num_files[COUNTER]).ndim == 3:
            stack[:,:,n] = io.imread(dir + '/' + num_files[COUNTER])[::,::,0]
        COUNTER += 1

    Stackk  = stack / 255.0      #white is 255 as color code
    stack   = (Stackk-1)**2      #electrode is zero, rest is one
    
    #Take a smaller teststack --> works fine :)
    #stack   = np.ones((50,50,50))
    #stack[0:50,0:50,0] = np.zeros((50,50))
    #START   = stack[:50:,:50:,:50:]
    
    START   = np.rint(stack)
    
    #==============================================================================
    OUTPATH  = os.path.join(Path, 'CT_DATA.npy')
    np.save(OUTPATH, START[::,::,::])
    
    if PLOTSTACK == 1:
        from mpl_toolkits.mplot3d import Axes3D
        def make_ax(grid=False):
            fig = plt.figure(figsize = (8,8))
            ax = fig.gca(projection='3d')
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            ax.grid(grid)
            return ax
        
        filled = (START-1)**2
        ax = make_ax(True)
        ax.grid(False)
        ax.voxels(filled, facecolors = 'k', edgecolors='gray')
        
        ax.view_init(30, 30)
        plt.show()
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


#Get_CT_Data()
