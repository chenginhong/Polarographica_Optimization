# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:30:20 2023

@author: timtic
"""

####################################################################################################
####################################################################################################
#                    Importing required modules
####################################################################################################
####################################################################################################
import numpy as np
import struct
"""
####################################################################################################
# tkinter for building GUIs 
####################################################################################################"""
import tkinter as tk 
from   tkinter import ttk                           # Python 3
from   tkinter import *                             # Python 3
from   tkinter import messagebox                    # Python 3
from   tkinter import Menu                          # Python 3
from   tkinter import filedialog                    # Python 3
from   tkinter.scrolledtext  import ScrolledText    # Python 3
####################################################################################################
# matplotlib related modules for graphics
####################################################################################################
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
####################################################################################################
# Implement the default Matplotlib key bindings.
####################################################################################################
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
'''
####################################################################################################
# time for timed events, serial for serial
# communication of Arduino with the PC. The 
# serial requires the installation of pyserial
# on the operating PC.
####################################################################################################'''
import time
import serial
import threading
import os
####################################################################################################
from datetime import datetime        # import datetome for checking recent cal file     
####################################################################################################
# Everything for threading the Data aqcuisition out of the GUI
####################################################################################################
from    threading                                import Thread
from    PolArStat_CV_Script_Acquisition          import ACQUISITION_BACKEND     as CVAQB
from    PolArStat_CV_Script_Write_Outputs        import PolArStat_CV_Write_Read as PCVWR
from    PolArStat_Serial_Communication_Functions import Serial_Communication    as SERCO
from    PolArStat_Serial_Communication_Functions import *
####################################################################################################
####################################################################################################  





def PotentiostatScript_CV():
    ####################################################################################################
    #    Define plotting options to make the plot look cool :)
    ####################################################################################################
    font = {'family': 'Times New Roman', 'color':  'black','weight': 'normal','size': 15,}
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    plt.rcParams['font.sans-serif'] = ['Times new Roman']
    ####################################################################################################
    # Initialize some global variables
    ####################################################################################################
    global PlotType      ;   PlotType      = 1
    global arduino       ;   arduino       = None      # initialize the device as None until it is defined
    global UpCounter     ;   UpCounter     = 0         # used for counting in the output window
    global CalData       ;   CalData       = 0
    global Zero_IDX_E    ;   Zero_IDX_E    = 13250
    global Zero_IDX_I    ;   Zero_IDX_I    = 13250
    global Initiate_Plot ;   Initiate_Plot = True      # For first setup of a plot 
    global thread_1
    ####################################################################################################
    # Check, if there is a recent calibration file. If so, load it and use it. If not, 
    # show a respective waring and start the interface
    ####################################################################################################
    try: 
        CalFile    = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CALIBRATION', datetime.today().strftime('%Y_%m_%d_Cal.txt'))
        CalData    = np.genfromtxt(CalFile , skip_header = 0, skip_footer = 1)
        Neg_IDX_E  = np.average(CalData[CalData[::,1]==1,3])
        Pos_IDX_E  = np.average(CalData[CalData[::,1]==2,3])
        Neg_IDX_I  = np.average(CalData[CalData[::,1]==1,4])
        Pos_IDX_I  = np.average(CalData[CalData[::,1]==2,4])
        Zero_IDX_E = int(0.5*(Neg_IDX_E + Pos_IDX_E))
        Zero_IDX_I = int(0.5*(Neg_IDX_I + Pos_IDX_I))
    except Exception:
            messagebox.showwarning(title="No Recent Calibration found!", message="You did not provide a recent calibration file!\n Your data might be incorrect.")
            txt_filename = ""
            
    ####################################################################################################
    # Change the variables in the class of the PolArStat_CV_Script_Write_Outputs script
    # according to the classical values or the values from the calibration file
    ####################################################################################################
    PCVWR.Zero_IDX_E = Zero_IDX_E
    PCVWR.Zero_IDX_I = Zero_IDX_I  
     
    ####################################################################################################
    # The following function asks for a filename, where the data will be stored
    ####################################################################################################
    def AskSaveFile_AndStart():
        try:
            txt_filename = filedialog.asksaveasfile(defaultextension='.txt').name
        except Exception:
            messagebox.showwarning(title="No File selected!", message="You did not select a file to save to!")
            txt_filename = ""
        else:
            outfile_CV_Data   = open(txt_filename, 'w')
            PCVWR.Output_File = outfile_CV_Data
            START_CV()
        
    ######################################################################################################
    # The following functionwill find available Arduino-ports which are able for 
    # serial communication. This function wil be embedded into a GUI later on.
   #######################################################################################################
    def StartSerial():
        global arduino
        SSer_btn.config(state='disabled')
        def _do_serial():
            global arduino
            try:
                arduino = serial.Serial(portList.get(), baudrate = 115200, timeout=.1)
                time.sleep(2)  # Wait for Arduino bootloader reset after DTR
                arduino.reset_input_buffer()  # Flush any bootloader garbage
                randfloat  = 11.01
                SEND_BYTES = b'\x44\x66' + struct.pack('f', randfloat)
                arduino.write(SEND_BYTES)
                time.sleep(2)  # Wait for Arduino to process and echo back
                SER_OUT_1  = arduino.readline()[:-2]
                SER_OUT_2  = arduino.readline()[:-2]
                if (SER_OUT_1 == SEND_BYTES) & (float(SER_OUT_2) == np.round(randfloat, decimals = 2)):
                    return ('ok',)
                else:
                    return ('corrupted',)
            except (serial.SerialException, OSError) as e:
                return ('error', str(e))
            except Exception as e:
                return ('error', str(e))
        def _on_serial_done(result):
            SSer_btn.config(state='normal')
            status = result[0]
            if status == 'ok':
                PCVWR.WRITE_TEXT_OUTPUT(Output_text, "Serial communication at Port %s established successfully" %portList.get())
            elif status == 'corrupted':
                messagebox.showwarning(title="Serial Connection Corrupted!",
                                       message="It seems like the data recieved over Serial is corrupted!")
            else:
                messagebox.showerror(title="Serial Communication failed!",
                                     message="Unable to communicate with selected Port or no Port selected")
        def _run():
            result = _do_serial()
            PotentiostatWindow.after(0, _on_serial_done, result)
        threading.Thread(target=_run, daemon=True).start()
                
    ###################################################################################################### 
    # The following function will set the inputs for a CV-measurement
    ###################################################################################################### 
    def Set_CV_Params():
        global ReadResist
        CVSet_btn.config(state='disabled')
        try:
            ReadResist   = float(R_read_Eingabe.get())
            try:
                Exper_Notes = str(Exper_Notes_Eingabe.get())
                if Exper_Notes != "":
                    PCVWR.Exper_Notes = Exper_Notes
                if Exper_Notes == "":
                    PCVWR.Exper_Notes = "No notes specified"
            except Exception:
                pass
            PCVWR.ReadResist = ReadResist
            E_initial    = float(E_in_Eingabe.get())
            E_vert1      = float(E_vert1_Eingabe.get())
            E_vert2      = float(E_vert2_Eingabe.get())
            E_fin        = float(E_fin_Eingabe.get())
            Scanrate     = float(Scanrate_Eingabe.get())
            Cycles       = float(Cycles_Eingabe.get())
            Conditime    = float(Condi_Eingabe.get())
        except (ValueError, TypeError) as e:
            CVSet_btn.config(state='normal')
            messagebox.showerror(title="Incomplete Parameters!",
                                 message="All inputs have to be filled out in order to run a CV!\n\nDetail: %s" % str(e))
            return
        if arduino is None or (hasattr(arduino, 'is_open') and not arduino.is_open):
            CVSet_btn.config(state='normal')
            messagebox.showerror(title="Serial Communication failed!",
                                 message="No serial connection. Please Start Serial first!")
            return
        def _do_work():
            SEND_E_in  = b'\x11' + b'\x10' + struct.pack('f', E_initial) + b'\x00\x00'
            SEND_E_v1  = b'\x11' + b'\x11' + struct.pack('f', E_vert1) + b'\x00\x00'
            SEND_E_v2  = b'\x11' + b'\x12' + struct.pack('f', E_vert2) + b'\x00\x00'
            SEND_E_fi  = b'\x11' + b'\x13' + struct.pack('f', E_fin) + b'\x00\x00'
            SEND_Cycl  = b'\x11' + b'\x14' + struct.pack('f', Cycles) + b'\x00\x00'
            SEND_ScaR  = b'\x11' + b'\x15' + struct.pack('f', Scanrate) + b'\x00\x00'
            SEND_Cond  = b'\x11' + b'\x16' + struct.pack('f', Conditime) + b'\x00\x00'
            for pkt in [SEND_E_in, SEND_E_v1, SEND_E_v2, SEND_E_fi, SEND_Cycl, SEND_ScaR, SEND_Cond]:
                arduino.write(pkt)
                time.sleep(0.05)  # Small delay for Arduino to process each parameter
        def _on_done():
            CVSet_btn.config(state='normal')
            PCVWR.WRITE_TEXT_OUTPUT(Output_text, "\n Data transmission complete :) \n---------------------------------------------")
        def _on_error(err):
            CVSet_btn.config(state='normal')
            messagebox.showerror(title="Send Parameters Failed!",
                                 message="Failed to send CV parameters to device.\n\nDetail: %s" % str(err))
        def _run():
            try:
                _do_work()
                PotentiostatWindow.after(0, _on_done)
            except (serial.SerialException, OSError) as e:
                print("[CV Set_CV_Params Error]", e)
                PotentiostatWindow.after(0, _on_error, e)
            except Exception as e:
                print("[CV Set_CV_Params Error]", e)
                PotentiostatWindow.after(0, _on_error, e)
        threading.Thread(target=_run, daemon=True).start()
        
    ###################################################################################################### 
    # The following function will retreive the transmitted parameters for a CV-measurement
    ###################################################################################################### 
    def Get_CV_Params():
        CVGet_btn.config(state='disabled')
        def _do_work():
            print('[DEBUG] arduino:', arduino, 'is_open:', arduino.is_open if arduino else 'None')
            print('[DEBUG] DEVICE:', CVAQB.DEVICE, 'running:', CVAQB.running)
            CVAQB.DEVICE = arduino
            CVAQB.running = False
            print('[DEBUG] Sending query command...')
            Readback_Array     = np.zeros(8)
            GETTING_COMMAND_BYTES  = b'\x22\x00\x00\x00\x00\x00\x00\x00'
            arduino.reset_input_buffer()
            arduino.write(GETTING_COMMAND_BYTES)
            time.sleep(2)  # Wait for Arduino to process and send back parameters
            print('[DEBUG] Reading response...')
            for i in range(7):
                raw = arduino.readline()[:-2]
                if not raw:
                    raw = arduino.readline()[:-2]
                Readback_Array[i] = float(raw)
            Readback_Array[7]  = ReadResist
            CVAQB.CV_Inputs    = Readback_Array
            return Readback_Array
        def _on_done(Readback_Array):
            CVGet_btn.config(state='normal')
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"The following Parameters were set! \n")
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"E_init. \t= \t %s V vs. Re"     %Readback_Array[0])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"E_vertex 1.\t = \t %s V vs. Re" %Readback_Array[1])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"E_vertex 2.\t = \t %s V vs. Re" %Readback_Array[2])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"E_final.\t = \t %s V vs. Re"    %Readback_Array[3])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"n-Cycles \t = \t %s "           %Readback_Array[4])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"Scanrate \t= \t %s mV/s"        %Readback_Array[5])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"Cond. time \t = \t %s s"        %Readback_Array[6])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"R Read \t = \t %s Ohm"          %Readback_Array[7])
        def _on_error(err):
            CVGet_btn.config(state='normal')
            messagebox.showerror(title="Read Parameters Failed!",
                                 message="Failed to read CV parameters from device.\n\nDetail: %s" % str(err))
        def _run():
            try:
                result = _do_work()
                PotentiostatWindow.after(0, _on_done, result)
            except Exception as e:
                PotentiostatWindow.after(0, _on_error, e)
        threading.Thread(target=_run, daemon=True).start()
        
    
    ####################################################################################################
    # The following function will start a CV measurment
    ####################################################################################################
    def START_CV():
        global UpCounter    ;   UpCounter= 0
        if arduino is None:
            messagebox.showwarning(title="No device in list!",
                                   message="No device selected. Initialize serial communication before running an experiment!")
            return
        CVRun_btn.config(state='disabled')
        def _do_work():
            CVAQB.DEVICE = arduino
            CVAQB.running = False
            Readback_Array     = np.zeros(8)
            GETTING_COMMAND_BYTES  = b'\x22\x00\x00\x00\x00\x00\x00\x00'
            arduino.reset_input_buffer()
            arduino.write(GETTING_COMMAND_BYTES)
            time.sleep(2)  # Wait for Arduino to process and send back parameters
            for i in range(7):
                raw = arduino.readline()[:-2]
                if not raw:
                    raw = arduino.readline()[:-2]
                Readback_Array[i] = float(raw)
            Readback_Array[7]  = ReadResist
            CVAQB.CV_Inputs    = Readback_Array
            SEND_BYTES = b'\x33\x01\x00\x00\x00\x00\x00\x00'
            SERCO.Send_And_Wait_For_ACK(arduino, SEND_BYTES)
            Starter = arduino.readline()[:-2]
            while Starter != b'10101010':
                time.sleep(0.01)
                Starter = arduino.readline()[:-2]
            return Readback_Array
        def _on_done(Readback_Array):
            CVRun_btn.config(state='normal')
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"The following Parameters were set! \n")
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"E_init. \t= \t %s V vs. Re"     %Readback_Array[0])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"E_vertex 1.\t = \t %s V vs. Re" %Readback_Array[1])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"E_vertex 2.\t = \t %s V vs. Re" %Readback_Array[2])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"E_final.\t = \t %s V vs. Re"    %Readback_Array[3])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"n-Cycles \t = \t %s "           %Readback_Array[4])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"Scanrate \t= \t %s mV/s"        %Readback_Array[5])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"Cond. time \t = \t %s s"        %Readback_Array[6])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"R Read \t = \t %s Ohm"          %Readback_Array[7])
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"\n Preparing CV completed.\n\n")
            PCVWR.WRITE_TEXT_OUTPUT(Output_text,"t/s\t E/V\t I/mA\t Cyc.No.\n")
            CVAQB.CLEAR_STORAGE()
            CVAQB.DEVICE  = arduino
            CVAQB.running = True
            CVAQB.Looper_ON = True
            Serial_Looper(plot1, canvas)
        def _on_error(err):
            CVRun_btn.config(state='normal')
            messagebox.showerror(title="Incomplete Parameters!",
                                 message="All inputs have to be filled out in order to run a CV!")
        def _run():
            try:
                result = _do_work()
                PotentiostatWindow.after(0, _on_done, result)
            except Exception as e:
                PotentiostatWindow.after(0, _on_error, e)
        threading.Thread(target=_run, daemon=True).start()

    ####################################################################################################
    # Function, to stop a CV measurement. The Portlist_Refresh() function is defined in the 
    # PolArStat_Serial_Communication_Functions script and re-imported above
    ####################################################################################################   
    def Stop_CV():
        if CVAQB.running == False:
            messagebox.showwarning(title="No active measurement!", 
                                   message="There is no measurement which can be stopped!")
        else:
            if CVAQB.running == True:
                PCVWR.WRITE_TEXT_OUTPUT(Output_text,"\n Status: Measurment was manually stopped by User.\n") 
            CVAQB.StopMeasure = True
            CVAQB.Looper_ON   = False
       
    ####################################################################################################
    # The following function will kill thread_1 and exit GUI cleanly.
    ####################################################################################################
    def _quit():
        CVAQB.ALIVE = False
        CVAQB.running = False
        CVAQB.Looper_ON = False
        thread_1.join(timeout = 3)
        PotentiostatWindow.quit()
        PotentiostatWindow.destroy()

    ####################################################################################################
    # Function, calling itself
    ####################################################################################################   
    def Serial_Looper(plot1, canvas):          # This function will be called in a loop back and forth with the main_loop of the GUI to update the output
        if CVAQB.Looper_ON == True:
            #===========================================================================================
            # Everything important for updating the plot
            #===========================================================================================
            sliced = CVAQB.Storage_Array[1::,::]
            if len(sliced) < 2:
                PotentiostatWindow.after(100, Serial_Looper, plot1, canvas)
                return
            InputArray = CVAQB.ARRAYCONDENSER(ARRAY = sliced)
            Redraw_Plot()
            #===========================================================================================
            # Everything important for updating the serial Output window
            #===========================================================================================
            global UpCounter
            PCVWR.WRITE_DATA_OUTPUT(WHERE_TO_WRITE = Output_text, DATA_TO_WRITE = InputArray[UpCounter::,::])
            UpCounter = len(InputArray[::,0])
            Output_text.see(tk.END)
            #===========================================================================================
            # Callback to looper function - only if the respective bool is True
            #===========================================================================================
            PotentiostatWindow.after(100, Serial_Looper, plot1, canvas)
       
    ####################################################################################################
    # Shared function to redraw the plot based on current PlotType and available data
    ####################################################################################################
    def Redraw_Plot():
        sliced = CVAQB.Storage_Array[1::,::]
        if len(sliced) < 2:
            return
        InputArray = CVAQB.ARRAYCONDENSER(ARRAY = sliced)
        if len(InputArray) < 3:
            return
        plot1.clear()
        if PlotType == 1:
            plot1.plot(-0.000249*(InputArray[2::,1]-Zero_IDX_E) ,  -0.12452*(InputArray[2::,2]-Zero_IDX_I)/ReadResist   )
            plot1.set_xlabel("$E$ vs. RE in V", fontsize = 13)
            plot1.set_ylabel("$I$ in mA", fontsize = 13)
        if PlotType == 2:
            plot1.plot((1e-6)*InputArray[2::,0], -0.12452*(InputArray[2::,2]-Zero_IDX_I)/ReadResist )
            plot1.set_xlabel("$t$ in s", fontsize = 13)
            plot1.set_ylabel("$I$ in mA", fontsize = 13)
        if PlotType == 3:
            plot1.plot((1e-6)*InputArray[2::,0], -0.000249*(InputArray[2::,1]-Zero_IDX_E) )
            plot1.set_xlabel("$t$ in s", fontsize = 13)
            plot1.set_ylabel("$E$ vs. RE in V", fontsize = 13)
        plot1.grid(which='both', linestyle = '--')
        plot1.tick_params(direction = 'in', length=4, width=0.5, colors='k', labelsize = 13)
        canvas.draw()

    ####################################################################################################
    # The following three functions are used to define the way in which the data is displayed
    ####################################################################################################
    def Show_E_vs_I():
        global PlotType
        PlotType = 1
        Redraw_Plot()

    def Show_t_vs_I():
        global PlotType
        PlotType = 2
        Redraw_Plot()

    def Show_t_vs_E():
        global PlotType
        PlotType = 3
        Redraw_Plot()
    
    # Stop any leftover thread from a previous window session
    CVAQB.ALIVE = False
    CVAQB.running = False
    CVAQB.Looper_ON = False
    time.sleep(0.5)
    ####################################################################################################
    # Start the conditional data-acquisition on its own thread, that the GUI does not freeze
    # if data collection is busy
    ####################################################################################################
    thread_1 = Thread(target=CVAQB.LINE_READER_FOR_THREAD, daemon=True)
    thread_1.start()
    ####################################################################################################
    # Run the PolArStat_CV_Script in a mainloop
    ####################################################################################################
    if __name__ == 'PolArStat_CV_Script_GUI':
        #==============================================================================================
        # The following part is the main loop which calls the GUI for interfacing with 
        # the Arduino. The Arduino has to contain the firmware before this script works.
        # First, the size of the GUI-window will be defined. Subsequently, all input field are set
        #==============================================================================================
        PotentiostatWindow = tk.Tk()
        PotentiostatWindow.title("Main Interface for CV-measurement with PolArStat Potentiostat")
        PotentiostatWindow.geometry('800x460')
        menu = Menu(PotentiostatWindow)                                               
        PotentiostatWindow.config(menu=menu)
        #==============================================================================================
        # Set all buttons and the dropdown menu for initializing serial communication
        #==============================================================================================
        labelTop = tk.Label(PotentiostatWindow, text = "Choose COM Port")
        labelTop.place(x = 25, y = 10)
        #-----------------------------------------------------------------------------------------------
        portsl = SERCO.portsl
        #-----------------------------------------------------------------------------------------------
        portList = ttk.Combobox(PotentiostatWindow, postcommand=lambda: portList.configure(values=SERCO.portsl))
        portList.place(x = 140, y = 10, width = 135, height = 19)
         #-----------------------------------------------------------------------------------------------       
        Refresh_btn = ttk.Button(PotentiostatWindow, text="Refresh Portlist", command=SERCO.FindPorts)   
        Refresh_btn.place(x = 25, y = 35, width = 120, height = 25)
        #-----------------------------------------------------------------------------------------------
        SSer_btn = ttk.Button(PotentiostatWindow, text="Start Serial", command=StartSerial)   
        SSer_btn.place(x = 155, y = 35, width = 120, height = 25)
        #==============================================================================================
        # Set all the input buttons for CV measurement and place them in the GUI
        #==============================================================================================
        E_in_Label = Label(PotentiostatWindow,text="E initial vs. Ref [V]*")
        E_in_Label.place(x = 25, y = 70)
        E_in_Eingabe = Entry(PotentiostatWindow)
        E_in_Eingabe.insert(END, -0.5) 
        E_in_Eingabe.place(x = 150, y = 70)
        #-----------------------------------------------------------------------------------------------
        E_vert1_Label = Label(PotentiostatWindow,text="E vertex 1 vs. Ref [V]*")
        E_vert1_Label.place(x = 25, y = 95)
        E_vert1_Eingabe = Entry(PotentiostatWindow)
        E_vert1_Eingabe.insert(END, 0.5) 
        E_vert1_Eingabe.place(x = 150, y = 95)
        #-----------------------------------------------------------------------------------------------
        E_vert2_Label = Label(PotentiostatWindow,text="E vertex 2 vs. Ref [V]*")
        E_vert2_Label.place(x = 25, y = 120)
        E_vert2_Eingabe = Entry(PotentiostatWindow)
        E_vert2_Eingabe.insert(END, -0.5)
        E_vert2_Eingabe.place(x = 150, y = 120)
        #-----------------------------------------------------------------------------------------------
        E_fin_Label = Label(PotentiostatWindow,text="E final vs. Ref [V]*")
        E_fin_Label.place(x = 25, y = 145)
        E_fin_Eingabe = Entry(PotentiostatWindow)
        E_fin_Eingabe.insert(END, -0.5)
        E_fin_Eingabe.place(x = 150, y = 145)
        #-----------------------------------------------------------------------------------------------
        Scanrate_Label = Label(PotentiostatWindow,text="Scanrate [mV/s]*")
        Scanrate_Label.place(x = 25, y = 170)
        Scanrate_Eingabe = Entry(PotentiostatWindow)
        Scanrate_Eingabe.insert(END, 20)
        Scanrate_Eingabe.place(x = 150, y = 170)
        #-----------------------------------------------------------------------------------------------
        Cycles_Label = Label(PotentiostatWindow,text="Num. of Cycles*")
        Cycles_Label.place(x = 25, y = 195)
        Cycles_Eingabe = Entry(PotentiostatWindow)
        Cycles_Eingabe.insert(END, 1)
        Cycles_Eingabe.place(x = 150, y = 195)
        #-----------------------------------------------------------------------------------------------
        Condi_Label = Label(PotentiostatWindow,text="Condit. time [s]*")
        Condi_Label.place(x = 25, y = 220)
        Condi_Eingabe = Entry(PotentiostatWindow)
        Condi_Eingabe.insert(END, 5)
        Condi_Eingabe.place(x = 150, y = 220)
        #-----------------------------------------------------------------------------------------------
        R_read_Label = Label(PotentiostatWindow,text="R read [Ohm]*")
        R_read_Label.place(x = 25, y = 245)
        R_read_Eingabe = Entry(PotentiostatWindow)
        R_read_Eingabe.insert(END, 120)
        R_read_Eingabe.place(x = 150, y = 245)
        #==============================================================================================
        # Define all buttons in the PolArStat CV-GUI
        #==============================================================================================
        CVSet_btn = ttk.Button(PotentiostatWindow, text="Send CV inputs", command=Set_CV_Params)   
        CVSet_btn.place(x = 25, y = 370, width = 120, height = 25)
        #-----------------------------------------------------------------------------------------------
        CVGet_btn = ttk.Button(PotentiostatWindow, text="Check CV inputs", command=Get_CV_Params)   
        CVGet_btn.place(x = 155, y = 370, width = 120, height = 25)
        #-----------------------------------------------------------------------------------------------
        CVRun_btn = ttk.Button(PotentiostatWindow, text="Run CV", command=AskSaveFile_AndStart)   
        CVRun_btn.place(x = 25, y = 405, width = 120, height = 50)
        #-----------------------------------------------------------------------------------------------
        CVStop_btn = ttk.Button(PotentiostatWindow, text="Stop CV", command=Stop_CV)   
        CVStop_btn.place(x = 155, y = 405, width = 120, height = 50)
        #-----------------------------------------------------------------------------------------------
        plot_button = ttk.Button(PotentiostatWindow, text = "Display E vs. I", command=Show_E_vs_I)
        plot_button.place(x = 325, y = 8, width = 140, height = 25)
        #-----------------------------------------------------------------------------------------------
        plot_button = ttk.Button(PotentiostatWindow, text = "Display t vs. I", command=Show_t_vs_I)
        plot_button.place(x = 475, y = 8, width = 140, height = 25)
        #-----------------------------------------------------------------------------------------------
        plot_button = ttk.Button(PotentiostatWindow, text = "Display t vs. E", command=Show_t_vs_E)
        plot_button.place(x = 625, y = 8, width = 140, height = 25)
        #==============================================================================================
        # Define all text-containing windows in the GUI
        #==============================================================================================
        Exper_Notes_Label = Label(PotentiostatWindow,text="Experimental Notes")
        Exper_Notes_Label.place(x = 25, y = 270)
        Exper_Notes_Eingabe = Entry(PotentiostatWindow) 
        Exper_Notes_Eingabe.place(x = 25, y = 295, width = 250, height = 60)
        #-----------------------------------------------------------------------------------------------
        Output_text_Label = Label(PotentiostatWindow,text="Monitor serial output data")
        Output_text_Label.place(x = 325, y = 350)
        Output_text = tk.scrolledtext.ScrolledText(PotentiostatWindow,  wrap = tk.WORD,  font = ("Times New Roman", 9)) 
        Output_text.place(x = 325, y = 370,  width = 450,  height = 85) 
        Output_text.focus() 
        #==============================================================================================
        # Initiate the Plot in the main GUI of the PolArStat CV-Script
        #==============================================================================================
        if Initiate_Plot == True:
            global plot1
            global canvas
            fig = Figure(figsize = (6.25,4), dpi = 70)
            x = 0  ;  y = 0
            plot1 = fig.add_subplot(111)
            if PlotType == 1:
                plot1.plot(x,y)
                plot1.set_xlabel("$E$ vs. RE in V", fontsize = 13)
                plot1.set_ylabel("$I$ in mA", fontsize = 13)
            if PlotType == 2:
                plot1.plot(x,y)
                plot1.set_xlabel("$t$ in s", fontsize = 13)
                plot1.set_ylabel("$I$ in mA", fontsize = 13)
            if PlotType == 3:
                plot1.plot(x,y)
                plot1.set_xlabel("$t$ in s", fontsize = 13)
                plot1.set_ylabel("$E$ vs. RE in V", fontsize = 13)
            plot1.tick_params(direction = 'in', length=4, width=0.5, colors='k', labelsize = 13)
            plot1.grid(which='both', linestyle = '--')
            canvas = FigureCanvasTkAgg(fig, PotentiostatWindow)  
            canvas.draw()
            canvas.get_tk_widget().place(x = 325, y = 50)
            Initiate_Plot = False
        #==============================================================================================
        # Run the Potentiostat window in a mainloop
        #==============================================================================================
        PotentiostatWindow.mainloop()
        #==============================================================================================
        # If mainloop is killed, kill port (here named arduino)
        #==============================================================================================
        arduino.close()   
    

    

        
        

    