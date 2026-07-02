# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:30:20 2023

@author: timtic
"""

########################################################################################################
########################################################################################################
#                    Importing required modules
########################################################################################################
########################################################################################################
import numpy as np
import struct
"""
########################################################################################################
# tkinter for building GUIs 
########################################################################################################"""
import tkinter as tk 
from   tkinter import ttk                           # Python 3
from   tkinter import *                             # Python 3
from   tkinter import messagebox                    # Python 3
from   tkinter import Menu                          # Python 3
from   tkinter import filedialog                    # Python 3
from   tkinter.scrolledtext  import ScrolledText    # Python 3
########################################################################################################
# matplotlib related modules for graphics
########################################################################################################
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
########################################################################################################
# Implement the default Matplotlib key bindings.
########################################################################################################
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
"""
########################################################################################################
# time for timed events, serial for serial
# communication of Arduino with the PC. The 
# serial requires the installation of pyserial
# on the operating PC.
########################################################################################################"""   
import time
import serial
import threading
import os
########################################################################################################
from datetime import datetime        # import datetome for checking recent cal file     
########################################################################################################
# Everything for threading the Data aqcuisition out of the GUI
########################################################################################################
from    threading                                import Thread
from    PolArStat_CA_Script_Acquisition          import CA_ACQUISITION_BACKEND  as CAAQB
from    PolArStat_CA_Script_Write_Outputs        import PolArStat_CA_Write_Read as PCAWR
from    PolArStat_Serial_Communication_Functions import Serial_Communication    as SERCO
from    PolArStat_Serial_Communication_Functions import *    



def PotentiostatScript_CA():
    ####################################################################################################
    #    Define plotting options to make the plot look cool :)
    ####################################################################################################
    font = {'family': 'Times New Roman', 'color':  'black','weight': 'normal','size': 15,}
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    plt.rcParams['font.sans-serif'] = ['Times new Roman']
    ####################################################################################################
    # Initialize some global variables
    ####################################################################################################
    global arduino       ;   arduino       = None      # initialize the device as None until it is defined
    global UpCounter     ;   UpCounter     = 0
    global CalData       ;   CalData       = 0
    global Zero_IDX_E    ;   Zero_IDX_E    = 13250
    global Zero_IDX_I    ;   Zero_IDX_I    = 13250
    global Initiate_Plot ;   Initiate_Plot = True      # For first setup of a plot 
    global thread_2
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
    # Change the variables in the class of the PolArStat_CA_Script_Write_Outputs script
    # according to the classical values or the values from the calibration file
    ####################################################################################################
    PCAWR.Zero_IDX_E = Zero_IDX_E
    PCAWR.Zero_IDX_I = Zero_IDX_I         
         
    ####################################################################################################
    # The following function asks for a filename, where the data will be stored
    ####################################################################################################
    def AskSaveFile_AndStart_CA():
        try:
            txt_filename = filedialog.asksaveasfile(defaultextension='.txt').name #öffnet fenster um nach speicherort zu fragen
        except Exception:
            messagebox.showwarning(title="No File selected!", message="You did not select a file to save to!")
            txt_filename = ""
        else:
            outfile_CA_Data   = open(txt_filename, 'w')
            PCAWR.Output_File = outfile_CA_Data
            START_CA()
         
    ######################################################################################################
    # The following functionwill find available Arduino-ports which are able for 
    # serial communication. This function wil be embedded into a GUI later on.
    ###################################################################################################### 
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
            except Exception as e:
                return ('error', str(e))
        def _on_serial_done(result):
            SSer_btn.config(state='normal')
            status = result[0]
            if status == 'ok':
                PCAWR.WRITE_TEXT_OUTPUT(Output_text, "Serial communication at Port %s established successfully" %portList.get())
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
    # The following function will set the inputs for a CA-measurement
    ###################################################################################################### 
    def Set_CA_Params():
        global ReadResist
        CASet_btn.config(state='disabled')
        try:
            E_1          = float(E_1_Eingabe.get())
            E_2          = float(E_2_Eingabe.get())
            E_3          = float(E_3_Eingabe.get())
            E_4          = float(E_4_Eingabe.get())
            E_5          = float(E_5_Eingabe.get())
            t_1          = (1e3)*float(t_1_Eingabe.get())
            t_2          = (1e3)*float(t_2_Eingabe.get())
            t_3          = (1e3)*float(t_3_Eingabe.get())
            t_4          = (1e3)*float(t_4_Eingabe.get())
            t_5          = (1e3)*float(t_5_Eingabe.get())
            Repetitions  = float(Repetitions_Eingabe.get())
            ReadResist   = float(R_read_Eingabe.get())
            try:
                Exper_Notes = str(Exper_Notes_Eingabe.get())
                if Exper_Notes != "":
                    PCAWR.Exper_Notes = Exper_Notes
                if Exper_Notes == "":
                    PCAWR.Exper_Notes = "No notes specified"
            except Exception:
                pass
            PCAWR.ReadResist = ReadResist
        except (ValueError, TypeError) as e:
            CASet_btn.config(state='normal')
            messagebox.showerror(title="Incomplete Parameters!",
                                 message="All inputs have to be filled out in order to run a CA!\n\nDetail: %s" % str(e))
            return
        if arduino is None or (hasattr(arduino, 'is_open') and not arduino.is_open):
            CASet_btn.config(state='normal')
            messagebox.showerror(title="Serial Communication failed!",
                                 message="No serial connection. Please Start Serial first!")
            return
        def _do_work():
            print('[Set_CA] Building packets...')
            SEND_E_1   = b'\x12' + b'\x17' + struct.pack('f', E_1) + b'\x00\x00'
            SEND_E_2   = b'\x12' + b'\x18' + struct.pack('f', E_2) + b'\x00\x00'
            SEND_E_3   = b'\x12' + b'\x19' + struct.pack('f', E_3) + b'\x00\x00'
            SEND_E_4   = b'\x12' + b'\x20' + struct.pack('f', E_4) + b'\x00\x00'
            SEND_E_5   = b'\x12' + b'\x21' + struct.pack('f', E_5) + b'\x00\x00'
            SEND_t_1   = b'\x12' + b'\x22' + struct.pack('f', t_1) + b'\x00\x00'
            SEND_t_2   = b'\x12' + b'\x23' + struct.pack('f', t_2) + b'\x00\x00'
            SEND_t_3   = b'\x12' + b'\x24' + struct.pack('f', t_3) + b'\x00\x00'
            SEND_t_4   = b'\x12' + b'\x25' + struct.pack('f', t_4) + b'\x00\x00'
            SEND_t_5   = b'\x12' + b'\x26' + struct.pack('f', t_5) + b'\x00\x00'
            SEND_REP   = b'\x12' + b'\x27' + struct.pack('f', Repetitions) + b'\x00\x00'
            print('[Set_CA] E_1=', E_1, 'E_2=', E_2, 't_1=', t_1, 't_2=', t_2, 'Reps=', Repetitions)
            pkts = [SEND_E_1, SEND_E_2, SEND_E_3, SEND_E_4, SEND_E_5,
                    SEND_t_1, SEND_t_2, SEND_t_3, SEND_t_4, SEND_t_5, SEND_REP]
            for i, pkt in enumerate(pkts):
                arduino.write(pkt)
                time.sleep(0.05)  # Small delay for Arduino to process each parameter
                print('[Set_CA] pkt', i, 'sent')
        def _on_done():
            CASet_btn.config(state='normal')
            PCAWR.WRITE_TEXT_OUTPUT(Output_text, "\n Data transmission complete :) \n---------------------------------------------")
        def _on_error(err):
            CASet_btn.config(state='normal')
            messagebox.showerror(title="Send Parameters Failed!",
                                 message="Failed to send CA parameters to device.\n\nDetail: %s" % str(err))
        def _run():
            try:
                _do_work()
                PotentiostatWindow.after(0, _on_done)
            except (serial.SerialException, OSError) as e:
                print("[CA Set_CA_Params Error]", e)
                PotentiostatWindow.after(0, _on_error, e)
            except Exception as e:
                print("[CA Set_CA_Params Error]", e)
                PotentiostatWindow.after(0, _on_error, e)
        threading.Thread(target=_run, daemon=True).start()
        
    ###################################################################################################### 
    # The following function will retreive the transmitted parameters for a CA-measurement
    ###################################################################################################### 
    def Get_CA_Params():
        CAGet_btn.config(state='disabled')
        def _do_work():
            print('[DEBUG] arduino:', arduino, 'is_open:', arduino.is_open if arduino else 'None')
            print('[DEBUG] DEVICE:', CAAQB.DEVICE, 'running:', CAAQB.running)
            CAAQB.DEVICE = arduino
            CAAQB.running = False
            print('[DEBUG] Sending query command...')
            Readback_Array     = np.zeros(12)
            GETTING_COMMAND_BYTES  = b'\x13\x00\x00\x00\x00\x00\x00\x00'
            arduino.reset_input_buffer()
            arduino.write(GETTING_COMMAND_BYTES)
            time.sleep(2)  # Wait for Arduino to process and send back parameters
            print('[DEBUG] Reading response...')
            for i in range(5):
                raw = arduino.readline()[:-2]
                if not raw:
                    raw = arduino.readline()[:-2]
                Readback_Array[i] = float(raw)
            for i in range(5, 10):
                raw = arduino.readline()[:-2]
                if not raw:
                    raw = arduino.readline()[:-2]
                Readback_Array[i] = (1e-3)*float(raw)
            raw = arduino.readline()[:-2]
            if not raw:
                raw = arduino.readline()[:-2]
            Readback_Array[10] = float(raw)
            Readback_Array[11] = ReadResist
            CAAQB.CV_Inputs     = Readback_Array
            return Readback_Array
        def _on_done(Readback_Array):
            CAGet_btn.config(state='normal')
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"The following Parameters were set! \n")
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_1 \t= \t %s V vs. Re"         %Readback_Array[0] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_2 \t= \t %s V vs. Re"         %Readback_Array[1] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_3 \t= \t %s V vs. Re"         %Readback_Array[2] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_4 \t= \t %s V vs. Re"         %Readback_Array[3] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_5 \t= \t %s V vs. Re"         %Readback_Array[4] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t1 \t= \t %s s"                 %Readback_Array[5] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t2 \t= \t %s s"                 %Readback_Array[6] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t3 \t= \t %s s"                 %Readback_Array[7] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t4 \t= \t %s s"                 %Readback_Array[8] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t5 \t= \t %s s"                 %Readback_Array[9] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"n-Reps \t = \t %s "             %Readback_Array[10])
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"R_read \t = \t %s Ohm"          %Readback_Array[11])
        def _on_error(err):
            CAGet_btn.config(state='normal')
            messagebox.showerror(title="Read Parameters Failed!",
                                 message="Failed to read CA parameters from device.\n\nDetail: %s" % str(err))
        def _run():
            try:
                result = _do_work()
                PotentiostatWindow.after(0, _on_done, result)
            except Exception as e:
                PotentiostatWindow.after(0, _on_error, e)
        threading.Thread(target=_run, daemon=True).start()
       

    ####################################################################################################
    # The following function will start a CA measurment
    #################################################################################################### 
    def START_CA():
        global UpCounter    ;   UpCounter= 0
        if arduino is None:
            messagebox.showwarning(title="No device in list!",
                                   message="No device selected. Initialize serial communication before running an experiment!")
            return
        CARun_btn.config(state='disabled')
        def _do_work():
            print('[START_CA] arduino:', arduino, 'is_open:', arduino.is_open if arduino else 'None')
            print('[START_CA] DEVICE:', CAAQB.DEVICE, 'running:', CAAQB.running)
            CAAQB.DEVICE = arduino
            CAAQB.running = False
            print('[START_CA] Sending query 0x13...')
            Readback_Array     = np.zeros(12)
            GETTING_COMMAND_BYTES  = b'\x13\x00\x00\x00\x00\x00\x00\x00'
            arduino.reset_input_buffer()
            arduino.write(GETTING_COMMAND_BYTES)
            time.sleep(2)  # Wait for Arduino to process and send back parameters
            print('[START_CA] Reading parameters...')
            for i in range(5):
                raw = arduino.readline()[:-2]
                print('[START_CA] E_' + str(i+1) + ' raw:', repr(raw))
                if not raw:
                    raw = arduino.readline()[:-2]
                    print('[START_CA] E_' + str(i+1) + ' retry:', repr(raw))
                Readback_Array[i] = float(raw)
            for i in range(5, 10):
                raw = arduino.readline()[:-2]
                print('[START_CA] t_' + str(i-4) + ' raw:', repr(raw))
                if not raw:
                    raw = arduino.readline()[:-2]
                Readback_Array[i] = (1e-3)*float(raw)
            raw = arduino.readline()[:-2]
            print('[START_CA] Reps raw:', repr(raw))
            if not raw:
                raw = arduino.readline()[:-2]
            Readback_Array[10] = float(raw)
            Readback_Array[11] = ReadResist
            CAAQB.CV_Inputs     = Readback_Array
            print('[START_CA] Params read OK. Sending start 0x14...')
            SEND_BYTES = b'\x14\x00\x00\x00\x00\x00\x00\x00'
            SERCO.Send_And_Wait_For_ACK(arduino, SEND_BYTES)
            print('[START_CA] Waiting for 10101010...')
            Starter = arduino.readline()[:-2]
            while Starter != b'10101010':
                time.sleep(0.01)
                Starter = arduino.readline()[:-2]
            print('[START_CA] Got 10101010! Experiment started.')
            return Readback_Array
        def _on_done(Readback_Array):
            CARun_btn.config(state='normal')
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"The following Parameters were set! \n")
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_1 \t= \t %s V vs. Re"         %Readback_Array[0] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_2 \t= \t %s V vs. Re"         %Readback_Array[1] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_3 \t= \t %s V vs. Re"         %Readback_Array[2] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_4 \t= \t %s V vs. Re"         %Readback_Array[3] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"E_5 \t= \t %s V vs. Re"         %Readback_Array[4] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t1 \t= \t %s s"                 %Readback_Array[5] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t2 \t= \t %s s"                 %Readback_Array[6] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t3 \t= \t %s s"                 %Readback_Array[7] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t4 \t= \t %s s"                 %Readback_Array[8] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"t5 \t= \t %s s"                 %Readback_Array[9] )
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"n-Reps \t = \t %s "             %Readback_Array[10])
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"R_read \t = \t %s Ohm"          %Readback_Array[11])
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"\n Preparing CA completed.\n\n")
            PCAWR.WRITE_TEXT_OUTPUT(Output_text,"Loop.No.\t Step.No.\t t/s\t E/V\t I/mA\n")
            CAAQB.CLEAR_STORAGE()
            CAAQB.DEVICE  = arduino
            CAAQB.running = True
            CAAQB.Looper_ON = True
            Serial_Looper(plot1, plot2, canvas)
        def _on_error(err):
            CARun_btn.config(state='normal')
            messagebox.showerror(title="Start CA Failed!",
                                 message="Failed to start CA experiment.\n\nDetail: %s" % str(err))
        def _run():
            try:
                result = _do_work()
                PotentiostatWindow.after(0, _on_done, result)
            except Exception as e:
                PotentiostatWindow.after(0, _on_error, e)
        threading.Thread(target=_run, daemon=True).start()
        
    ####################################################################################################
    # Function, to stop a CA measurement. The Portlist_Refresh() function is defined in the 
    # PolArStat_Serial_Communication_Functions script and re-imported above
    ####################################################################################################   
    def Stop_CA():
        if CAAQB.running == False:
            messagebox.showwarning(title="No active measurement!", 
                                   message="There is no measurement which can be stopped!")
        else:
            if CAAQB.running == True:
                PCAWR.WRITE_TEXT_OUTPUT(Output_text,"\n Status: Measurment was manually stopped by User.\n") 
            CAAQB.StopMeasure = True
            CAAQB.Looper_ON   = False
     
    ####################################################################################################
    # The following function will kill thread_1 and exit GUI cleanly.
    ####################################################################################################
    def _quit():
        CAAQB.ALIVE = False
        CAAQB.running = False
        CAAQB.Looper_ON = False
        thread_2.join(timeout = 3)
        PotentiostatWindow.quit()
        PotentiostatWindow.destroy()
        
    ####################################################################################################
    # Function, calling itself
    ####################################################################################################   
    def Serial_Looper(plot1, plot2, canvas):          # This function will be called in a loop back and forth with the main_loop of the GUI to update the output
        if CAAQB.Looper_ON == True:
            #===========================================================================================
            # Everything important for updating the plot
            #===========================================================================================
            InputArray = CAAQB.Storage_Array
            if len(InputArray[::,0]) < 3:
                PotentiostatWindow.after(100, Serial_Looper, plot1, plot2, canvas)
                return
            # Remove old data lines without resetting axis config
            while plot1.lines:
                plot1.lines[0].remove()
            while plot2.lines:
                plot2.lines[0].remove()
            if len(InputArray[::,0]) <= 10000:
                plot1.plot((1e-3)*InputArray[2::,2], -0.12452*(InputArray[2::,4]-Zero_IDX_I)/ReadResist, color = 'blue' )
                plot2.plot((1e-3)*InputArray[2::,2], -0.000249*(InputArray[2::,3]-Zero_IDX_E), color = 'red' )
            if len(InputArray[::,0]) > 10000:
                plot1.plot((1e-3)*InputArray[2:10:,2], -0.12452*(InputArray[2:10:,4]-Zero_IDX_I)/ReadResist, color = 'blue' )
                plot2.plot((1e-3)*InputArray[2:10:,2], -0.000249*(InputArray[2:10:,3]-Zero_IDX_E), color = 'red' )
            canvas.draw()
            #===========================================================================================
            # Everything important for updating the serial Output window
            #===========================================================================================
            global UpCounter
            PCAWR.WRITE_DATA_OUTPUT(WHERE_TO_WRITE = Output_text, DATA_TO_WRITE = InputArray[UpCounter+1::,::])
            UpCounter = len(InputArray[::,0])
            Output_text.see(tk.END)
            #===========================================================================================
            # Callback to looper function - only if the respective bool is True
            #===========================================================================================
            PotentiostatWindow.after(100, Serial_Looper, plot1, plot2, canvas)
        
    # Stop any leftover thread from a previous window session
    CAAQB.ALIVE = False
    CAAQB.running = False
    CAAQB.Looper_ON = False
    time.sleep(0.2)
    ####################################################################################################
    # Start the conditional data-acquisition on its own thread, that the GUI does not freeze
    # if data collection is busy
    ####################################################################################################
    thread_2 = Thread(target=CAAQB.LINE_READER_FOR_THREAD, daemon=True)
    thread_2.start()
    ####################################################################################################
    # Run the PolArStat_CA_Script in a mainloop
    ####################################################################################################
    if __name__ == 'PolArStat_CA_Script_GUI':
        ################################################################################################
        # The following part is the main loop which calls the GUI for interfacing with 
        #the Arduino. The Arduino has to contain the firmware before this script works.
        ################################################################################################
        #===============================================================================================
        # The following part is the main loop which calls the GUI for interfacing with 
        # the Arduino. The Arduino has to contain the firmware before this script works.
        # First, the size of the GUI-window will be defined. Subsequently, all input field are set
        #===============================================================================================
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
         #----------------------------------------------------------------------------------------------       
        Refresh_btn = ttk.Button(PotentiostatWindow, text="Refresh Portlist", command=SERCO.FindPorts)   
        Refresh_btn.place(x = 25, y = 35, width = 120, height = 25)
        #-----------------------------------------------------------------------------------------------
        SSer_btn = ttk.Button(PotentiostatWindow, text="Start Serial", command=StartSerial)   
        SSer_btn.place(x = 155, y = 35, width = 120, height = 25)
        #==============================================================================================
        # Set all the input buttons for CV measurement and place them in the GUI
        #==============================================================================================
        E_1_Label = Label(PotentiostatWindow,text="E1 vs. Ref [V]*")
        E_1_Label.place(x = 25, y = 70)
        E_1_Eingabe = Entry(PotentiostatWindow)
        E_1_Eingabe.insert(END, -0.5) 
        E_1_Eingabe.place(x = 110, y = 70, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        E_2_Label = Label(PotentiostatWindow,text="E2 vs. Ref [V]*")
        E_2_Label.place(x = 25, y = 95)
        E_2_Eingabe = Entry(PotentiostatWindow)
        E_2_Eingabe.insert(END, 0.5) 
        E_2_Eingabe.place(x = 110, y = 95, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        E_3_Label = Label(PotentiostatWindow,text="E3 vs. Ref [V]*")
        E_3_Label.place(x = 25, y = 120)
        E_3_Eingabe = Entry(PotentiostatWindow)
        E_3_Eingabe.insert(END, 0.0) 
        E_3_Eingabe.place(x = 110, y = 120, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        E_4_Label = Label(PotentiostatWindow,text="E4 vs. Ref [V]*")
        E_4_Label.place(x = 25, y = 145)
        E_4_Eingabe = Entry(PotentiostatWindow)
        E_4_Eingabe.insert(END, 0.0) 
        E_4_Eingabe.place(x = 110, y = 145, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        E_5_Label = Label(PotentiostatWindow,text="E5 vs. Ref [V]*")
        E_5_Label.place(x = 25, y = 170)
        E_5_Eingabe = Entry(PotentiostatWindow)
        E_5_Eingabe.insert(END, 0.0) 
        E_5_Eingabe.place(x = 110, y = 170, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        t_1_Label = Label(PotentiostatWindow,text="t1 [s]*")
        t_1_Label.place(x = 175, y = 70)
        t_1_Eingabe = Entry(PotentiostatWindow)
        t_1_Eingabe.insert(END, 10) 
        t_1_Eingabe.place(x = 222, y = 70, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        t_2_Label = Label(PotentiostatWindow,text="t2 [s]*")
        t_2_Label.place(x = 175, y = 95)
        t_2_Eingabe = Entry(PotentiostatWindow)
        t_2_Eingabe.insert(END, 10) 
        t_2_Eingabe.place(x = 222, y = 95, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        t_3_Label = Label(PotentiostatWindow,text="t3 [s]*")
        t_3_Label.place(x = 175, y = 120)
        t_3_Eingabe = Entry(PotentiostatWindow)
        t_3_Eingabe.insert(END, 0) 
        t_3_Eingabe.place(x = 222, y = 120, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        t_4_Label = Label(PotentiostatWindow,text="t4 [s]*")
        t_4_Label.place(x = 175, y = 145)
        t_4_Eingabe = Entry(PotentiostatWindow)
        t_4_Eingabe.insert(END, 0) 
        t_4_Eingabe.place(x = 222, y = 145, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        t_5_Label = Label(PotentiostatWindow,text="t5 [s]*")
        t_5_Label.place(x = 175, y = 170)
        t_5_Eingabe = Entry(PotentiostatWindow)
        t_5_Eingabe.insert(END, 0) 
        t_5_Eingabe.place(x = 222, y = 170, width = 50, height = 20)
        #-----------------------------------------------------------------------------------------------
        Repetitions_Label = Label(PotentiostatWindow,text="Repetitions*")
        Repetitions_Label.place(x = 25, y = 215)
        Repetitions_Eingabe = Entry(PotentiostatWindow)
        Repetitions_Eingabe.insert(END, 1)
        Repetitions_Eingabe.place(x = 150, y = 215)
        #-----------------------------------------------------------------------------------------------
        R_read_Label = Label(PotentiostatWindow,text="R read [Ohm]*")
        R_read_Label.place(x = 25, y = 240)
        R_read_Eingabe = Entry(PotentiostatWindow)
        R_read_Eingabe.insert(END, 120)
        R_read_Eingabe.place(x = 150, y = 240)
        #==============================================================================================
        # Define all buttons in the PolArStat CA-GUI
        #==============================================================================================
        CASet_btn = ttk.Button(PotentiostatWindow, text="Send CA inputs", command=Set_CA_Params)   
        CASet_btn.place(x = 25, y = 370, width = 120, height = 25)
        #-----------------------------------------------------------------------------------------------
        CAGet_btn = ttk.Button(PotentiostatWindow, text="Check CA inputs", command=Get_CA_Params)   
        CAGet_btn.place(x = 155, y = 370, width = 120, height = 25)
        #-----------------------------------------------------------------------------------------------
        CARun_btn = ttk.Button(PotentiostatWindow, text="Run CA", command=AskSaveFile_AndStart_CA)   
        CARun_btn.place(x = 25, y = 405, width = 120, height = 50)
        #-----------------------------------------------------------------------------------------------
        CAStop_btn = ttk.Button(PotentiostatWindow, text="Stop CA", command=Stop_CA)   
        CAStop_btn.place(x = 155, y = 405, width = 120, height = 50)
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
        # Initiate the Plot in the main GUI of the PolArStat CA-Script
        #==============================================================================================
        if Initiate_Plot == True:
            global plot1
            global plot2
            global canvas
            fig = Figure(figsize = (6.5,4.5), dpi = 68)
            plot1 = fig.add_subplot(111)
            x = np.array([0])
            y = np.array([0])
            plot1.plot((x,y), linewidth = 0)
            plot1.set_xlabel("$t$ in s", fontsize = 13)
            plot1.set_ylabel("$I$ in mA", fontsize = 13, color = 'blue', labelpad = 10)
            plot1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
            plot1.grid(which='both', linestyle = '--')
            plot1.tick_params(direction = 'in', length=4, width=0.5, colors='k', labelsize = 13)
            plot2 = plot1.twinx()
            plot2.set_ylabel("$E$ vs. RE in V", fontsize = 13, color = 'red', labelpad = 10)
            plot2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
            plot2.tick_params(direction = 'in', length=4, width=0.5, colors='k', labelsize = 13)
            fig.subplots_adjust(left = 0.15, right = 0.85)
            canvas = FigureCanvasTkAgg(fig, PotentiostatWindow)
            canvas.draw()
            canvas.get_tk_widget().place(x = 325, y = 10)
            Initiate_Plot = False
        #==============================================================================================
        # Run the Potentiostat window in a mainloop
        #==============================================================================================
        PotentiostatWindow.mainloop()
        #==============================================================================================
        # If mainloop is killed, kill port (here named arduino)
        #==============================================================================================
        arduino.close()   
    

        
        























    