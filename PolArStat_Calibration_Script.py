# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:30:37 2022

@author: gisel
"""


"""
Created on Mon Jul 04 13:15:00 2022

@author: timtic
"""

###############################################################################
###############################################################################
#                    Importing required modules
###############################################################################
###############################################################################
import numpy as np
import struct
import sys
import glob
import os
"""
#------------------------------------------------------------------------------
# tkinter for building GUIs 
#------------------------------------------------------------------------------"""
import tkinter as tk 
from   tkinter import ttk                           # Python 3
from   tkinter import *                             # Python 3
from   tkinter import messagebox                    # Python 3
from   tkinter import Menu                          # Python 3
from   tkinter import filedialog                    # Python 3
from   tkinter.scrolledtext  import ScrolledText    # Python 3
from   PIL     import ImageTk, Image
#------------------------------------------------------------------------------
# matplotlib related modules for graphics
#------------------------------------------------------------------------------   
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
#------------------------------------------------------------------------------
# Implement the default Matplotlib key bindings.
#------------------------------------------------------------------------------
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
"""
#------------------------------------------------------------------------------
# time for timed events, serial for serial
# communication of Arduino with the PC. The 
# serial requires the installation of pyserial
# on the operating PC.
#------------------------------------------------------------------------------"""   
import time
import serial
import threading
#------------------------------------------------------------------------------
from datetime import datetime        # import datetome for generating recent cal file     
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------     



def PotentiostatScript_Calib():
    global running
    running       = False              # boolean for running the script deactivated
    Initiate_Plot = True               # For first setup of a plot
    Refresh_Count = 0                  # Count to a certain number before refreshing plot
    global PlotType
    PlotType      = 1
    global Data_Collected
    Data_Collected = False
    ###############################################################################
    #    Define plotting options to make the plot look cool :)
    ###############################################################################
    font = {'family': 'Times New Roman', 'color':  'black','weight': 'normal','size': 15,}
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    plt.rcParams['font.sans-serif'] = ['Times new Roman']
    ###############################################################################
    ###############################################################################
    ###############################################################################
    # The following function asks for a filename, where the data will be stored
    ###############################################################################
    def Create_Cal_File():
        cal_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CALIBRATION')
        os.makedirs(cal_dir, exist_ok=True)
        txt_filename = os.path.join(cal_dir, datetime.today().strftime('%Y_%m_%d_Cal.txt'))
        global outfile_Calib_Data
        outfile_Calib_Data = open(txt_filename, 'w')
    
    
    ###############################################################################
    # The following function  will find available Arduino-ports which are able for 
    # serial communication. This function wil be embedded into a GUI later on.
    ###############################################################################
    def FindPorts():
        if sys.platform.startswith('win'):
            ports = ['COM%s' % (i+1) for i in range(256)]
        elif sys.platform.startswith('linux'):
            ports = glob.glob('/dev/ttyUSB*') + glob.glob('/dev/ttyACM*')
        elif sys.platform.startswith('darwin'):
            ports = glob.glob('/dev/tty.usbserial*') + glob.glob('/dev/tty.usbmodem*')
        else:
            ports = []
        result = []
        for port in ports:
            try:
                s = serial.Serial(port)
                s.close()
                result.append(port)
            except (OSError, serial.SerialException):
                pass
        return result
    ###############################################################################
    # The following function is used to refresh the portlist (function callback)
    ###############################################################################
    def Portlist_Refresh():
        global portsl
        portsl = FindPorts()
    ###############################################################################
    # The following functionwill find available Arduino-ports which are able for 
    # serial communication. This function wil be embedded into a GUI later on.
    ###############################################################################
    
    def StartSerial():
        global arduino
        SSer_btn.config(state='disabled')
        def _do_serial():
            global arduino
            try:
                arduino = serial.Serial(portList.get(), baudrate = 115200, timeout=.1)
                time.sleep(2)
                arduino.reset_input_buffer()  # Flush any bootloader garbage
                randfloat  = 11.01
                SEND_BYTES = b'\x44\x66' + struct.pack('f', randfloat)
                arduino.write(SEND_BYTES)
                time.sleep(2)
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
                Write_Output_Calib("Serial communication at Port %s established successfully" %portList.get())
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
    ###############################################################################
    ###############################################################################
    
    
    def Set_Calib_Params():
        global ReadResist
        try:
            E_1          = -0.5        # Potential of first step of calibration
            E_2          =  0.5        # Potential of second step of calibration
            E_3          =  0.0        # Reqired to include to mis-use the CA script in the Arduino
            E_4          =  0.0        # Reqired to include to mis-use the CA script in the Arduino
            E_5          =  0.0        # Reqired to include to mis-use the CA script in the Arduino
            t_1          =  5000.0     # time of step in ms (CA script convention)
            t_2          =  5000.0     # time of step in ms (CA script convention)
            t_3          =  0.0        # Reqired to include to mis-use the CA script in the Arduino
            t_4          =  0.0        # Reqired to include to mis-use the CA script in the Arduino
            t_5          =  0.0        # Reqired to include to mis-use the CA script in the Arduino
            Repetitions  =  1.0        # Reqired to include to mis-use the CA script in the Arduino
            ReadResist   =  float(R_read_Eingabe.get())
            #---------------------------------------------------
            # Prepare sending the Info to Arduino
            #---------------------------------------------------
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
            #---------------------------------------------------
            # Send packets to Arduino (8 bytes each, 50ms delay)
            #---------------------------------------------------
            pkts = [SEND_E_1, SEND_E_2, SEND_E_3, SEND_E_4, SEND_E_5,
                        SEND_t_1, SEND_t_2, SEND_t_3, SEND_t_4, SEND_t_5, SEND_REP]
            for i, pkt in enumerate(pkts):
                arduino.write(pkt)
                time.sleep(0.05)

            #-----------------------------------------------------
            # In case of successful transmission, print statement   
            #-----------------------------------------------------
            Write_Output_Calib("\n Data transmission complete :) \n---------------------------------------------")
        
        except Exception:
            messagebox.showerror(title="Incomplete Parameters!", 
                                 message="All inputs have to be filled out in order to run a CV!")
        
    ###############################################################################
    ###############################################################################
    
    def StartRead_Calib():
        global running
        global file_error
        global Storage_Array
        global UpCounter
        global Data_Collected
        Create_Cal_File()
        messagebox.showwarning(title="Calibration takes time!", message="You have successfully started a calibration!\nThis takes some time so be patient. Approx 20 s for data check\nand 10 s for measuring. Green LED lights up and indicates calibration in progress.\nGreen LED flashes 3 times if calibation is complete.\n\nHit okay to start.")
        CalRun_btn.config(state='disabled')
        def _do_work():
            Set_Calib_Params()
            Storage_Array_local = np.zeros((1,5))
            SEND_BYTES = b'\x14\x00\x00\x00\x00\x00\x00\x00'
            arduino.write(SEND_BYTES)
            # Read raw bytes until 999999 arrives
            raw = b''
            while b'999999' not in raw:
                n = arduino.in_waiting
                if n > 0:
                    raw += arduino.read(n)
                else:
                    time.sleep(0.05)
            # Parse raw bytes
            text = raw.decode('utf-8', errors='ignore')
            # Find data between 10101010 and 999999
            start = text.find('10101010')
            end = text.find('999999')
            if start >= 0 and end > start:
                data_section = text[start+8:end]
                for line in data_section.split('\r\n'):
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) == 5:
                        try:
                            row = [float(p) for p in parts]
                            Storage_Array_local = np.vstack([Storage_Array_local, row])
                        except ValueError:
                            pass
            return Storage_Array_local
        def _on_work_done(storage):
            global running, file_error, Data_Collected, Storage_Array
            CalRun_btn.config(state='normal')
            Storage_Array  = storage
            file_error     = True
            Data_Collected = True
            # Write calibration data to file
            for i in range(len(storage[::,0])-1):
                outfile_Calib_Data.write(str(storage[i+1,0]))
                outfile_Calib_Data.write("\t")
                outfile_Calib_Data.write(str(storage[i+1,1]))
                outfile_Calib_Data.write("\t")
                outfile_Calib_Data.write(str(0.001*storage[i+1,2]))
                outfile_Calib_Data.write("\t")
                outfile_Calib_Data.write(str(storage[i+1,3]))
                outfile_Calib_Data.write("\t")
                outfile_Calib_Data.write(str(storage[i+1,4]))
                outfile_Calib_Data.write("\n")
            outfile_Calib_Data.close()
            Write_Output_Calib('\n Calibration completed successfully! File saved.')
            running = False
        def _on_work_error(err):
            CalRun_btn.config(state='normal')
            messagebox.showerror(title="Incomplete Parameters!",
                                 message="All inputs have to be filled out in order to run a CV!")
        def _run():
            try:
                result = _do_work()
                PotentiostatWindow.after(0, _on_work_done, result)
            except Exception as e:
                PotentiostatWindow.after(0, _on_work_error, e)
        threading.Thread(target=_run, daemon=True).start()
        
        
        
        
    
    ###############################################################################
    ###############################################################################
      
    def StopReadSuccess_Calib():
        global running
        global file_error
        file_error = False
        running = False
        #==================================================================
        #   Write output in the output-txt file
        #==================================================================
        for i in range(len(Storage_Array[::,0])-1):
            outfile_Calib_Data.write(str(Storage_Array[i+1,0]))
            outfile_Calib_Data.write("\t")
            outfile_Calib_Data.write(str(Storage_Array[i+1,1]))
            outfile_Calib_Data.write("\t")
            outfile_Calib_Data.write(str(0.001*Storage_Array[i+1,2]))
            outfile_Calib_Data.write("\t")
            outfile_Calib_Data.write(str(Storage_Array[i+1,3]))
            outfile_Calib_Data.write("\t")
            outfile_Calib_Data.write(str(Storage_Array[i+1,4]))
            outfile_Calib_Data.write("\n")
        outfile_Calib_Data.close()
        #==================================================================
        # write in output window that the measurement has terminated 
        #==================================================================  
        Write_Output_Calib("\n Status: Calibration terminated successfully.\n") 
        try:
            arduino.close() 
            Write_Output("\n Status: Serial port closed...\n") 
        except Exception:
            pass
        
        
    
    ###############################################################################
    # The following function will exit GUI cleanly.
    ###############################################################################
    def _quit():
        PotentiostatWindow.quit()
        PotentiostatWindow.destroy()
    ###############################################################################
    ############################################################################### 
    def Write_Output_Calib(Text): 
        Output_text.insert(INSERT, "%s \n" %Text )

    ###############################################################################
    #    Marius-functions reduced
    ###############################################################################   
    def Serial_Looper():          # This function will be called in a loop back and forth with the main_loop of the GUI to update the output
        global running
        global file_error
        global data 
        global Storage_Array
        global UpCounter
        if running == True:
            try:
                data = arduino.readline()[:-2]
                if data:
                    if data != b'999999':   # As soon as the Arduino sends this, the measurement is done
                        DECODED       = np.array(data.decode("utf-8").split('\t'))
                        DECODED_FLOAT = DECODED.astype(float)
                        Storage_Array = np.vstack([Storage_Array, DECODED_FLOAT])
                    if data == b'999999':
                        StopReadSuccess_Calib()
            except Exception:
                messagebox.showerror(title="Communication NOT started or broken!", message="It seems like communication was not started\nor is broken. Please Start Serial Communication first!")
                running = False
        PotentiostatWindow.after(1, Serial_Looper)
        
    
  
    
    if __name__ == 'PolArStat_Calibration_Script':
        
        ###############################################################################
        # The following part is the main loop which calls the GUI for interfacing with 
        #the Arduino. The Arduino has to contain the firmware before this script works.
        ###############################################################################
    
        PotentiostatWindow = tk.Tk()
        PotentiostatWindow.title("Main Interface for Calibrating the PolArStat")
        PotentiostatWindow.geometry('325x300')
        
        
        menu = Menu(PotentiostatWindow)                                               
        PotentiostatWindow.config(menu=menu)
        #=================================================================================
        #=================================================================================
           
        labelTop = tk.Label(PotentiostatWindow, text = "Choose COM Port")
        labelTop.place(x = 25, y = 10)
    
        portsl = FindPorts()
    
        portList = ttk.Combobox(PotentiostatWindow, postcommand=lambda: portList.configure(values=portsl))
        portList.place(x = 140, y = 10, width = 135, height = 19)
        
        Refresh_btn = ttk.Button(PotentiostatWindow, text="Refresh Portlist", command=Portlist_Refresh)   
        Refresh_btn.place(x = 25, y = 35, width = 120, height = 25)
        
        SSer_btn = ttk.Button(PotentiostatWindow, text="Start Serial", command=StartSerial)   
        SSer_btn.place(x = 155, y = 35, width = 120, height = 25)
        
        R_read_Label = Label(PotentiostatWindow,text="R read [Ohm]*")
        R_read_Label.place(x = 25, y = 75)
        R_read_Eingabe = Entry(PotentiostatWindow)
        R_read_Eingabe.insert(END, 120)
        R_read_Eingabe.place(x = 150, y = 75)
        
        Output_text_Label = Label(PotentiostatWindow,text="Monitor serial output data")
        Output_text_Label.place(x = 25, y = 125)
        Output_text = tk.scrolledtext.ScrolledText(PotentiostatWindow,  wrap = tk.WORD,  font = ("Times New Roman", 9)) 
        Output_text.place(x = 25, y = 150,  width = 250,  height = 85) 
        # Placing cursor in the text area 
        Output_text.focus() 
        
        CalRun_btn = ttk.Button(PotentiostatWindow, text="Run Calibration", command=StartRead_Calib)   
        CalRun_btn.place(x = 25, y = 250, width = 250, height = 50)
        
        
        
        PotentiostatWindow.after(50, Serial_Looper)     # this command pulls back and forth with the Serial_Looper function
        PotentiostatWindow.mainloop()
        arduino.close()   
    

        
        
