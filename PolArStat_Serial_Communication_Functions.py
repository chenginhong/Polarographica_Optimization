# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 17:30:20 2023

@author: timtic
"""

import serial
import sys
import glob
import time

class Serial_Communication():
    
    portsl = None
    ACK_TIMEOUT = 0.2  # seconds to wait for ACK
    
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
        Serial_Communication.portsl = result

    ###############################################################################
    # Send a command and wait for ACK byte from Arduino.
    # If ACK is received, return immediately (fast path).
    # If ACK times out, fall back to a short sleep (legacy firmware compatibility).
    ###############################################################################
    def Send_And_Wait_For_ACK(arduino, cmd_bytes, ack_byte=None, timeout=None):
        """
        Send cmd_bytes to arduino and wait for ack_byte echo.
        If ack_byte is None, uses the first byte of cmd_bytes as expected ACK.
        Returns True if ACK received, False if timed out.
        """
        if ack_byte is None:
            ack_byte = cmd_bytes[0:1]
        if timeout is None:
            timeout = Serial_Communication.ACK_TIMEOUT
        # Flush stale input
        arduino.reset_input_buffer()
        # Send command
        arduino.write(cmd_bytes)
        # Wait for ACK — Arduino replies within microseconds
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if arduino.in_waiting > 0:
                byte = arduino.read(1)
                if byte == ack_byte:
                    return True
            else:
                time.sleep(0.001)
        # ACK not received within timeout — proceed anyway
        return False

