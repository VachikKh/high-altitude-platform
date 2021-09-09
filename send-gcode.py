#!/usr/bin/env python3

import serial
import sys
import time
s = serial.Serial('/dev/ttyUSB0', baudrate=115200, bytesize=8, parity=serial.PARITY_NONE, stopbits=1)

time.sleep(1)
m = b''
while s.in_waiting != 0:
    m = m + s.read()

print(m.decode('ascii'))

for line in sys.stdin.readlines():
    print('>', line)
    s.write(line.encode('ascii'))
    a = s.readline()
    print('<', a.decode('ascii'))

s.close()
