#! /usr/bin/env python
import sys
sys.path.insert(0, sys.path[0]+'/source')
sys.path.insert(0, sys.path[0]+'/Objects')
sys.path.insert(0, sys.path[1]+'/Calculator')
sys.path.insert(0, sys.path[2]+'/Solver')
sys.path.insert(0, sys.path[3]+'/Models')
sys.path.insert(0, sys.path[4]+'/Models/AGN_Models')
sys.path.insert(0, sys.path[5]+'/Models/Surface_Brightness_Models')
sys.path.insert(0, sys.path[6]+'/Models/Event_Map_Models')
sys.path.insert(0, sys.path[7]+'/XTRA_pkg')
sys.path.insert(0, sys.path[8]+'/Utilities')


FLUX_MODE = False
if FLUX_MODE:
   import matplotlib
   matplotlib.use('Agg')

from XTRA_pkg import __logo__
from mainPipeline import makingXRayCatalogs,\
                         makingXRayRealization,\
                         makingEventMap

# print(__logo__)
print "XTRA (c) Arya Farahi" 
if len(sys.argv) < 1: 
  print "(1) HALOS MODE"
  print "(2) SB MAP MODE"
  print "(3) EVENT MAP MODE"

try: 
  ans = int(sys.argv[1]) # int(raw_input("Please enter the mode number : "))
except IndexError: 
  print "Usage: >> [MODE NUMBER] [INPUT]"
  raise SystemExit

if (ans == 1):
    makingXRayCatalogs()
if (ans == 2):
    makingXRayRealization()
if (ans == 3):
    makingEventMap()
else:
    pass


