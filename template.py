#!/usr/bin/python

# jiri kvita, 19.5.2014

from myAll import *


filename = ''
rfile = ROOT.TFile(filename, 'read')

can = nextCan.nextTCanvas()
can.Divide(2,2)
can.cd(1)
ROOT.gPad.SetGridx() ; ROOT.gPad.SetGridy()

# ROOT.gApplication.Run()

