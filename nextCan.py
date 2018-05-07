import commands
import os
import tempfile
import ROOT

# code by Jan Veverka
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/JanVeverka/FWLite/Tools/python/canvases.py?view=markup

canvases = []
legends = []
graphs = []
hists1 = []
hists2 = []
anything = []
tf1s = []
rfiles = []

xperiod = 30
yperiod = 5
wheight = 500
wwidth = 1000



def nextTFile(name, opt = 'read'):
    rfile = ROOT.TFile(name, opt)
    rfiles.append(rfile)
    return rfile

def nextTF1(name = 'tmp', title = 'tmp', x1=0, x2=10):
    tf1 = ROOT.TF1(name, title, x1, x2)
    tf1s.append(tf1)
    return tf1

def nextH1(name = 'tmp', title = 'tmp', n=10, x1=0, x2=10):
    hist = ROOT.TH1D(name, title, n, x1, x2)
    hists1.append(hist)
    return hist

def nextH2(name = 'tmp', title = 'tmp', nx=10, x1=0, x2=10, ny = 10, y1=0, y2=10):
    hist = ROOT.TH2D(name, title, nx, x1, x2, ny, y1, y2)
    hists2.append(hist)
    return hist


nextTH1F = nextH1
nextTH1D = nextH1
nextTH2F = nextH2
nextTH2D = nextH2

def nextGraph(marker = 20,color = 1):
    gr = ROOT.TGraph()
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(1)
    gr.SetMarkerColor(color)
    gr.SetLineColor(color)
    gr.SetLineStyle(1)
    graphs.append(gr)
    return gr

def nextGraphErrors(marker = 20,color = 1):
    gr = ROOT.TGraphErrors()
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(1)
    gr.SetMarkerColor(color)
    gr.SetLineColor(color)
    gr.SetLineStyle(1)
    graphs.append(gr)
    return gr

nextTGraphErrors = nextGraphErrors
nextTGraph = nextGraph

def nextLeg(lx1=0.65,ly1=0.65,lx2=0.87,ly2=0.87):
    leg1 = ROOT.TLegend(lx1,ly1,lx2,ly2)
    leg1.SetBorderSize(0)
    leg1.SetFillColor(10)
    legends.append(leg1)
    return leg1

nextTLegend = nextLeg

#______________________________________________________________________________
def update():
    for c in canvases:
        if c:
            c.Update()
    ## end of loop over canvases
## end of update()

#______________________________________________________________________________

def nextCan(name=None, title=None, xoff = 0, yoff = 0, w=wwidth, h=wheight):
    update()
    i = len(ROOT.gROOT.GetListOfCanvases())
    wtopx = 20 * (i % xperiod)
    wtopy = 20 * (i % yperiod)
    
    if not title:
        title = name
        
    if name:
        if ROOT.gROOT.GetListOfCanvases().FindObject(name):
            i = 0
            while ROOT.gROOT.GetListOfCanvases().FindObject(name + '_%d' % i):
                i += 1
            name = name + '_%d' % i
            if title:
                title = title + ' %d' % i
        c1 = ROOT.TCanvas(name, title)
    else:
        c1 = ROOT.TCanvas()
        if title:
            c1.SetTitle(title)

    c1.SetWindowPosition(wtopx+xoff, wtopy+yoff)
    c1.SetWindowSize(w, h)

    canvases.append(c1)
    return c1
## end of next()

next = nextCan
nextTCanvas = nextCan
