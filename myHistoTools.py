#!/usr/bin/python


import math
import ROOT

# 26.5.2014
def GetBinContents(histo):
    nx = histo.GetXaxis().GetNbins()
    content = []
    for binx in range(1,nx+1):
        content.append(histo.GetBinContent(binx))
    return content


def GetBinEdges(histo):
    nx = histo.GetXaxis().GetNbins()
    edges = []
    for binx in range(1,nx+1):
        edges.append(histo.GetBinLowEdge(binx))
    edges.append( (histo.GetBinLowEdge(nx) + histo.GetBinWidth(nx)) )
    return edges



# 26.5.2014
def GetXYarray(histo):
    nx = histo.GetXaxis().GetNbins()
    content = []
    for binx in range(1,nx+1):
        content.append( [ histo.GetBinCenter(binx), histo.GetBinContent(binx) ])
    return content



def PrintBinContent(histo):
    nx = histo.GetXaxis().GetNbins()
    line=''
    for binx in range(0,nx+2):
        line = '%s, %4.7f' % (line, histo.GetBinContent(binx),)
    print line

def Flip2D(histo):
    flip = histo.Clone(histo.GetName() + "_flip")
    flip.Reset()
    nx = histo.GetXaxis().GetNbins()
    ny = histo.GetYaxis().GetNbins()
    for binx in range(0,nx+2):
        for biny in range(0,ny+2):
            flip.SetBinContent(binx,biny, histo.GetBinContent(nx - binx + 1,biny))
            flip.SetBinError(binx,biny, histo.GetBinError(nx - binx + 1,biny))
    flip.SetEntries(histo.GetEntries())
    return flip



minNsig = 15. # 15

def MakeSignificantDiff2D(h1,h2, doAbs = False, doDiff = False):
    histo = h1.Clone(h1.GetName() + "_sig")
    histo.Reset()
    nx = histo.GetXaxis().GetNbins()
    ny = histo.GetYaxis().GetNbins()
    for binx in range(0,nx+2):
        for biny in range(0,ny+2):
            val1 = h1.GetBinContent(binx,biny)
            val2 = h2.GetBinContent(binx,biny)
            if not doDiff:
                if (val1+val2 > 0. and val1 > minNsig and val2 > minNsig):
                    if doAbs:
                        histo.SetBinContent(binx,biny, math.fabs(val1-val2) / math.sqrt(val1+val2))
                    else:
                        histo.SetBinContent(binx,biny, (val1-val2) / math.sqrt(val1+val2))
            else:
                histo.SetBinContent(binx,biny, (val1-val2) )
    return histo



def Flip1D(histo):
    flip = histo.Clone(histo.GetName() + "_flip")
    flip.Reset()    
    nx = histo.GetXaxis().GetNbins()
    for binx in range(0,nx+2):
        flip.SetBinContent(binx, histo.GetBinContent(nx - binx + 1))
        flip.SetBinError(binx, histo.GetBinError(nx - binx + 1))
    flip.SetEntries(histo.GetEntries())
    #for binx in range(0,nx+2):
    #    print histo.GetBinContent(binx), flip.GetBinContent(binx)
    return flip

def MakeSignificantDiff1D(h1,h2, doAbs = False, doDiff = False):
    histo = h1.Clone(h1.GetName() + "_sig")
    histo.Reset()
    nx = histo.GetXaxis().GetNbins()
    for binx in range(0,nx+2):
            val1 = h1.GetBinContent(binx)
            val2 = h2.GetBinContent(binx)
            # print binx, val1, val2, math.sqrt(val1+val2)
            if not doDiff:
                if (val1+val2 > 0. and val1 > minNsig and val2 > minNsig):
                    sig = (val1-val2) / math.sqrt(val1+val2)
                    if doAbs:
                        histo.SetBinContent(binx, math.fabs(sig))
                    else:
                        histo.SetBinContent(binx, sig)
            else:
                histo.SetBinContent(binx, val1-val2)
    histo.Scale(1.)
    if not doDiff:
        histo.SetYTitle('Significance')
    else:
        histo.SetYTitle('Difference')

    if 0:
        for binx in range(1,nx+1):
            print '(%f,%f)=%f' % (histo.GetBinLowEdge(binx), histo.GetBinLowEdge(binx) + histo.GetBinWidth(binx), histo.GetBinContent(binx) )
    return histo
