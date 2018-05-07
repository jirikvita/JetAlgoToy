#!/usr/bin/python
# Jiri Kvita, Oct 30 2015
# jiri.kvita@gmail.com
# jet clustering algorithm simple and inefficient implementation

from myAll import *

Pi = math.pi
MinPhi = 0
MaxPhi = 2*Pi
SF = 1000
debug = 0

###############################################
# make empty eta phi histo
def MakeEmpty(hname = 'etaphi', htitle = 'etaphi;#eta;#phi', neta = 32, eta = 2.4, nphi = 32 ):
    hist = nextCan.nextH2(hname, htitle, neta, -eta, eta, nphi, MinPhi, MaxPhi)
    hist.SetStats(0)
    return hist

###############################################

class pitem:
    # massless initializer:
    def __init__(self, pt, eta, phi):
        if debug > 1: print 'MassLess eta,phi,pt ', eta, phi, pt
        self.rap = eta
        self.phi = phi
        # massless:
        self.E = pt*math.cosh(eta)
        self.pt = pt
        self.px = pt*math.cos(phi)
        self.py = pt*math.sin(phi)
        self.pz = pt*math.sinh(eta)
        self.p = self.E  # math.sqrt( pow(self.pt,2) + pow(self.pz,2) )
        self.m = 0. # math.sqrt( pow(self.E,2) + pow(self.p,2) )
        if debug > 1: print '  px,py,pz,E', self.px, self.py, self.pz, self.E
        if self.E < self.pz:
            print '    WARNING, E<pz!'
        # self.m = 0.

    # mass full initializer:
    def SetMassFull(self, E, pt, rap, phi):
        if debug > 1: print 'MassFull rap,phi,pt,E ', rap, phi, pt, E
        self.rap = rap
        self.phi = phi
        self.E = E
        self.pt = pt
        self.px = pt*math.cos(phi)
        self.py = pt*math.sin(phi)
        self.pz = E*math.tanh(rap)
        self.p = math.sqrt( pow(self.pt,2) + pow(self.pz,2) )
        self.m = math.sqrt( pow(self.E,2) + pow(self.p,2) )
        if debug > 1: print '  px,py,pz,E', self.px, self.py, self.pz, self.E
        if self.E < self.pz:
            print '   WARNING, E<pz!'
        # self.m = 

    def Print(self):
        print '  px,py,pz,E || rap,phi,pt', self.px, self.py, self.pz, self.E, ' || ', self.rap, self.phi, self.pt

###############################################
def Adjust(phi):
    if phi > MaxPhi:
        phi = MinPhi + phi - MaxPhi
    if phi < MinPhi:
        phi = MaxPhi - math.fabs(MinPhi - phi)
    return phi

###############################################
def combine(it1, it2):
    E = it1.E + it2.E
    px = it1.px + it2.px
    py = it1.py + it2.py
    if debug:
        print 'combining rap,phi,pt of [%f,%f,%f], [%f,%f,%f]' % (it1.rap,it1.phi,it1.pt,it2.rap,it2.phi,it2.pt,)
    pt = math.sqrt( math.pow(px,2) + math.pow(py,2)) 
    pz = it1.pz + it2.pz
    rap = 0.5*math.log( (E+pz) / (E-pz)) 
    phi = math.atan2( it1.py + it2.py, it1.px + it2.px)
    phi = Adjust(phi)
    if debug > 0: print '   ...combined rap, phi, E, pt,pz: ', rap, phi, E, pt, pz
    it = pitem(0., 0., 0.)
    it.SetMassFull(E, pt, rap, phi)
    return it

###############################################

# generate same jets in event:
def AddJetParticles(plist, hist, jeteta, jetphi, wEta, wPhi, N = 100):
    eta0 = hist.GetXaxis().GetXmin()
    eta1 = hist.GetXaxis().GetXmax()
    phi0 = hist.GetYaxis().GetXmin()
    phi1 = hist.GetYaxis().GetXmax()
    fun = ROOT.TF2('gen', '[0]*exp(-(x-[1])^2/(2*[2]^2))*exp(-(y-[3])^2/(2*[4]^3))', eta0, eta1, phi0 - 2*Pi, phi1 + 2*Pi)
    fun.SetParameters(1/(2*Pi*wEta*wPhi), jeteta, wEta, jetphi, wPhi)
    for i in xrange(0, N):
        eta = ROOT.Double(0) 
        phi = ROOT.Double(0) 
        fun.GetRandom2(eta,phi)
        # care for the cyclic phi!       
        phi = Adjust(phi)
        #ieta = hist.GetXaxis().FindBin(eta)
        #iphi = hist.GetYaxis().FindBin(phi)
        pt = SF*fun.Eval(eta, phi)
        hist.Fill(eta, phi)
        plist.append(pitem(pt, eta, phi))

###############################################

# generate same jets in event:
def AddPileUpParticles(plist, hist, N = 1000):
    eta0 = hist.GetXaxis().GetXmin()
    eta1 = hist.GetXaxis().GetXmax()
    phi0 = hist.GetYaxis().GetXmin()
    phi1 = hist.GetYaxis().GetXmax()
    fun = ROOT.TF2('gen', '[0]', eta0, eta1, phi0 - 2*Pi, phi1 + 2*Pi)
    fun.SetParameter(0, 0.05)
    for i in xrange(0, N):
        eta = ROOT.Double(0) 
        phi = ROOT.Double(0) 
        fun.GetRandom2(eta,phi)
        # care for the cyclic phi!       
        phi = Adjust(phi)
        #ieta = hist.GetXaxis().FindBin(eta)
        #iphi = hist.GetYaxis().FindBin(phi)
        pt = SF*fun.Eval(eta, phi)
        hist.Fill(eta, phi)
        plist.append(pitem(pt, eta, phi))


###############################################
def DeltaPhi(phi1, phi2):
    dphi = math.fabs(phi1 - phi2)
    if dphi > 2*Pi:
        dphi = 2*Pi - dphi 
    if debug > 2: print 'dphi = %f ' % (dphi,)
    return dphi

###############################################
def DeltaR2(i1, i2):
    return  pow(DeltaPhi(i1.phi,i2.phi), 2) + pow(i1.rap - i2.rap, 2)


###############################################
def MakeDistances(pList, R, p):
    distances = []
    for i in range(0, len(pList)):
        for j in range(0, len(pList)):
            distance = -1.
            if i != j:
                minpt = min( pow(pList[i].pt, 2*p), pow(pList[j].pt, 2*p) )
                distance = DeltaR2(pList[i], pList[j])*minpt
            else:
                distance = pow(pList[i].pt, 2*p) * pow(R, 2)
            distances.append( [i, j, distance] )
    return distances
            

###############################################
def GetMin(distances):
    imin = -1
    dmin = 1.e32
    i = 0
    for d in distances:
        if d[2] < dmin:
            dmin = d[2]
            imin = i
        i = i + 1
    return imin

###############################################
def RunJetAlgo(pList, R = 0.5, p = -1.):

    Jets = []

    ipass = -1
    while len(pList) > 0:
        ipass = ipass+1
        if ipass % 10 == 0: print '=== iteration %i ===' % (ipass,)
        if debug > 1:
            for particle in pList:
                particle.Print()

    # index of min distance
        distances = MakeDistances(pList, R, p)
        imind = GetMin(distances) 
        mind = distances[imind][2]

        # is the min distance a self-distance?
        # compare indices of the objects, w.r.t. the pList:
        if distances[imind][0] == distances[imind][1]:
            # we have a jet!
            Jets.append(pList[distances[imind][0]])
            pList.pop(distances[imind][0])
        else:
            # combine the particles:
            newitem = combine(pList[distances[imind][0]], pList[distances[imind][1]])
            indices = [distances[imind][0], distances[imind][1] ]
            for index in sorted(indices, reverse=True):
                if debug: print 'deleting item ', index
                del pList[index]
            pList.append(newitem)
    
    return Jets


###############################################
###############################################
###############################################



# STEERING

addPU = True

# the R parameter of the jet algorithm
Rs = [0.6]
#Rs = [ 0.3 + 0.1*x for x in range(0, 10) ]

# the p parameter of the jet algorithm
# -1 ... anti-kt
#  0 ... cone
#  1 ... kt
#ps = [-1.]
ps = [1, 0, -1.]
ps = [1, 0]


for p in ps:
    print '--- processig p=%1.1f ---' % (p,)
    for R in Rs:
        print '--- processig R=%1.1f ---' % (R,)
        # generate same jets in event:
        hist = MakeEmpty()
        ROOT.gRandom.SetSeed(1)

        pList = []
        # still crashes on cyclic phi issues...
        #AddJetParticles(pList, hist, 0, 0, 0.3, 0.35)

        AddJetParticles(pList, hist, 1, 1.5, 0.3, 0.35)
        AddJetParticles(pList, hist, -1.2, 4.3, 0.3, 0.35)
        AddJetParticles(pList, hist, 1.5, 5.1, 0.3, 0.35)
        if addPU:
            # Pile-up!
            AddPileUpParticles(pList, hist, 300)
        #print pList


        ###################################
        # Draw bare:
        canname = 'GeneratedParticles'
        if addPU: canname = canname + '_PU'
        can = nextCan.nextTCanvas(canname, canname, 100, 100, 1000, 800)
        hcopy = hist.DrawCopy('lego2')
        #can.Print(canname + '.eps')
        can.Print(canname + '.png')
        hcopy = hcopy.Draw('colz')
        #can.Print(canname + '_colz.eps')
        can.Print(canname + '_colz.png')


        # run jet algorithm on event

        Jets = RunJetAlgo(pList, R, p)
        # keep associated particles in jet somehow?


        # plots jets constituents with different color, mark jet axis and draw a circle



        ###################################
        # Draw jets:
        canname = 'Jets_R%1.1f_p%1.1f' % (R,p,)
        if addPU: canname = canname + '_PU'
        can = nextCan.nextTCanvas(canname, canname, 0, 0, 1000, 800)


        hist.Draw('colz')
        print '============================================='
        print '                  RESULTS                    '
        print '============================================='
        print '===> %i jets found:' % (len(Jets),)

        cirs = []
        for jet in Jets:
            jet.Print()
            # draw a circle:
            # TODO
            circ = ROOT.TEllipse(jet.rap, jet.phi, R, R)
            delta = ROOT.Double(0)
            circ.SetLineStyle(3)
            if jet.pt > 10.:
                circ.SetLineWidth(2)
                circ.SetLineStyle(2)
            if jet.pt > 1000.:
                circ.SetLineWidth(3)
                circ.SetLineStyle(1)
                #delta = max(0,math.log(jet.pt))
                #circ.SetLineWidth(1.+delta/10.)
            circ.SetFillStyle(-1)
            circ.Draw()
            # one more circle to account for cyclic phi?

            cirs.append(circ)


        ROOT.gPad.SetGridx() ; ROOT.gPad.SetGridy()


        #can.Print(canname + '.eps')
        can.Print(canname + '.png')


ROOT.gApplication.Run()

