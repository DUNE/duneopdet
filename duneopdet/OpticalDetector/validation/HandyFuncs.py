from ROOT import *
from array import array
import os, math
import progressbar as pb

def PrintError(val, err):
    precision = max(int(-math.log10(err)), 0)+1
    form = "{0:.%if} #pm {1:.%if}" % (precision, precision)
    return form.format(val, err)


def printlist(mylist):
    rows, columns = os.popen('stty size', 'r').read().split()
    maxw = max(map(len,mylist))
    ncols = int(columns)/(maxw+1)

    col = 0
    for word in mylist:
        print word.ljust(maxw), " ",
        col += 1
        if col >= ncols:
            print ""
            col = 0
    print ""


def GetHists(c1):
    prims = list(c1.GetListOfPrimitives())
    return filter(lambda x: x.InheritsFrom("TH1"), prims)        


def VerticalRange(hists, xrange=(0,0), ratio=False, forceOne=True, ignoreError=False, maxerr=0.25, absMax=-999, absMin=-999, buffer=0.05):
    if xrange == (0,0):
        minbin = 1
        maxbin = hists[0].GetNbinsX()
    else:
        minbin = hists[0].FindBin(xrange[0])
        maxbin = hists[0].FindBin(xrange[1])
    
    min = 999999.
    max = -999999.

    for hist in hists:
        for b in range(minbin, maxbin+1):
            denom = hist.GetBinContent(b)
            if denom == 0: denom = 0.00001
            errratio = abs(hist.GetBinError(b) / denom)

            if ignoreError:
                lowval = hist.GetBinContent(b)
                highval = hist.GetBinContent(b)
            else:
                lowval = hist.GetBinContent(b) - hist.GetBinError(b)
                highval = hist.GetBinContent(b) + hist.GetBinError(b)
            #print "Bin = %i, Cont = %f, Low = %f, High = %f, ER = %f" %(b, hist.GetBinContent(b), lowval, highval, errratio)
        
        
            if (errratio < maxerr) or ignoreError:
                if lowval < min and not (hist.GetBinContent(b) == 0 and ratio):
                    min = lowval
                if highval > max:
                    max = highval
                    
    if ratio:
        if forceOne:
            if min > 1: min = 1
            if max < 1: max = 1
        adjust = (max - min)*buffer/2.
        max = max + adjust
        min = min - adjust
    else:
        min = 0
        max = (1+buffer)*max

    if absMax != -999:
        if max > absMax:
            max = absMax
    if absMin != -999:
        if min < absMin:
            min = absMin
                    

    #print "Min = %f, Max = %f" %(min, max)
    for hist in hists:
        hist.SetMinimum(min)
        hist.SetMaximum(max)


def HorizontalRange(hists, mustinclude=[]):
    low = 1e90
    high = -1e90

    for hist in hists:
        mean = hist.GetMean()
        rms  = hist.GetRMS()
        low  = min(low, mean-2*rms)
        high = max(high, mean+2*rms)

    for m in mustinclude:
        low  = min(low, m)
        high = max(high, m)


    adjust = (high-low)*0.05
    low  -= adjust
    high += adjust

    for hist in hists:
        hist.GetXaxis().SetRangeUser(low,high)


def ChisqCount(h1, h2):
    sum1 = 0
    sum2 = 0
    w1 = 0
    w2 = 0

    for bin in range(1, h1.GetNbinsX()+1):
        sum1 += h1.GetBinContent(bin)
        sum2 += h2.GetBinContent(bin)
        ew1   = h1.GetBinError(bin)
        ew2   = h2.GetBinError(bin)
        w1   += ew1*ew1
        w2   += ew2*ew2

    # My simple chisquared method
    delta = sum1 - sum2
    sigmasq = w1 + w2
    chi2 = delta*delta / sigmasq
    prob = TMath.Prob(chi2, 1)

    return prob

def MakeGradient(nsteps, start, end):
    from ROOT import TColor
    r1, g1, b1 = start
    r2, g2, b2 = end
    
#    print start, "to", end

    
    gradient = []

    rstep = float(r2 - r1)/float(nsteps-1)
    gstep = float(g2 - g1)/float(nsteps-1)
    bstep = float(b2 - b1)/float(nsteps-1)
    
#    print "steps", rstep, gstep, bstep

    for i in range(nsteps):
        r = r1 + rstep*i
        g = g1 + gstep*i
        b = b1 + bstep*i
        color = TColor.GetColor(r, g, b)
#        print "(%.2f, %.2f, %.2f) = %i" %(r,g,b, color)
        gradient.append(color)

    return gradient

def Intersections(hist, value, interp = False):
    intersections = []
    for xb in range(1,hist.GetNbinsX()):
        x1 = hist.GetXaxis().GetBinLowEdge(xb)
        y1 = hist.GetBinContent(xb)
        x2 = hist.GetXaxis().GetBinLowEdge(xb+1)
        y2 = hist.GetBinContent(xb+1)
        x3 = hist.GetXaxis().GetBinLowEdge(xb+2)
        lowleft = (y1 < value and y2 > value)
        lowright = (y1 > value and y2 < value)

        if lowleft or lowright:
            if interp:
                intersections += [ Interpolate(x1, y1, x2, y2, value) ]
            elif lowleft:
                intersections += [ x3 ]
            elif lowright:
                intersections += [ x2 ]

    return intersections

def gIntersections(graph, value, xrng, nsteps=1000):
    intersections = []

    xmin, xmax = xrng
    xstp = (xmax - xmin)/(nsteps - 1)

    for xi in range(nsteps-1):
        x1 = xmin + (xi+0)*xstp
        x2 = xmin + (xi+1)*xstp
        y1 = graph.Eval(x1)
        y2 = graph.Eval(x2)
        lowleft = (y1 < value and y2 > value)
        lowright = (y1 > value and y2 < value)

        if lowleft or lowright:
            intersections += [ Interpolate(x1, y1, x2, y2, value) ]
    return intersections


def Interpolate(x1, y1, x2, y2, yvalue):
    m = (y2-y1)/(x2-x1)
    return x1 + (yvalue-y1)/m


def pbloop(iterable, name = "Entries"):
    widgets = [ name+': ',pb.Value(), '/', pb.Total(), ' ', pb.Percentage(), ' ',
                pb.Bar(marker='=',left='[',right=']'),
                ' ',pb.ETA() ]
    pbar = pb.ProgressBar(widgets=widgets, maxval=len(iterable), term_width=100)
    pbar.start()
    for i, val in enumerate(iterable):
        pbar.update(i)
        yield val
    pbar.finish()
    print ""


def BinWidthNormalize(h, width = -1):
    if width < 0: 
        width = h.GetBinWidth(1)

    for i in range(1, h.GetNbinsX()+1):
        wi = h.GetBinWidth(i)
        w = wi/width
        c = h.GetBinContent(i)
        e = h.GetBinError(i)
        h.SetBinContent(i, c/w)
        h.SetBinError(i, e/w)



def MakeBins(bmin, bmax, nbins, log=False):
    bins = []
    if log:
        bmin = log10(bmin)
        bmax = log10(bmax)
    bins = [ x*(bmax-bmin)/(nbins - 1.) + bmin for x in range(nbins) ]
    if log:
        bins = map(lambda x: 10.**x, bins)
    return array('d',bins)
        
        

def DrawDiv(x, c1):
    c1.cd()
    c1.Update()
    ldiv = TLine(x, c1.GetUymin(), x, c1.GetUymax())
    ldiv.SetLineStyle(7)
    ldiv.SetLineWidth(2)
    ldiv.Draw()
    return ldiv
    
def DrawLine(c1, val = 0):
    c1.cd()
    c1.Update()
    ldiv = TLine(c1.GetUxmin(), val, c1.GetUxmax(), val)
    ldiv.SetLineStyle(1)
    ldiv.SetLineWidth(1)
    ldiv.Draw()
    return ldiv
    
def mean(lst):
    return sum(lst)/len(lst)

def stdev(lst):
    m = mean(lst)
    return sqrt( sum(map(lambda x: (x-m)**2, lst)) / (len(lst)-1.) )

def rms(lst):
    return sqrt( sum(map(lambda x: (x)**2, lst)) / len(lst) )


def ProfileX(hist):
    axis = hist.GetXaxis()
    axisother = hist.GetYaxis()

    prof = TH1D(hist.GetName()+"_profX", hist.GetTitle(), axis.GetNbins(), axis.GetBinLowEdge(1), axis.GetBinLowEdge(axis.GetNbins()+1))
    prof.SetXTitle(axis.GetTitle())

    for nb in range(1, axis.GetNbins()+1):
        val = 99999999.
        for nbo in range(1, axisother.GetNbins()+1):
            val = min(val, hist.GetBinContent(nb, nbo))
        prof.SetBinContent(nb, val)

    return prof


def ProfileY(hist):
    axis = hist.GetYaxis()
    axisother = hist.GetXaxis()

    prof = TH1D(hist.GetName()+"_profX", hist.GetTitle(), axis.GetNbins(), axis.GetBinLowEdge(1), axis.GetBinLowEdge(axis.GetNbins()+1))
    prof.SetXTitle(axis.GetTitle())

    for nb in range(1, axis.GetNbins()+1):
        val = 99999999.
        for nbo in range(1, axisother.GetNbins()+1):
            val = min(val, hist.GetBinContent(nbo, nb))
        prof.SetBinContent(nb, val)

    return prof

