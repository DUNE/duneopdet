#!/usr/bin/python
import sys, os




###################
# Parse arguments #
###################

startup = False
if sys.argv[0] and "PYTHONSTARTUP" in os.environ:
    startup = True
    
from optparse import OptionParser

usage = "usage: %prog -n <nominal file> -t <test file> [-o <dir>] [-s] [-d] [-e 0,1] [--opdets=0,1] "
parser = OptionParser(usage=usage)
parser.add_option("-n", "--nominal", dest="nominalfile", help="Nominal flashan file", metavar="F")
parser.add_option("-t", "--test",    dest="testfile",    help="Test flashan file",    metavar="F")
parser.add_option("-o", "--outdir",  dest="dir",         help="Directory to output plots", default="plots")
parser.add_option(      "--flash",   dest="flash",       help="Flash-by-flash compare", default=False, action="store_true")
parser.add_option(      "--hit",     dest="hit",         help="Hit-by-hit compare",  default=False, action="store_true")
parser.add_option(      "--assoc",   dest="assoc",       help="Assosciation-by-assosciation compare",  default=False, action="store_true")
parser.add_option(      "--plots",   dest="plots",       help="Hit and flash comparison plots",  default=False, action="store_true")
parser.add_option(      "--rootdir", dest="rootdir",     help="TDirectory in root file from analyzer (%default)", default="flashana")

(options, args) = parser.parse_args()

if not (options.nominalfile and options.testfile):
    print "Both nominal file (-n) and test file (-t) required."
    sys.exit(1)

if not (options.flash or options.hit or options.assoc or options.plots):
    print "No plots to make.  Specify at least one of --flash, --hit, --assoc, --plots."
    sys.exit(1)

# Load libraries after parsing arguments

if startup:
    execfile(os.environ["PYTHONSTARTUP"])

from ROOT import *
from HandyFuncs import VerticalRange, GetHists, pbloop
from collections import defaultdict
import copy
import progressbar as pb


############################
# Define utility functions #
############################

def ApproxCompare(left, right, tolerance = 0.001):
    try:
        return ( abs(right / left - 1.) < tolerance )
    except ZeroDivisionError:
        return (right == 0)

def CountFailures(left, right):
    return len(ListFailures(left,right))

def ListFailures(left, right):
    if type(left) is Flash: return ListFailuresFlash(left, right)
    if type(left) is Assoc: return ListFailuresAssoc(left, right)
    if type(left) is Hit:   return ListFailuresHit(left, right)

def ListFailuresFlash(left, right):
    failures = []
    tolerance = 0.001
    if left.InBeamFrame != right.InBeamFrame: failures += [ 'InBeamFrame']
    if left.OnBeamTime  != right.OnBeamTime:  failures += [ 'OnBeamTime' ]
    if not ApproxCompare(left.YCenter,   right.YCenter,   tolerance): failures += [ 'YCenter' ]
    if not ApproxCompare(left.ZCenter,   right.ZCenter,   tolerance): failures += [ 'ZCenter' ]
    if not ApproxCompare(left.YWidth,    right.YWidth,    tolerance): failures += [ 'YWidth' ]
    if not ApproxCompare(left.ZWidth,    right.ZWidth,    tolerance): failures += [ 'YWidth' ]
    if not ApproxCompare(left.FlashTime, right.FlashTime, tolerance): failures += [ 'FlashTime' ]
    if not ApproxCompare(left.FlashFrame,right.FlashFrame,tolerance): failures += [ 'Frame' ]
    if not ApproxCompare(left.AbsTime,   right.AbsTime,   tolerance): failures += [ 'AbsTime' ] 
    if not ApproxCompare(left.TotalPE,   right.TotalPE,   0.5): failures += [ 'TotalPE' ]
    return failures

def ListFailuresAssoc(left, right):
    failures = []
    tolerance = 0.000001
    if left.EventID != right.EventID: failures += [ 'EventID' ]
    if left.HitFrame != right.HitFrame: failures += [ 'HitFrame' ]
    if left.FlashFrame != right.FlashFrame: failures += [ 'FlashFrame' ]
    #if left.FlashID != right.FlashID: failures += [ 'FlashID' ]
    if not ApproxCompare(left.HitPeakTime,      right.HitPeakTime,    tolerance): failures += [ 'HitPeakTime' ]
    if not ApproxCompare(left.HitPE,            right.HitPE,          tolerance): failures += [ 'HitPE' ]
    if not ApproxCompare(left.FlashTime,        right.FlashTime,      tolerance): failures += [ 'FlashTime' ]
    if not ApproxCompare(left.FlashPE,          right.FlashPE,        tolerance): failures += [ 'FlashPE' ]
    return failures

def ListFailuresHit(left, right):
    failures = []
    tolerance = 0.000000001
    if left.Frame != right.Frame: failures += [ 'Frame' ]
    if not ApproxCompare(left.PeakTimeAbs,   right.PeakTimeAbs, tolerance): failures += [ 'PeakTimeAbs' ] 
    if not ApproxCompare(left.PeakTime,      right.PeakTime,    tolerance): failures += [ 'PeakTime' ]
    if not ApproxCompare(left.Width,         right.Width,       tolerance): failures += [ 'Width' ]
    if not ApproxCompare(left.Area,          right.Area,        tolerance): failures += [ 'Area' ]
    if not ApproxCompare(left.Amplitude,     right.Amplitude,   tolerance): failures += [ 'Amplitude' ]
    if not ApproxCompare(left.PE,            right.PE,          tolerance): failures += [ 'PE' ]
    if not ApproxCompare(left.FastToTotal,   right.FastToTotal, tolerance): failures += [ 'FastToTotal' ]
    return failures



##################
# Define objects #
##################

class Flash(object):

    def __init__(self, tree):
        self.EventID     = tree.EventID
        self.FlashID     = tree.FlashID
        self.YCenter     = tree.YCenter
        self.ZCenter     = tree.ZCenter
        self.YWidth      = tree.YWidth
        self.ZWidth      = tree.ZWidth
        self.FlashFrame  = 0 #tree.FlashFrame
        self.FlashTime   = tree.FlashTime
        self.AbsTime     = tree.AbsTime
        self.InBeamFrame = tree.InBeamFrame
        self.OnBeamTime  = tree.OnBeamTime
        self.TotalPE     = tree.TotalPE

    def __lt__(self, other):
        if not isinstance(other, Flash):
            return NotImplemented
        return self.AbsTime < other.AbsTime

    def __le__(self, other):
        if not isinstance(other, Flash):
            return NotImplemented
        return self.AbsTime <= other.AbsTime

    def __gt__(self, other):
        if not isinstance(other, Flash):
            return NotImplemented
        return self.AbsTime > other.AbsTime

    def __ge__(self, other):
        if not isinstance(other, Flash):
            return NotImplemented
        return self.AbsTime >= other.AbsTime

        
    def __eq__(self, other):
        if not isinstance(other, Flash):
            return NotImplemented
        
        if self.EventID     != other.EventID:     return False
            
        failures = ListFailures(self, other)
        if "AbsTime" in failures or "FlashTime" in failures: return False

        if len(failures) > 3:
            return False

        return True

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def PrintStr(self):
        return "t={0.AbsTime:9.3f} q={0.TotalPE:7.2f} FlashID={0.FlashID:2} Frame={0.FlashFrame:1}".format(self)


class Assoc(object):

    def __init__(self, tree):
        self.EventID        = tree.EventID
        self.FlashID        = tree.FlashID
        self.HitID          = tree.HitID
        self.OpChannel      = tree.OpChannel
        self.HitPeakTimeAbs = tree.HitPeakTimeAbs
        self.HitPeakTime    = tree.HitPeakTime
        self.HitPE          = tree.HitPE
        self.FlashPE        = tree.FlashPE
        self.FlashTimeAbs   = tree.FlashTimeAbs
        self.FlashTime      = tree.FlashTime
        self.HitFrame       = tree.HitFrame
        self.FlashFrame     = tree.FlashFrame

    def __lt__(self, other):
        if not isinstance(other, Assoc):
            return NotImplemented
        if self.HitPeakTimeAbs < other.HitPeakTimeAbs:
            return true
        if self.HitPeakTimeAbs == other.HitPeakTimeAbs and self.OpChannel < other.OpChannel:
            return true
        return False

    def __le__(self, other):
        if not isinstance(other, Assoc):
            return NotImplemented
        if self < other:
            return True
        if self.HitPeakTimeAbs == other.HitPeakTimeAbs and self.OpChannel == other.OpChannel:
            return true
        return False

    def __gt__(self, other):
        if not isinstance(other, Assoc):
            return NotImplemented
        return not (self <= other)

    def __ge__(self, other):
        if not isinstance(other, Assoc):
            return NotImplemented
        return not (self < other)
        
    def __eq__(self, other):
        if not isinstance(other, Assoc):
            return NotImplemented
        if self.EventID     != other.EventID:   return False
        if self.OpChannel   != other.OpChannel: return False
        failures = ListFailures(self, other)
        if "HitPeakTimeAbs" in failures \
          or "HitPeakTime"  in failures \
          or "PE"           in failures :       return False
        return True

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def PrintStr(self):
        return "{0.OpChannel:2d} t={0.HitPeakTimeAbs:8.3f}/{0.FlashTimeAbs:8.3f}  q={0.HitPE:7.2f}/{0.FlashPE:7.2f} frame={0.HitFrame:1}/{0.FlashFrame:1}".format(self)
    #.OpChannel, self.HitPeakTimeAbs, self.HitPE, self.HitFrame, self.FlashFrame, self.FlashID)



class Hit(object):

    def __init__(self, tree):
        self.EventID     = tree.EventID
        self.OpChannel   = tree.OpChannel
        self.PeakTimeAbs = tree.PeakTimeAbs
        self.PeakTime    = tree.PeakTime
        self.Frame       = tree.Frame
        self.Width       = tree.Width
        self.Area        = tree.Area
        self.Amplitude   = tree.Amplitude
        self.PE          = tree.PE
        self.FastToTotal = tree.FastToTotal


    def __lt__(self, other):
        if not isinstance(other, Hit):
            return NotImplemented
        return self.PeakTimeAbs < other.PeakTimeAbs

    def __le__(self, other):
        if not isinstance(other, Hit):
            return NotImplemented
        return self.PeakTimeAbs <= other.PeakTimeAbs

    def __gt__(self, other):
        if not isinstance(other, Hit):
            return NotImplemented
        return self.PeakTimeAbs > other.PeakTimeAbs

    def __ge__(self, other):
        if not isinstance(other, Hit):
            return NotImplemented
        return self.PeakTimeAbs >= other.PeakTimeAbs
        
    def __eq__(self, other):
        if not isinstance(other, Hit):
            return NotImplemented
        if self.EventID     != other.EventID:   return False
        if self.OpChannel   != other.OpChannel: return False
        failures = ListFailures(self, other)
        if "PeakTimeAbs" in failures \
          or "PeakTime" in failures:            return False
        if len(failures) > 0:                   return False
        return True

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def PrintStr(self):
        return "{0:2d} t={1:9.3f}  q={2:7.2f}".format(self.OpChannel, self.PeakTimeAbs, self.PE)

    


################
# Shared setup #
################

versions = [ "nominal", "test" ]
S = versions[0]
R = versions[1]

files = {}
files[S] = TFile(options.nominalfile)
files[R] = TFile(options.testfile)

hittrees = {}
flashtrees = {}
assoctrees = {}
for v in versions:
    hittrees[v]   = files[v].Get(options.rootdir+"/PerOpHitTree")
    flashtrees[v] = files[v].Get(options.rootdir+"/PerFlashTree")
    assoctrees[v] = files[v].Get(options.rootdir+"/FlashHitMatchTree")





####################
# Flash Comparison #
####################

if options.flash:

    all = {}
    for v in versions:
        for e in range(flashtrees[v].GetEntries()):
            flashtrees[v].GetEntry(e)
            if flashtrees[v].EventID not in all:
                all[flashtrees[v].EventID] = defaultdict(list)
            all[flashtrees[v].EventID][v].append(Flash(flashtrees[v]))

    for EventID in sorted(all.keys()):
        for v in versions:
            all[EventID][v].sort()
        print EventID, ":", len(all[EventID][S]), len(all[EventID][R])

    onlies = copy.deepcopy(all)

    for E in sorted(all.keys()):
        iS = 0

        matches = []
        failures = [0]*8

        while True:
            #print iS, "/", len(onlies[E][S])
            found = False
            for iR in range(len(onlies[E][R])):
                if onlies[E][S][iS] == onlies[E][R][iR]:
                    matches.append( { S: onlies[E][S].pop( iS ),
                                      R: onlies[E][R].pop( iR ) } )
                    found = True
                    break
            if not found:
                iS += 1
            if iS >= len(onlies[E][S]):
                break

        print "Event", E
        print "  matches", len(matches)
        for v in versions:
            print "  only", v, len(onlies[E][v])

        for v in versions:
            if len(onlies[E][v]):
                print "  Only in",v
                for flsh in onlies[E][v]:
                    index = all[E][v].index(flsh)
                    print "    ",flsh.PrintStr()#, " preceeded by ", all[E][v][index-1].PrintStr()


        print "  Matches"
        for flsh in matches:
            failures_here = ListFailures(flsh[S],flsh[R])
            if failures_here:
                print "    ", S+":", flsh[S].PrintStr(), "  ", R+":", flsh[R].PrintStr(), " ".join(failures_here)
        




##################
# Hit Comparison #
##################

if options.hit:

    onlies = {}
    for v in versions:
        for e in pbloop(range(hittrees[v].GetEntries())):
            hittrees[v].GetEntry(e)
            if hittrees[v].EventID not in onlies:
                onlies[hittrees[v].EventID]= defaultdict(list)
            onlies[hittrees[v].EventID][v].append(Hit(hittrees[v]))

    for EventID in sorted(onlies.keys()):
        for v in versions:
            onlies[EventID][v].sort()
        print EventID, ":", len(onlies[EventID][S]), len(onlies[EventID][R])


    for E in sorted(onlies.keys()):
        iS = 0

        matches = []
        failures = [0]*8


        maxval = len(onlies[E][S])
        widgets = [ 'Entries: ',pb.Value(), '/', pb.Total(), ' ', pb.Percentage(), ' ',
                    pb.Bar(marker='=',left='[',right=']'),
                    ' ',pb.ETA() ]
        pbar = pb.ProgressBar(widgets=widgets, maxval=maxval, term_width=100)
        pbar.start()
        while True:
            pbar.update(maxval-len(onlies[E][S]))
            #print iS, "/", len(onlies[E][S])
            found = False
            for iR in range(len(onlies[E][R])):
                if onlies[E][S][iS] == onlies[E][R][iR]:
                    matches.append( { S: onlies[E][S].pop( iS ),
                                      R: onlies[E][R].pop( iR ) } )
                    found = True
                    break
            if not found:
                iS += 1
            if iS >= len(onlies[E][S]):
                break
        pbar.finish()
        print ""

        print "Event", E
        print "  matches", len(matches)
        for v in versions:
            print "  only", v, len(onlies[E][v])

        for v in versions:
            if len(onlies[E][v]):
                print "  Only in",v
                for flsh in onlies[E][v]:
                    print "    ",flsh.PrintStr()

        #print "  Matches"
        #for flsh in matches:
        #    print "    ", S+":", flsh[S].PrintStr(), "  ", R+":", flsh[R].PrintStr(), " ".join(ListFailures(flsh[S],flsh[R]))




    
####################
# Assoc Comparison #
####################

if options.assoc:


    onlies = {}
    for v in versions:
        for e in pbloop(range(assoctrees[v].GetEntries())):
            assoctrees[v].GetEntry(e)
            if assoctrees[v].EventID not in onlies:
                onlies[assoctrees[v].EventID] = defaultdict(list)
            onlies[assoctrees[v].EventID][v].append(Assoc(assoctrees[v]))

    for EventID in sorted(onlies.keys()):
        for v in versions:
            onlies[EventID][v].sort()
        print EventID, ":", len(onlies[EventID][S]), len(onlies[EventID][R])


    for E in sorted(onlies.keys()):
        iS = 0

        matches = []
        failures = [0]*8
        print "Event", E
        print "Unmatched lengths", S, len(onlies[E][S]), R, len(onlies[E][R])
        maxval = len(onlies[E][S])
        widgets = [ 'Event %i: '%E,pb.Value(), '/', pb.Total(), ' ', pb.Percentage(), ' ',
                    pb.Bar(marker='=',left='[',right=']'),
                    ' ',pb.ETA() ]
        pbar = pb.ProgressBar(widgets=widgets, maxval=maxval, term_width=100)
        pbar.start()
        while onlies[E][S] and onlies[E][R]:
            pbar.update(maxval-len(onlies[E][S]))
            if onlies[E][S][0] == onlies[E][R][0]:
                matches.append( { S: onlies[E][S].pop( 0 ),
                                  R: onlies[E][R].pop( 0 ) } )
            elif onlies[E][S][0] < onlies[E][R][0]:
                matches.append( { S: onlies[E][S].pop( 0 ) } )
            elif onlies[E][S][0] > onlies[E][R][0]:
                matches.append( { R: onlies[E][R].pop( 0 ) } )
            else: # Same time, unmatched
                matches.append( { S: onlies[E][S].pop( 0 ) } )
                matches.append( { R: onlies[E][R].pop( 0 ) } )
        pbar.finish()
        print ""
        print "matched lengths  ", S, len(onlies[E][S]), R,  len(onlies[E][R])

        while onlies[E][S]:
            matches.append( { S: onlies[E][S].pop( 0 ) } )
        while onlies[E][R]:
            matches.append( { R: onlies[E][R].pop( 0 ) } )

        ## while True:
        ##     pbar.update(maxval-len(onlies[E][S]))
        ##     #print iS, "/", len(onlies[E][S])
        ##     found = False
        ##     for iR in range(len(onlies[E][R])):
        ##         if onlies[E][S][iS] == onlies[E][R][iR]:
        ##             matches.append( { S: onlies[E][S].pop( iS ),
        ##                               R: onlies[E][R].pop( iR ) } )
        ##             found = True
        ##             break
        ##     if not found:
        ##         matches.append( { S: 
        ##         iS += 1
        ##     if iS >= len(onlies[E][S]):
        ##         break
        ## pbar.finish()
        ## print ""

        #print "Event", E
        #print "  matches", len(matches)
        #for v in versions:
        #    print "  only", v, len(onlies[E][v])

        #for v in versions:
        #    if len(onlies[E][v]):
        #        print "  Only in",v
        #        for flsh in onlies[E][v]:
        #            print "    ",flsh.PrintStr()

        print "  Matches"
        for hit in matches:
            print "",
            if S in hit:   print  S+":", hit[S].PrintStr(),
            else:          print "{0:65}".format(""),
            print "",
            if R in hit:   print "    ", R+":", hit[R].PrintStr(),
            else:          print "{0:65}".format(""),
            print "  ",
            if S in hit and R in hit:
                failures_here = ListFailures(hit[S],hit[R])
                print " ".join(failures_here)
            else:
                print






    



#######################
# Flash and Hit Plots #
#######################

if options.plots:
    if not os.path.exists(options.dir):
        os.makedirs(options.dir)

    color = { S:kBlack, R:kRed }
    hstyle = "hist"
    mstyle = 20
    lstyle = 2

    nbins = 100


    c1 = TCanvas("c1","c1")

    for var, name, select in [ ("EventID","EventID",""),
                               ("OpChannel","OpChannel",""),
                               ("PeakTimeAbs","PeakTimeAbs",""),
                               ("PeakTime","PeakTime",""),
                               ("Frame","Frame",""),
                               ("Width","Width","Width<1"),
                               ("TMath::Log10(Area)","Area",""),
                               ("TMath::Log10(Amplitude)","Amplitude",""),
                               ("TMath::Log10(PE)","PE",""),
                               ("FastToTotal","FastToTotal","")]:

        hittrees[S].Draw("{0}>>h{1}{2:d}({3:d})".format(var,name,1,nbins), select, "hist")
        hittrees[R].Draw("{0}>>h{1}{2:d}({3:d})".format(var,name,2,nbins), select, hstyle+"same")
        h1 = gROOT.FindObject("h"+name+"1")
        h2 = gROOT.FindObject("h"+name+"2")
        h1.SetLineColor(color[S])
        h2.SetLineColor(color[R])
        h2.SetLineStyle(2)
        h2.SetMarkerColor(color[R])
        h2.SetMarkerStyle(mstyle)
        c1.Print(os.path.join(options.dir,"perhit_"+name+".png"))

        hRatio = h2.Clone("hRatio")
        hRatio.Divide(h1)
        hRatio.Draw("p")
        #hRatio.GetYaxis().SetRangeUser(0.5,1.1)
        VerticalRange([hRatio], ratio=True, forceOne=True, ignoreError=True)
        gPad.Update()
        l1 = TLine(gPad.GetUxmin(), 1, gPad.GetUxmax(), 1)
        l1.SetLineWidth(1)
        l1.SetLineColor(kBlack)
        l1.Draw()
        hRatio.Draw("psame")
        c1.Print(os.path.join(options.dir,"perhit_"+name+"_ratio.png"))


    for var, name, select in [("EventID","EventID",""),
                              ("FlashID","FlashID",""),
                              ("YCenter","YCenter",""),
                              ("ZCenter","ZCenter",""),
                              ("YWidth","YWidth",""),
                              ("ZWidth","ZWidth",""),
                              ("FlashTime","FlashTime",""),
                              ("AbsTime","AbsTime",""),
                              ("InBeamFrame","InBeamFrame",""),
                              ("OnBeamTime","OnBeamTime",""),
                              ("TMath::Log10(TotalPE)","TotalPE","") ]:

        flashtrees[S].Draw("{0}>>h{1}{2:d}({3:d})".format(var,name,1,nbins), select, "hist")
        flashtrees[R].Draw("{0}>>h{1}{2:d}({3:d})".format(var,name,2,nbins), select, hstyle+"same")
        h1 = gROOT.FindObject("h"+name+"1")
        h2 = gROOT.FindObject("h"+name+"2")
        h1.SetLineColor(color[S])
        h2.SetLineColor(color[R])
        h2.SetLineStyle(lstyle)
        h2.SetMarkerColor(color[R])
        h2.SetMarkerStyle(mstyle)
        c1.Print(os.path.join(options.dir,"perflash_"+name+".png"))

        hRatio = h2.Clone("hRatio")
        hRatio.Divide(h1)
        hRatio.Draw("p")
        #hRatio.GetYaxis().SetRangeUser(0.5,1.1)
        VerticalRange([hRatio], ratio=True, forceOne=True, ignoreError=True)
        gPad.Update()
        l1 = TLine(gPad.GetUxmin(), 1, gPad.GetUxmax(), 1)
        l1.Draw()
        l1.SetLineWidth(1)
        l1.SetLineColor(kBlack)
        hRatio.Draw("psame")
        c1.Print(os.path.join(options.dir,"perflash_"+name+"_ratio.png"))

