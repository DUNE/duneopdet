#!/usr/bin/python
import sys, os

startup = False
if sys.argv[0] and "PYTHONSTARTUP" in os.environ:
    startup = True
    
from optparse import OptionParser

usage = "usage: %prog -n <nominal file> -t <test file> [-o <dir>] [-s] [-d] [-e 0,1] [--opdets=0,1] "
parser = OptionParser(usage=usage)
parser.add_option("-n", "--nominal", dest="nominalfile", help="Nominal hist file", metavar="F")
parser.add_option("-t", "--test",    dest="testfile",    help="Test hist file",    metavar="F")
parser.add_option("-o", "--outdir",  dest="dir",         help="Directory to output plots", default="plots")
parser.add_option("-s", "--sim",     dest="sim",         help="Make simulation plots", default=False, action="store_true")
parser.add_option("-d", "--digi",    dest="digi",        help="Make digitizer plots",  default=False, action="store_true")
parser.add_option("-e", "--events",  dest="events",      help="List of events to draw an individual opdet, comma-separated")
parser.add_option(      "--opdets",  dest="opdets",      help="List of opdets to draw, comma-separated", default="0")
parser.add_option(      "--wide",    dest="wide",        help="Make widely binned plots (for uBooNE)")

(options, args) = parser.parse_args()

if not (options.nominalfile and options.testfile):
    print "Both nominal file (-n) and test file (-t) required."
    sys.exit(1)

if not (options.sim or options.digi or options.event):
    print "No plots to make.  Specify at least one of --sim, --digi, or --events=M,N,O."
    sys.exit(1)

if not os.path.exists(options.dir):
    os.makedirs(options.dir)

# Load libraries after parsing arguments

if startup:
    execfile(os.environ["PYTHONSTARTUP"])

from ROOT import *
from HandyFuncs import VerticalRange, GetHists, pbloop
from collections import defaultdict

# Load the files
files = {}
versions = [ "nominal", "test" ]
colors = { "nominal":kBlue, "test":kRed }
N = versions[0]
T = versions[1]

files[N] = TFile(options.nominalfile)
files[T] = TFile(options.testfile)

c = {}
        

############################
# Compare simulation plots #
############################

if options.sim:
    if options.wide:
        plotinfo = [ ("CountAll", "OpDetEvents", 100, 0, 100000),
                     ("CountAll", "OpDets",      100, 0, 5000) ]
    else:
        plotinfo = [ ("CountAll", "OpDetEvents", 100, 0, 4500),
                     ("CountAll", "OpDets",      100, 0, 500) ]
    for var, treename, nbins, xmin, xmax in plotinfo:
        trees = {}
        for v in versions:
            trees[v] = files[v].Get("pmtresponse/"+treename)

        c1 = TCanvas("_".join(["c",var,treename]),"_".join(["c",var,treename]))
        dopt = "hist"
        hists = {}
        for v in versions:
            hists[v] = TH1D("h"+v,v, nbins, xmin, xmax)
            hists[v].SetLineColor(colors[v])
            if "same" in dopt: hists[v].SetLineColor(2)
            trees[v].Draw(var+">>h"+v, "", dopt)
            if "same" not in dopt: dopt+="same"
            print v, var, hists[v].Integral()

        VerticalRange(GetHists(c1), ignoreError=True)
        c1.Print(os.path.join(options.dir, "g4_{0}_by_{1}.png".format(treename, var)))
        c[var,treename] = c1
        

###########################
# Compare digitizer plots #
###########################

nom_channel_map = {  0:0,  1:0, 2:0,
                     3:1,  4:1, 5:1,
                     6:2,  7:2,
                     8:3,  9:3, 10:3,
                    11:4, 12:4, 13:4,
                    14:5, 15:5, 16:5,
                    17:6, 18:6, 19:6,
                    20:7, 21:7, 22:7 }


if options.digi:

    hists = { "per_Event":{},
              "per_Channel":{}, "by_Channel":{}, "by_Channel_cnt":{}, 
              "per_OpDet":{},   "by_OpDet":{},   "by_OpDet_cnt":{},
            }

    for name in hists:
        c[name] = TCanvas("c"+name, "c"+name)

    dopt = "hist"
    
    for v in versions:
        eventintegral = defaultdict(int)
        opdetintegral = defaultdict(int)
        
        hists["per_Event"][v]      = TH1D("hperevent"     +v, "ADCs Per Event",           100, 0,  4e4)
        hists["per_Channel"][v]    = TH1D("hperchannel"   +v, "ADCs Per Channel",         100, 0,  3e3)
        hists["by_Channel"][v]     = TH1D("hbychannel"    +v, "ADCs by Channel Number",    96, 0,   96)
        hists["by_Channel_cnt"][v] = TH1D("hbychannelcnt" +v, "Count by Channel Number",   96, 0,   96)
        hists["per_OpDet"][v]      = TH1D("hperopdet"     +v, "ADCs Per OpDet",           100, 0, 5000)
        hists["by_OpDet"][v]       = TH1D("hbyopdet"      +v, "ADCs by OpDet Number",       8, 0,    8)
        hists["by_OpDet_cnt"][v]   = TH1D("hbyopdetcnt"   +v, "Count by OpDet Number",      8, 0,    8)

        for name in hists:
            hists[name][v].SetLineColor(colors[v])

                
        tdir = files[v].Get("opdigiana")
        keys = list(tdir.GetListOfKeys())
        for rootkey in pbloop(keys, v+":"):
            hist = tdir.Get(rootkey.GetName())
            
            event = int(rootkey.GetName().split("_")[1])
            channel = int(rootkey.GetName().split("_")[3])
            integral = hist.Integral()
            if v == T: opdet = int(channel/12)
            else:      opdet = nom_channel_map[channel]

            eventintegral[event] += integral
            hists["per_Channel"][v].Fill(integral)
            hists["by_Channel"][v].Fill(channel, integral)
            if integral > 0: hists["by_Channel_cnt"][v].Fill(channel)
            opdetintegral[opdet,event] += integral
            hists["by_OpDet"][v].Fill(opdet, integral)
            if integral > 0: hists["by_OpDet_cnt"][v].Fill(opdet)
            
        for val in eventintegral.values():
            if val > 4e4: print "Overflow: ",val
            hists["per_Event"][v].Fill(val)
        for val in opdetintegral.values():
            hists["per_OpDet"][v].Fill(val)

        for name in hists:
            c[name].cd()
            hists[name][v].Draw(dopt)

        if "same" not in dopt: dopt+="same"

    for name in hists:
        c[name].cd()
        if name == "per_Channel":
            VerticalRange(GetHists(c[name]), ignoreError = True, absMin = 0.5)
            c[name].SetLogy()
        else:
            VerticalRange(GetHists(c[name]), ignoreError = True)
        c[name].Update()
        c[name].Print(os.path.join(options.dir, "ADC_"+name+".png"))




if options.events:
    events = map(int, options.events.split(","))
    opdets = map(int, options.opdets.split(","))
    for event, opdet in [ (e,o) for e in events for o in opdets ]:
        channels = { N:[], T:[] }
        for c,o in nom_channel_map.items():
            if o == opdet:
                channels[N].append(c)
        channels[T] = range(o*12,(o+1)*12)

        hists = {}
        trees = {}
        for v in versions:
            trees[v] = files[v].Get("pmtresponse/OpDets")
            trees[v].GetEntry(8*(event-1)+opdet)

            tdir = files[v].Get("opdigiana")
            for c in channels[v]:
                htemp = tdir.Get("Event_{0:d}_OpDet_{1:d}_Pulse_0".format(event, c))
                print v, event, c, htemp.Integral()
                if v not in hists:
                    hists[v] = htemp.Clone("h"+v)
                    hists[v].SetLineColor(colors[v])
                    hists[v].GetXaxis().SetRangeUser(0, 8000)
                else:
                    hists[v].Add(htemp)

        ratio = {}
        for v in versions:
            ratio[v] = hists[v].Clone("hratio"+v)
            ratio[v].Divide(hists[N])

        c1 = TCanvas("c1_event{0:d}_opdet{1:d}".format(event, opdet),"c_event{0:d}_opdet{1:d}")
        dopt = ""
        for v in versions:
            hists[v].Draw(dopt)
            print v, hists[v].Integral(), hists[v].Integral()/trees[v].CountAll
            if "same" not in dopt: dopt += "same"
        c1.Print(os.path.join(options.dir,"waveform_event{0:d}_opdet{1:d}.png".format(event, opdet)))

        c1 = TCanvas("c2_event{0:d}_opdet{1:d}".format(event, opdet),"c_event{0:d}_opdet{1:d}")
        dopt = ""
        for v in versions:
            ratio[v].Draw(dopt)
            if "same" not in dopt: dopt += "same"
        c2.Print(os.path.join(options.dir,"waveform_event{0:d}_opdet{1:d}_ratio.png".format(event, opdet)))

        c[c1.GetName()] = c1
        c[c2.GetName()] = c2

