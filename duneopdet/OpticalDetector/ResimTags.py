#!/usr/bin/env python


def CreateTags():
    preqes  = [ 15, 25, 45 ]
    prethrs = [ 2, 3 ]
    prerefs  = [ "Opt", "Pes" ]

    defqes  = [ 35 ]
    defthrs = [ 1 ]
    defrefs   = [ "Non" ]

    areas = {}
    noises = {}
    snrs = {}
    reflected = {}


    sets = [ ("DEF", defqes, defrefs, defthrs),
             ("EFF", preqes, defrefs, defthrs),
             ("REF", defqes, prerefs, defthrs),
             ("THR", defqes, defrefs, prethrs) ]


    params = {}

    # Fixed noise value
    noise = 100
    snr = 5
    peakadc = 9
    
    for ( name, iqes, irefs, ithrs ) in sets:
        for qe in iqes:
            for ref in irefs:
                for thr in ithrs:
                    tag = "{0}{1:02d}QE{2}Refl{3:d}PE".format(name, qe, ref, thr)

                    # Reduce Eff by 30% for shadowing
                    qeval = float(qe) / 1000.
                    qeval *= 0.7

                    # Produce the parameter tables
                    params[tag] = {}
                    params[tag]["qe"] = qeval
                    params[tag]["darkrate"] = noise
                    params[tag]["linerms"] = float(peakadc)/snr
                    params[tag]["thresh"] = int( (float(thr)-0.3) * peakadc )
                    
                    # Reflections
                    if ref == "Opt":
                        params[tag]["refqe"] = 1.44*qeval
                    elif ref == "Pes":
                        params[tag]["refqe"] = 0.84*qeval
                        params[tag]["qe"] /= 2
                    else:
                        params[tag]["refqe"] = 0
                    

    tags = sorted(params.keys())
    return tags, params


def print_seq(tags, form):
    for tag in tags[:-1]:
        print form.format(tag)+","
    else:
        print form.format(tags[-1])


tags, params = CreateTags()



print """
#include "largeantmodules_dune.fcl"
#include "detsimmodules_dune.fcl"
#include "opticaldetectormodules_dune.fcl"
#include "OpSlicer.fcl"
#include "FlashMatchAna.fcl"
#include "SNAna.fcl"
"""
print "BEGIN_PROLOG"
print

print 
print "############################################################################"
print "pd_detsim_modules: {"
for tag in tags:
    print "      opdigi{0}:    @local::dunefd_opdigi_threegang".format(tag)
print "}"
print

for tag in tags:
    pm = params[tag]
    print "pd_detsim_modules.opdigi{0}.QEOverride:                  {1:.6f}".format(tag, pm["qe"])
    print "pd_detsim_modules.opdigi{0}.QERefOverride:               {1:.6f}".format(tag, pm["refqe"])
    print "pd_detsim_modules.opdigi{0}.DarkNoiseRate:               {1:d} #Hz".format(tag, pm["darkrate"])
    print "pd_detsim_modules.opdigi{0}.LineNoiseRMS:                {1:.3f}".format(tag, pm["linerms"])
    print "pd_detsim_modules.opdigi{0}.algo_threshold.ADCThreshold: {1:.3f}".format(tag, pm["thresh"])
    print 

print "pd_detsim_path: ["
print_seq(tags, "                 opdigi{0}")
print "                ]"


    

print ""
print "############################################################################"
print "pd_reco_modules: {"
for tag in tags:
    print "      ophit{0}:     @local::dunefd_ophit".format(tag)
for tag in tags:
    print "      opflash{0}:   @local::dunefd_opflash".format(tag)
for tag in tags:
    print "      opslicer{0}:   @local::standard_opslicer".format(tag)
print "}"
print

for tag in tags:
    pm = params[tag]
    print "pd_reco_modules.ophit{0}.InputModule:             opdigi{0}".format(tag)
    print "pd_reco_modules.opflash{0}.InputModule:           ophit{0}".format(tag)
    print "pd_reco_modules.opslicer{0}.OpHitModuleLabel:     ophit{0}".format(tag)
    print 

recoparts = []
print "pd_reco_path: ["
print_seq(tags, "                 ophit{0}, opflash{0}, opslicer{0}")
print "              ]"




    

print ""
print "############################################################################"
print "pd_ana_modules: {"
for tag in tags:
    print "      flashmatch{0}:  @local::marley_flashmatchana".format(tag)
for tag in tags:
    print "      slicematch{0}:  @local::marley_flashmatchana".format(tag)
for tag in tags:
    print "      snana{0}:       @local::standard_snana".format(tag)

print "}"
print

for tag in tags:
    pm = params[tag]
    print "pd_ana_modules.flashmatch{0}.OpDetWaveformLabel: opdigi{0}".format(tag)
    print "pd_ana_modules.flashmatch{0}.OpHitModuleLabel:   ophit{0}".format(tag)
    print "pd_ana_modules.flashmatch{0}.OpFlashModuleLabel: opflash{0}".format(tag)
    print "pd_ana_modules.slicematch{0}.OpDetWaveformLabel: opdigi{0}".format(tag)
    print "pd_ana_modules.slicematch{0}.OpHitModuleLabel:   ophit{0}".format(tag)
    print "pd_ana_modules.slicematch{0}.OpFlashModuleLabel: opslice{0}".format(tag)
    print "pd_ana_modules.snana{0}.OpHitModuleLabel:        ophit{0}".format(tag)
    print
    
print "pd_ana_path: ["
print_seq(tags, "                 flashmatch{0}, slicematch{0}, snana{0}")
print "             ]"


print "END_PROLOG"
