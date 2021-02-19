////////////////////////////////////////////////////////////////////////
// Class:       ScintTimeXeDoping
// Plugin Type: tool
// File:        ScintTimeXeDoping_tool.cc ScintTimeXeDoping.h
// Description:
// Generate a random number for timing of LAr scintillation
// Oct. 8 by Mu Wei
////////////////////////////////////////////////////////////////////////
#ifndef ScintTimeXeDoping_H
#define ScintTimeXeDoping_H

#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

// Random number engine
#include "CLHEP/Random/RandFlat.h"

namespace phot
{
    class ScintTimeXeDoping : public ScintTime
    {
    public:
        explicit ScintTimeXeDoping(fhicl::ParameterSet const& pset);
        void GenScintTime(bool is_fast, CLHEP::HepRandomEngine& engine)  ;
        
    private:
        int           fLogLevel;

        // From fhicl configuration
        double fXeConcentration; // ppm
        double fArSingletTime;   // ns
        double fArTripletTime;   // ns
        double fXe150nmTime;     // ns ArXe* 150 nm decay time
        double fTauAX;           // ns ArXe* formation time
        double fTauXX;           // ns XeXe* formation time
        double fTauN2;           // ns N2 quenching time

        // Calculated during setup
        double fTauTAs;
        double fTauTAt;
        double fTauTX;
        double fMaxProbs;
        double fMaxProbt;
        double fMaxTs;
        double fMaxTt;
        
        // utility functions
        double exp_diff(double t, double tau1, double tau2) const;
        double singlet_distro(double t) const { return exp_diff(t, fTauTAs, fTauTX); };
        double triplet_distro(double t) const { return exp_diff(t, fTauTAt, fTauTX); };
    };
}
#endif
