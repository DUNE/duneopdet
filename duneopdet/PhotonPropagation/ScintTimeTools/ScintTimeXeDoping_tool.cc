////////////////////////////////////////////////////////////////////////
// Class:       ScintTimeXeDoping
// Plugin Type: tool
// File:        ScintTimeXeDoping_tool.cc ScintTimeXeDoping.h
// Description:
// Generate a random number for timing of LAr scintillation
// Oct. 8 by Mu Wei
////////////////////////////////////////////////////////////////////////
#include "dune/PhotonPropagation/ScintTimeTools/ScintTimeXeDoping.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace phot
{

    //......................................................................    
    ScintTimeXeDoping::ScintTimeXeDoping(fhicl::ParameterSet const& pset)
    : ScintTime()
    , fXeConcentration{pset.get<double>("XeConcentration")}
    , fArSingletTime  {pset.get<double>("ArSingletTime")}
    , fArTripletTime  {pset.get<double>("ArTripletTime")}
    , fXe150nmTime    {pset.get<double>("Xe150nmTime")}
    , fTauAX          {pset.get<double>("TauAX")}
    , fTauXX          {pset.get<double>("TauXX")}
    , fTauN2          {pset.get<double>("TauN2", 0)}
    {
        double invTauN2 = 0;
        if (fTauN2 > 0) invTauN2 = 1./fTauN2;

        // Calculations from D. Totani, Feb. 2021
        fTauTAs = 1./(1./fArSingletTime + 4*fXeConcentration/fTauAX);
        fTauTAt = 1./(1./fArTripletTime +   fXeConcentration/fTauAX + invTauN2);
        fTauTX  = 1./(1./fXe150nmTime   +   fXeConcentration/fTauXX + invTauN2);

        double step = 0.1; //ns

        fMaxProbs = 0;
        fMaxTs = 0;
        double integral = 0;
        while (integral < 0.999) {
            double val = singlet_distro(fMaxTs);
            if (val > fMaxProbs) fMaxProbs = val;
            integral += val*step;
            fMaxTs += step;
        }

        fMaxProbt = 0;
        fMaxTt = 0;
        integral = 0;
        while (integral < 0.999) {
            double val = triplet_distro(fMaxTt);
            if (val > fMaxProbt) fMaxProbt = val;
            integral += val*step;
            fMaxTt += step;
        }

        mf::LogInfo("ScintTimeXeDoping")
            << "Configured:" << "\n"
            << "  XeConcentration: " << fXeConcentration << " ppm\n"
            << "  ArTripletTime:   " << fArTripletTime << " ns\n"
            << "  Xe150nmTime:     " << fXe150nmTime << " ns\n"
            << "  TauAX:           " << fTauAX << " ns\n"
            << "  TauXX:           " << fTauXX << " ns\n"
            << "Calculated:" << "\n"
            << "  TauTA singlet:   " << fTauTAs << " ns\n"
            << "  TauTA triplet:   " << fTauTAt << " ns\n"
            << "  TauTX:           " << fTauTX  << " ns\n"
            << "  MaxTime Singlet: " << fMaxTs  << " ns\n"
            << "  MaxTime Triplet: " << fMaxTt  << " ns\n";
    }
            
    //......................................................................    
    double ScintTimeXeDoping::exp_diff(double t, double tau1, double tau2) const
    {
        using std::exp;
        return ( exp(-t/tau1) - exp(-t/tau2) ) / (tau1 - tau2);
    }
    
    //......................................................................    
    // Returns the time within the time distribution of the scintillation process, when the photon was created.
    // Rejection sampling is used to get a random time from the difference-of-exponentials distribution
    void ScintTimeXeDoping::GenScintTime(bool is_fast, CLHEP::HepRandomEngine& engine)
    {
        CLHEP::RandFlat randflatscinttime{engine};

        //ran1, ran2 = random numbers for the algorithm        
        while (1)
        {
            double ran1 = randflatscinttime();
            double ran2 = randflatscinttime();

            // Use ran1 to pick a proposal time
            // uniformly within the allowed range
            double t = ran1 * (is_fast ? fMaxTs 
                                       : fMaxTt);

            // Rejection prob for this point
            auto p = (is_fast ? singlet_distro(t) / fMaxProbs
                              : triplet_distro(t) / fMaxProbt);

            // Keep the point if ran2 below p
            // (larger value of scint_distro -> more likely to sample)
            if (ran2 < p) {
                timing = t;
                return;
            }
        }
    }
}
    
DEFINE_ART_CLASS_TOOL(phot::ScintTimeXeDoping)
