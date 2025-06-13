#ifndef ScintTimeXeDoping2_H
#define ScintTimeXeDoping2_H

#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"
#include "canvas/Utilities/Exception.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

// Random number engine
#include "CLHEP/Random/RandFlat.h"

namespace phot
{
    class ScintTimeXeDoping2 : public ScintTime
    {
    public:
        explicit ScintTimeXeDoping2(fhicl::ParameterSet const& pset);
        void GenScintTime(bool is_fast, CLHEP::HepRandomEngine& engine) override;
        double fastScintTime() override;
        double slowScintTime() override;

        
    private:
        //dla int           fLogLevel;
        std::unique_ptr<CLHEP::RandFlat> fUniformGen;

        // From fhicl configuration
        double fXeConcentration; // ppm
        double fN2Concentration; // ppm
        double fArSingletTime;   // ns
        double fArTripletTime;   // ns
        double fXe150nmTime;     // ns ArXe* 150 nm decay time
        double fTauAX;           // ns ArXe* formation time
        double fTauXX;           // ns XeXe* formation time
        double fTauN2ArAr;           // ns N2 quenching time
        double fTauN2ArXe;           // ns N2 quenching time
        bool Is128nm; //To simulate light at 128nm (fast and slow)

        // Calculated during setup
        double fTauTA;
        double fTauTX;
        double fMaxProbs;
        double fMaxProbt;
        double fMaxTs;
        double fMaxTt;

        double FRTime=0.0;
        double FDTime;
        double SRTime=0.0;
        double SDTime;
        
        // utility functions
        double exp_diff(double t, double tau1, double tau2) const;
        double singlet_distro(double t) const { return (exp(-t/fArSingletTime) / fArSingletTime); };
        double triplet_distro(double t) const { return exp_diff(t, fTauTA, fTauTX); };

        void initRand(CLHEP::HepRandomEngine& engine) override;
        
        double single_exp(double t, double tau2);
        double bi_exp(double t, double tau1, double tau2);
    };
}
#endif
