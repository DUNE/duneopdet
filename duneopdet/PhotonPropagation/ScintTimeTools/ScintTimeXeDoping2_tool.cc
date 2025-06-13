////////////////////////////////////////////////////////////////////////
// Class:       ScintTimeXeDoping
// Plugin Type: tool
// File:        ScintTimeXeDoping_tool.cc ScintTimeXeDoping.h
// Description:
// Tool to simulate the xenon doping model, as described in https://arxiv.org/abs/2203.16134 and https://repository.cern/records/qjx7b-18v64
////////////////////////////////////////////////////////////////////////

#include "duneopdet/PhotonPropagation/ScintTimeTools/ScintTimeXeDoping2.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace phot
{

    //......................................................................    
    ScintTimeXeDoping2::ScintTimeXeDoping2(fhicl::ParameterSet const& pset)
    : ScintTime()
    , fXeConcentration{pset.get<double>("XeConcentration")}
    , fN2Concentration{pset.get<double>("N2Concentration")}
    , fArSingletTime  {pset.get<double>("ArSingletTime",6.0)}
    , fArTripletTime  {pset.get<double>("ArTripletTime")}
    , fXe150nmTime    {pset.get<double>("Xe150nmTime")}
    , fTauAX          {pset.get<double>("TauAX")}
    , fTauXX          {pset.get<double>("TauXX")}
    , fTauN2ArAr      {pset.get<double>("TauN2ArAr", 0)}
    , fTauN2ArXe      {pset.get<double>("TauN2ArXe", 0)}
    , Is128nm         {pset.get<bool>("ArLight",false)}
    {

        FDTime=fArSingletTime;

        fTauTA = 1./(1./fArTripletTime +   fXeConcentration/fTauAX + fN2Concentration/fTauN2ArAr);
        fTauTX  = 1./(1./fXe150nmTime   +   fXeConcentration/fTauXX + fN2Concentration/fTauN2ArXe);

        SDTime=fTauTA;

        mf::LogInfo("ScintTimeXeDoping2")
            << "Configured:" << "\n"
            << "  ArLight: " << Is128nm << "\n"
            << "  XeConcentration: " << fXeConcentration << " ppm\n"
            << "  N2Concentration: " << fN2Concentration << " ppm\n"
            << "  ArSingletTime:   " << fArSingletTime << " ns\n"
            << "  ArTripletTime:   " << fArTripletTime << " ns\n"
            << "  Xe150nmTime:     " << fXe150nmTime << " ns\n"
            << "  TauAX:           " << fTauAX << " ns\n"
            << "  TauXX:           " << fTauXX << " ns\n"
            << "  TauN2ArAr:           " << fTauN2ArAr << " ns\n"
            << "  TauN2ArXe:           " << fTauN2ArXe << " ns\n"
            << "Calculated:" << "\n"
            << "  TauTA triplet:   " << fTauTA << " ns\n"
            << "  TauTX:           " << fTauTX  << " ns\n"
            << "  fMaxProbs:       " << fMaxProbs << " \n"
            << "  fMaxProbt:       " << fMaxProbt << " \n"
            << "  fMaxTs:          " << fMaxTs << " ns\n"
            << "  fMaxTt:          " << fMaxTt << " ns\n"
            << "  FRTime:          " << FRTime << " ns\n"
            << "  FDTime:          " << FDTime << " ns\n"
            << "  SRTime:          " << SRTime << " ns\n"
            << "  SDTime:          " << SDTime << " ns\n";
    }
            
    //......................................................................    
    double ScintTimeXeDoping2::exp_diff(double t, double tau1, double tau2) const
    {
        using std::exp;
        return (- exp(-t/tau1) + exp(-t/tau2) ) / (tau2 - tau1);
    }

    double ScintTimeXeDoping2::single_exp(double t, double tau2)
    {
        return std::exp((-1.0 * t) / tau2) / tau2;
    }
    double ScintTimeXeDoping2::bi_exp(double t, double tau1, double tau2)
    {
        return (((std::exp((-1.0 * t) / tau2) * (1.0 - std::exp((-1.0 * t) / tau1))) / tau2) / tau2) * (tau1 + tau2);
    }
    
    //......................................................................    
    // Returns the time within the time distribution of the scintillation process, when the photon was created.
    // There is no fast component for the xenon light! All photons are treated as slow light.
    void ScintTimeXeDoping2::GenScintTime(bool is_fast, CLHEP::HepRandomEngine& engine)
    {
      CLHEP::RandFlat randflatscinttime{engine};
      double tau1, tau2;

      if(Is128nm) //Ar light
      {
         double tau2;
         if(is_fast == 1)  // fast scinitllation
         {
          tau2 = FDTime;
         }
         else
         {
          tau2 = SDTime;
         }
         CLHEP::RandFlat randflatscinttime{engine};
         timing = -tau2 * std::log(randflatscinttime());
         return;
       }
       else //Xe light -> Double exponential, no fast considered for Xe light
       {
         tau1=fTauTA;//rise
         tau2=fTauTX;//decay
         while (1)
         {
            auto ran1 = randflatscinttime();
            auto ran2 = randflatscinttime();
            auto d = (tau1 + tau2) / tau2;
            auto t = -tau2 * std::log(1 - ran1);
            auto g = d * single_exp(t, tau2);
            if (ran2 <= bi_exp(t, tau1, tau2) / g)
            {
               timing = t;
               return;
            }
         }
      }
    }
    
    
    void ScintTimeXeDoping2::initRand(CLHEP::HepRandomEngine& engine)
    {
        fUniformGen = std::make_unique<CLHEP::RandFlat>(engine);
    }
    
  double ScintTimeXeDoping2::fastScintTime()
  {
    if(Is128nm) //Ar light
    {
      return  -FDTime * std::log(fUniformGen->fire());
    }
    else //Xe light -> Double exponential, no fast considered for Xe light
    {

        throw art::Exception(art::errors::Configuration)
        << "ScintTimeXeDoping2_tool bad configured, it is requested to simulate fast component for Xe light "
        <<  "but the current model do not include that. Check scintYieldRatio parameter."
        << "\n";
    }
  }

  double ScintTimeXeDoping2::slowScintTime()
  {
    if(Is128nm) //Ar light
    {
      return -SDTime * std::log(fUniformGen->fire());
    }
    else //Xe light -> Double exponential, no fast considered for Xe light
    {
      double tau1, tau2;
      tau1=fTauTA;//rise
      tau2=fTauTX;//decay
      while (1)
      {
        auto ran1 = fUniformGen->fire();
        auto ran2 = fUniformGen->fire();
        auto d = (tau1 + tau2) / tau2;
        auto t = -tau2 * std::log(1 - ran1);
        auto g = d * single_exp(t, tau2);
        if (ran2 <= bi_exp(t, tau1, tau2) / g)
        {
           return t;
        }
      }
    }
  }

}
    
DEFINE_ART_CLASS_TOOL(phot::ScintTimeXeDoping2)
