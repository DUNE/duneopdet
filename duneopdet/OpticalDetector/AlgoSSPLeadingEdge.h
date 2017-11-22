//=========================================================
// AlgoSSPLeadingEdge.h
// This class provides the digital output for the 
// Leading Edge Discrimintator, used in SSP module.
//
// Vitor Luzio, UFABC, 2017
// Based on AlgoSiPM.cxx
//=========================================================



#ifndef AlgoSSPLeadingEdge_H
#define AlgoSSPLeadingEdge_H 1

#include "fhiclcpp/ParameterSet.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include <vector>

namespace pmtana {

  class AlgoSSPLeadingEdge : public PMTPulseRecoBase {

  public:

    // Default constructor
    AlgoSSPLeadingEdge(const fhicl::ParameterSet &pset,const std::string name="AlgoSSPLeadingEdge");
   
    // Default destructor
    ~AlgoSSPLeadingEdge();
   
    // Implementation of PMTPulseRecoBase::Reset() method
    void Reset();
  protected:

    bool RecoPulse( const pmtana::Waveform_t&,
                    const pmtana::PedestalMean_t&,
                    const pmtana::PedestalSigma_t&     );
  
    // A variable holder for a user-defined absolute ADC threshold value
    double _adc_thres;
    
    // Minimum width for a hit to be recorded
    int _min_width;
    
    // Start recording hit information after this threshold is reached
    double _2nd_thres;
    
    // Use this pedestal instead of the one given by the pedestal algorithm
    double _pedestal;

  };

}

#endif

