/////////////////////////////////////////////////////////////////////////
// SaturationSmoothBreak_tool.cc
//
// Created: Oct 14 by Laura Paulucci (lpaulucc@fnal.gov)
//
// Saturation correction tool for a "smooth break" fit of the ToS method.
// Fit introduced by Jacob Boza
////////////////////////////////////////////////////////////////////////


#include "larana/OpticalDetector/OpHitFinder/SaturationCorrectionTool.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include <map>
#include <vector>
#include <cmath>

namespace opdet {

  class SaturationSmoothBreak : public pmtana::SaturationCorrectionTool {
  public:
    explicit SaturationSmoothBreak(fhicl::ParameterSet const& pset) {
      fCalibThrs        = pset.get<double>("calibThrs");
      fMinPlateauSize   = pset.get<int>("MinPlateauSize", 3);

      auto channelVec  = pset.get<std::vector<size_t>>("channelVec");
      auto calibParam1 = pset.get<std::vector<double>>("calibParam1");
      auto calibParam2 = pset.get<std::vector<double>>("calibParam2");
      auto calibParam3 = pset.get<std::vector<double>>("calibParam3");
      auto calibParam4 = pset.get<std::vector<double> >("calibParam4");
      auto calibParam5 = pset.get<std::vector<double> >("calibParam5");
      for (size_t i = 0; i < channelVec.size(); ++i) {
        fCalibMap[channelVec[i]] = {calibParam1[i], calibParam2[i], calibParam3[i], calibParam4[i], calibParam5[i]};
      }
    }
  private:
    double fCalibThrs;
    int fMinPlateauSize;
    std::map<size_t, std::vector<double>> fCalibMap;
    //---------------------------------------------------------------------------
    void CorrectSaturation(pmtana::Waveform_t const& wf, pmtana::pulse_param& _pulse, int ch) override
    {
      bool in_saturation = false;
      int satcounter=0, plateau_count=0;
      std::vector<int> plateau_lengths;
      double tolerance = 0.1, max_margin = 100.0; // Derivative and max allowed drop (in ADC) from the pulse peak to count as saturation; //to be optimized

      // Ensure valid indices
      if (_pulse.t_start >= _pulse.t_end || _pulse.t_end >= static_cast<int>(wf.size())) return;
      // Start loop at max(1, t_start) to prevent wf[i-1] out-of-bounds access
      size_t start_tick = std::max(static_cast<size_t>(1), static_cast<size_t>(_pulse.t_start));

      for (size_t i = start_tick; i <= static_cast<size_t>(_pulse.t_end); ++i) {//go over the full lenght of a pulse looking for saturation (slope ~ 0 near peak)
        if (std::abs((static_cast<double>(wf[i]) - static_cast<double>(wf[i - 1])) <= tolerance) && (static_cast<double>(wf[i]) >= fCalibThrs) && (static_cast<double>(wf[i]) >= (_pulse.peak - max_margin))){//the last one is just a safegard against spurious fluctuations 
          if (!in_saturation) {
            in_saturation = true;
            satcounter = 1;
          }else{ satcounter++; }
        }else{ // End of a saturated region
          if(in_saturation) {
            if(satcounter >= fMinPlateauSize) {
              plateau_count++;
              plateau_lengths.push_back(satcounter);
            }
            in_saturation = false;
            satcounter = 0;
          }
        }
      }

      // Catch trailing saturation at the end of pulse window
      if (in_saturation && satcounter >= fMinPlateauSize) { plateau_lengths.push_back(satcounter);}  

      if (!plateau_lengths.empty()) {
        //For the moment, the correction applied here is based on the "time-over-saturation (ToS)" method, where the missing area is added to the original area calculation
        //parameters of the calib curve, channel dependent, are from the "smoothbreak", in order: a, b, c, d, x0
        std::vector<double> calibvec = fCalibMap[ch];
        double power = (calibvec[2]-calibvec[1])/calibvec[3];
        for(size_t i =0; i<plateau_lengths.size(); i++){
          double ToT = plateau_lengths[i];
          double missingArea = 1 + pow(ToT/calibvec[4], calibvec[3]);
          missingArea = calibvec[0] * pow(ToT,calibvec[1]) * pow(missingArea, power);
          //for debugging 
          std::cout << "Saturated peak " << i << " ch " << ch << " old area " << _pulse.area << " sat ticks " << plateau_lengths[i] <<  " missing area " << missingArea*_pulse.peak << std::endl;
          _pulse.area += missingArea*_pulse.peak;
        }
      }
    }

  };
}

DEFINE_ART_CLASS_TOOL(opdet::SaturationSmoothBreak)

