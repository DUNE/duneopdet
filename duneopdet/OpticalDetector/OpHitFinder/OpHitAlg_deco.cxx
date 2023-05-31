// ========================================================================================
// OpHitAlg_deco.cxx
// This module is based on the larana/OpHitAlg.cxx. It has been updated to deal with   
// deconvolved signals. These are the algorithms used by OpHitFinderDeco to produce optical 
// hits. recob::OpWaveform object has been included inside the RunHitFinder_deco function.
// Added the scaling factor inside the new RunHitFinder_deco function. It scales the values 
// of the deconvolved signals before the hit finder.
// 
// @authors     : Daniele Guffanti, Maritza Delgado, Sergio Manthey Corchado
// @created     : Oct, 2022 
//=========================================================================================

#include "OpHitAlg_deco.h"

#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>

namespace opdet {

 
  void RunHitFinder(std::vector<raw::OpDetWaveform> const& opDetWaveformVector,
                    std::vector<recob::OpHit>& hitVector,
                    pmtana::PulseRecoManager const& pulseRecoMgr,
                    pmtana::PMTPulseRecoBase const& threshAlg,
                    geo::GeometryCore const& geometry,
                    float hitThreshold,
                    detinfo::DetectorClocksData const& clocksData,
                    calib::IPhotonCalibrator const& calibrator,
                    bool use_start_time)
  {

    for (auto const& waveform : opDetWaveformVector) {

      const int channel = static_cast<int>(waveform.ChannelNumber());

      if (!geometry.IsValidOpChannel(channel)) {
        mf::LogError("OpHitFinder")
          << "Error! unrecognized channel number " << channel << ". Ignoring pulse";
        continue;
      }

      pulseRecoMgr.Reconstruct(waveform);

      // Get the result
      auto const& pulses = threshAlg.GetPulses();

      const double timeStamp = waveform.TimeStamp();

      for (auto const& pulse : pulses)
        ConstructHit(hitThreshold,
                     channel,
                     timeStamp,
                     pulse,
                     hitVector,
                     clocksData,
                     calibrator,
                     use_start_time);
    }
  }

  //----------------------------------------------------------------------------
  void RunHitFinder_deco(std::vector<recob::OpWaveform>const& opWaveformVector,
                    std::vector<raw::OpDetWaveform> const& opDetWaveformVector,
                    std::vector<recob::OpHit>& hitVector,
                    pmtana::PulseRecoManager const& pulseRecoMgr,
                    pmtana::PMTPulseRecoBase const& threshAlg,
                    geo::GeometryCore const& geometry,
                    float hitThreshold,
                    float scale,  //It scales the values of the deconvolved signals.
                    detinfo::DetectorClocksData const& clocksData,
                    calib::IPhotonCalibrator const& calibrator,
                    bool use_start_time)
  {
      
     for (int i=0; i< int (opWaveformVector.size()); i++){
        recob::OpWaveform deco_waveform=opWaveformVector.at(i);
        int channel = static_cast<int>(deco_waveform.Channel());
        //The raw timestamp is used
        raw::OpDetWaveform waveform=opDetWaveformVector.at(i);
        const double timeStamp = waveform.TimeStamp();
       
        if (!geometry.IsValidOpChannel(channel)) {
          mf::LogError("OpHitFinder")
          << "Error! unrecognized channel number " << channel << ". Ignoring pulse";
         continue;
        }
      
       // Loop to convert from float to short
       std::vector<short int> short_deco_waveform; //vector used to convert from float to short.
       for (unsigned int i_tick=0; i_tick < deco_waveform.Signal().size(); ++i_tick)
       {
       short_deco_waveform.emplace_back(static_cast<short int>(scale*deco_waveform.Signal().at(i_tick)));
       }
         
       pulseRecoMgr.Reconstruct(short_deco_waveform);
      
      // Get the result
      auto const& pulses = threshAlg.GetPulses();
      for (auto const& pulse : pulses)
        ConstructHit(hitThreshold,
                     channel,
                     timeStamp,
                     pulse,
                     hitVector,
                     clocksData,
                     calibrator,
                     use_start_time);
     }
 }
  
  //----------------------------------------------------------------------------
  void ConstructHit(float hitThreshold,
                    int channel,
                    double timeStamp,
                    pmtana::pulse_param const& pulse,
                    std::vector<recob::OpHit>& hitVector,
                    detinfo::DetectorClocksData const& clocksData,
                    calib::IPhotonCalibrator const& calibrator,
                    bool use_start_time)
  {

    if (pulse.peak < hitThreshold) return;

     double absTime = timeStamp + clocksData.OpticalClock().TickPeriod() * (use_start_time ? pulse.t_start : pulse.t_max);
     double relTime = absTime - clocksData.TriggerTime();
     double startTime = timeStamp + clocksData.OpticalClock().TickPeriod() * pulse.t_start - clocksData.TriggerTime();
     double riseTime = clocksData.OpticalClock().TickPeriod() * pulse.t_rise;
     int frame = clocksData.OpticalClock().Frame(timeStamp);
     double PE = 0.0;

     if (calibrator.UseArea())   
       PE = calibrator.PE(pulse.area, channel);
     else    
      PE = calibrator.PE(pulse.peak, channel);   

     double width = (pulse.t_end - pulse.t_start) * clocksData.OpticalClock().TickPeriod();

     hitVector.emplace_back(channel,
                           relTime,
                           absTime,
                           startTime,
                           riseTime,
                           frame,
                           width,
                           pulse.area,
                           pulse.peak,
                           PE,
                           0.0);
  }
  

} // End namespace opdet
