// ========================================================================================
// OpHitAlg_deco.h
// This module is based on the larana/OpHitAlg.h. It has been updated to deal with   
// deconvolved signals. These are the algorithms used by OpHitFinderDeco to produce optical 
// hits. recob::OpWaveform object has been included inside the RunHitFinder_deco function.
// Added the scaling factor inside the new RunHitFinder_deco function. It scales the values 
// of the deconvolved signals before the hit finder.
// 
// @authors     : Daniele Guffanti, Maritza Delgado, Sergio Manthey Corchado
// @created     : Oct, 2022 
//=========================================================================================

#ifndef OPHITALG_DECO_H
#define OPHITALG_DECO_H

#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"

#include <vector>

namespace calib {
  class IPhotonCalibrator;
}
namespace detinfo {
  class DetectorClocksData;
}
namespace geo {
  class GeometryCore;
}
namespace pmtana {
  class PulseRecoManager;
}

namespace opdet {
  
  void RunHitFinder(std::vector<raw::OpDetWaveform> const&,
                    std::vector<recob::OpHit>&,
                    pmtana::PulseRecoManager const&,
                    pmtana::PMTPulseRecoBase const&,
                    geo::GeometryCore const&,
                    float,
                    detinfo::DetectorClocksData const&,
                    calib::IPhotonCalibrator const&,
                    bool use_start_time = false);

  void RunHitFinder_deco(std::vector<recob::OpWaveform> const&,
                    std::vector<raw::OpDetWaveform> const&,                 
                    std::vector<recob::OpHit>&,
                    pmtana::PulseRecoManager const&,
                    pmtana::PMTPulseRecoBase const&,
                    geo::GeometryCore const&,
                    float,
                    float,
                    detinfo::DetectorClocksData const&,
                    calib::IPhotonCalibrator const&,
                    bool use_start_time = false);

  void ConstructHit(float,                                   
                    int,
                    double,
                    pmtana::pulse_param const&,
                    std::vector<recob::OpHit>&,
                    detinfo::DetectorClocksData const&,
                    calib::IPhotonCalibrator const&,
                    bool use_start_time = false);

} // End opdet namespace

#endif
