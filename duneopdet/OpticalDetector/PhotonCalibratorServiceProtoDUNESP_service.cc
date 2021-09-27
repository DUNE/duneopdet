/**
 * Required minimal implementation file for calibrator service
 * which only returns a provider.
 */
#include "dune/OpticalDetector/PhotonCalibratorServiceProtoDUNESP.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

DEFINE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorServiceProtoDUNESP, 
                                  calib::IPhotonCalibratorService)
