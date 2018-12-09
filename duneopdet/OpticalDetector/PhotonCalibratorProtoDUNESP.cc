// dunetpc includes
#include "dune/OpticalDetector/PhotonCalibratorProtoDUNESP.h"

// LArSoft Includes
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
   
#include <vector>


namespace calib {

  PhotonCalibratorProtoDUNESP::PhotonCalibratorProtoDUNESP(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  {
    // Get the geometry service for getting max number of channels
    //auto const& geometry(*lar::providerFrom< geo::Geometry >());

    // Initialize the SPE vectors with default values
    //for (unsigned int channel = 0; channel < geometry.MaxOpChannel(); channel++) {
    // fSPESizes.push_back(1.);
    // fSPEShifts.push_back(0.);
    //}
    fSPESizes={ 
      //APA1 
      //SSP101
      1766., 1860., 1852., 1875.,  //Channels 0-3, SensL-A, DS
      1859., 1845., 1846., 1838.,  //Channels 4-6, SensL-A, DC
      1827., 1860., 1831., 1853.,  //Channels 8-11, SensL-A, DS
      
      //SSP102
      1845., 1803., 1829., 1815., //Channels 12-15, SensL-A, DC
      1829., 1835., 1859., 1919., //Channels 16-19, SensL-A, DS
      1827., 1816., 1818., 1811., //Channels 20-23, SensL-A, DC

      //SSP103
      1836., 1., 1819., 1807., //Channels 24-27, SensL-A, DS
      1851., 1859., 1827., 1862., //Channels 28-31, SensL-A, DC
      1846., 1866., 1814., 1811., //Channels 32-35, SensL-A, DS
      
      //SSP104
      1., 1785., 1805., 1., //Channels 36-39, SensL-A, DC
      1., 1., 1., 1., //Channels 40-43 --EMPTY
      1., 1., 1., 1., //Channels 44-47 --EMPTY

      //APA2
      //SSP201
      1818., 1675., 1813., 2215., //Channels 48-51, SensL-A, DS
      1825., 1857., 1827., 1848., //Channels 52-55, SensL-A, DC
      1627., 1826., 1., 1819., //Channels 56-59, SensL-A, DS
      
      //SSP202
      1849., 1828., 1., 1889., //Channels 60-63, SensL-A, DC
      1840., 1., 1817., 1853., //Channels 64-67, SensL-A, DS
      1826., 1851., 1809., 1860., //Channels 68-71, SensL-A, DC
      
      //SSP203
      1840., 1., 1848., 1., //Channels 72-75, SensL-A, DS
      1846., 2100., 1848., 1825., //Channels 76-79, SensL-A, DC
      1830., 1823., 1., 1812., //Channels 80-83, SensL-A, DS
      
      //SSP204
      1832., 1846., 1846., 1854., //Channels 84-87, SensL-A, DC
      1., 1., 1., 1., //Channels 88-91 --EMPTY
      1., 1., 1., 1., //Channels 92-95 --EMPTY

      //APA3
      //SSP301
      1810., 1859., 1829., 1838., //Channels 96-99, SensL-A, DS
      1817., 1., 1840., 1840., //Channels 100-103, SensL-A, DC
      1838., 1857., 1993., 1929., //Channels 104-107, SensL-A, DS
      
      //SSP304
      //ARAPUCA
      1., 1., 1., 1., //Channels 132-135, MPPC, ARAPUCA
      1., 1., 1., 1., //Channels 136-139, MPPC, ARAPUCA
      1., 1., 1., 1., //Channels 140-143, MPPC, ARAPUCA
      
      //SSP302
      1824., 1844., 1., 1828., //Channels 108-111, SensL-A, DC
      1843., 1853., 1828., 1842., //Channels 112-115, SensL-A, DS
      1831., 1855., 1841., 1., //Channels 116-119, SensL-A, DC

      //SSP303
      2108., 1820., 1859., 1828., //Channels 120-123, SensL-A, DC
      2023., 2195., 2112., 2064., //Channels 124-127, SensL-A, DS
      1838., 1827., 1823., 1832., //Channels 128-131, SensL-C, DC
      
      
      //APA4
      //SSP401
      1845., 1819., 1848., 1819., //Channels 144-147, SensL-C, DC
      1830., 1837., 1842., 1840., //Channels 148-151, SensL-C, DS
      1806., 1842., 1823., 1825., //Channels 152-155, SensL-A, DC
      
      //SSP402
      1., 1839., 1860., 1837., //Channels 156-159, SensL-C, DS
      1., 1841., 1835., 1812., //Channels 160-163, SensL-C, DC
      1827., 1816., 1842., 1837., //Channels 164-167, SensL-C, DS
      
      //SSP403
      1831., 1844., 1835., 1825., //Channels 168-171, SensL-C, DC
      2112., 1876., 1848., 2102., //Channels 172-175, SensL-C, DS
      1., 1838., 1808., 1817., //Channels 176-179, SensL-C, DS
      
      //SSP404
      2049., 1841., 1828., 1798., //Channels 180-183, SensL-C, DC
      1., 1., 1., 1., //Channels 184-187 --EMPTY
      1., 1., 1., 1., //Channels 188-191 --EMPTY

      //APA5
      //SSP503
      1., 1., 1., 1., //Channels 216-219, MPPC, DC
      1., 1., 1., 1., //Channels 220-223, MPPC, DS
      1., 1., 1., 1., //Channels 224-227, MPPC, DC
      
      //SSP504
      1., 1., 1., 1., //Channels 228-231, MPPC, DS
      
      //SSP501
      1835., 1838., 1827., 1822., //Channels 192-195, SensL-C, DC
      
      //SSP504
      1., 1., 1., 1., //Channels 232-235, MPPC, DS
      
      //SSP501
      1852., 1825., 1831., 1836., //Channels 196-199, SensL-C, DC
      1848., 2058., 1862., 1841., //Channels 200-203, SensL-C, DS
      
      //SSP504
      1., 1., 1., 1., //Channels 236-239, MPPC, DC
      
      //SSP502
      1858., 1868., 1863., 1813., //Channels 204-207, SensL-C DS
      
      //APA6
      //SSP601
      1., 1., 1., 1., //Channels 240-243, MPPC, DC
      1., 1., 1., 1., //Channels 244-247, MPPC, DS
      1., 1., 1., 1., //Channels 248-251, MPPC, DC
      
      //SSP602
      1., 1., 1., 1., //Channels 252-255, MPPC, DS
      1., 1., 1., 1., //Channels 256-259, MPPC, DC
      
      //SSP603
      //ARAPUCA
      1., 1., 1., 1., //Channels 264-267, MPPC, ARAPUCA
      1., 1., 1., 1., //Channels 268-271, MPPC, ARAPUCA
      1., 1., 1., 1., //Channels 272-275, MPPC, ARAPUCA
      
      //SSP602
      1., 1., 1., 1., //Channels 260-263, MPPC, DC
      
      //SSP604
      1., 1., 1., 1., //Channels 276-279, MPPC, DS
      1., 1., 1., 1., //Channels 280-283, MPPC, DC
      1., 1., 1., 1. //Channels 284-287, MPPC, DS
    };
    
    fSPEShifts = {
      //APA1 
      //SSP101
      0.067, 0.048, 0.050, 0.054,  //Channels 0-3, SensL-A, DS
      0.057, 0.058, 0.058, 0.059,  //Channels 4-6, SensL-A, DC
      0.051, 0.042, 0.052, 0.055,  //Channels 8-11, SensL-A, DS
      
      //SSP102
      0.051, 0.052, 0.056, 0.059, //Channels 12-15, SensL-A, DC
      0.048, 0.051, 0.043, 0.165, //Channels 16-19, SensL-A, DS
      0.031, 0.071, 0.058, 0.060, //Channels 20-23, SensL-A, DC

      //SSP103
      0.053, 0., 0.062, 0.050, //Channels 24-27, SensL-A, DS
      0.068, 0.039, 0.057, 0.071, //Channels 28-31, SensL-A, DC
      0.057, 0.194, 0.053, 0.052, //Channels 32-35, SensL-A, DS
      
      //SSP104
      0., 0.056, 0.052, 0., //Channels 36-39, SensL-A, DC
      0., 0., 0., 0., //Channels 40-43 --EMPTY
      0., 0., 0., 0., //Channels 44-47 --EMPTY

      //APA2
      //SSP201
      0.062, 0.327, 0.065, 0.135, //Channels 48-51, SensL-A, DS
      0.036, 0.044, 0.067, 0.062, //Channels 52-55, SensL-A, DC
      0.149, 0.043, 0., 0.041, //Channels 56-59, SensL-A, DS
      
      //SSP202
      0.041, 0.045, 0., 0.033, //Channels 60-63, SensL-A, DC
      0.058, 0., 0.060, 0.054, //Channels 64-67, SensL-A, DS
      0.085, 0.047, 0.050, 0.051, //Channels 68-71, SensL-A, DC
      
      //SSP203
      0.045, 0., 0.040, 0., //Channels 72-75, SensL-A, DS
      0.061, 0.067, 0.049, 0.068, //Channels 76-79, SensL-A, DC
      0.063, 0.050, 0., 0.049, //Channels 80-83, SensL-A, DS
      
      //SSP204
      0.057, 0.051, 0.042, 0.042, //Channels 84-87, SensL-A, DC
      0., 0., 0., 0., //Channels 88-91 --EMPTY
      0., 0., 0., 0., //Channels 92-95 --EMPTY

      //APA3
      //SSP301
      0.049, 0.046, 0.049, 0.039, //Channels 96-99, SensL-A, DS
      0.056, 0., 0.059, 0.055, //Channels 100-103, SensL-A, DC
      0.217, 0.045, 0.125, 0.158, //Channels 104-107, SensL-A, DS
      
      //SSP304
      //ARAPUCA
      0., 0., 0., 0., //Channels 132-135, MPPC, ARAPUCA
      0., 0., 0., 0., //Channels 136-139, MPPC, ARAPUCA
      0., 0., 0., 0., //Channels 140-143, MPPC, ARAPUCA
      
      //SSP302
      0.048, 0.061, 0., 0.033, //Channels 108-111, SensL-A, DC
      0.048, 0.041, 0.056, 0.056, //Channels 112-115, SensL-A, DS
      0.052, 0.037, 0.053, 0., //Channels 116-119, SensL-A, DC

      //SSP303
      0.052, 0.048, 0.052, 0.064, //Channels 120-123, SensL-A, DC
      0.085, 0.149, 0.062, 0.065, //Channels 124-127, SensL-A, DS
      0.045, 0.050, 0.054, 0.065, //Channels 128-131, SensL-C, DC
      
      
      //APA4
      //SSP401
      0.060, 0.058, 0.060, 0.053, //Channels 144-147, SensL-C, DC
      0.042, 0.051, 0.058, 0.031, //Channels 148-151, SensL-C, DS
      0.071, 0.058, 0.057, 0.040, //Channels 152-155, SensL-A, DC
      
      //SSP402
      0., 0.053, 0.052, 0.053, //Channels 156-159, SensL-C, DS
      0., 0.033, 0.051, 0.058, //Channels 160-163, SensL-C, DC
      0.049, 0.061, 0.049, 0.055, //Channels 164-167, SensL-C, DS
      
      //SSP403
      0.037, 0.053, 0.049, 0.059, //Channels 168-171, SensL-C, DC
      0.061, 0.038, 0.039, 0.069, //Channels 172-175, SensL-C, DS
      0., 0.054, 0.054, 0.072, //Channels 176-179, SensL-C, DS
      
      //SSP404
      0.083, 0.061, 0.041, 0.083, //Channels 180-183, SensL-C, DC
      0., 0., 0., 0., //Channels 184-187 --EMPTY
      0., 0., 0., 0., //Channels 188-191 --EMPTY

      //APA5
      //SSP503
      0., 0., 0., 0., //Channels 216-219, MPPC, DC
      0., 0., 0., 0., //Channels 220-223, MPPC, DS
      0., 0., 0., 0., //Channels 224-227, MPPC, DC
      
      //SSP504
      0., 0., 0., 0., //Channels 228-231, MPPC, DS
      
      //SSP501
      0.050, 0.062, 0.055, 0.044, //Channels 192-195, SensL-C, DC
      
      //SSP504
      0., 0., 0., 0., //Channels 232-235, MPPC, DS
      
      //SSP501
      0.067, 0.082, 0.041, 0.063, //Channels 196-199, SensL-C, DC
      0.051, 0.087, 0.048, 0.035, //Channels 200-203, SensL-C, DS
      
      //SSP504
      0., 0., 0., 0., //Channels 236-239, MPPC, DC
      
      //SSP502
      0.050, 0.061, 0.060, 0.079, //Channels 204-207, SensL-C DS
      
      //APA6
      //SSP601
      0., 0., 0., 0., //Channels 240-243, MPPC, DC
      0., 0., 0., 0., //Channels 244-247, MPPC, DS
      0., 0., 0., 0., //Channels 248-251, MPPC, DC
      
      //SSP602
      0., 0., 0., 0., //Channels 252-255, MPPC, DS
      0., 0., 0., 0., //Channels 256-259, MPPC, DC
      
      //SSP603
      //ARAPUCA
      0., 0., 0., 0., //Channels 264-267, MPPC, ARAPUCA
      0., 0., 0., 0., //Channels 268-271, MPPC, ARAPUCA
      0., 0., 0., 0., //Channels 272-275, MPPC, ARAPUCA
      
      //SSP602
      0., 0., 0., 0., //Channels 260-263, MPPC, DC
      
      //SSP604
      0., 0., 0., 0., //Channels 276-279, MPPC, DS
      0., 0., 0., 0., //Channels 280-283, MPPC, DC
      0., 0., 0., 0. //Channels 284-287, MPPC, DS
    };


  }
  double PhotonCalibratorProtoDUNESP::PE(double adcs, int opchannel) const
  {
    return adcs/fSPESizes[opchannel] + fSPEShifts[opchannel];
  }

}

//DEFINE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorProtoDUNESP, calib::IPhotonCalibrator)
