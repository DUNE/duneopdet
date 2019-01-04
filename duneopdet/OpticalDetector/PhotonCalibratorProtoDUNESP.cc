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
      1411., 1416., 1405., 1436.,  //Channels 0-3, SensL-A, DS
      1437., 1442., 1432., 1436.,  //Channels 4-6, SensL-A, DC
      1409., 1421., 1423., 1445.,  //Channels 8-11, SensL-A, DS
      
      //SSP102
      1427., 1400., 1419., 1421., //Channels 12-15, SensL-A, DC
      1410., 1402., 1451., 1637., //Channels 16-19, SensL-A, DS
      1419., 1419., 1412., 1406., //Channels 20-23, SensL-A, DC

      //SSP103
      1421., 1., 1407., 1387., //Channels 24-27, SensL-A, DS
      1439., 1434., 1423., 1463., //Channels 28-31, SensL-A, DC
      1434., 1604., 1407., 1416., //Channels 32-35, SensL-A, DS
      
      //SSP104
      1., 1400., 1392., 1., //Channels 36-39, SensL-A, DC
      1., 1., 1., 1., //Channels 40-43 --EMPTY
      1., 1., 1., 1., //Channels 44-47 --EMPTY

      //APA2
      //SSP201
      1410., 1824., 1386., 1841., //Channels 48-51, SensL-A, DS
      1411., 1439., 1408., 1435., //Channels 52-55, SensL-A, DC
      1389., 1392., 1., 1368., //Channels 56-59, SensL-A, DS
      
      //SSP202
      1427., 1409., 1., 1459., //Channels 60-63, SensL-A, DC
      1430., 1., 1401., 1456., //Channels 64-67, SensL-A, DS
      1419., 1431., 1405., 1437., //Channels 68-71, SensL-A, DC
      
      //SSP203
      1433., 1., 1421., 1., //Channels 72-75, SensL-A, DS
      1435., 1643., 1425., 1428., //Channels 76-79, SensL-A, DC
      1378., 1372., 1., 1364., //Channels 80-83, SensL-A, DS
      
      //SSP204
      1417., 1430., 1434., 1432., //Channels 84-87, SensL-A, DC
      1., 1., 1., 1., //Channels 88-91 --EMPTY
      1., 1., 1., 1., //Channels 92-95 --EMPTY

      //APA3
      //SSP301
      1396., 1432., 1392., 1394., //Channels 96-99, SensL-A, DS
      1392., 1., 1414., 1413., //Channels 100-103, SensL-A, DC
      1624., 1399., 1624., 1621., //Channels 104-107, SensL-A, DS
      
      //SSP304
      //ARAPUCA
      745., 668., 633., 632., //Channels 132-135, MPPC, ARAPUCA
      714., 531., 649., 582., //Channels 136-139, MPPC, ARAPUCA
      1., 1., 1., 1., //Channels 140-143, MPPC, ARAPUCA
      
      //SSP302
      1386., 1418., 1., 1398., //Channels 108-111, SensL-A, DC
      1425., 1446., 1421., 1436., //Channels 112-115, SensL-A, DS
      1446., 1470., 1436., 1., //Channels 116-119, SensL-A, DC

      //SSP303
      1643., 1416., 1438., 1414., //Channels 120-123, SensL-A, DC
      1577., 1848., 1634., 1603., //Channels 124-127, SensL-A, DS
      1413., 1406., 1418., 1405., //Channels 128-131, SensL-C, DC
      
      
      //APA4
      //SSP401
      1398., 1395., 1416., 1397., //Channels 144-147, SensL-C, DC
      1396., 1415., 1421., 1397., //Channels 148-151, SensL-C, DS
      1403., 1434., 1417., 1418., //Channels 152-155, SensL-A, DC
      
      //SSP402
      1., 1419., 1443., 1410., //Channels 156-159, SensL-C, DS
      1., 1430., 1433., 1419., //Channels 160-163, SensL-C, DC
      1401., 1409., 1424., 1423., //Channels 164-167, SensL-C, DS
      
      //SSP403
      1417., 1445., 1437., 1418., //Channels 168-171, SensL-C, DC
      1662., 1450., 1424., 1629., //Channels 172-175, SensL-C, DS
      1815., 1424., 1402., 1428., //Channels 176-179, SensL-C, DS
      
      //SSP404
      1626., 1442., 1413., 1417., //Channels 180-183, SensL-C, DC
      1., 1., 1., 1., //Channels 184-187 --EMPTY
      1., 1., 1., 1., //Channels 188-191 --EMPTY

      //APA5
      //SSP503
      908., 950., 891., 876., //Channels 216-219, MPPC, DC
      890., 871., 876., 908., //Channels 220-223, MPPC, DS
      1107., 1., 1., 946., //Channels 224-227, MPPC, DC
      
      //SSP504
      986., 1., 1017., 1046., //Channels 228-231, MPPC, DS
      
      //SSP501
      1433., 1426., 1418., 1418., //Channels 192-195, SensL-C, DC
      
      //SSP504
      1., 1., 1., 1., //Channels 232-235, MPPC, DS
      
      //SSP501
      1436., 1427., 1419., 1430., //Channels 196-199, SensL-C, DC
      1440., 1634., 1436., 1431., //Channels 200-203, SensL-C, DS
      
      //SSP504
      1., 1., 1., 1., //Channels 236-239, MPPC, DC
      
      //SSP502
      1439., 1469., 1455., 1429., //Channels 204-207, SensL-C DS
      
      //APA6
      //SSP601
      945., 977., 990., 924., //Channels 240-243, MPPC, DC
      963., 983., 980., 974., //Channels 244-247, MPPC, DS
      1., 1., 1., 1., //Channels 248-251, MPPC, DC
      
      //SSP602
      1., 1., 1., 1., //Channels 252-255, MPPC, DS
      1., 1., 1., 1., //Channels 256-259, MPPC, DC
      
      //SSP603
      //ARAPUCA
      720., 766., 739., 540., //Channels 264-267, MPPC, ARAPUCA
      872., 570., 614., 703., //Channels 268-271, MPPC, ARAPUCA
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
      0.007, 0.003, 0.014, 0.003,  //Channels 0-3, SensL-A, DS
      -0.003, -0.017, -0.010, -0.012,  //Channels 4-6, SensL-A, DC
      -0.021, -0.015, -0.021, -0.023,  //Channels 8-11, SensL-A, DS
      
      //SSP102
      -0.020, -0.022, -0.014, -0.015, //Channels 12-15, SensL-A, DC
      -0.017, 0.000, -0.025, -0.007, //Channels 16-19, SensL-A, DS
      -0.040, -0.004, -0.012, -0.011, //Channels 20-23, SensL-A, DC

      //SSP103
      -0.023, 0., -0.009, -0.008, //Channels 24-27, SensL-A, DS
      -0.005, -0.029, -0.022, -0.012, //Channels 28-31, SensL-A, DC
      -0.015, 0.005, -0.022, -0.040, //Channels 32-35, SensL-A, DS
      
      //SSP104
      0., -0.036, -0.005, 0., //Channels 36-39, SensL-A, DC
      0., 0., 0., 0., //Channels 40-43 --EMPTY
      0., 0., 0., 0., //Channels 44-47 --EMPTY

      //APA2
      //SSP201
      -0.009, 0.001, 0.001, 0.005, //Channels 48-51, SensL-A, DS
      -0.046, -0.013, 0.008, -0.012, //Channels 52-55, SensL-A, DC
      -0.028, -0.002, 0., -0.001, //Channels 56-59, SensL-A, DS
      
      //SSP202
      -0.018, -0.007, 0., -0.033, //Channels 60-63, SensL-A, DC
      -0.019, 0., -0.004, -0.028, //Channels 64-67, SensL-A, DS
      0.023, -0.016, -0.026, -0.010, //Channels 68-71, SensL-A, DC
      
      //SSP203
      -0.030, 0., -0.027, 0., //Channels 72-75, SensL-A, DS
      -0.009, -0.012, -0.009, -0.010, //Channels 76-79, SensL-A, DC
      0.030, 0.010, 0., 0.008, //Channels 80-83, SensL-A, DS
      
      //SSP204
      -0.009, -0.016, -0.031, -0.028, //Channels 84-87, SensL-A, DC
      0., 0., 0., 0., //Channels 88-91 --EMPTY
      0., 0., 0., 0., //Channels 92-95 --EMPTY

      //APA3
      //SSP301
      -0.036, -0.026, -0.003, -0.008, //Channels 96-99, SensL-A, DS
      -0.003, 0., 0.000, -0.006, //Channels 100-103, SensL-A, DC
      0.012, -0.004, 0.004, 0.005, //Channels 104-107, SensL-A, DS
      
      //SSP304
      //ARAPUCA
      -0.009, 0.199, 0.342, 0.490, //Channels 132-135, MPPC, ARAPUCA
      0.224, 0.439, 0.263, 0.232, //Channels 136-139, MPPC, ARAPUCA
      0., 0., 0., 0., //Channels 140-143, MPPC, ARAPUCA
      
      //SSP302
      0.005, 0.002, 0., -0.024, //Channels 108-111, SensL-A, DC
      -0.024, -0.036, -0.018, -0.013, //Channels 112-115, SensL-A, DS
      -0.047, -0.071, -0.032, 0., //Channels 116-119, SensL-A, DC

      //SSP303
      -0.013, -0.028, -0.008, -0.004, //Channels 120-123, SensL-A, DC
      0.020, 0.011, -0.001, 0.009, //Channels 124-127, SensL-A, DS
      -0.009, -0.009, -0.016, 0.021, //Channels 128-131, SensL-C, DC
      
      
      //APA4
      //SSP401
      0.010, 0.001, 0.008, 0.001, //Channels 144-147, SensL-C, DC
      -0.027, -0.002, -0.033, -0.019, //Channels 148-151, SensL-C, DS
      -0.007, -0.013, -0.013, -0.032, //Channels 152-155, SensL-A, DC
      
      //SSP402
      0., -0.005, -0.016, -0.004, //Channels 156-159, SensL-C, DS
      0., -0.036, -0.018, -0.021, //Channels 160-163, SensL-C, DC
      -0.011, -0.001, -0.011, -0.007, //Channels 164-167, SensL-C, DS
      
      //SSP403
      -0.025, -0.024, -0.025, -0.014, //Channels 168-171, SensL-C, DC
      -0.011, -0.028, -0.025, 0.024, //Channels 172-175, SensL-C, DS
      0.005, -0.010, -0.021, -0.013, //Channels 176-179, SensL-C, DS
      
      //SSP404
      0.004, -0.011, -0.028, -0.003, //Channels 180-183, SensL-C, DC
      0., 0., 0., 0., //Channels 184-187 --EMPTY
      0., 0., 0., 0., //Channels 188-191 --EMPTY

      //APA5
      //SSP503
      0.152, 0.115, 0.189, 0.185, //Channels 216-219, MPPC, DC
      0.167, 0.225, 0.196, 0.144, //Channels 220-223, MPPC, DS
      0.017, 0., 0., 0.230, //Channels 224-227, MPPC, DC
      
      //SSP504
      0.113, 0., 0.062, 0.035, //Channels 228-231, MPPC, DS
      
      //SSP501
      -0.018, -0.003, -0.011, -0.025, //Channels 192-195, SensL-C, DC
      
      //SSP504
      0., 0., 0., 0., //Channels 232-235, MPPC, DS
      
      //SSP501
      0.010, 0.014, -0.029, -0.011, //Channels 196-199, SensL-C, DC
      -0.015, 0.003, -0.006, -0.034, //Channels 200-203, SensL-C, DS
      
      //SSP504
      0., 0., 0., 0., //Channels 236-239, MPPC, DC
      
      //SSP502
      -0.019, -0.013, -0.013, -0.008, //Channels 204-207, SensL-C DS
      
      //APA6
      //SSP601
      0.117, 0.110, 0.136, 0.191, //Channels 240-243, MPPC, DC
      0.082, 0.081, 0.129, 0.115, //Channels 244-247, MPPC, DS
      0., 0., 0., 0., //Channels 248-251, MPPC, DC
      
      //SSP602
      0., 0., 0., 0., //Channels 252-255, MPPC, DS
      0., 0., 0., 0., //Channels 256-259, MPPC, DC
      
      //SSP603
      //ARAPUCA
      0.070, 0.055, 0.042, 0.657, //Channels 264-267, MPPC, ARAPUCA
      -0.092, 0.233, 0.248, 0.140, //Channels 268-271, MPPC, ARAPUCA
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
