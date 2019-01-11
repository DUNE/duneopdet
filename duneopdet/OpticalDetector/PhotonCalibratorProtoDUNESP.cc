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
      1418., 1419., 1412., 1406., //Channels 20-23, SensL-A, DC

      //SSP103
      1421., 1., 1407., 1387., //Channels 24-27, SensL-A, DS
      1439., 1434., 1423., 1463., //Channels 28-31, SensL-A, DC
      1434., 1604., 1407., 1416., //Channels 32-35, SensL-A, DS
      
      //SSP104
      1., 1400., 1392., 1423., //Channels 36-39, SensL-A, DC
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
      782., 750., 780., 783., //Channels 132-135, MPPC, ARAPUCA
      812., 665., 737., 631., //Channels 136-139, MPPC, ARAPUCA
      673., 628., 721., 744., //Channels 140-143, MPPC, ARAPUCA
      
      //SSP302
      1387., 1418., 1., 1398., //Channels 108-111, SensL-A, DC
      1425., 1446., 1421., 1436., //Channels 112-115, SensL-A, DS
      1446., 1470., 1436., 1., //Channels 116-119, SensL-A, DC

      //SSP303
      1643., 1416., 1438., 1414., //Channels 120-123, SensL-A, DC
      1577., 1848., 1634., 1603., //Channels 124-127, SensL-A, DS
      1413., 1406., 1418., 1405., //Channels 128-131, SensL-C, DC
      
      
      //APA4
      //SSP401
      1398., 1395., 1416., 1397., //Channels 144-147, SensL-C, DC
      1404., 1432., 1428., 1418., //Channels 148-151, SensL-C, DS
      1403., 1434., 1417., 1418., //Channels 152-155, SensL-A, DC
      
      //SSP402
      1., 1419., 1443., 1410., //Channels 156-159, SensL-C, DS
      1., 1430., 1433., 1419., //Channels 160-163, SensL-C, DC
      1419., 1431., 1405., 1437., //Channels 164-167, SensL-C, DS
      
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
      1107., 1098., 1131., 946., //Channels 224-227, MPPC, DC
      
      //SSP504
      986., 1267., 1017., 1046., //Channels 228-231, MPPC, DS
      
      //SSP501
      1433., 1426., 1418., 1418., //Channels 192-195, SensL-C, DC
      
      //SSP504
      1012., 1014., 997., 1020., //Channels 232-235, MPPC, DS
      
      //SSP501
      1436., 1427., 1419., 1430., //Channels 196-199, SensL-C, DC
      1440., 1634., 1436., 1431., //Channels 200-203, SensL-C, DS
      
      //SSP504
      1029., 996., 992., 1019., //Channels 236-239, MPPC, DC
      
      //SSP502
      1439., 1469., 1455., 1429., //Channels 204-207, SensL-C DS
      
      //APA6
      //SSP601
      945., 977., 990., 924., //Channels 240-243, MPPC, DC
      963., 983., 980., 974., //Channels 244-247, MPPC, DS
      1013., 1038., 1023., 1083., //Channels 248-251, MPPC, DC
      
      //SSP602
      1017., 1023., 987., 973., //Channels 252-255, MPPC, DS
      960., 993., 1021., 992., //Channels 256-259, MPPC, DC
      
      //SSP603
      //ARAPUCA
      758., 780., 772., 628., //Channels 264-267, MPPC, ARAPUCA
      804., 637., 668., 738., //Channels 268-271, MPPC, ARAPUCA
      638., 755., 709., 831., //Channels 272-275, MPPC, ARAPUCA
      
      //SSP602
      1038., 1010., 1015., 1068., //Channels 260-263, MPPC, DC
      
      //SSP604
      1003., 987., 1009., 1022., //Channels 276-279, MPPC, DS
      1008., 1032., 981., 1045., //Channels 280-283, MPPC, DC
      987., 992., 1009., 989. //Channels 284-287, MPPC, DS
    };
    
    fSPEShifts = {
      //APA1 
      //SSP101
      0.01, 0.00, 0.01, 0.00,  //Channels 0-3, SensL-A, DS
      0.00, -0.02, -0.01, -0.01,  //Channels 4-6, SensL-A, DC
      -0.02, 0.00, 0.01, 0.00,  //Channels 8-11, SensL-A, DS
      
      //SSP102
      -0.02, -0.02, -0.02, -0.02, //Channels 12-15, SensL-A, DC
      -0.02, 0.00, -0.03, -0.01, //Channels 16-19, SensL-A, DS
      -0.04, -0.00, -0.01, -0.01, //Channels 20-23, SensL-A, DC

      //SSP103
      -0.02, 0., -0.03, -0.01, //Channels 24-27, SensL-A, DS
      0.00, -0.03, -0.02, -0.01, //Channels 28-31, SensL-A, DC
      -0.01, 0.01, -0.02, -0.04, //Channels 32-35, SensL-A, DS
      
      //SSP104
      0., -0.04, -0.01, -0.03, //Channels 36-39, SensL-A, DC
      0., 0., 0., 0., //Channels 40-43 --EMPTY
      0., 0., 0., 0., //Channels 44-47 --EMPTY

      //APA2
      //SSP201
      -0.01, 0.00, 0.01, 0.01, //Channels 48-51, SensL-A, DS
      0.00, -0.02, 0.01, -0.01, //Channels 52-55, SensL-A, DC
      -0.03, 0.00, 0., 0.00, //Channels 56-59, SensL-A, DS
      
      //SSP202
      -0.02, -0.01, 0., -0.03, //Channels 60-63, SensL-A, DC
      -0.02, 0., 0.00, -0.03, //Channels 64-67, SensL-A, DS
      0.02, -0.02, -0.03, -0.01, //Channels 68-71, SensL-A, DC
      
      //SSP203
      -0.03, 0., -0.03, 0., //Channels 72-75, SensL-A, DS
      -0.01, -0.01, -0.01, -0.01, //Channels 76-79, SensL-A, DC
      0.03, 0.01, 0., 0.01, //Channels 80-83, SensL-A, DS
      
      //SSP204
      -0.01, -0.02, -0.03, -0.03, //Channels 84-87, SensL-A, DC
      0., 0., 0., 0., //Channels 88-91 --EMPTY
      0., 0., 0., 0., //Channels 92-95 --EMPTY

      //APA3
      //SSP301
      -0.04, -0.03, 0.00, -0.01, //Channels 96-99, SensL-A, DS
      -0.00, 0., 0.00, -0.01, //Channels 100-103, SensL-A, DC
      0.01, 0.00, 0.00, 0.01, //Channels 104-107, SensL-A, DS
      
      //SSP304
      //ARAPUCA
      -0.05, -0.05, -0.05, -0.03, //Channels 132-135, MPPC, ARAPUCA
      -0.07, -0.06, -0.09, -0.02, //Channels 136-139, MPPC, ARAPUCA
      -0.08, 0.00, 0., 0., //Channels 140-143, MPPC, ARAPUCA
      
      //SSP302
      0.00, 0.00, 0., -0.02, //Channels 108-111, SensL-A, DC
      -0.02, -0.04, -0.02, -0.01, //Channels 112-115, SensL-A, DS
      -0.05, -0.07, -0.03, 0., //Channels 116-119, SensL-A, DC

      //SSP303
      -0.01, -0.03, -0.01, 0.00, //Channels 120-123, SensL-A, DC
      0.02, 0.01, 0.00, 0.01, //Channels 124-127, SensL-A, DS
      -0.01, -0.01, -0.02, 0.02, //Channels 128-131, SensL-C, DC
      
      
      //APA4
      //SSP401
      0.01, 0.00, 0.01, 0.00, //Channels 144-147, SensL-C, DC
      -0.02, -0.02, -0.01, -0.04, //Channels 148-151, SensL-C, DS
      -0.01, -0.01, -0.01, -0.03, //Channels 152-155, SensL-A, DC
      
      //SSP402
      0., 0.00, -0.02, -0.00, //Channels 156-159, SensL-C, DS
      0., -0.04, -0.02, -0.02, //Channels 160-163, SensL-C, DC
      0.02, -0.02, -0.03, -0.01, //Channels 164-167, SensL-C, DS
      
      //SSP403
      -0.02, -0.02, -0.03, -0.01, //Channels 168-171, SensL-C, DC
      -0.01, -0.03, -0.02, 0.02, //Channels 172-175, SensL-C, DS
      0.00, -0.01, -0.02, -0.01, //Channels 176-179, SensL-C, DS
      
      //SSP404
      0.00, -0.01, -0.03, 0.00, //Channels 180-183, SensL-C, DC
      0., 0., 0., 0., //Channels 184-187 --EMPTY
      0., 0., 0., 0., //Channels 188-191 --EMPTY

      //APA5
      //SSP503
      0.15, 0.12, 0.19, 0.19, //Channels 216-219, MPPC, DC
      0.17, 0.22, 0.20, 0.14, //Channels 220-223, MPPC, DS
      0.02, 0., 0., 0.23, //Channels 224-227, MPPC, DC
      
      //SSP504
      0.11, 0., 0.06, 0.04, //Channels 228-231, MPPC, DS
      
      //SSP501
      -0.02, 0.00, -0.01, -0.03, //Channels 192-195, SensL-C, DC
      
      //SSP504
      0., 0., 0., 0., //Channels 232-235, MPPC, DS
      
      //SSP501
      0.01, 0.01, -0.03, -0.01, //Channels 196-199, SensL-C, DC
      -0.01, 0.00, -0.01, -0.03, //Channels 200-203, SensL-C, DS
      
      //SSP504
      0., 0., 0., 0., //Channels 236-239, MPPC, DC
      
      //SSP502
      -0.02, -0.01, -0.01, -0.01, //Channels 204-207, SensL-C DS
      
      //APA6
      //SSP601
      0.12, 0.11, 0.14, 0.19, //Channels 240-243, MPPC, DC
      0.08, 0.08, 0.13, 0.11, //Channels 244-247, MPPC, DS
      0., 0., 0., 0., //Channels 248-251, MPPC, DC
      
      //SSP602
      0., 0., 0., 0., //Channels 252-255, MPPC, DS
      0., 0., 0., 0., //Channels 256-259, MPPC, DC
      
      //SSP603
      //ARAPUCA
      -0.04, -0.03, -0.05, -0.01, //Channels 264-267, MPPC, ARAPUCA
      -0.06, -0.06, -0.04, -0.01, //Channels 268-271, MPPC, ARAPUCA
      -0.08, -0.05, 0., -0.08, //Channels 272-275, MPPC, ARAPUCA
      
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
