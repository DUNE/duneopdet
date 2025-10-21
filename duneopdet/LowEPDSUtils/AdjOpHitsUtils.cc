#include "AdjOpHitsUtils.h"

using namespace producer;

namespace solar
{
  AdjOpHitsUtils::AdjOpHitsUtils(fhicl::ParameterSet const &p)
      : fOpHitTimeVariable(p.get<std::string>("OpHitTimeVariable", "PeakTime")), // Variable to use for time sorting ("StartTime" or "PeakTime")
        fOpFlashAlgoNHit(p.get<int>("OpFlashAlgoNHit", 3)),       // Minimum number of OpHits in a cluster to consider it for flash creation.
        fOpFlashAlgoMinTime(p.get<float>("OpFlashAlgoMinTime", 0.010)), // Negative time window to look for adj. OpHits. Default for HD 10 ns [0.6 tick]
        fOpFlashAlgoMaxTime(p.get<float>("OpFlashAlgoMaxTime", 0.016)), // Positive time window to look for adj. OpHits. Default for HD 16 ns [1 tick]
        fOpFlashAlgoRad(p.get<float>("OpFlashAlgoRad", 300.0)),     // Radius to look for adj. OpHits in [cm] units.
        fOpFlashAlgoPE(p.get<float>("OpFlashAlgoPE", 1.5)),         // Minimum PE of OpHit to consider it for flash creation.
        fOpFlashAlgoTriggerPE(p.get<float>("OpFlashAlgoTriggerPE", 1.5)), // Minimum PE of OpHit to consider it as a trigger for flash creation.
        fOpFlashAlgoHotVertexThld(p.get<float>("OpFlashAlgoHotVertexThld", 0.3)), // Fraction of MaxPE to consider an OpHit as part of the flash center calculation.
        fOpFlashAlgoHitDuplicates(p.get<bool>("OpFlashAlgoHitDuplicates", true)), // Whether to allow OpHits to be part of multiple clusters.
        fXACathodeX(p.get<float>("XACathodeX", -327.5)),
        fXAMembraneY(p.get<float>("XAMembraneY", 743.302)),
        fXAFinalCapZ(p.get<float>("XAFinalCapZ", 2188.38)),
        fXAStartCapZ(p.get<float>("XAStartCapZ", -96.5))
  {
  }
  void AdjOpHitsUtils::MakeFlashVector(std::vector<FlashInfo> &FlashVec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &OpHitClusters, art::Event const &evt)
  {
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto TickPeriod = clockData.OpticalClock().TickPeriod();
    if (fOpFlashAlgoHotVertexThld > 1)
    {
      ProducerUtils::PrintInColor("Hot vertex threshold must be between 0 and 1", ProducerUtils::GetColor("red"), "Error");
      return;
    }
    for (std::vector<art::Ptr<recob::OpHit>> Cluster : OpHitClusters)
    {
      if (!Cluster.empty())
      {
        if (fOpHitTimeVariable == "StartTime")
          std::stable_sort(Cluster.begin(), Cluster.end(), [](art::Ptr<recob::OpHit> a, art::Ptr<recob::OpHit> b)
                           { return a->StartTime() < b->StartTime(); });
        else // Default to PeakTime
          std::stable_sort(Cluster.begin(), Cluster.end(), [](art::Ptr<recob::OpHit> a, art::Ptr<recob::OpHit> b)
                          { return a->PeakTime() < b->PeakTime(); });
      }
      int Plane = GetOpHitPlane(Cluster[0], 0.1);
      int NHit = 0;
      double Time = -1e6;
      double TimeWidth = 0;
      double TimeSum = 0;
      double PE = 0;
      double MaxPE = 0;
      std::vector<double> PEperOpDet;
      double FastToTotal = 1;
      double X = 0;
      double Y = 0;
      double Z = 0;
      double XWidth = 0;
      double YWidth = 0;
      double ZWidth = 0;
      double XSum = 0;
      double YSum = 0;
      double ZSum = 0;
      double STD = 0;

      // Compute total number of PE and MaxPE.
      for (art::Ptr<recob::OpHit> PDSHit : Cluster)
      {
        NHit++;
        auto thisPlane = GetOpHitPlane(PDSHit, 0.1);
        if (thisPlane != Plane) {
          ProducerUtils::PrintInColor("OpHit in cluster not in same plane: CH " + ProducerUtils::str(PDSHit->OpChannel()) + " Plane " + ProducerUtils::str(thisPlane) + " Expected Plane " + ProducerUtils::str(Plane), ProducerUtils::GetColor("red"), "Error");
          Plane = -1; // Set plane to -1 if hits in cluster are not in the same plane
        }

        PE += PDSHit->PE();
        if (PDSHit->PE() > MaxPE)
          MaxPE = PDSHit->PE();

        PEperOpDet.push_back(PDSHit->PE());
        if (fOpHitTimeVariable == "StartTime")
          TimeSum += PDSHit->StartTime() * TickPeriod * PDSHit->PE();
        else // Default to PeakTime
          TimeSum += PDSHit->PeakTime() * TickPeriod * PDSHit->PE();
      }

      float HotPE = 0;
      Time = TimeSum / PE;
      // Compute flash center from weighted average of "hottest" ophits.
      for (art::Ptr<recob::OpHit> PDSHit : Cluster)
      {
        auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(PDSHit->OpChannel()).GetCenter();
        if (PDSHit->PE() >= fOpFlashAlgoHotVertexThld * MaxPE) {
          XSum += OpHitXYZ.X() * PDSHit->PE();
          YSum += OpHitXYZ.Y() * PDSHit->PE();
          ZSum += OpHitXYZ.Z() * PDSHit->PE();
          HotPE += PDSHit->PE();
        }
      }
      X = XSum / HotPE;
      Y = YSum / HotPE;
      Z = ZSum / HotPE;

      // Compute the flash width and STD from divergence of 1/r² signal decay.
      std::vector<float> varXY;
      std::vector<float> varYZ;
      std::vector<float> varXZ;
      for (art::Ptr<recob::OpHit> PDSHit : Cluster)
      {
        auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(PDSHit->OpChannel()).GetCenter();
        if (fOpHitTimeVariable == "StartTime")
          TimeWidth += (PDSHit->StartTime() * TickPeriod - Time) * (PDSHit->StartTime() * TickPeriod - Time);
        else // Default to PeakTime
          TimeWidth += (PDSHit->PeakTime() * TickPeriod - Time) * (PDSHit->PeakTime() * TickPeriod - Time);
        
        XWidth += (OpHitXYZ.X() - X) * (OpHitXYZ.X() - X);
        YWidth += (OpHitXYZ.Y() - Y) * (OpHitXYZ.Y() - Y);
        ZWidth += (OpHitXYZ.Z() - Z) * (OpHitXYZ.Z() - Z);
        varXY.push_back(sqrt(pow(X - OpHitXYZ.X(), 2) + pow(Y - OpHitXYZ.Y(), 2)) * PDSHit->PE());
        varYZ.push_back(sqrt(pow(Y - OpHitXYZ.Y(), 2) + pow(Z - OpHitXYZ.Z(), 2)) * PDSHit->PE());
        varXZ.push_back(sqrt(pow(X - OpHitXYZ.X(), 2) + pow(Z - OpHitXYZ.Z(), 2)) * PDSHit->PE());
      }

      TimeWidth = sqrt(TimeWidth / Cluster.size());
      XWidth = sqrt(XWidth / Cluster.size());
      YWidth = sqrt(YWidth / Cluster.size());
      ZWidth = sqrt(ZWidth / Cluster.size());

      // Compute STD
      STD = GetOpFlashPlaneSTD(Plane, varXY, varYZ, varXZ);
      // Compute FastToTotal according to the #PEs arriving within the first 10% of the time window wrt the total #PEs
      for (art::Ptr<recob::OpHit> PDSHit : Cluster)
      {
        auto thisTime = 0.0;
        if (fOpHitTimeVariable == "StartTime")
          thisTime = PDSHit->StartTime() * TickPeriod;
        else // Default to PeakTime
          thisTime = PDSHit->PeakTime() * TickPeriod;
        
        if (thisTime < Time + TimeWidth / 10)
          FastToTotal += PDSHit->PE();
      }
      FastToTotal /= PE;
      FlashVec.push_back(FlashInfo{Plane, NHit, Time, TimeWidth, PE, MaxPE, PEperOpDet, FastToTotal, X, Y, Z, XWidth, YWidth, ZWidth, STD});
    }
    return;
  }

  float AdjOpHitsUtils::GetOpFlashPlaneSTD(const int Plane, const std::vector<float> varXY, const std::vector<float> varYZ, const std::vector<float> varXZ)
  {
    std::vector<float> var;
    std::string geoName = geom->DetectorName();
    if (geoName.find("dune10kt") != std::string::npos) {
      var = varYZ; // In HD geometry, only 1 plane
    }
    else if (geoName.find("dunevd10kt") != std::string::npos) {
      if (Plane == 0) {
        var = varYZ;
      }
      else if (Plane == 1 || Plane == 2) {
        var = varXZ;
      }
      else if (Plane == 3 || Plane == 4) {
        var = varXY;
      }
      else {
        ProducerUtils::PrintInColor("Plane not recognized: Must be between 0 and 4", ProducerUtils::GetColor("red"), "Error");
        var = varYZ; // Default to varYZ
      }
    }
    else {
      ProducerUtils::PrintInColor("Geometry not recognized: Must be 'HD' or 'VD'", ProducerUtils::GetColor("red"), "Error");
      var = varYZ; // Default to varYZ
    }

    float varmean = 0;
    for (float v : var)
    {
      varmean += v;
    }
    varmean /= var.size();
    float varstd = 0;
    for (float v : var)
    {
      varstd += pow(v - varmean, 2);
    }
    varstd = sqrt(varstd / var.size());
    return varstd;
  }

  void AdjOpHitsUtils::CalcAdjOpHits(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector, std::vector<std::vector<art::Ptr<recob::OpHit>>> &OpHitClusters, std::vector<std::vector<int>> &OpHitClusterIdx, art::Event const &evt)
  {
    // This is the low energy flash (hit clustering) algorithm aka flip-flop. It is based on the idea that the hits in the same flash are close in time and space and follow a 1/r² signal decay.
    auto clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto TickPeriod = clockData.OpticalClock().TickPeriod();
    // Initialize the vector of OpHitClusters and the vector of indices
    OpHitClusters.clear();
    OpHitClusterIdx.clear();

    // Create a sorting index based on the input vector
    std::vector<int> sortingIndex(OpHitVector.size());
    std::iota(sortingIndex.begin(), sortingIndex.end(), 0);
    
    if (fOpHitTimeVariable == "StartTime")
      std::stable_sort(sortingIndex.begin(), sortingIndex.end(), [&](int a, int b)
                       { return OpHitVector[a]->StartTime() < OpHitVector[b]->StartTime(); });
    else // Default to PeakTime
      std::stable_sort(sortingIndex.begin(), sortingIndex.end(), [&](int a, int b)
                      { return OpHitVector[a]->PeakTime() < OpHitVector[b]->PeakTime(); });

    // Create a vector to track if a hit has been clustered or not
    std::vector<bool> ClusteredHits(OpHitVector.size(), false);

    // Create the OpHitPlane map
    std::map<int, int> OpHitPlane = GetOpHitPlaneMap(OpHitVector);
    std::string sOpHitClustering = "AdjOpHitsUtils::CalcAdjOpHits " + ProducerUtils::str(OpHitVector.size()) + " OpHits found in the event\n";
    for (auto it = sortingIndex.begin(); it != sortingIndex.end(); ++it)
    {
      if (ClusteredHits[*it])
        continue; // Skip if the hit has already been clustered

      std::vector<int> AdjHitIdx = {};
      std::vector<art::Ptr<recob::OpHit>> AdjHitVec = {};

      const auto &hit = OpHitVector[*it];
      if (hit->PE() < fOpFlashAlgoTriggerPE)
        continue;

      bool main_hit = true;

      // If a trigger hit is found, start a new cluster with the hits around it that are within the time and radius range
      ClusteredHits[*it] = true;
      AdjHitIdx.push_back(*it);
      AdjHitVec.push_back(hit);
      sOpHitClustering += "Trigger hit found: PE " + ProducerUtils::str(hit->PE()) + " CH " + ProducerUtils::str(hit->OpChannel()) + " Time " + ProducerUtils::str(hit->PeakTime() * TickPeriod) + "\n";

      int refHit1 = hit->OpChannel();
      auto ref1 = wireReadout.OpDetGeoFromOpChannel(refHit1).GetCenter();
      
      // Make use of the fact that the hits are sorted in time to only consider the hits that are adjacent in the vector up to a certain time range
      for (auto it2 = it; it2 != sortingIndex.end(); ++it2)
      {
        // make sure we don't go out of bounds and the pointer is valid
        if (it == sortingIndex.end())
          break;
        if (it2 == it)
          continue;

        auto &adjHit = OpHitVector[*it2]; // Update adjHit here

        auto thisHitTime = 0.0;
        auto thisAdjHitTime = 0.0;
        if (fOpHitTimeVariable == "StartTime") {
          thisHitTime = hit->StartTime();
          thisAdjHitTime = adjHit->StartTime();
        }
        else { // Default to PeakTime 
          thisHitTime = hit->PeakTime();
          thisAdjHitTime = adjHit->PeakTime();
        }

        if (std::abs(thisAdjHitTime - thisHitTime) * TickPeriod > fOpFlashAlgoMaxTime) {
          sOpHitClustering += "Breaking time loop at dT " + ProducerUtils::str(std::abs(thisAdjHitTime - thisHitTime) * TickPeriod) + " us\n";
          break;
        }
        
        if (adjHit->PE() < fOpFlashAlgoPE)
          continue;

        int refHit2 = adjHit->OpChannel();
        // Check if the hits are in the same plane and continue if they are not
        if (CheckOpHitPlane(OpHitPlane, refHit1, refHit2) == false)
          continue;

        // If hit has already been clustered, skip
        if (ClusteredHits[*it2] == true && !fOpFlashAlgoHitDuplicates) {
          sOpHitClustering += "Skipping already clustered hit: CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(adjHit->PeakTime() * TickPeriod) + "\n";
          continue;
        }

        auto ref2 = wireReadout.OpDetGeoFromOpChannel(refHit2).GetCenter();
        auto ref3 = TVector3(ref1.X(), ref1.Y(), ref1.Z()) - TVector3(ref2.X(), ref2.Y(), ref2.Z());
        if (ref3.Mag() < fOpFlashAlgoRad)
        {
          if (adjHit->PE() > hit->PE())
          {
            sOpHitClustering += "¡¡¡Hit with PE > TriggerPE found: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(adjHit->PeakTime() * TickPeriod) + "\n";
            ClusteredHits[*it] = false;
            main_hit = false;

            // Reset the ClusteredHits values for the hits that have been added to the cluster
            for (auto it3 = AdjHitIdx.begin(); it3 != AdjHitIdx.end(); ++it3)
            {
              sOpHitClustering += "---Removing hit: CH " + ProducerUtils::str(OpHitVector[*it3]->OpChannel()) + " Time " + ProducerUtils::str(OpHitVector[*it3]->PeakTime() * TickPeriod) + "\n";
              ClusteredHits[*it3] = false;
            }
            break;
          }
          else {
            sOpHitClustering += "+++Adding hit: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(adjHit->PeakTime() * TickPeriod) + "\n";
            AdjHitVec.push_back(adjHit);
            AdjHitIdx.push_back(*it2);
            ClusteredHits[*it2] = true;
          }
        }
      }

      if (main_hit == false) {
        sOpHitClustering += "Skipping non-main hit\n";
        ClusteredHits[*it] = false;
        continue;
      }

      for (auto it4 = it; it4 != sortingIndex.begin(); --it4)
      {
        // make sure we don't go out of bounds and the pointer is valid
        if (it == sortingIndex.begin()) {
          sOpHitClustering += "Breaking time loop at beginning of vector\n";
          break;
        }

        if (it4 == it)
          continue;
        
        auto thisHitTime = 0.0;
        auto thisAdjHitTime = 0.0;
        auto &adjHit = OpHitVector[*it4];
        if (fOpHitTimeVariable == "StartTime") {
          thisHitTime = hit->StartTime();
          thisAdjHitTime = adjHit->StartTime();
        }
        else { // Default to PeakTime 
          thisHitTime = hit->PeakTime();
          thisAdjHitTime = adjHit->PeakTime();
        }

        if (std::abs(thisHitTime - thisAdjHitTime) * TickPeriod > fOpFlashAlgoMinTime)
          break;
        
        if (adjHit->PE() < fOpFlashAlgoPE)
          continue;

        int refHit4 = adjHit->OpChannel();
        // Check if the hits are in the same plane and continue if they are not
        if (CheckOpHitPlane(OpHitPlane, refHit1, refHit4) == false)
          continue;

        // if hit has already been clustered, skip
        if (ClusteredHits[*it4] == true && !fOpFlashAlgoHitDuplicates) {
          sOpHitClustering += "Skipping already clustered hit: CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(adjHit->PeakTime() * TickPeriod) + "\n";
          continue;
        }

        auto ref4 = wireReadout.OpDetGeoFromOpChannel(refHit4).GetCenter();
        auto ref5 = TVector3(ref1.X(), ref1.Y(), ref1.Z()) - TVector3(ref4.X(), ref4.Y(), ref4.Z());
        if (ref5.Mag() < fOpFlashAlgoRad)
        {
          if (adjHit->PE() > hit->PE())
          {
            sOpHitClustering += "¡¡¡Hit with PE > TriggerPE found: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(adjHit->PeakTime() * TickPeriod) + "\n";
            ClusteredHits[*it] = false;
            main_hit = false;

            // Reset the ClusteredHits values for the hits that have been added to the cluster
            for (auto it5 = AdjHitIdx.begin(); it5 != AdjHitIdx.end(); ++it5)
            {
              sOpHitClustering += "---Removing hit: CH " + ProducerUtils::str(OpHitVector[*it5]->OpChannel()) + " Time " + ProducerUtils::str(OpHitVector[*it5]->PeakTime() * TickPeriod) + "\n";
              ClusteredHits[*it5] = false;
            }
            break;
          }
          else {
            sOpHitClustering += "Adding hit: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(adjHit->PeakTime() * TickPeriod) + "\n";
            AdjHitVec.push_back(adjHit);
            AdjHitIdx.push_back(*it4);
            ClusteredHits[*it4] = true;
          }
        }
      }

      if (main_hit && int(AdjHitVec.size()) >= fOpFlashAlgoNHit) {
        // Store the original indices of the clustered hits
        OpHitClusters.push_back(std::move(AdjHitVec));
        OpHitClusterIdx.push_back(std::move(AdjHitIdx));
        sOpHitClustering += "***Cluster size: " + ProducerUtils::str(int(OpHitClusters.back().size())) + "\n";
        ProducerUtils::PrintInColor(sOpHitClustering, ProducerUtils::GetColor("green"), "Debug");
      }
      else {
        sOpHitClustering += "Cluster rejected: Size " + ProducerUtils::str(int(AdjHitVec.size())) + " < MinSize " + ProducerUtils::str(fOpFlashAlgoNHit) + "\n";
        ProducerUtils::PrintInColor(sOpHitClustering, ProducerUtils::GetColor("red"), "Debug");
      }
    }
    return;
  }

  void AdjOpHitsUtils::FlashMatchResidual(float &Residual, std::vector<art::Ptr<recob::OpHit>> Hits, double x, double y, double z)
  {
    if (Hits.size() == 0)
    {
      Residual = 1e6;
      ProducerUtils::PrintInColor("Failed Residual Evaluation: Empty Flash!", ProducerUtils::GetColor("yellow"), "Error");
      return;
    }
    // Initialize variables
    Residual = 0;
    float PE = 0;

    // Find index of the hit with the highest PE
    int maxPEIdx = 0;
    for (unsigned int i = 0; i < Hits.size(); i++)
    {
      if (Hits[i]->PE() > Hits[maxPEIdx]->PE())
        maxPEIdx = i;
    }

    // Start with the first hit in the flash as reference point
    double firstHitY = wireReadout.OpDetGeoFromOpChannel(Hits[maxPEIdx]->OpChannel()).GetCenter().Y();
    double firstHitZ = wireReadout.OpDetGeoFromOpChannel(Hits[maxPEIdx]->OpChannel()).GetCenter().Z();

    // Get the first hit PE and calculate the squared distance and angle to the reference point
    float firstHitPE = Hits[maxPEIdx]->PE();
    float firstHitDistSq = pow(firstHitY - y, 2) + pow(firstHitZ - z, 2);

    // Calculate the expected PE value for the reference point based on the first hit PE and the squared distance + angle
    // float firstHitAngle = atan2(sqrt(firstHitDistSq), abs(x));
    // float refHitPE = firstHitPE * (pow(x, 2) + firstHitDistSq) / pow(x, 2) / cos(firstHitAngle);
    float refHitPE = firstHitPE * (pow(x, 2) + firstHitDistSq) / pow(x, 2);

    // Loop over all OpHits in the flash and compute the squared distance to the reference point
    for (const auto &hit : Hits)
    {
      double hitY = wireReadout.OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Y();
      double hitZ = wireReadout.OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Z();

      // The expected distribution of PE corresponds to a decrease of 1/r² with the distance from the flash center. Between adjacent OpHits, the expected decrease in charge has the form r²/(r²+d²)
      float hitDistSq = pow(hitY - y, 2) + pow(hitZ - z, 2);

      // float hitAngle = atan2(sqrt(hitDistSq), abs(x));
      // float predPE = refHitPE * cos(hitAngle) * pow(x, 2) / (pow(x, 2) + hitDistSq);
      float predPE = refHitPE * pow(x, 2) / (pow(x, 2) + hitDistSq);

      Residual += pow(hit->PE() - predPE, 2);
      PE += hit->PE();
    }

    Residual /= float(Hits.size());
    std::string debug = "PE: " + ProducerUtils::str(PE) +
      " X: " + ProducerUtils::str(x) +
      " RefPE: " + ProducerUtils::str(refHitPE) +
      " NHits: " + ProducerUtils::str(int(Hits.size())) +
      " Residual: " + ProducerUtils::str(Residual);

    ProducerUtils::PrintInColor(debug, ProducerUtils::GetColor("yellow"), "Debug");
    return;
  }

  int AdjOpHitsUtils::GetOpHitPlane(const art::Ptr<recob::OpHit> &hit, float buffer)
  {
    std::string geoName = geom->DetectorName();
    auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter();
    if (geoName.find("dune10kt") != std::string::npos) {
      if (abs(OpHitXYZ.X()) < buffer)
        return 0;
      else
        return -1; // Return -1 if the hit not at the center
    }
    else if(geoName.find("dunevd10kt") != std::string::npos) {
      if (OpHitXYZ.X() > fXACathodeX - buffer && OpHitXYZ.X() < fXACathodeX + buffer)
        return 0;
      else if (OpHitXYZ.Y() > fXAMembraneY - buffer && OpHitXYZ.Y() < fXAMembraneY + buffer)
        return 1;
      else if (OpHitXYZ.Y() > -fXAMembraneY - buffer && OpHitXYZ.Y() < -fXAMembraneY + buffer)
        return 2;
      else if (OpHitXYZ.Z() > fXAFinalCapZ - buffer && OpHitXYZ.Z() < fXAFinalCapZ + buffer)
        return 3;
      else if (OpHitXYZ.Z() > fXAStartCapZ - buffer && OpHitXYZ.Z() < fXAStartCapZ + buffer)
        return 4;
      else {
        ProducerUtils::PrintInColor("Hit not close to any plane: CH " + ProducerUtils::str(hit->OpChannel()) + " X " + ProducerUtils::str(OpHitXYZ.X()) + " Y " + ProducerUtils::str(OpHitXYZ.Y()) + " Z " + ProducerUtils::str(OpHitXYZ.Z()), ProducerUtils::GetColor("red"), "Error");
        return -1; // Return -1 if the hit is not close to any plane
      }
    }
    else {
      ProducerUtils::PrintInColor("Unknown geometry: " + geoName, ProducerUtils::GetColor("red"), "Error");
      return -1;
    }
  }

  // Define a function and map to get the plane of the hit. If geometry is VD the number of planes is 5 (cathode, leftmembrane, rightmembrane, startcap, finalcap). If geometry is HD the number of planes is 2 (leftdrift, rightdrift).
  std::map<int, int> AdjOpHitsUtils::GetOpHitPlaneMap(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector)
  {
    std::map<int, int> OpHitPlane;
    for (const auto &hit : OpHitVector)
    {
      OpHitPlane[hit->OpChannel()] = GetOpHitPlane(hit, 0.1);
    }
    return OpHitPlane;
  }

  // Define a function to check if the adjhit is in the same plane as the reference hit using the OpHitPlane map
  bool AdjOpHitsUtils::CheckOpHitPlane(std::map<int, int> OpHitPlane, int refHit, int adjHit)
  {
    // Check that both hits are in the map
    if (OpHitPlane.find(refHit) == OpHitPlane.end() || OpHitPlane.find(adjHit) == OpHitPlane.end()) {
      ProducerUtils::PrintInColor("Hit not found in OpHitPlane map: refHit " + ProducerUtils::str(refHit) + " adjHit " + ProducerUtils::str(adjHit), ProducerUtils::GetColor("red"), "Error");
      return false;
    }
    
    if (OpHitPlane[refHit] == OpHitPlane[adjHit]) 
      return true;
    else
      return false;
  }

} // namespace solar