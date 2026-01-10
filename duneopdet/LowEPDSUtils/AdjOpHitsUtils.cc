#include "AdjOpHitsUtils.h"

using namespace producer;

namespace solar
{
  AdjOpHitsUtils::AdjOpHitsUtils(fhicl::ParameterSet const &p)
      : fOpWaveformLabel(p.get<std::string>("OpWaveformLabel", "opdetwaveform")), // Label for OpDetWaveform collection.
        fOpHitLabel(p.get<std::string>("OpHitLabel", "ophit")), // Label for OpHit collection.
        fOpHitTimeVariable(p.get<std::string>("OpHitTimeVariable", "PeakTime")), // Variable to use for time sorting ("StartTime" or "PeakTime").
        fOpHitTime2us(p.get<bool>("OpHitTime2us", false)), // Conversion factor from OpHit time units to microseconds. Default factor is TickPeriod() from DetectorClocksService.
        fOpFlashAlgoNHit(p.get<int>("OpFlashAlgoNHit", 3)),       // Minimum number of OpHits in a cluster to consider it for flash creation.
        fOpFlashAlgoMinTime(p.get<float>("OpFlashAlgoMinTime", 0.32)), // Negative time window to look for adj. OpHits in [us].
        fOpFlashAlgoMaxTime(p.get<float>("OpFlashAlgoMaxTime", 0.96)), // Positive time window to look for adj. OpHits in [us].
        fOpFlashAlgoWeightedTime(p.get<bool>("OpFlashAlgoWeightedTime", false)), // Whether to use weighted time of trigger time for flash time reference.
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
    if (fOpFlashAlgoHotVertexThld > 1) {
      ProducerUtils::PrintInColor("Hot vertex threshold must be between 0 and 1", ProducerUtils::GetColor("red"), "Error");
      return;
    }
    
    for (std::vector<art::Ptr<recob::OpHit>> Cluster : OpHitClusters)
    {
      if (Cluster.empty())
      {
        continue;
      }

      int Plane = GetOpHitPlane(Cluster[0], 0.1);
      int MaxIdx = 0;
      int Idx = 0;
      int NHit = 0;
      double TimeMax = -1e6;
      double TimeWidth = 0;
      double TimeWeighted = -1e6;
      double TimeSum = 0;
      double Amplitude = 0;
      double PE = 0;
      double MaxPE = 0;
      std::vector<double> PEperOpDet = {};
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
      std::vector<int> MainOpWaveform = {}; // Make vector for main OpWaveform with 1000 entries (max waveform size)
      bool MainOpWaveformValid = false;
      float MainOpWaveformTime = -1e6;

      std::vector<bool> OpHitWvfValid = {};
      std::vector<float> OpHitWvfTime = {};
      std::vector<std::vector<int>> OpHitWvfIntVector = {};
      GetOpHitSignal(Cluster, OpHitWvfIntVector, OpHitWvfTime, OpHitWvfValid, evt); // Get OpHit waveforms

      // Compute total number of PE and MaxPE.
      for (art::Ptr<recob::OpHit> PDSHit : Cluster)
      {
        NHit++;
        double thisTime = -1e6;
        double thisPE = PDSHit->PE();
        double thisAmp = PDSHit->Amplitude();
        auto thisPlane = GetOpHitPlane(PDSHit, 0.1);
        
        if (thisPlane != Plane) {
          ProducerUtils::PrintInColor("OpHit in cluster not in same plane: CH " + ProducerUtils::str(PDSHit->OpChannel()) + " Plane " + ProducerUtils::str(thisPlane) + " Expected Plane " + ProducerUtils::str(Plane), ProducerUtils::GetColor("red"), "Error");
          Plane = -1; // Set plane to -1 if hits in cluster are not in the same plane
        }

        if (fOpHitTimeVariable == "StartTime") {
          thisTime = PDSHit->StartTime(); // Use StartTime
        }
        else {
          thisTime = PDSHit->PeakTime(); // Default to PeakTime
        }
        if (fOpHitTime2us) {
          thisTime *= TickPeriod; // Convert to microseconds
        }
        
        PE += thisPE;
        TimeSum += thisTime * thisPE;
        
        if (thisPE > MaxPE) {
          Amplitude = thisAmp;
          MaxPE = thisPE;
          MaxIdx = Idx;
          TimeMax = thisTime;
        }

        PEperOpDet.push_back(thisPE);
        Idx++;
      }

      if (OpHitWvfValid[MaxIdx]) {
        MainOpWaveform = OpHitWvfIntVector[MaxIdx];
        MainOpWaveformTime = OpHitWvfTime[MaxIdx];
        MainOpWaveformValid = true;
      }

      float HotPE = 0;
      TimeWeighted = TimeSum / PE;

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
        float ThisOpHitTime = -1e6;
        auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(PDSHit->OpChannel()).GetCenter();
        
        if (fOpHitTimeVariable == "StartTime")
          ThisOpHitTime = PDSHit->StartTime();
        else // Default to PeakTime
          ThisOpHitTime = PDSHit->PeakTime();
        if (fOpHitTime2us) {
          ThisOpHitTime *= TickPeriod;
        }
        
        TimeWidth += (TimeMax - ThisOpHitTime) * (TimeMax - ThisOpHitTime);
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
        auto thisTime = -1e6;
        if (fOpHitTimeVariable == "StartTime")
          thisTime = PDSHit->StartTime();
        else // Default to PeakTime
          thisTime = PDSHit->PeakTime();
        if (fOpHitTime2us) {
          thisTime *= TickPeriod;
        }
        
        // Check if thisTime is within 10% of the time window in positive and negative direction
        if (std::abs(thisTime - TimeMax) <= 0.05 * TimeWidth)
          FastToTotal += PDSHit->PE();
      }
      FastToTotal /= PE;
      FlashVec.push_back(FlashInfo{
        Plane, 
        NHit, 
        MaxPE, 
        TimeMax, 
        Amplitude, 
        TimeWidth, 
        TimeWeighted, 
        PE, 
        PEperOpDet, 
        FastToTotal, 
        X, 
        Y, 
        Z, 
        XWidth, 
        YWidth, 
        ZWidth, 
        STD, 
        MainOpWaveform, 
        MainOpWaveformTime, 
        MainOpWaveformValid
        }
      );
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
      float hitTime = -1e6;
      if (fOpHitTimeVariable == "StartTime")
        hitTime = hit->StartTime();
      else // Default to PeakTime
        hitTime = hit->PeakTime();
      if (fOpHitTime2us) {
        hitTime *= TickPeriod;
      }
      bool main_hit = true;

      // If a trigger hit is found, start a new cluster with the hits around it that are within the time and radius range
      ClusteredHits[*it] = true;
      AdjHitIdx.push_back(*it);
      AdjHitVec.push_back(hit);
      sOpHitClustering += "Trigger hit found: PE " + ProducerUtils::str(hit->PE()) + " CH " + ProducerUtils::str(hit->OpChannel()) + " Time " + ProducerUtils::str(hitTime) + "\n";

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
        if (fOpHitTime2us) {
          thisHitTime *= TickPeriod;
          thisAdjHitTime *= TickPeriod;
        }

        if (std::abs(thisAdjHitTime - thisHitTime) > fOpFlashAlgoMaxTime) {
          sOpHitClustering += "Breaking time loop at dT " + ProducerUtils::str(std::abs(thisAdjHitTime - thisHitTime)) + " us\n";
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
          sOpHitClustering += "Skipping already clustered hit: CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(thisAdjHitTime) + "\n";
          continue;
        }

        auto ref2 = wireReadout.OpDetGeoFromOpChannel(refHit2).GetCenter();
        auto ref3 = TVector3(ref1.X(), ref1.Y(), ref1.Z()) - TVector3(ref2.X(), ref2.Y(), ref2.Z());
        if (ref3.Mag() < fOpFlashAlgoRad)
        {
          if (adjHit->PE() > hit->PE())
          {
            sOpHitClustering += "¡¡¡Hit with PE > TriggerPE found: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(thisAdjHitTime) + "\n";
            ClusteredHits[*it] = false;
            main_hit = false;

            // Reset the ClusteredHits values for the hits that have been added to the cluster
            for (auto it3 = AdjHitIdx.begin(); it3 != AdjHitIdx.end(); ++it3)
            {
              auto refTime = OpHitVector[*it3]->PeakTime();
              auto refChannel = OpHitVector[*it3]->OpChannel();
              if (fOpHitTimeVariable == "StartTime")
                refTime = OpHitVector[*it3]->StartTime();
              if (fOpHitTime2us) {
                refTime *= TickPeriod;
              }
              sOpHitClustering += "---Removing hit: CH " + ProducerUtils::str(refChannel) + " Time " + ProducerUtils::str(refTime) + "\n";
              ClusteredHits[*it3] = false;
            }
            break;
          }
          else {
            sOpHitClustering += "+++Adding hit: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(thisAdjHitTime) + "\n";
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
        if (fOpHitTime2us) {
          thisHitTime *= TickPeriod;
          thisAdjHitTime *= TickPeriod;
        }

        if (std::abs(thisHitTime - thisAdjHitTime) > fOpFlashAlgoMinTime)
          break;
        
        if (adjHit->PE() < fOpFlashAlgoPE)
          continue;

        int refHit4 = adjHit->OpChannel();
        // Check if the hits are in the same plane and continue if they are not
        if (CheckOpHitPlane(OpHitPlane, refHit1, refHit4) == false)
          continue;

        // if hit has already been clustered, skip
        if (ClusteredHits[*it4] == true && !fOpFlashAlgoHitDuplicates) {
          sOpHitClustering += "Skipping already clustered hit: CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(thisAdjHitTime) + "\n";
          continue;
        }

        auto ref4 = wireReadout.OpDetGeoFromOpChannel(refHit4).GetCenter();
        auto ref5 = TVector3(ref1.X(), ref1.Y(), ref1.Z()) - TVector3(ref4.X(), ref4.Y(), ref4.Z());
        if (ref5.Mag() < fOpFlashAlgoRad)
        {
          if (adjHit->PE() > hit->PE())
          {
            sOpHitClustering += "¡¡¡Hit with PE > TriggerPE found: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(thisAdjHitTime) + "\n";
            ClusteredHits[*it] = false;
            main_hit = false;

            // Reset the ClusteredHits values for the hits that have been added to the cluster
            for (auto it5 = AdjHitIdx.begin(); it5 != AdjHitIdx.end(); ++it5)
            {
              auto refTime = OpHitVector[*it5]->PeakTime();
              auto refChannel = OpHitVector[*it5]->OpChannel();
              if (fOpHitTimeVariable == "StartTime")
                refTime = OpHitVector[*it5]->StartTime();
              if (fOpHitTime2us) {
                refTime *= TickPeriod;
              }
              sOpHitClustering += "---Removing hit: CH " + ProducerUtils::str(refChannel) + " Time " + ProducerUtils::str(refTime) + "\n";
              ClusteredHits[*it5] = false;
            }
            break;
          }
          else {
            sOpHitClustering += "Adding hit: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(thisAdjHitTime) + "\n";
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


  std::map<int, int> AdjOpHitsUtils::GetOpHitPlaneMap(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector)
  { // Define a function and map to get the plane of the hit. If geometry is VD the number of planes is 5 (cathode, leftmembrane, rightmembrane, startcap, finalcap). If geometry is HD the number of planes is 2 (leftdrift, rightdrift).
    std::map<int, int> OpHitPlane;
    for (const auto &hit : OpHitVector)
    {
      OpHitPlane[hit->OpChannel()] = GetOpHitPlane(hit, 0.1);
    }
    return OpHitPlane;
  }


  bool AdjOpHitsUtils::CheckOpHitPlane(std::map<int, int> OpHitPlane, int refHit, int adjHit)
  { // Define a function to check if the adjhit is in the same plane as the reference hit using the OpHitPlane map
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


  void AdjOpHitsUtils::GetOpHitWaveforms(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector, std::vector<art::Ptr<raw::OpDetWaveform>> &OpHitWvfVector, std::vector<bool> &OpHitWvfValid, art::Event const &evt)
  { // Define a function to get the OpHit waveforms that correspond to each OpHit in the input vector
    // Get the OpDetWaveform handle
    art::Handle<std::vector<raw::OpDetWaveform>> opDetWvfHandle;
    evt.getByLabel(fOpWaveformLabel, opDetWvfHandle);
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    if (!opDetWvfHandle.isValid())
    {
      ProducerUtils::PrintInColor("Invalid OpDetWaveform handle", ProducerUtils::GetColor("red"), "Debug");
      for (size_t i = 0; i < OpHitVector.size(); i++)
      {
        OpHitWvfVector.push_back(art::Ptr<raw::OpDetWaveform>());  // Fill with a null pointer to keep the indices consistent
        OpHitWvfValid.push_back(false);
      }
      return;
    }

    // Loop over the OpHit vector and get the corresponding OpDetWaveform
    for (const auto &hit : OpHitVector)
    {
      bool found = false;
      float hitTime = -1e6;
      unsigned int opChannel = hit->OpChannel();
      auto TickPeriod = clockData.OpticalClock().TickPeriod();
      
      if (fOpHitTimeVariable == "StartTime")
        hitTime = hit->StartTime();
      else // Default to PeakTime
        hitTime = hit->PeakTime(); 
      if (fOpHitTime2us) {
        hitTime *= TickPeriod;
      }
      
      for (size_t i = 0; i < opDetWvfHandle->size(); i++)
      { // Match by channel number and amplitude
        if ((*opDetWvfHandle)[i].ChannelNumber() == opChannel && hitTime >= (*opDetWvfHandle)[i].TimeStamp() && hitTime <= (*opDetWvfHandle)[i].TimeStamp() + (*opDetWvfHandle)[i].Waveform().size() * TickPeriod)
        {
          ProducerUtils::PrintInColor("Matching OpDetWaveform found for OpHit channel " + ProducerUtils::str(int(opChannel)), ProducerUtils::GetColor("green"), "Debug");
          art::Ptr<raw::OpDetWaveform> wvfPtr(opDetWvfHandle, i);
          OpHitWvfVector.push_back(wvfPtr);
          OpHitWvfValid.push_back(true);
          found = true;
          break;
        }
      }
      if (!found)
      {
        ProducerUtils::PrintInColor("No matching OpDetWaveform found for OpHit channel " + ProducerUtils::str(int(opChannel)), ProducerUtils::GetColor("red"), "Debug");
        OpHitWvfVector.push_back(art::Ptr<raw::OpDetWaveform>());  // Fill with a null pointer to keep the indices consistent
        OpHitWvfValid.push_back(false);
      }
    }
    return;
  }


  void AdjOpHitsUtils::GetOpHitWaveforms(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector, std::vector<art::Ptr<recob::OpWaveform>> &OpHitWvfVector, std::vector<bool> &OpHitWvfValid, art::Event const &evt)
  { // Repeat the function but for waveforms of type std::vector<recob::OpWaveform>
    // Get the OpWaveform handle
    art::Handle<std::vector<recob::OpWaveform>> opWvfHandle;
    evt.getByLabel(fOpWaveformLabel, opWvfHandle);
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    if (!opWvfHandle.isValid())
    {
      ProducerUtils::PrintInColor("Invalid OpWaveform handle", ProducerUtils::GetColor("red"), "Debug");
      for (size_t i = 0; i < OpHitVector.size(); i++)
      {
        OpHitWvfVector.push_back(art::Ptr<recob::OpWaveform>());  // Fill with a null pointer to keep the indices consistent
        OpHitWvfValid.push_back(false);
      }
      return;
    }

    // Loop over the OpHit vector and get the corresponding OpWaveform
    for (const auto &hit : OpHitVector)
    {
      bool found = false;
      float hitTime = -1e6;
      unsigned int opChannel = hit->OpChannel();
      auto TickPeriod = clockData.OpticalClock().TickPeriod();

      if (fOpHitTimeVariable == "StartTime")
        hitTime = hit->StartTime();
      else // Default to PeakTime
        hitTime = hit->PeakTime();
      if (fOpHitTime2us) {
        hitTime *= TickPeriod;
      }

      for (size_t i = 0; i < opWvfHandle->size(); i++)
      { // Match by channel number and amplitude
        if ((*opWvfHandle)[i].Channel() == opChannel && hitTime >= (*opWvfHandle)[i].TimeStamp() && hitTime <= (*opWvfHandle)[i].TimeStamp() + (*opWvfHandle)[i].Signal().size() * TickPeriod)
        {
          ProducerUtils::PrintInColor("Matching OpWaveform found for OpHit channel " + ProducerUtils::str(int(opChannel)), ProducerUtils::GetColor("green"), "Debug");
          art::Ptr<recob::OpWaveform> wvfPtr(opWvfHandle, i);
          OpHitWvfVector.push_back(wvfPtr);
          OpHitWvfValid.push_back(true);
          found = true;
          break;
        }
      }

      if (!found) {
        ProducerUtils::PrintInColor("No matching OpWaveform found for OpHit channel " + ProducerUtils::str(int(opChannel)), ProducerUtils::GetColor("red"), "Debug");
        OpHitWvfVector.push_back(art::Ptr<recob::OpWaveform>());
        OpHitWvfValid.push_back(false);
      }
    }
    return;
  }


  void AdjOpHitsUtils::GetOpHitSignal(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector, std::vector<std::vector<int>> &OpHitWvfIntVector, std::vector<float> &OpHitWvfTime, std::vector<bool> &OpHitWvfValid, art::Event const &evt)
  {
    auto geoName = geom->DetectorName();
    if (geoName.find("dune10kt") != std::string::npos) {
      std::vector<art::Ptr<recob::OpWaveform>> OpHitWvfVector;
      GetOpHitWaveforms(OpHitVector, OpHitWvfVector, OpHitWvfValid, evt);
      // Get entries in OpHitWvfvector and save to OpHitWvfIntVector as integers
      for (int i = 0; i < int(OpHitWvfVector.size()); i++)
      {
        art::Ptr<recob::OpWaveform> wvf = OpHitWvfVector[i];
        if (OpHitWvfValid[i]) {
          std::vector<int> wvfInt;
          for (auto adc : wvf->Signal()) {
            wvfInt.push_back(int(adc));
          }
          OpHitWvfIntVector.push_back(wvfInt);
          OpHitWvfTime.push_back(wvf->TimeStamp());
        }
        else {
          OpHitWvfIntVector.push_back(std::vector<int>{});
          OpHitWvfTime.push_back(-1e6);
        }
      } 
    }
    else if (geoName.find("dunevd10kt") != std::string::npos) {
      std::vector<art::Ptr<raw::OpDetWaveform>> OpHitWvfVector;
      GetOpHitWaveforms(OpHitVector, OpHitWvfVector, OpHitWvfValid, evt);
      for (int i = 0; i < int(OpHitWvfVector.size()); i++)
      {
        art::Ptr<raw::OpDetWaveform> wvf = OpHitWvfVector[i];
        if (OpHitWvfValid[i]) {
          std::vector<int> wvfInt;
          for (auto adc : wvf->Waveform()) {
            wvfInt.push_back(int(adc));
          }
          OpHitWvfIntVector.push_back(wvfInt);
          OpHitWvfTime.push_back(wvf->TimeStamp());
        }
        else {
          OpHitWvfIntVector.push_back(std::vector<int>{});
          OpHitWvfTime.push_back(-1e6);
        }
      }
    }
  }
} // namespace solar