#include "AdjOpHitsUtils.h"

using namespace producer;

namespace solar
{
  AdjOpHitsUtils::AdjOpHitsUtils(fhicl::ParameterSet const &p)
      : fGeometry(p.get<std::string>("Geometry")),
        fOpFlashAlgoNHit(p.get<int>("OpFlashAlgoNHit")),
        fOpFlashAlgoMinTime(p.get<float>("OpFlashAlgoMinTime")),
        fOpFlashAlgoMaxTime(p.get<float>("OpFlashAlgoMaxTime")),
        fOpFlashAlgoRad(p.get<float>("OpFlashAlgoRad")),
        fOpFlashAlgoPE(p.get<float>("OpFlashAlgoPE")),
        fOpFlashAlgoTriggerPE(p.get<float>("OpFlashAlgoTriggerPE")),
        fOpFlashAlgoHotVertexThld(p.get<float>("OpFlashAlgoHotVertexThld")),
        fXACathodeX(p.get<float>("XACathodeX")),
        fXAMembraneY(p.get<float>("XAMembraneY")),
        fXAFinalCapZ(p.get<float>("XAFinalCapZ")),
        fXAStartCapZ(p.get<float>("XAStartCapZ"))
        // fOpFlashAlgoCentroid(p.get<bool>("OpFlashAlgoCentroid"))
  {
  }
  void AdjOpHitsUtils::MakeFlashVector(std::vector<FlashInfo> &FlashVec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &OpHitClusters, art::Event const &evt)
  {
    if (fOpFlashAlgoHotVertexThld > 1)
    {
      ProducerUtils::PrintInColor("Hot vertex threshold must be between 0 and 1", ProducerUtils::GetColor("red"), "Error");
      return;
    }
    for (std::vector<art::Ptr<recob::OpHit>> Cluster : OpHitClusters)
    {
      if (!Cluster.empty())
      {
        std::stable_sort(Cluster.begin(), Cluster.end(), [](art::Ptr<recob::OpHit> a, art::Ptr<recob::OpHit> b)
                         { return a->PeakTime() < b->PeakTime(); });
      }
      int Plane = -1;
      int NHit = 0;
      double Time = Cluster[0]->PeakTime();
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
        Plane = GetOpHitPlane(PDSHit);

        NHit++;
        PE += PDSHit->PE();
        if (PDSHit->PE() > MaxPE)
        {
          MaxPE = PDSHit->PE();
        }

        PEperOpDet.push_back(PDSHit->PE());
        TimeSum += PDSHit->PeakTime() * PDSHit->PE();
      }
      Time = TimeSum / PE;

      // Compute flash center from weighted average of "hottest" ophits.
      float HotPE = 0;
      for (art::Ptr<recob::OpHit> PDSHit : Cluster)
      {
        auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(PDSHit->OpChannel()).GetCenter();
        if (PDSHit->PE() >= fOpFlashAlgoHotVertexThld * MaxPE)
        {
          XSum += OpHitXYZ.X() * PDSHit->PE();
          YSum += OpHitXYZ.Y() * PDSHit->PE();
          ZSum += OpHitXYZ.Z() * PDSHit->PE();
          HotPE += PDSHit->PE();
        }
      }
      X = XSum / HotPE;
      Y = YSum / HotPE;
      Z = ZSum / HotPE;

      // Alternatively compute the centroid of the flash in 3D space. NEEDS TO BE IMPLEMENTED!
      // if (fOpFlashAlgoCentroid)
      // {
      //   CalcCentroid(Cluster, X, Y, Z);
      // }

      // Compute the flash width and STD from divergence of 1/r² signal decay.
      std::vector<float> varXY;
      std::vector<float> varYZ;
      std::vector<float> varXZ;
      for (art::Ptr<recob::OpHit> PDSHit : Cluster)
      {
        auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(PDSHit->OpChannel()).GetCenter();
        TimeWidth += (PDSHit->PeakTime() - Time) * (PDSHit->PeakTime() - Time);
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
        if (PDSHit->PeakTime() < Time + TimeWidth / 10)
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
    if (Plane == 0)
    {
      var = varYZ;
    }
    if (Plane == 1 || Plane == 2)
    {
      var = varXZ;
    }
    if (Plane == 3 || Plane == 4)
    {
      var = varXY;
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

  void AdjOpHitsUtils::CalcAdjOpHits(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector, std::vector<std::vector<art::Ptr<recob::OpHit>>> &OpHitClusters, std::vector<std::vector<int>> &OpHitClusterIdx)
  {
    // This is the low energy flash (hit clustering) algorithm aka flip-flop. It is based on the idea that the hits in the same flash are close in time and space and follow a 1/r² signal decay.

    // Initialize the vector of OpHitClusters and the vector of indices
    OpHitClusters.clear();
    OpHitClusterIdx.clear();

    // Create a sorting index based on the input vector
    std::vector<int> sortingIndex(OpHitVector.size());
    std::iota(sortingIndex.begin(), sortingIndex.end(), 0);
    std::stable_sort(sortingIndex.begin(), sortingIndex.end(), [&](int a, int b)
                     { return OpHitVector[a]->PeakTime() < OpHitVector[b]->PeakTime(); });

    // Create a vector to track if a hit has been clustered or not
    std::vector<bool> ClusteredHits(OpHitVector.size(), false);

    // Create the OpHitPlane map
    std::map<int, int> OpHitPlane = GetOpHitPlaneMap(OpHitVector);

    for (auto it = sortingIndex.begin(); it != sortingIndex.end(); ++it)
    {
      std::string sOpHitClustering = "";
      std::vector<art::Ptr<recob::OpHit>> AdjHitVec = {};
      std::vector<int> AdjHitIdx = {};

      const auto &hit = OpHitVector[*it];
      if (hit->PE() < fOpFlashAlgoTriggerPE)
        continue;

      bool main_hit = true;

      // If a trigger hit is found, start a new cluster with the hits around it that are within the time and radius range
      ClusteredHits[*it] = true;
      AdjHitVec.push_back(hit);
      AdjHitIdx.push_back(*it);
      sOpHitClustering += "Trigger hit found: PE " + ProducerUtils::str(hit->PE()) + " CH " + ProducerUtils::str(hit->OpChannel()) + " Time " + ProducerUtils::str(hit->PeakTime()) + "\n";

      int refHit1 = hit->OpChannel();
      auto ref1 = wireReadout.OpDetGeoFromOpChannel(refHit1).GetCenter();
      
      // Make use of the fact that the hits are sorted in time to only consider the hits that are adjacent in the vector up to a certain time range
      for (auto it2 = it + 1; it2 != sortingIndex.end(); ++it2)
      {
        // make sure we don't go out of bounds and the pointer is valid
        if (it == sortingIndex.end())
          break;
        if (it2 == sortingIndex.end())
          break;

        auto &adjHit = OpHitVector[*it2]; // Update adjHit here

        if (std::abs(adjHit->PeakTime() - hit->PeakTime()) > fOpFlashAlgoMaxTime)
          break;
        if (adjHit->PE() < fOpFlashAlgoPE)
          continue;

        int refHit2 = adjHit->OpChannel();
        auto ref2 = wireReadout.OpDetGeoFromOpChannel(refHit2).GetCenter();

        // Check if the hits are in the same plane and continue if they are not
        if (!CheckOpHitPlane(OpHitPlane, refHit1, refHit2))
          continue;

        // If hit has already been clustered, skip
        if (ClusteredHits[*it2])
          continue;

        auto ref4 = TVector3(ref1.X(), ref1.Y(), ref1.Z()) - TVector3(ref2.X(), ref2.Y(), ref2.Z());
        if (ref4.Mag() < fOpFlashAlgoRad)
        {
          if (adjHit->PE() > hit->PE())
          {
            main_hit = false;
            sOpHitClustering += "Hit with PE > TriggerPE found: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(adjHit->PeakTime()) + "\n";

            // Reset the ClusteredHits values for the hits that have been added to the cluster
            for (auto it3 = AdjHitIdx.begin(); it3 != AdjHitIdx.end(); ++it3)
            {
              ClusteredHits[*it3] = false;
              sOpHitClustering += "Removing hit: CH " + ProducerUtils::str(OpHitVector[*it3]->OpChannel()) + " Time " + ProducerUtils::str(OpHitVector[*it3]->PeakTime()) + "\n";
            }
            break;
          }
          AdjHitVec.push_back(adjHit);
          AdjHitIdx.push_back(*it2);
          ClusteredHits[*it2] = true;
          sOpHitClustering += "Adding hit: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(adjHit->PeakTime()) + "\n";
        }
      }

      for (auto it4 = it - 1; it4 != sortingIndex.begin(); --it4)
      {
        // make sure we don't go out of bounds and the pointer is valid
        if (it == sortingIndex.begin())
          break;
        if (it4 == sortingIndex.begin())
          break;
        
        auto &adjHit = OpHitVector[*it4];

        if (std::abs(hit->PeakTime() - adjHit->PeakTime()) > fOpFlashAlgoMinTime)
          break;
        if (adjHit->PE() < fOpFlashAlgoPE)
          continue;

        int refHit2 = adjHit->OpChannel();
        auto ref2 = wireReadout.OpDetGeoFromOpChannel(refHit2).GetCenter();

        // Check if the hits are in the same plane and continue if they are not
        if (!CheckOpHitPlane(OpHitPlane, refHit1, refHit2))
          continue;

        // if hit has already been clustered, skip
        if (ClusteredHits[*it4])
          continue;

        auto ref4 = TVector3(ref1.X(), ref1.Y(), ref1.Z()) - TVector3(ref2.X(), ref2.Y(), ref2.Z());
        if (ref4.Mag() < fOpFlashAlgoRad)
        {
          if (adjHit->PE() > hit->PE())
          {
            main_hit = false;
            sOpHitClustering += "*** Hit with PE > TriggerPE found: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(adjHit->PeakTime()) + "\n";

            // Reset the ClusteredHits values for the hits that have been added to the cluster
            for (auto it5 = AdjHitIdx.begin(); it5 != AdjHitIdx.end(); ++it5)
            {
              ClusteredHits[*it5] = false;
              sOpHitClustering += "Removing hit: CH " + ProducerUtils::str(OpHitVector[*it5]->OpChannel()) + " Time " + ProducerUtils::str(OpHitVector[*it5]->PeakTime()) + "\n";
            }
            break;
          }
          AdjHitVec.push_back(adjHit);
          AdjHitIdx.push_back(*it4);
          ClusteredHits[*it4] = true;
          sOpHitClustering += "Adding hit: PE " + ProducerUtils::str(adjHit->PE()) + " CH " + ProducerUtils::str(adjHit->OpChannel()) + " Time " + ProducerUtils::str(adjHit->PeakTime()) + "\n";
        }
      }

      if (main_hit && int(AdjHitVec.size()) >= fOpFlashAlgoNHit)
      {
        // Store the original indices of the clustered hits
        OpHitClusters.push_back(std::move(AdjHitVec));
        OpHitClusterIdx.push_back(std::move(AdjHitIdx));
        sOpHitClustering += "Cluster size: " + ProducerUtils::str(int(OpHitClusters.back().size())) + "\n";
        ProducerUtils::PrintInColor(sOpHitClustering, ProducerUtils::GetColor("green"), "Debug");
      }
      else
      {
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

  int AdjOpHitsUtils::GetOpHitPlane(const art::Ptr<recob::OpHit> &hit)
  {
    auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter();
    if (fGeometry == "VD")
    {
      if (OpHitXYZ.X() > fXACathodeX -1 && OpHitXYZ.X() < fXACathodeX + 1)
        return 0;
      else if (OpHitXYZ.Y() > fXAMembraneY - 1 && OpHitXYZ.Y() < fXAMembraneY + 1)
        return 1;
      else if (OpHitXYZ.Y() > -fXAMembraneY - 1 && OpHitXYZ.Y() < -fXAMembraneY + 1)
        return 2;
      else if (OpHitXYZ.Z() > fXAFinalCapZ - 1 && OpHitXYZ.Z() < fXAFinalCapZ + 1)
        return 3;
      else if (OpHitXYZ.Z() > fXAStartCapZ - 1 && OpHitXYZ.Z() < fXAStartCapZ + 1)
        return 4;
    }
    else if (fGeometry == "HD")
    {
      if (OpHitXYZ.X() > 0)
        return 0;
      else if (OpHitXYZ.X() < 0)
        return 1;
    }
    else
    {
      ProducerUtils::PrintInColor("Geometry not recognized: Must be 'HD' or 'VD'", ProducerUtils::GetColor("red"), "Error");
    }
    return -1;
  }

  // Define a function and map to get the plane of the hit. If geometry is VD the number of planes is 5 (cathode, leftmembrane, rightmembrane, startcap, finalcap). If geometry is HD the number of planes is 2 (leftdrift, rightdrift).
  std::map<int, int> AdjOpHitsUtils::GetOpHitPlaneMap(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector)
  {
    std::map<int, int> OpHitPlane;
    for (const auto &hit : OpHitVector)
    {
      OpHitPlane[hit->OpChannel()] = GetOpHitPlane(hit);
    }
    return OpHitPlane;
  }

  // Define a function to check if the adjhit is in the same plane as the reference hit using the OpHitPlane map
  bool AdjOpHitsUtils::CheckOpHitPlane(std::map<int, int> OpHitPlane, int refHit, int adjHit)
  {
    if (OpHitPlane[refHit] == OpHitPlane[adjHit])
      return true;
    else
      return false;
  }

  // Function to calculate the Gaussian probability density function
  // double AdjOpHitsUtils::GaussianPDF(double x, double mean, double sigma)
  // {
  //   return exp(-0.5 * pow((x - mean) / sigma, 2)) / (sqrt(2 * M_PI) * sigma);
  // }

  // void AdjOpHitsUtils::CalcCentroid(std::vector<art::Ptr<recob::OpHit>> Hits, double x, double y, double z)
  // {
  //   const double sigma = fOpFlashAlgoRad; // Gaussian sigma (range in cm)

  //   // Initialize variables
  //   double maxLikelihood = 0.0;
  //   double bestY = 0.0;
  //   double bestZ = 0.0;

  //   // Loop over possible x positions
  //   double firstHitX = wireReadout.OpDetGeoFromOpChannel(Hits[0]->OpChannel()).GetCenter().X();
  //   double firstHitY = wireReadout.OpDetGeoFromOpChannel(Hits[0]->OpChannel()).GetCenter().Y();
  //   double firstHitZ = wireReadout.OpDetGeoFromOpChannel(Hits[0]->OpChannel()).GetCenter().Z();

  //   for (double yPos = firstHitY - fOpFlashAlgoRad; yPos <= firstHitY + fOpFlashAlgoRad; yPos += 5)
  //   {
  //     for (double zPos = firstHitZ - fOpFlashAlgoRad; zPos <= firstHitZ + fOpFlashAlgoRad; zPos += 5)
  //     {
  //       // Skipt the yPos and zPos that are outside the circle of radius sigma around the first hit
  //       if (pow(yPos - firstHitY, 2) + pow(zPos - firstHitZ, 2) > pow(sigma, 2))
  //         continue;
  //       double likelihood = 0.0;
  //       double sumY = 0.0;
  //       double sumZ = 0.0;
  //       int count = 0;

  //       // Loop over hits
  //       for (const auto &hit : Hits)
  //       {
  //         double hitYPos = wireReadout.OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Y();
  //         double hitZPos = wireReadout.OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Z();

  //         // Calculate the likelihood for the hit
  //         double hitLikelihood = GaussianPDF(hitYPos, yPos, sigma) * GaussianPDF(hitZPos, zPos, sigma);

  //         // Accumulate the likelihood and calculate the weighted sum of Y and Z coordinateslarsoft_v09_91_02/work/prodmarley_nue_cc_flat_radiological_decay0_dune10kt_1x2x6_centralAPA/solar_ana_flash_dune10kt_1x2x6.fcl
  //         likelihood += log(hitLikelihood);
  //         sumY += hitLikelihood * hitYPos;
  //         sumZ += hitLikelihood * hitZPos;
  //         count++;
  //       }

  //       // Check if the current likelihood is the maximum
  //       if (likelihood > maxLikelihood)
  //       {
  //         maxLikelihood = likelihood;
  //         bestY = sumY / count;
  //         bestZ = sumZ / count;
  //       }
  //     }
  //   }

  //   // Set the best 3D spacepoint
  //   x = firstHitX;
  //   y = bestY;
  //   z = bestZ;
  // }

} // namespace solar