#include "AdjOpHitsUtils.h"
#include "SolarAuxUtils.h"

namespace solar
{
  AdjOpHitsUtils::AdjOpHitsUtils(fhicl::ParameterSet const &p)
      : fGeometry(p.get<std::string>("Geometry")),
        fOpFlashAlgoNHit(p.get<int>("OpFlashAlgoNHit")),
        fOpFlashAlgoTime(p.get<float>("OpFlashAlgoTime")),
        fOpFlashAlgoRad(p.get<float>("OpFlashAlgoRad")),
        fOpFlashAlgoPE(p.get<float>("OpFlashAlgoPE")),
        fOpFlashAlgoTriggerPE(p.get<float>("OpFlashAlgoTriggerPE")),
        fOpFlashAlgoHotVertexThld(p.get<float>("OpFlashAlgoHotVertexThld")),
        fDetectorSizeX(p.get<double>("DetectorSizeX")) // Changed type to double
                                                       // fOpFlashAlgoCentroid(p.get<bool>("OpFlashAlgoCentroid"))
  {
  }
  void AdjOpHitsUtils::MakeFlashVector(std::vector<FlashInfo> &FlashVec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &OpHitClusters, art::Event const &evt)
  {
    if (fOpFlashAlgoHotVertexThld > 1)
    {
      SolarAuxUtils::PrintInColor("Hot vertex threshold must be between 0 and 1", SolarAuxUtils::GetColor("red"), "Error");
      return;
    }
    for (std::vector<art::Ptr<recob::OpHit>> Cluster : OpHitClusters)
    {
      if (!Cluster.empty())
      {
        std::stable_sort(Cluster.begin(), Cluster.end(), [](art::Ptr<recob::OpHit> a, art::Ptr<recob::OpHit> b)
                         { return a->PeakTime() < b->PeakTime(); });
      }
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
        auto OpHitXYZ = geo->OpDetGeoFromOpChannel(PDSHit->OpChannel()).GetCenter();
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
      std::vector<float> varYZ;
      for (art::Ptr<recob::OpHit> PDSHit : Cluster)
      {
        auto OpHitXYZ = geo->OpDetGeoFromOpChannel(PDSHit->OpChannel()).GetCenter();
        TimeWidth += (PDSHit->PeakTime() - Time) * (PDSHit->PeakTime() - Time);
        YWidth += (OpHitXYZ.Y() - Y) * (OpHitXYZ.Y() - Y);
        ZWidth += (OpHitXYZ.Z() - Z) * (OpHitXYZ.Z() - Z);
        varYZ.push_back(sqrt(pow(Y - OpHitXYZ.Y(), 2) + pow(Z - OpHitXYZ.Z(), 2)) * PDSHit->PE());
      }

      TimeWidth = sqrt(TimeWidth / Cluster.size());
      YWidth = sqrt(YWidth / Cluster.size());
      ZWidth = sqrt(ZWidth / Cluster.size());

      // Compute STD of varYZ
      float varYZmean = 0;
      for (float var : varYZ)
      {
        varYZmean += var;
      }
      varYZmean /= varYZ.size();
      float varYZstd = 0;
      for (float var : varYZ)
      {
        varYZstd += pow(var - varYZmean, 2);
      }
      varYZstd = sqrt(varYZstd / varYZ.size());
      STD = varYZstd;

      // Compute FastToTotal according to the #PEs arriving within the first 10% of the time window wrt the total #PEs
      for (art::Ptr<recob::OpHit> PDSHit : Cluster)
      {
        if (PDSHit->PeakTime() < Time + TimeWidth / 10)
          FastToTotal += PDSHit->PE();
      }
      FastToTotal /= PE;
      FlashVec.push_back(FlashInfo{NHit, Time, TimeWidth, PE, MaxPE, PEperOpDet, FastToTotal, X, Y, Z, YWidth, ZWidth, STD});
    }
    return;
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
      sOpHitClustering += "Trigger hit found: PE " + SolarAuxUtils::str(hit->PE()) + " CH " + SolarAuxUtils::str(hit->OpChannel()) + " Time " + SolarAuxUtils::str(hit->PeakTime()) + "\n";

      // Make use of the fact that the hits are sorted in time to only consider the hits that are adjacent in the vector up to a certain time range
      for (auto it2 = it + 1; it2 != sortingIndex.end(); ++it2)
      {
        // make sure we don't go out of bounds and the pointer is valid
        if (it2 == sortingIndex.end())
          break;

        auto &adjHit = OpHitVector[*it2]; // Update adjHit here

        if (std::abs(adjHit->PeakTime() - hit->PeakTime()) > fOpFlashAlgoTime)
          break;
        if (adjHit->PE() < fOpFlashAlgoPE)
          continue;

        int refHit1 = hit->OpChannel();
        int refHit2 = adjHit->OpChannel();
        auto ref1 = geo->OpDetGeoFromOpChannel(refHit1).GetCenter();
        auto ref2 = geo->OpDetGeoFromOpChannel(refHit2).GetCenter();

        // If sign of x is the same (HD), then the two hits are in the same drift volume and can be clustered, else skip.
        // If x is smaller than drift (VD), then one hit is in membrane XAs and we skip (awaiting better implementation!!).
        if (fGeometry == "HD")
        {
          if (ref1.X() * ref2.X() < 0)
            continue;
        }
        else if (fGeometry == "VD")
        {
          // Only use cathode hits in VD for now
          if (ref1.X() > -fDetectorSizeX || ref2.X() > -fDetectorSizeX)
            continue;
        }
        else // If the geometry is not HD or VD, skip
        {
          SolarAuxUtils::PrintInColor("Geometry not recognized: Must be 'HD' or 'VD'", SolarAuxUtils::GetColor("red"), "Error");
          continue;
        }
        // If hit has already been clustered, skip
        if (ClusteredHits[*it2])
          continue;

        auto ref4 = TVector3(ref1.X(), ref1.Y(), ref1.Z()) - TVector3(ref2.X(), ref2.Y(), ref2.Z());
        if (ref4.Mag() < fOpFlashAlgoRad)
        {
          if (adjHit->PE() > hit->PE())
          {
            main_hit = false;
            sOpHitClustering += "Hit with PE > TriggerPE found: PE " + SolarAuxUtils::str(adjHit->PE()) + " CH " + SolarAuxUtils::str(adjHit->OpChannel()) + " Time " + SolarAuxUtils::str(adjHit->PeakTime()) + "\n";

            // Reset the ClusteredHits values for the hits that have been added to the cluster
            for (auto it3 = AdjHitIdx.begin(); it3 != AdjHitIdx.end(); ++it3)
            {
              ClusteredHits[*it3] = false;
              sOpHitClustering += "Removing hit: CH " + SolarAuxUtils::str(OpHitVector[*it3]->OpChannel()) + " Time " + SolarAuxUtils::str(OpHitVector[*it3]->PeakTime()) + "\n";
            }
            break;
          }
          AdjHitVec.push_back(adjHit);
          AdjHitIdx.push_back(*it2);
          ClusteredHits[*it2] = true;
          sOpHitClustering += "Adding hit: PE " + SolarAuxUtils::str(adjHit->PE()) + " CH " + SolarAuxUtils::str(adjHit->OpChannel()) + " Time " + SolarAuxUtils::str(adjHit->PeakTime()) + "\n";
        }
      }

      for (auto it4 = it - 1; it4 != sortingIndex.begin(); --it4)
      {
        // make sure we don't go out of bounds and the pointer is valid
        if (it4 == sortingIndex.begin())
          break;
        auto &adjHit = OpHitVector[*it4];

        if (std::abs(adjHit->PeakTime() - hit->PeakTime()) > fOpFlashAlgoTime)
          break;
        if (adjHit->PE() < fOpFlashAlgoPE)
          continue;

        int refHit1 = hit->OpChannel();
        int refHit2 = adjHit->OpChannel();
        auto ref1 = geo->OpDetGeoFromOpChannel(refHit1).GetCenter();
        auto ref2 = geo->OpDetGeoFromOpChannel(refHit2).GetCenter();

        if (fGeometry == "HD")
        {
          if (ref1.X() * ref2.X() < 0)
            continue;
        }
        else if (fGeometry == "VD")
        {
          // Only use cothode hits in VD for now
          if (ref1.X() > -fDetectorSizeX || ref2.X() > -fDetectorSizeX)
            continue;
        }
        else // If the geometry is not HD or VD, skip
        {
          SolarAuxUtils::PrintInColor("Geometry not recognized: Must be 'HD' or 'VD'", SolarAuxUtils::GetColor("red"), "Error");
          continue;
        }

        // if hit has already been clustered, skip
        if (ClusteredHits[*it4])
          continue;

        auto ref4 = TVector3(ref1.X(), ref1.Y(), ref1.Z()) - TVector3(ref2.X(), ref2.Y(), ref2.Z());
        if (ref4.Mag() < fOpFlashAlgoRad)
        {
          if (adjHit->PE() > hit->PE())
          {
            main_hit = false;
            sOpHitClustering += "*** Hit with PE > TriggerPE found: PE " + SolarAuxUtils::str(adjHit->PE()) + " CH " + SolarAuxUtils::str(adjHit->OpChannel()) + " Time " + SolarAuxUtils::str(adjHit->PeakTime()) + "\n";

            // Reset the ClusteredHits values for the hits that have been added to the cluster
            for (auto it5 = AdjHitIdx.begin(); it5 != AdjHitIdx.end(); ++it5)
            {
              ClusteredHits[*it5] = false;
              sOpHitClustering += "Removing hit: CH " + SolarAuxUtils::str(OpHitVector[*it5]->OpChannel()) + " Time " + SolarAuxUtils::str(OpHitVector[*it5]->PeakTime()) + "\n";
            }
            break;
          }
          AdjHitVec.push_back(adjHit);
          AdjHitIdx.push_back(*it4);
          ClusteredHits[*it4] = true;
          sOpHitClustering += "Adding hit: PE " + SolarAuxUtils::str(adjHit->PE()) + " CH " + SolarAuxUtils::str(adjHit->OpChannel()) + " Time " + SolarAuxUtils::str(adjHit->PeakTime()) + "\n";
        }
      }

      if (main_hit && int(AdjHitVec.size()) >= fOpFlashAlgoNHit)
      {
        // Store the original indices of the clustered hits
        OpHitClusters.push_back(std::move(AdjHitVec));
        OpHitClusterIdx.push_back(std::move(AdjHitIdx));
        sOpHitClustering += "Cluster size: " + SolarAuxUtils::str(int(OpHitClusters.back().size())) + "\n";
        SolarAuxUtils::PrintInColor(sOpHitClustering, SolarAuxUtils::GetColor("green"), "Debug");
      }
      else
      {
        SolarAuxUtils::PrintInColor(sOpHitClustering, SolarAuxUtils::GetColor("red"), "Debug");
      }
    }
    return;
  }

  void AdjOpHitsUtils::FlashMatchResidual(float &Residual, std::vector<art::Ptr<recob::OpHit>> Hits, double x, double y, double z)
  {
    if (Hits.size() == 0)
    {
      Residual = 1e6;
      SolarAuxUtils::PrintInColor("Failed Residual Evaluation: Empty Flash!", SolarAuxUtils::GetColor("yellow"), "Error");
      return;
    }
    // Initialize variables
    Residual = 0;
    float PE = 0;

    // Find index of the hit with the highest PE
    int maxPEIdx = 0;
    float refHitPE = 0;
    for (unsigned int i = 0; i < Hits.size(); i++)
    {
      if (Hits[i]->PE() > Hits[maxPEIdx]->PE())
        maxPEIdx = i;
    }

    int hotHits = 0;
    for (const auto &hit : Hits)
    {
      if (hit->PE() >= fOpFlashAlgoHotVertexThld * Hits[maxPEIdx]->PE())
      {
        // Start with the first hit in the flash as reference point
        double firstHitY = geo->OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Y();
        double firstHitZ = geo->OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Z();

        // Get the first hit PE and calculate the squared distance and angle to the reference point
        float firstHitPE = hit->PE();
        float firstHitDistSq = pow(firstHitY - y, 2) + pow(firstHitZ - z, 2);
        float firstHitAngle = atan2(sqrt(firstHitDistSq), abs(x));

        // Calculate the expected PE value for the reference point based on the first hit PE and the squared distance + angle
        refHitPE += firstHitPE * (pow(x, 2) + firstHitDistSq) / pow(x, 2) / cos(firstHitAngle);
        // float refHitPE = firstHitPE * (pow(x, 2) + firstHitDistSq) / pow(x, 2);
        hotHits++;
      }
    }
    refHitPE /= hotHits;

    // Loop over all OpHits in the flash and compute the squared distance to the reference point
    for (const auto &hit : Hits)
    {
      double hitY = geo->OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Y();
      double hitZ = geo->OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Z();

      // The expected distribution of PE corresponds to a decrease of 1/r² with the distance from the flash center. Between adjacent OpHits, the expected decrease in charge has the form r²/(r²+d²)
      float hitDistSq = pow(hitY - y, 2) + pow(hitZ - z, 2);
      float hitAngle = atan2(sqrt(hitDistSq), abs(x));
      float predPE = refHitPE * cos(hitAngle) * pow(x, 2) / (pow(x, 2) + hitDistSq);
      // float predPE = refHitPE * pow(x, 2) / (pow(x, 2) + hitDistSq);

      Residual += pow(hit->PE() - predPE, 2);
      PE += hit->PE();
    }

    Residual /= PE;
    Residual /= float(Hits.size());
    std::string debug = "PE: " + SolarAuxUtils::str(PE) +
                        " X: " + SolarAuxUtils::str(x) +
                        " RefPE: " + SolarAuxUtils::str(refHitPE) +
                        " NHits: " + SolarAuxUtils::str(int(Hits.size())) +
                        " Residual: " + SolarAuxUtils::str(Residual);

    SolarAuxUtils::PrintInColor(debug, SolarAuxUtils::GetColor("yellow"), "Debug");
    return;
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
  //   double firstHitX = geo->OpDetGeoFromOpChannel(Hits[0]->OpChannel()).GetCenter().X();
  //   double firstHitY = geo->OpDetGeoFromOpChannel(Hits[0]->OpChannel()).GetCenter().Y();
  //   double firstHitZ = geo->OpDetGeoFromOpChannel(Hits[0]->OpChannel()).GetCenter().Z();

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
  //         double hitYPos = geo->OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Y();
  //         double hitZPos = geo->OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Z();

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