// ========================================================================================
// AdjOpHits.h
// It is used to find adjacent ophits in time and space to create flashes in the context of
// the SolarNuAna module and DUNE's solar neutrino analysis.
//
// @authors     : Sergio Manthey Corchado
// @created     : Apr, 2024
//=========================================================================================

#ifndef AdjOpHitsTool_h
#define AdjOpHitsTool_h

#include <cmath>
#include <numeric>
#include <iostream>
#include <vector>
#include <fcntl.h>

#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1I.h"
#include "TH1F.h"

namespace solar
{
    class AdjOpHitsUtils
    {
    public:
        struct FlashInfo
        {
            int NHit;
            double Time;
            double TimeWidth;
            double PE;
            double MaxPE;
            std::vector<double> PEperOpDet;
            double FastToTotal;
            double X;
            double Y;
            double Z;
            double YWidth;
            double ZWidth;
            double STD;
        };
        explicit AdjOpHitsUtils(fhicl::ParameterSet const &p);
        void CalcAdjOpHits(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector, std::vector<std::vector<art::Ptr<recob::OpHit>>> &OpHitClusters, std::vector<std::vector<int>> &OpHitClusterIdx);
        void MakeFlashVector(std::vector<FlashInfo> &FlashVec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &OpHitClusters, art::Event const &evt);
        void FlashMatchResidual(float &Residual, std::vector<art::Ptr<recob::OpHit>> Hits, double x, double y, double z);
        // void CalcCentroid(std::vector<art::Ptr<recob::OpHit>> Hits, double x, double y, double z);
        // double GaussianPDF(double x, double mean, double sigma);

    private:
        art::ServiceHandle<geo::Geometry> geo;
        // From fhicl configuration
        const std::string fGeometry;
        const int fOpFlashAlgoNHit;
        const float fOpFlashAlgoTime;
        const float fOpFlashAlgoRad;
        const float fOpFlashAlgoPE;
        const float fOpFlashAlgoTriggerPE;
        const float fOpFlashAlgoHotVertexThld;
        const double fDetectorSizeX;
        // const bool fOpFlashAlgoCentroid;
    };
}
#endif