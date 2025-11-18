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
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/WireReadout.h"
#include "fhiclcpp/ParameterSet.h"

#include "dunecore/ProducerUtils/ProducerUtils.h"

#include "TH1I.h"
#include "TH1F.h"

namespace solar
{
    class AdjOpHitsUtils
    {
    public:
        struct FlashInfo
        {
            int Plane;
            int NHit;
            double Time;
            double TimeWidth;
            double TimeWeighted;
            double PE;
            double MaxPE;
            std::vector<double> PEperOpDet;
            double FastToTotal;
            double X;
            double Y;
            double Z;
            double XWidth;
            double YWidth;
            double ZWidth;
            double STD;
            std::vector<int> MainOpWaveform;
        };
        struct OpFlashes
        {
            std::vector<double> OpFlashPE;
            std::vector<double> OpFlashMaxPE;
            std::vector<double> OpFlashX;
            std::vector<double> OpFlashY;
            std::vector<double> OpFlashZ;
            std::vector<double> OpFlashT;
            std::vector<int> OpFlashNHit;
            std::vector<int> OpFlashGen;
            std::vector<double> OpFlashPur;
        };
        explicit AdjOpHitsUtils(fhicl::ParameterSet const &p);
        void CalcAdjOpHits(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector, std::vector<std::vector<art::Ptr<recob::OpHit>>> &OpHitClusters, std::vector<std::vector<int>> &OpHitClusterIdx, art::Event const &evt);
        void MakeFlashVector(std::vector<FlashInfo> &FlashVec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &OpHitClusters, art::Event const &evt);
        void FlashMatchResidual(float &Residual, std::vector<art::Ptr<recob::OpHit>> Hits, double x, double y, double z);
        float GetOpFlashPlaneSTD(const int Plane, const std::vector<float> varXY, const std::vector<float> varYZ, const std::vector<float> varXZ);
        int GetOpHitPlane(const art::Ptr<recob::OpHit> &hit, float buffer);
        std::map<int, int> GetOpHitPlaneMap(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector);
        bool CheckOpHitPlane(std::map<int, int> OpHitPlane, int refHit1, int refHit2);
        void GetOpHitWaveforms(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector, std::vector<art::Ptr<raw::OpDetWaveform>> &OpDetWvfVector, std::vector<bool> &OpHitWvfValid, art::Event const &evt);
        void GetOpHitWaveforms(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector, std::vector<art::Ptr<recob::OpWaveform>> &OpDetWvfVector, std::vector<bool> &OpHitWvfValid, art::Event const &evt);
        void GetOpHitSignal(const std::vector<art::Ptr<recob::OpHit>> &OpHitVector, std::vector<std::vector<int>> &OpDetWvfIntVector, std::vector<bool> &OpHitWvfValid, art::Event const &evt);

    private:
        art::ServiceHandle<geo::Geometry> geom;
        geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
        // From fhicl configuration
        const std::string fOpWaveformLabel;
        const std::string fOpHitLabel;
        const std::string fOpHitTimeVariable;
        const int fOpFlashAlgoNHit;
        const float fOpFlashAlgoMinTime;
        const float fOpFlashAlgoMaxTime;
        const float fOpFlashAlgoWeightedTime;
        const float fOpFlashAlgoRad;
        const float fOpFlashAlgoPE;
        const float fOpFlashAlgoTriggerPE;
        const float fOpFlashAlgoHotVertexThld;
        const bool  fOpFlashAlgoHitDuplicates;
        const float fXACathodeX;
        const float fXAMembraneY;
        const float fXAFinalCapZ;
        const float fXAStartCapZ;
    };
}
#endif
