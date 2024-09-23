// ========================================================================================
// AdjHits.h
// This library is based on Michael Baird's DAQSimAna_module.
// It is used to find adjacent hits in time and space to create clusters in the context of
// the SolarNuAna module and DUNE's solar neutrino analysis.
// 
// @authors     : Sergio Manthey Corchado
// @created     : Apr, 2024 
//=========================================================================================

#ifndef AdjHitsTool_h
#define AdjHitsTool_h

#include <iostream>
#include <vector>
#include <fcntl.h>

#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "TH1I.h"
#include "TH1F.h"

namespace solar
{
    class AdjHitsUtils
    {
        public:
            explicit AdjHitsUtils( fhicl::ParameterSet const& p);
            void CalcAdjHits(std::vector<recob::Hit> MyVec, std::vector<std::vector<recob::Hit>> &Clusters, TH1I *MyHist, TH1F *ADCIntHist, bool HeavDebug);
            // void CalcAdjOpHits(std::vector<recob::OpHit> MyVec, std::vector<std::vector<recob::OpHit>> &Clusters, TH1I *MyHist, TH1F *ADCIntHist, bool HeavDebug);
        
        private:
            // From fhicl configuration
            const double fClusterAlgoTime;
            const int fClusterAlgoAdjChannel;
    };
}
#endif