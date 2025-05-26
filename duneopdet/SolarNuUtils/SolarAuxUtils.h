// ========================================================================================
// SolarAuxUtils.h
// This library is an auxiliary library for the SolarNuAna module.
// It contains utility functions that are used in the SolarNuAna module.
//
// @authors     : Sergio Manthey Corchado
// @created     : Jul, 2024
//=========================================================================================

#ifndef SolarAuxTool_h
#define SolarAuxTool_h

#include <cmath>
#include <iostream>
#include <vector>
#include <fcntl.h>
#include <iomanip>
#include <limits>

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

namespace solar
{
    class SolarAuxUtils
    {
    public:
        explicit SolarAuxUtils(fhicl::ParameterSet const &p);
        // void ComputeDistanceX(double &ClusterDistance, double t1, double t2);
        // void ComputeDistance3D(double &ClusterDistance, double t1, double y1, double z1, double t2, double y2, double z2);
        void FillMyMaps(std::map<int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle<std::vector<simb::MCTruth>> Hand);
        static long unsigned int WhichGeneratorType(std::vector<std::map<int, simb::MCParticle>> GeneratorParticles, int TrID);
        static bool InMyMap(int TrID, std::map<int, float> TrackIDMap);
        static bool InMyMap(int TrID, std::map<int, simb::MCParticle> ParMap);
        static void PrintInColor(std::string MyString, int Color, std::string Type = "Info");
        static void resume_stdout(int fd);
        static int supress_stdout();
        static int GetColor(std::string ColorName);
        static std::string str(int MyInt);
        static std::string str(bool MyBool);
        static std::string str(float MyFloat, int MyPrecision = 2);
        static std::string str(double MyDouble, int MyPrecision = 2);
        static std::string str(std::vector<int> MyVec);
        static std::string str(std::vector<float> MyVec, int MyPrecision = 2);
        static std::string str(std::vector<double> MyVec, int MyPrecision = 2);

    private:
        // From fhicl configuration
        const std::string fGeometry;
        const double fDetectorSizeX;
        const double fDetectorDriftTime;
    };

}
#endif