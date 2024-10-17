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

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

namespace solar
{
    class SolarAuxUtils
    {
    public:
        explicit SolarAuxUtils(fhicl::ParameterSet const &p);
        void ComputeDistanceX(double &ClusterDistance, double t1, double t2);
        void ComputeDistance3D(double &ClusterDistance, double t1, double y1, double z1, double t2, double y2, double z2);
        static int GetColor(std::string ColorName);
        static void PrintInColor(std::string MyString, int Color, std::string Type = "Info");
        static void resume_stdout(int fd);
        static int supress_stdout();
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