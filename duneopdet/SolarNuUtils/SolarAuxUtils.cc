#include "SolarAuxUtils.h"

namespace solar
{
    SolarAuxUtils::SolarAuxUtils(fhicl::ParameterSet const &p)
        : fGeometry(p.get<std::string>("Geometry")),
          fDetectorSizeX(p.get<double>("DetectorSizeX")), // Changed type to double
          fDetectorDriftTime(p.get<double>("DetectorDriftTime"))
    {
    }

    void SolarAuxUtils::ComputeDistanceX(double &ClusterDistance, double t1, double t2)
    {
        ClusterDistance = 0;
        if (fGeometry == "HD")
        {
            ClusterDistance = abs(t1 - t2) * fDetectorSizeX / fDetectorDriftTime;
        }
        else if (fGeometry == "VD")
        {
            ClusterDistance = abs(t1 - t2) * fDetectorSizeX / (fDetectorDriftTime / 2);
        }
    }

    void SolarAuxUtils::ComputeDistance3D(double &ClusterDistance, double t1, double y1, double z1, double t2, double y2, double z2)
    {
        ClusterDistance = 0;
        double x12 = 0;
        ComputeDistanceX(x12, t1, t2);
        ClusterDistance = sqrt(pow(x12, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
    }

    // This function creates a terminal color printout
    void SolarAuxUtils::PrintInColor(std::string MyString, int Color, std::string Type)
    {
        if (Type == "Info")
        {
            mf::LogInfo("SolarNuAna") << "\033[" << Color << "m" << MyString << "\033[0m";
        }
        if (Type == "Debug")
        {
            mf::LogDebug("SolarNuAna") << "\033[" << Color << "m" << MyString << "\033[0m";
        }
        if (Type == "Error")
        {
            mf::LogError("SolarNuAna") << "\033[" << Color << "m" << MyString << "\033[0m";
        }
        return;
    }

    // ......................................................
    // This function returns an integer that corresponds to a given color name
    int SolarAuxUtils::GetColor(std::string ColorName)
    {
        if (ColorName == "black")
            return 30;
        else if (ColorName == "red")
            return 31;
        else if (ColorName == "green")
            return 32;
        else if (ColorName == "yellow")
            return 33;
        else if (ColorName == "blue")
            return 34;
        else if (ColorName == "magenta")
            return 35;
        else if (ColorName == "cyan")
            return 36;
        else if (ColorName == "white")
            return 37;
        else if (ColorName == "bright_black")
            return 90;
        else if (ColorName == "bright_red")
            return 91;
        else if (ColorName == "bright_green")
            return 92;
        else if (ColorName == "bright_yellow")
            return 93;
        else if (ColorName == "bright_blue")
            return 94;
        else if (ColorName == "bright_magenta")
            return 95;
        else if (ColorName == "bright_cyan")
            return 96;
        else if (ColorName == "bright_white")
            return 97;
        else
        {
            mf::LogError("SolarNuAna") << "Color " << ColorName << " not recognized. Returning white.";
            return 37;
        }
        return 0;
    }

    std::string SolarAuxUtils::str(int i)
    {
        std::stringstream ss;
        ss << i;
        return ss.str();
    }
    std::string SolarAuxUtils::str(bool i)
    {
        std::string j = "";
        if (i)
            j = "true";
        else
            j = "false";
        std::stringstream ss;
        ss << j;
        return ss.str();
    }
    std::string SolarAuxUtils::str(double i, int prec)
    {
        // Use prec to define the precision of the output in terms of number of decimal places
        std::stringstream ss;
        ss << std::fixed << std::setprecision(prec) << i;
        return ss.str();
    }
    std::string SolarAuxUtils::str(float i, int prec)
    {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(prec) << i;
        return ss.str();
    }
    std::string SolarAuxUtils::str(std::vector<int> i)
    {
        std::stringstream ss;
        for (int j = 0; j < int(i.size()); j++)
        {
            ss << i[j] << " ";
        }
        return ss.str();
    }
    std::string SolarAuxUtils::str(std::vector<double> i, int prec)
    {
        std::stringstream ss;
        for (int j = 0; j < int(i.size()); j++)
        {
            ss << std::fixed << std::setprecision(prec) << i[j] << " ";
        }
        return ss.str();
    }
    std::string SolarAuxUtils::str(std::vector<float> i, int prec)
    {
        std::stringstream ss;
        for (int j = 0; j < int(i.size()); j++)
        {
            ss << std::fixed << std::setprecision(prec) << i[j] << " ";
        }
        return ss.str();
    }

    int SolarAuxUtils::supress_stdout()
    {
        std::fflush(stdout);

        int ret = dup(1);
        int nullfd = open("/dev/null", O_WRONLY);
        // check nullfd for error omitted
        dup2(nullfd, 1);
        close(nullfd);

        return ret;
    }

    void SolarAuxUtils::resume_stdout(int fd)
    {
        std::fflush(stdout);
        dup2(fd, 1);
        close(fd);
    }
}