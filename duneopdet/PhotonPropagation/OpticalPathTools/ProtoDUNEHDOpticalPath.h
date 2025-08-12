// ProtoDUNE HD optical path tool

#ifndef ProtoDUNEHDOpticalPath_H
#define ProtoDUNEHDOpticalPath_H

#include "art/Utilities/ToolMacros.h" 
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

namespace phot {
    class ProtoDUNEHDOpticalPath : public phot::OpticalPath {
    public:
        explicit ProtoDUNEHDOpticalPath(fhicl::ParameterSet const& ps) {};
        ~ProtoDUNEHDOpticalPath() noexcept override = default;

        const bool isOpDetVisible(geo::Point_t const& ScintPoint, geo::Point_t const& OpDetPoint) override {           
            // special case for ProtoDUNE-HD
            // check x coordinate has same sign
            if ((ScintPoint.X() < 0.) != (OpDetPoint.X() < 0.)) return false;
            else return true; 
        }
    };
}

#endif