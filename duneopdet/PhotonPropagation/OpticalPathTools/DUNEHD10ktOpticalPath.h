// DUNE HD 10kt optical path tool

#ifndef DUNEHD10ktOpticalPath_H
#define DUNEHD10ktOpticalPath_H

#include "art/Utilities/ToolMacros.h" 
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include <iostream>

namespace phot {
    class DUNEHD10ktOpticalPath : public phot::OpticalPath {
    public:
        explicit DUNEHD10ktOpticalPath(fhicl::ParameterSet const& ps) {};
        ~DUNEHD10ktOpticalPath() noexcept override = default;

        const bool isOpDetVisible(geo::Point_t const& ScintPoint, geo::Point_t const& OpDetPoint) override {           
            // special case for DUNE-HD 10kt
            // check whether distance in drift direction > 1 drift distance (360.375 cm)
            if (std::abs(ScintPoint.X() - OpDetPoint.X()) > 360.375) return false;
            else return true; 
        }
    };
}

#endif