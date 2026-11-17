// DUNE VD 10kt optical path tool

#ifndef DUNEVD10ktOpticalPath_H
#define DUNEVD10ktOpticalPath_H

#include "art/Utilities/ToolMacros.h" 
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include <iostream>

namespace phot {
    class DUNEVD10ktOpticalPath : public phot::OpticalPath {
    public:
        explicit DUNEVD10ktOpticalPath(fhicl::ParameterSet const& ps) {};
        ~DUNEVD10ktOpticalPath() noexcept override = default;

        const bool isOpDetVisible(geo::Point_t const& ScintPoint, geo::Point_t const& OpDetPoint) override {           
            // special case for DUNE-VD 10kt
            // Cathode PDs 
            // at X ~ 0cm, double-sided, visible from both drift volumes  
            if (std::abs(OpDetPoint.X()) < 10.0) return true;

            // Lateral/Membrane PDs [Long and Short]
            // at X > 350 and X < -350, only visible from within same drift volume
            else if (OpDetPoint.X() > 300.0 && ScintPoint.X() > 0.0) return true;
            else if (OpDetPoint.X() < -300.0 && ScintPoint.X() < 0.0) return true;

            else return false; 
        }
    };
}

#endif