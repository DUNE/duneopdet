/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : PhotonLibraryVisModel.hh
 * @created     : Thursday Apr 16, 2026 06:29:45 CDT
 */

#ifndef PHOTONLIBRARYVISMODEL_HH

#define PHOTONLIBRARYVISMODEL_HH

#include "duneopdet/PhotonPropagation/NeedsTensorflow/VisibilityModelTools/VisibilityModel.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"

namespace phot { 
  namespace vismodel {
    /**
     * @class PhotonLibraryVisModel
     * @brief Implementation of Photon-Library visibility model
     */
    class PhotonLibraryVisModel : public VisibilityModel<phot::PhotonVisibilityService> {
      public:
        inline explicit PhotonLibraryVisModel(const fhicl::ParameterSet& config) 
          : VisibilityModel<phot::PhotonVisibilityService>() {
            const fhicl::ParameterSet& pset = config.get<fhicl::ParameterSet>("settings");
            BuildVisModel(pset); 
          }

        inline EVisModel GetModel() const override{return kPhotonLibrary;}

        inline void BuildVisModel(const fhicl::ParameterSet& settings) override {
          printf("Creating PhotonLibrary visibility model\n"); 
          fEngineHolder = art::ServiceHandle<phot::PhotonVisibilityService>();
          fEngine = fEngineHolder.get();
        }

        inline void GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) override
        {
          auto mapped_vis = fEngine->GetAllVisibilities(xspot); 
          size_t iopdet = 0; 
          for (auto &vis : mapped_vis) {
            vis_vec.at(iopdet) = vis;
            iopdet++; 
          }
          return;
        }

      private: 
        art::ServiceHandle<phot::PhotonVisibilityService> fEngineHolder;
    };


  }
}



#endif /* end of include guard PHOTONLIBRARYVISMODEL_HH */

