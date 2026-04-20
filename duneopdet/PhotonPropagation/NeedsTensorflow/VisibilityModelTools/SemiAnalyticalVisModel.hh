/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : SemiAnalyticalVisModel
 * @created     : Thursday Apr 16, 2026 04:58:16 CDT
 */

#ifndef SEMIANALYTICALVISMODEL_HH

#define SEMIANALYTICALVISMODEL_HH

#include "duneopdet/PhotonPropagation/NeedsTensorflow/VisibilityModelTools/VisibilityModel.hh"

#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace phot {
  namespace vismodel {
    /**
     * @class SemiAnalyticalVisModel
     * @brief Implementation of Semi-analytical visibility model
     */
    class SemiAnalyticalVisModel : public VisibilityModel<phot::SemiAnalyticalModel> {
      public: 
        inline explicit SemiAnalyticalVisModel(const fhicl::ParameterSet& config) 
          : VisibilityModel<phot::SemiAnalyticalModel>() { 
            const fhicl::ParameterSet& pset = config.get<fhicl::ParameterSet>("settings");
            BuildVisModel( pset ); 
          }

        inline ~SemiAnalyticalVisModel() { delete fEngine; };

        inline EVisModel GetModel() const override {return kSemiAnalytical;}

        inline void BuildVisModel(const fhicl::ParameterSet& pset) override {
          printf("Creating Semi-analytical visibility model\n");
          (pset.get<bool>("do_refl", false) ) ? 
            printf("Reflections included\n") : 
            printf("Reflections NOT included\n"); 
          fEngine = new phot::SemiAnalyticalModel(
              pset.get<fhicl::ParameterSet>("vuvhitspars"), 
              pset.get<fhicl::ParameterSet>("vishitspars"),
              std::shared_ptr<phot::OpticalPath>(art::make_tool<phot::OpticalPath>(pset.get<fhicl::ParameterSet>("OpticalPathTool"))), 
              pset.get<bool>("do_refl", false), 
              pset.get<bool>("do_include_anode_refl", false),
              pset.get<bool>("do_include_xe_absorption", false)
              ); 
          fIncludeAnodeReflections = pset.get<bool>("do_include_anode_refl", false);
        }
        inline bool DoReflections() const override {return fIncludeReflections;}
        inline bool DoIncludeAnodeReflections() const {return fIncludeAnodeReflections;}
        inline void SetAnodeReflections(const bool& anode_refl) {fIncludeAnodeReflections = anode_refl;}
        inline void EnableReflections(const bool& refl) {fIncludeReflections = refl;}
        inline void GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) override
        {
          fEngine->detectedDirectVisibilities(vis_vec, xspot);
          return;
        }
        inline void GetReflectedVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) override
        {
          fEngine->detectedReflectedVisibilities(vis_vec, xspot, fIncludeAnodeReflections);
          return;
        }
      private: 
        bool fIncludeReflections = false;
        bool fIncludeAnodeReflections = false;
    };
  }
}

#endif /* end of include guard SEMIANALYTICALVISMODEL_HH */

