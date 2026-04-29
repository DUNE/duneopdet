/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : ComputableGraphVisModel.hh
 * @created     : Thursday Apr 16, 2026 06:26:44 CDT
 */

#ifndef COMPUTABLEGRAPHVISMODEL_HH

#define COMPUTABLEGRAPHVISMODEL_HH

#include "duneopdet/PhotonPropagation/NeedsTensorflow/VisibilityModelTools/VisibilityModel.hh"
#include "larsimdnn/PhotonPropagation/TFLoaderTools/TFLoader.h"

namespace phot {
  namespace vismodel {
    /**
     * @class CompGraphVisModel
     * @brief Implementation of Computable Graph visibility model
     */
    class ComputableGraphVisModel : public VisibilityModel<phot::TFLoader> {
      public: 
        inline explicit ComputableGraphVisModel(const fhicl::ParameterSet& config) 
          : VisibilityModel<phot::TFLoader>() { 
            const fhicl::ParameterSet& pset = config.get<fhicl::ParameterSet>("settings");
            BuildVisModel( pset ); 
          }

        inline EVisModel GetModel() const override {return kCompGraph;}

        inline void BuildVisModel(const fhicl::ParameterSet& settings) override {
          printf("Creating Computable-Graph visibility model\n"); 
          fEngineHolder = 
            art::make_tool<phot::TFLoader>(settings.get<fhicl::ParameterSet>("tfloaderpars"));
          fEngine = fEngineHolder.get();
          fEngine->Initialization();
        }

        inline void GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) override 
        {
          std::vector<double> xx(3, 0.0); 
          xx[0] = xspot.X(); xx[1] = xspot.Y(); xx[2] = xspot.Z();
          fEngine->Predict( xx ); 
          vis_vec = fEngine->GetPrediction(); 
          return;
        }

      private: 
        std::unique_ptr<phot::TFLoader> fEngineHolder;
    };

  }
}


#endif /* end of include guard COMPUTABLEGRAPHVISMODEL_HH */

