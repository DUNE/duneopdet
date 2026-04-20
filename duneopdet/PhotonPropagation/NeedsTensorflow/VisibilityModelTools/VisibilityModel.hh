/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : VisibilityModel.hh
 * @created     : Thursday Apr 16, 2026 04:54:16 CDT
 */

#ifndef VISIBILITYMODEL_HH

#define VISIBILITYMODEL_HH

#include <vector>
#include <map>
#include <string>

#include "larcorealg/CoreUtils/counter.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcore/Geometry/Geometry.h"

#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

namespace phot {
  namespace vismodel {
    /**
     * @brief Descriptor of visibility model
     */
    enum EVisModel {kSemiAnalytical = 0, kCompGraph = 1, kPhotonLibrary = 2};
    static const std::map<std::string, EVisModel> 
      vismodel_label = { 
        {"semianalytical", EVisModel::kSemiAnalytical},
        {"compgraph", EVisModel::kCompGraph}, 
        {"photonlibrary", EVisModel::kPhotonLibrary}
      };

    /**
     * @brief Converts visibility model name into the corresponding enum 
     *
     * @param label visibility model name 
     * @return visibility model enum
     */
    static EVisModel get_vismodel_code(const std::string& label) {
      const auto& v = vismodel_label.find(label);
      if ( v == vismodel_label.end() ) {
        fprintf(stderr, "PhotonVisibilityExport::get_vismodel_code() ERROR: "); 
        fprintf(stderr, "unknown label %s. Quit.\n", label.data()); 
        fprintf(stderr, "Available options are:\n"); 
        for (const auto &line : vismodel_label) {
          fprintf(stderr, "[%i] - %s\n", line.second, line.first.data()); 
        }
        exit(EXIT_FAILURE);
      }

      return v->second;
    }

    /**
     * @class VisModelSettings_t
     * @brief Visibility model configuration data structure
     */
    struct VisModelSettings_t {
      std::string fLabel = {};
      EVisModel fModelCode = {};
      fhicl::ParameterSet fFHiCLSettings = {};
      VisModelSettings_t() {}; 
      VisModelSettings_t(const fhicl::ParameterSet& pset) {
        fLabel = pset.get<std::string>("model_type"); 
        fModelCode = get_vismodel_code( fLabel ); 
        fFHiCLSettings = pset.get<fhicl::ParameterSet>("settings");
      }
    };

    /**
     * @brief checks if visibilities are valid
     *
     * This function loops over the optical detector visibilities and check
     * whether the values are smaller than 1. If `hard_mode` is `true`, the
     * method returns `false` at the first instance of `vis > 1`, otherwise
     * the method will "sanitize" the visibility capping it at 1.0.
     *
     * @param visibilities vector of optical detector visibilities
     * @param hard_mode return false if vis > 1
     */
    inline bool is_valid(std::vector<double>& visibilities, const bool hard_mode, const geo::Point_t& xspot, const geo::Geometry* geom)
    {
      int iopdet = 0; 

      for (auto& v : visibilities) {
        if (v > 1.0) {
          if (v > 1.1) {
            geo::OpDetGeo const& opDet = geom->OpDetGeoFromOpDet(iopdet);
            auto& center = opDet.GetCenter();
            const double d = sqrt((xspot - center).Mag2()); 
            printf("duneopdet::vismodel::is_valid(%.1f, %.1f, %.1f) WARNING: OpDet%i visibility %g (distance = %.1f cm).\n", 
                xspot.X(), xspot.Y(), xspot.Z(), iopdet, v, d); 
          }

          if (hard_mode) return false;
          else {
            if (v > 1.1) printf("vis OpDet %i %g -> 1.0\n", iopdet, v);
            v = 1.0; 
          }
        }
        iopdet++;
      }
      return true;
    }

    /**
     * @class IVisibilityModel
     * @brief Abstract visibility model interface class
     */
    class IVisibilityModel {
      public: 
        inline virtual ~IVisibilityModel() = default;
        virtual void BuildVisModel (const fhicl::ParameterSet& pset) = 0; 
        virtual void GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) = 0; 
        inline virtual void GetReflectedVisibilities(const geo::Point_t& xsport, std::vector<double>& vis_vec) {
          std::fill(vis_vec.begin(), vis_vec.end(), 0.0);
        }
        inline virtual bool DoReflections() const {return false;}
        virtual EVisModel GetModel() const = 0; 
    };

    /**
     * @brief Template visibility model class 
     *
     * @tparam T Visibility model engine
     */
    template<typename T> 
      class VisibilityModel : public IVisibilityModel {
        protected: 
          T* fEngine = nullptr; 
      };
  }
}

#endif /* end of include guard VISIBILITYMODEL_HH */

