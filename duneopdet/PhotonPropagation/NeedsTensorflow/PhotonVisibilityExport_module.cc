/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : PhotonVisibilityExport_module.cc
 * @created     : 2024-10-22 07:50
 */

#ifndef PHOTONVISIBILITYEXPORT_MODULE_CC

#define PHOTONVISIBILITYEXPORT_MODULE_CC

// ROOT includes
#include "TH1D.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TTree.h"
#include "TVectorF.h"
#include "TRandom3.h"
#include "TParameter.h"

// C++ includes
#include <cstdio>
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#include "larsimdnn/PhotonPropagation/TFLoaderTools/TFLoader.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larcorealg/CoreUtils/counter.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/type_traits.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace opdet {
  /**
   * @class PhotonVisibilityExport
   * @brief Export the visibility map of the detector volume and of each individual OpDet
   *
   * This module produces a visibility map of the detector volume and of each individual 
   * optical detector, both inside the TPC volume (photoVisMap) and in the buffer region
   * (photoVisMapBuffer). The TTree products can be converted into a 3D THn histogram with the 
   * macro in the VisibilityMapTools/ folder. 
   * In addition to the visibility tree, the module exports a tree containing the placement
   * of the individual optical detector and a tree with the size and placementes of the 
   * detector's TPCs. Finally, the tDimensions tree can be used to retrieve the overall size 
   * of the detector. 
   *
   * The voxel grid used for sampling the detector visibility can be recovered on the basis
   * of the three histograms hgrid[0-2].
   */
  class PhotonVisibilityExport : public art::EDAnalyzer 
  {
    public:
      enum EVisModel {kSemiAnalytical = 0, kCompGraph = 1, kPhotonLibrary = 2};
      static const std::map<TString, EVisModel> vismodel_label;      
      static EVisModel get_vismodel_code(const TString& label) {
        const auto& v = vismodel_label.find(label);
        if ( v == vismodel_label.end() ) {
          fprintf(stderr, "PhotonVisibilityExport::get_vismodel_code() ERROR: "); 
          fprintf(stderr, "unknown label %s. Quit.\n", label.Data()); 
          fprintf(stderr, "Available options are:\n"); 
          for (const auto &line : vismodel_label) {
            fprintf(stderr, "[%i] - %s\n", line.second, line.first.Data()); 
          }
          exit(EXIT_FAILURE);
        }

        return v->second;
      }

      struct VisModelSettings_t {
        TString fLabel = {};
        EVisModel fModelCode = {};
        fhicl::ParameterSet fFHiCLSettings = {};
        VisModelSettings_t() {}; 
        VisModelSettings_t(const fhicl::ParameterSet& pset) {
          fLabel = pset.get<std::string>("model_type"); 
          fModelCode = get_vismodel_code( fLabel ); 
          fFHiCLSettings = pset.get<fhicl::ParameterSet>("settings");
        }
      };
      
      static bool is_valid(const std::vector<double>& visibilities) ;

      class IVisibilityModel {
        public: 
          inline IVisibilityModel() {};
          inline virtual ~IVisibilityModel() = default;
          virtual void BuildVisModel (const fhicl::ParameterSet& pset) = 0; 
          virtual void GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) = 0; 
      };

      template<typename T> 
        class VisibilityModel : public IVisibilityModel {
          public:
            inline VisibilityModel() {}
            inline ~VisibilityModel() {}
            inline virtual void BuildVisMode(const fhicl::ParameterSet& pset) {}
            virtual void GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) {}; 

          protected: 
            T* fEngine = nullptr; 
        };

      class SemiAnalyticalVisModel : public VisibilityModel<phot::SemiAnalyticalModel> {
        public: 
          inline SemiAnalyticalVisModel() : VisibilityModel<phot::SemiAnalyticalModel>() {}
          inline ~SemiAnalyticalVisModel() {};
          inline void BuildVisModel(const fhicl::ParameterSet& pset) override {
            printf("Creating Semi-analytical visibility model\n");
            (pset.get<bool>("do_refl", false) ) ? 
              printf("Reflections included\n") : 
              printf("Reflections NOT included\n"); 
            fEngine = new phot::SemiAnalyticalModel(
                pset.get<fhicl::ParameterSet>("vuvhitspars"), 
                pset.get<fhicl::ParameterSet>("vishitspars"), 
                pset.get<bool>("do_refl", false), 
                pset.get<bool>("do_include_anode_refl", false),
                pset.get<bool>("do_include_xe_absorption", false)
                ); 
            fIncludeAnodeReflections = pset.get<bool>("do_include_anode_refl", false);

          }
          inline bool DoReflections() const {return fIncludeReflections;}
          inline bool DoIncludeAnodeReflections() const {return fIncludeAnodeReflections;}
          inline void SetAnodeReflections(const bool& anode_refl) {fIncludeAnodeReflections = anode_refl;}
          inline void EnableReflections(const bool& refl) {fIncludeReflections = refl;}
          inline void GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) override;
          inline void GetReflectedVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec)
          {
            fEngine->detectedReflectedVisibilities(vis_vec, xspot, fIncludeAnodeReflections);
            return;
          }
        private: 
          bool fIncludeReflections = false;
          bool fIncludeAnodeReflections = false;
      };

      class CompGraphVisModel : public VisibilityModel<phot::TFLoader> {
        public: 
          inline CompGraphVisModel() : VisibilityModel<phot::TFLoader>() {}
          inline void BuildVisModel(const fhicl::ParameterSet& settings) override {
            printf("Creating Computable-Graph visibility model\n"); 
            fEngineHolder = 
              art::make_tool<phot::TFLoader>(settings.get<fhicl::ParameterSet>("tfloaderpars"));
            fEngine = fEngineHolder.get();
            fEngine->Initialization();
          }
          inline void GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) override;

        private: 
          std::unique_ptr<phot::TFLoader> fEngineHolder;
      };
      
      class PhotonLibraryVisModel : public VisibilityModel<phot::PhotonVisibilityService> {
        public:
          inline PhotonLibraryVisModel() : VisibilityModel<phot::PhotonVisibilityService>() {}
          inline void BuildVisModel(const fhicl::ParameterSet& settings) override {
            printf("Creating PhotonLibrary visibility model\n"); 
            fEngine = art::ServiceHandle<phot::PhotonVisibilityService>().get(); 
          }
          inline void GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) override;

        private: 
          std::unique_ptr<phot::PhotonVisibilityService> fEngineHolder; 
      };

      explicit PhotonVisibilityExport(fhicl::ParameterSet const&);
      virtual ~PhotonVisibilityExport() {};

      void beginJob();
      void endJob() {}
      void analyze (const art::Event&);

    private: 
      size_t fNOpDets;
      std::vector<geo::Point_t> fOpDetCenter;

      double fVoxelSizeX = {};
      double fVoxelSizeY = {};
      double fVoxelSizeZ = {};

      int fNSamplings = 1;
      bool fIsDone = false;
      bool fIncludeBuffer = false;

      std::array<double, 3> fCryostatMin = {};
      std::array<double, 3> fCryostatMax = {};
      std::array<double, 3> fTPCMin = {}; 
      std::array<double, 3> fTPCMax = {};

      TH1D* fhGrid[3] = {}; 

      VisModelSettings_t fTPCVisModelConfig; 
      VisModelSettings_t fBufVisModelConfig;
      IVisibilityModel* fTPCVisModel;
      IVisibilityModel* fBufVisModel;

      void ExportTPCMap(); 
      void ExportOpDetMap(); 
      void ExportVoxelGrid(); 
      void ExportVisibility();

      IVisibilityModel* build_vismodel(const VisModelSettings_t&);

      const geo::Point_t SampleVoxel(const geo::Point_t& vxc) const; 

  };

} // close opdet namespace

namespace opdet {
  DEFINE_ART_MODULE(PhotonVisibilityExport)
}


#endif /* end of include guard VISMAPDUMP_MODULE_CC */

namespace opdet {
    inline void PhotonVisibilityExport::CompGraphVisModel::GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) {
      std::vector<double> xx(3, 0.0); 
      xx[0] = xspot.X(); xx[1] = xspot.Y(); xx[2] = xspot.Z();
      fEngine->Predict( xx ); 
      vis_vec = fEngine->GetPrediction(); 
      return;
    }

    inline void PhotonVisibilityExport::SemiAnalyticalVisModel::GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) {
      fEngine->detectedDirectVisibilities(vis_vec, xspot);
      return;
    }

    inline void PhotonVisibilityExport::PhotonLibraryVisModel::GetVisibilities(const geo::Point_t& xspot, std::vector<double>& vis_vec) {
      auto mapped_vis = fEngine->GetAllVisibilities(xspot); 

      size_t iopdet = 0; 
      for (auto &vis : mapped_vis) {
        vis_vec[iopdet] += vis;
        iopdet++; 
      }
      return;
    }

  const std::map<TString, PhotonVisibilityExport::EVisModel> 
    PhotonVisibilityExport::vismodel_label = { 
      {"semianalytical", EVisModel::kSemiAnalytical},
      {"compgraph", EVisModel::kCompGraph}, 
      {"photonlibrary", EVisModel::kPhotonLibrary}
    };

  PhotonVisibilityExport::PhotonVisibilityExport(fhicl::ParameterSet const& params) :
    art::EDAnalyzer{params},
    fVoxelSizeX{ params.get<double>("voxel_dx") }, 
    fVoxelSizeY{ params.get<double>("voxel_dy") }, 
    fVoxelSizeZ{ params.get<double>("voxel_dz") }, 
    fNSamplings{ params.get<int>("n_vis_samplings") },
    fIncludeBuffer{ params.get<bool>("do_include_buffer")}, 
    fTPCVisModelConfig{ params.get<fhicl::ParameterSet>("tpc_vis_model")}, 
    fBufVisModelConfig{ params.get<fhicl::ParameterSet>("buf_vis_model")} 
  { 
  }

  void PhotonVisibilityExport::beginJob() {
    printf("PhotonVisibilityExport::BeginJob() Building visibility models...\n"); 
    
    fTPCVisModel = build_vismodel( fTPCVisModelConfig );
    if (fIncludeBuffer) {
      fBufVisModel = build_vismodel( fBufVisModelConfig );
    }

    return;
  }

  PhotonVisibilityExport::IVisibilityModel* PhotonVisibilityExport::build_vismodel(const VisModelSettings_t& vmodel_config)
  {
    const auto& settings = vmodel_config.fFHiCLSettings;
    if (vmodel_config.fModelCode == kSemiAnalytical) {
       auto model = new SemiAnalyticalVisModel();
       model->BuildVisModel(settings);
       return model;
    }
    else if (vmodel_config.fModelCode == kCompGraph) {
      auto model = new CompGraphVisModel();
      model->BuildVisModel( settings ); 
      return model;
    }
    else if (vmodel_config.fModelCode == kPhotonLibrary) {
      auto model = new PhotonLibraryVisModel();
      model->BuildVisModel(settings); 
      return model; 
    }
    else {
      fprintf(stderr, "PhotonVisibilityExport::build_vismodel ERROR in reading visibility model config.\n"); 
      throw cet::exception("PhotonVisibilityExport") << "Wrong vis_model configuration " << vmodel_config.fLabel.Data();
      exit(EXIT_FAILURE);
    }
  }

  void PhotonVisibilityExport::analyze(const art::Event&) {

    if (fIsDone == false) {
      ExportOpDetMap();

      ExportTPCMap();

      ExportVisibility(); 
    }
    return;
  }

  void PhotonVisibilityExport::ExportVoxelGrid() {
    art::ServiceHandle<art::TFileService> tfs;

    const double voxelDim[3] = {fVoxelSizeX, fVoxelSizeY, fVoxelSizeZ};
    // define mesh points with an histogram helper to align tpcs and buffer
    // sampling coordinates
    for (int i=0; i<3; i++) {
      double   xc = 0.5*(fCryostatMax[i] + fCryostatMin[i]); 
      int    nbin = ceil( (fCryostatMax[i]-fCryostatMin[i]) / voxelDim[i] );
      double xmin = xc - nbin*0.5*voxelDim[i]; 
      double xmax = xc + nbin*0.5*voxelDim[i]; 
      fhGrid[i] = tfs->make<TH1D>(Form("hgrid%i", i), Form("mesh points axis %i", i), 
          nbin, xmin, xmax); 
    }
  }

  void PhotonVisibilityExport::ExportTPCMap() {
    printf("PhotonVisibilityExport::ExportTPCMap()...\n"); 
    // retrieve geometry
    const auto geom = art::ServiceHandle<geo::Geometry>();

    art::ServiceHandle< art::TFileService > tfs;

    float tpcH = 0;  
    float tpcW = 0;
    float tpcL = 0;
    float tpcPos[3]; 

    TTree* tTPC = tfs->make<TTree>("tpcMap", "tpcMap");
    tTPC->Branch("tpcH", &tpcH); 
    tTPC->Branch("tpcL", &tpcL); 
    tTPC->Branch("tpcW", &tpcW);
    tTPC->Branch("tpcPos", &tpcPos, "tpcPos[3]/F");

    // get cryostat dimensions 
    UInt_t icryo = 0; 
    for (geo::CryostatGeo const& cryo : geom->Iterate<geo::CryostatGeo>()) {
      printf("cryostat %u: [%g, %g, %g] - [%g, %g, %g]\n", icryo,
          cryo.MinX(), cryo.MinY(), cryo.MinZ(), cryo.MaxX(), cryo.MaxY(), cryo.MaxZ()); 
      fCryostatMin[0] = cryo.MinX();   fCryostatMax[0] = cryo.MaxX(); 
      fCryostatMin[1] = cryo.MinY();   fCryostatMax[1] = cryo.MaxY(); 
      fCryostatMin[2] = cryo.MinZ();   fCryostatMax[2] = cryo.MaxZ(); 
      icryo++; 
    }


    // loop over all TPCs to get the active volume dimensions
    for (geo::TPCGeo const& tpc : geom->Iterate<geo::TPCGeo>()) {
      auto point_center = tpc.GetCenter();
      double hlfW = tpc.ActiveHalfWidth (); tpcW = 2*hlfW;
      double hlfH = tpc.ActiveHalfHeight(); tpcH = 2*hlfH;
      double hlfL = tpc.ActiveHalfLength(); tpcL = 2*hlfL;

      const double hlfDim[3] = {hlfW, hlfH, hlfL}; 

      double x_min[3] = {0}; double x_max[3] = {0}; 
      point_center.GetCoordinates( x_min );
      point_center.GetCoordinates( x_max ); 
      for (int idim =0; idim<3; idim++) {
        x_min[idim] -= hlfDim[idim]; 
        x_max[idim] += hlfDim[idim]; 

        if (fTPCMin[idim] > x_min[idim]) fTPCMin[idim] = x_min[idim]; 
        if (fTPCMax[idim] < x_max[idim]) fTPCMax[idim] = x_max[idim];
      }
      tTPC->Fill();

      printf("TPC dimensions: hlfW %.2f - hlfH %.2f hlfL %.2f, ", 
          hlfW, hlfH, hlfL); 
      printf("TPC center:  %.2f -  %.2f - %.2f\n", 
          point_center.x(), point_center.y(), point_center.z()); 
    }

    printf("Cryostat min-max: [%g, %g, %g] - [%g, %g, %g]\n", 
        fCryostatMin[0], fCryostatMin[1], fCryostatMin[2], 
        fCryostatMax[0], fCryostatMax[1], fCryostatMax[2]);
    printf("TPC min-max: [%g, %g, %g] - [%g, %g, %g]\n", 
        fTPCMin[0], fTPCMin[1], fTPCMin[2], 
        fTPCMax[0], fTPCMax[1], fTPCMax[2]);

    TString descriptor = "";
    Double_t tmp[3]; 

    TTree* tDimensions = tfs->make<TTree>("tDimensions", "TPC and cryostat dimensions"); 
    tDimensions->Branch("descriptor", &descriptor); 
    tDimensions->Branch("dimension", &tmp, "dimension[3]/D"); 

    descriptor = "tpc_min"; 
    for (size_t i = 0; i < 3; i++) tmp[i] = fTPCMin[i]; 
    tDimensions->Fill(); 

    descriptor = "tpc_max"; 
    for (size_t i = 0; i < 3; i++) tmp[i] = fTPCMax[i]; 
    tDimensions->Fill(); 

    descriptor = "cryostat_min"; 
    for (size_t i = 0; i < 3; i++) tmp[i] = fCryostatMin[i]; 
    tDimensions->Fill(); 

    descriptor = "cryostat_max"; 
    for (size_t i = 0; i < 3; i++) tmp[i] = fCryostatMax[i]; 
    tDimensions->Fill(); 

    ExportVoxelGrid();
    
    return;
  }

  void PhotonVisibilityExport::ExportOpDetMap() {
    printf("PhotonVisibilityExport::ExportOpDetMap()...\n"); 
    // retrieve geometry
    const auto geom = art::ServiceHandle<geo::Geometry>();

    art::ServiceHandle< art::TFileService > tfs;

    float opDetH = 0; 
    float opDetW = 0;
    float opDetL = 0;
    float opDetPos[3]; 

    TTree* tOpDet = tfs->make<TTree>("opDetMap", "opDetMap");
    tOpDet->Branch("opDetH", &opDetH); 
    tOpDet->Branch("opDetL", &opDetL); 
    tOpDet->Branch("opDetW", &opDetW);
    tOpDet->Branch("opDetPos", &opDetPos, "opDetPos[3]/F");

    // store info from Geometry service
    fNOpDets = geom->NOpDets();
    for (size_t i : util::counter(fNOpDets)) {
      geo::OpDetGeo const& opDet = geom->OpDetGeoFromOpDet(i);
      auto center = opDet.GetCenter();
      center.GetCoordinates( opDetPos );
      opDetH = opDet.Height();
      opDetW = opDet.Width();
      opDetL = opDet.Length();

      tOpDet->Fill();
    }
    return;
  }

 const geo::Point_t PhotonVisibilityExport::SampleVoxel(const geo::Point_t& vxc) const
 {
   const geo::Vector_t delta(
       gRandom->Uniform(-0.5*fVoxelSizeX, 0.5*fVoxelSizeX), 
       gRandom->Uniform(-0.5*fVoxelSizeY, 0.5*fVoxelSizeY), 
       gRandom->Uniform(-0.5*fVoxelSizeZ, 0.5*fVoxelSizeZ) );  

   const geo::Point_t xspot = vxc + delta; 

   return xspot;
 }

 //bool PhotonVisibilityExport::is_valid(const phot::MappedCounts_t& visibilities) {
   //for (const auto& v : visibilities) {
     //if (v > 1.0) return false;
   //}
   //return true;
 //}

 bool PhotonVisibilityExport::is_valid(const std::vector<double>& visibilities) {
   for (const auto& v : visibilities) {
     if (v > 1.0) return false;
   }
   return true;
 }

 void PhotonVisibilityExport::ExportVisibility() {
   //const auto photonVisService = art::ServiceHandle<phot::PhotonVisibilityService>(); 
   // open file
   art::ServiceHandle< art::TFileService > tfs;

   Double_t point_[3];
   double total_visDirect = 0;
   double total_visReflct = 0;
   double* opDet_visDirect = new double[fNOpDets];
   double* opDet_visReflct = new double[fNOpDets]; 
   Double_t pointBuff_[3];
   double total_visDirectBuff = 0;
   double total_visReflctBuff = 0;
   double* opDet_visDirectBuff = new double[fNOpDets];
   double* opDet_visReflctBuff = new double[fNOpDets]; 

   // open tree
   TTree* tMap = tfs->make<TTree>("photoVisMap", "photoVisMap");
   tMap->SetNameTitle("photoVisMap", "photoVisMap"); 
   tMap->Branch("coords", &point_, "coords[3]/D");
   tMap->Branch("opDet_visDirect", opDet_visDirect, Form("opDet_visDirect[%ld]/D", fNOpDets));
   tMap->Branch("opDet_visReflct", opDet_visReflct, Form("opDet_visReflct[%ld]/D", fNOpDets));
   tMap->Branch("total_visDirect", &total_visDirect, "total_visDirect/D");
   tMap->Branch("total_visReflct", &total_visReflct, "total_visReflct/D");

   TTree* tMapBuffer = tfs->make<TTree>("photoVisMapBuffer", "photoVisMapBuffer"); 
   tMapBuffer->SetNameTitle("photoVisMapBuffer", "photoVisMapBuffer"); 
   tMapBuffer->Branch("coords", &pointBuff_, "coords[3]/D");
   tMapBuffer->Branch("opDet_visDirectBuff", opDet_visDirectBuff, Form("opDet_visDirectBuff[%ld]/D", fNOpDets));
   tMapBuffer->Branch("opDet_visReflctBuff", opDet_visReflctBuff, Form("opDet_visReflctBuff[%ld]/D", fNOpDets));
   tMapBuffer->Branch("total_visDirectBuff", &total_visDirectBuff, "total_visDirectBuff/D");
   tMapBuffer->Branch("total_visReflctBuff", &total_visReflctBuff, "total_visReflctBuff/D");

   std::vector<double> opdetvis_dir;
   std::vector<double> opdetvis_rfl;
   
   bool tpc_range_x = false, tpc_range_y = false, tpc_range_z = false; 

   // here we can loop over the points of the pre-defined grid
   double x_= 0, y_= 0, z_= 0;
   printf("starting loop: fHGrid histos: [%p, %p, %p]\n", 
       static_cast<void*>(fhGrid[0]), static_cast<void*>(fhGrid[1]), static_cast<void*>(fhGrid[2])); 

   for (int ix=1; ix<=fhGrid[0]->GetNbinsX(); ix++) {
       x_ = fhGrid[0]->GetBinCenter(ix); 
       if (x_> fTPCMin[0] && x_< fTPCMax[0]) tpc_range_x = true;
       else tpc_range_x = false; 

       for (int iy=1; iy<=fhGrid[1]->GetNbinsX(); iy++) {
         y_= fhGrid[1]->GetBinCenter(iy);            
         if (y_> fTPCMin[1] && y_< fTPCMax[1]) tpc_range_y = true; 
         else tpc_range_y = false; 

         for (int iz=1; iz<=fhGrid[2]->GetNbinsX(); iz++) {
           z_= fhGrid[2]->GetBinCenter(iz);            
           if (z_> fTPCMin[2] && z_< fTPCMax[2]) tpc_range_z = true; 
           else tpc_range_z = false; 

           opdetvis_dir.resize(fNOpDets, 0.0);
           opdetvis_rfl.resize(fNOpDets, 0.0); 

           if ( (tpc_range_x && tpc_range_y && tpc_range_z) == false) {
             if (fIncludeBuffer) {
               total_visDirectBuff = 0; 
               total_visReflctBuff = 0;
               for (size_t i = 0; i < fNOpDets; i++) {
                  opDet_visDirectBuff[i] = 0.;
                  opDet_visReflctBuff[i] = 0.; 
               }

               const geo::Point_t point( x_, y_, z_); 

               int nvalid = 0; 
               while (nvalid < fNSamplings) {
                 const auto xpoint = SampleVoxel( point ); 
                 fBufVisModel->GetVisibilities( xpoint, opdetvis_dir);
                 if ( is_valid(opdetvis_dir) == false ) {
                   continue; 
                 }

                 size_t iopdet = 0; 
                 for (auto &vis : opdetvis_dir) {
                   opDet_visDirectBuff[iopdet] += vis;
                   total_visDirectBuff += vis; 
                   iopdet++; 
                 }
                 nvalid++;
               }

               total_visDirectBuff = total_visDirectBuff / static_cast<double>(fNSamplings); 
               for (size_t i = 0; i < fNOpDets; i++) {
                 opDet_visDirectBuff[i] /= static_cast<double>(fNSamplings);
               }

               point.GetCoordinates( pointBuff_ ); 
               tMapBuffer->Fill(); 
             }
           } 
           else {
             const geo::Point_t vpoint( x_, y_, z_); 
             total_visDirect = 0.;
             total_visReflct = 0.;
             for (size_t i = 0; i < fNOpDets; i++) {
               opDet_visDirect[i] = 0.0; 
               opDet_visReflct[i] = 0.0; 
             }

             for (int i = 0; i < fNSamplings; i++) {
               const geo::Point_t xspot = SampleVoxel( vpoint ); 

               fTPCVisModel->GetVisibilities(xspot, opdetvis_dir);

               if (auto model = dynamic_cast<SemiAnalyticalVisModel*>(fTPCVisModel)) {
                 if (model->DoReflections()) {
                   model->GetReflectedVisibilities(xspot, opdetvis_rfl);
                 }
               }

               size_t iopdet = 0; 
               for (const auto &vis : opdetvis_dir) {
                 opDet_visDirect[iopdet] += vis;
                 total_visDirect += vis;
                 iopdet++;
               }

               iopdet = 0; 
               for (const auto &vis : opdetvis_rfl) {
                 opDet_visReflct[iopdet] += vis; 
                 total_visReflct += vis;
                 iopdet++;
               } 
             }

             total_visDirect = total_visDirect / static_cast<double>(fNSamplings); 
             total_visReflct = total_visReflct / static_cast<double>(fNSamplings); 

             for (size_t i = 0; i < fNOpDets; i++) {
               opDet_visDirect[i] /= static_cast<double>(fNSamplings);
               opDet_visReflct[i] /= static_cast<double>(fNSamplings);
             }

             // Fill tree
             vpoint.GetCoordinates( point_ ); 
             tMap->Fill();
           }
         }
       }
     }

   fIsDone = true;
   delete[] opDet_visReflct; 
   delete[] opDet_visDirect;
   delete[] opDet_visDirectBuff; 
   delete[] opDet_visReflctBuff;
   return;
 }
  
 }
 
