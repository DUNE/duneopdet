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
      enum EVisModel {kSemiAnalytical = 0, kCompGraph = 1};

      PhotonVisibilityExport(const fhicl::ParameterSet&);
      virtual ~PhotonVisibilityExport() {};

      void beginJob();
      void endJob() {}
      void analyze (const art::Event&);

    private: 
      std::unique_ptr<phot::SemiAnalyticalModel> fVisibilityModel;
      std::unique_ptr<phot::TFLoader> fTFGenerator; 
      size_t fNOpDets;
      std::vector<geo::Point_t> fOpDetCenter;

      EVisModel kVisModel; 

      bool fDoReflectedLight = {};
      bool fIncludeAnodeReflections = {};
      bool fIncludeBuffer = {}; 
      bool fUseXeAbsorption = {};

      double fVoxelSizeX = {};
      double fVoxelSizeY = {};
      double fVoxelSizeZ = {};

      fhicl::ParameterSet fVUVHitsParams;
      fhicl::ParameterSet fVISHitsParams;
      fhicl::ParameterSet fTFLoaderPars;

      bool fIsDone = false;

      std::array<double, 3> fCryostatMin = {};
      std::array<double, 3> fCryostatMax = {};
      std::array<double, 3> fTPCMin = {}; 
      std::array<double, 3> fTPCMax = {};

      TH1D* fhGrid[3] = {}; 

      void ExportTPCMap(); 
      void ExportOpDetMap(); 
      void ExportVoxelGrid(); 
      void ExportVisibility();

  };

} // close opdet namespace

namespace opdet {
  DEFINE_ART_MODULE(PhotonVisibilityExport)
}


#endif /* end of include guard VISMAPDUMP_MODULE_CC */

namespace opdet {
  PhotonVisibilityExport::PhotonVisibilityExport(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset)
  {
    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    // Read inputs from fcl file
    fVoxelSizeX = pset.get<double>("voxel_dx", 10.0);
    fVoxelSizeY = pset.get<double>("voxel_dy", 10.0);
    fVoxelSizeZ = pset.get<double>("voxel_dz", 10.0);

    fDoReflectedLight = pset.get<bool>("do_refl", false); 
    fIncludeAnodeReflections = pset.get<bool>("do_include_anode_refl", false); 
    fIncludeBuffer = pset.get<bool>("do_include_buffer", false); 
    fUseXeAbsorption = pset.get<bool>("do_include_xe_absorption", false);

    fVUVHitsParams = pset.get<fhicl::ParameterSet>("vuvhitspars"); 
    fVISHitsParams = pset.get<fhicl::ParameterSet>("vishitspars"); 

    fTFLoaderPars = pset.get<fhicl::ParameterSet>("tfloaderpars");

    TString vis_model_str = pset.get<std::string>("vis_model"); 

    if (vis_model_str.Contains("compgraph")) {
      kVisModel = kCompGraph;
    }
    else if (vis_model_str.Contains("semianalytical")) {
      kVisModel = kSemiAnalytical;
    }
  }

  void PhotonVisibilityExport::beginJob() {
    // create the photo-detector visibility model
    if (kVisModel == kSemiAnalytical) {
      printf("Creating Semi-analytical visibility model\n");
      (fDoReflectedLight) ? 
        printf("Reflections included\n") : 
        printf("Reflections NOT included\n"); 
      fVisibilityModel = std::make_unique<phot::SemiAnalyticalModel>(
          fVUVHitsParams, fVISHitsParams, 
          fDoReflectedLight, fIncludeAnodeReflections, fUseXeAbsorption
          ); 
    }
    else if (kVisModel == kCompGraph) {
      fTFGenerator = art::make_tool<phot::TFLoader>(fTFLoaderPars);
      fTFGenerator->Initialization();
    }

    return;
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
    for (geo::CryostatGeo const& cryo : geom->Iterate<geo::CryostatGeo>()) {
      printf("cryostat: %u [%g, %g, %g] - [%g, %g, %g]\n", cryo.ID().getIndex(),
          cryo.MinX(), cryo.MinY(), cryo.MinZ(), cryo.MaxX(), cryo.MaxY(), cryo.MaxZ()); 
      fCryostatMin[0] = cryo.MinX();   fCryostatMax[0] = cryo.MaxX(); 
      fCryostatMin[1] = cryo.MinY();   fCryostatMax[1] = cryo.MaxY(); 
      fCryostatMin[2] = cryo.MinZ();   fCryostatMax[2] = cryo.MaxZ(); 
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
    // retrieve geometry
    const auto geom = art::ServiceHandle<geo::Geometry>();

    art::ServiceHandle< art::TFileService > tfs;

    float opDetH = 0; 
    float opDetW = 0;
    float opDetL = 0;
    float opDetPos[3]; 
    std::vector<size_t> opChannel; 

    TTree* tOpDet = tfs->make<TTree>("opDetMap", "opDetMap");
    tOpDet->Branch("opDetH", &opDetH); 
    tOpDet->Branch("opDetL", &opDetL); 
    tOpDet->Branch("opDetW", &opDetW);
    tOpDet->Branch("opDetPos", &opDetPos, "opDetPos[3]/F");
    tOpDet->Branch("opDetCh", &opChannel); 

    // store info from Geometry service
    fNOpDets = geom->NOpDets();
    for (size_t i : util::counter(fNOpDets)) {
      opChannel.clear(); 
      geo::OpDetGeo const& opDet = geom->OpDetGeoFromOpDet(i);
      auto center = opDet.GetCenter();
      center.GetCoordinates( opDetPos );
      opDetH = opDet.Height();
      opDetW = opDet.Width();
      opDetL = opDet.Length();

      size_t n_ch = geom->NOpHardwareChannels(i); 
      opChannel.resize(n_ch, 0); 
      for (size_t ich = 0; ich < n_ch; ich++) {
        opChannel.at(ich) = geom->OpChannel(i, ich);  
      }

      tOpDet->Fill();
    }
    return;
  }


 void PhotonVisibilityExport::ExportVisibility() {
   const auto photonVisService = art::ServiceHandle<phot::PhotonVisibilityService>(); 
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
   const double voxelDim[3] = {fVoxelSizeX, fVoxelSizeY, fVoxelSizeZ}; 

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

           opdetvis_dir.clear(); 
           opdetvis_rfl.clear(); 

           if ( (tpc_range_x && tpc_range_y && tpc_range_z) == false) {
             if (fIncludeBuffer) {
               total_visDirectBuff = 0; 
               total_visReflctBuff = 0;
               for (size_t i = 0; i < fNOpDets; i++) {
                  opDet_visDirectBuff[i] = 0.;
                  opDet_visReflctBuff[i] = 0.; 
               }

               const geo::Point_t point( x_, y_, z_); 
               auto mapped_vis = photonVisService->GetAllVisibilities(point); 
               size_t iopdet = 0; 
               for (auto &vis : mapped_vis) {
                 opDet_visDirectBuff[iopdet] = vis;
                 total_visDirectBuff += vis; 
                 iopdet++; 
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

             const double n_samplings = 5.0; 

             for (int i = 0; i < n_samplings; i++) {
               const geo::Vector_t delta(
                   gRandom->Uniform(-0.5*voxelDim[0], 0.5*voxelDim[0]), 
                   gRandom->Uniform(-0.5*voxelDim[1], 0.5*voxelDim[1]), 
                   gRandom->Uniform(-0.5*voxelDim[2], 0.5*voxelDim[2]) );  

               const geo::Point_t xspot = vpoint + delta; 

               if (kVisModel == kSemiAnalytical ) {
                 fVisibilityModel->detectedDirectVisibilities   (opdetvis_dir, xspot);
                 if (fDoReflectedLight) {
                   fVisibilityModel->detectedReflectedVisibilities(opdetvis_rfl, xspot, fIncludeAnodeReflections);
                 }
               }
               else if (kVisModel == kCompGraph ) {
                 std::vector<Double_t> pos_tmp {xspot.x(), xspot.y(), xspot.z()}; 
                 fTFGenerator->Predict( pos_tmp ); 
                 opdetvis_dir = fTFGenerator->GetPrediction(); 
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

             total_visDirect = total_visDirect / n_samplings; 
             total_visReflct = total_visReflct / n_samplings; 

             for (size_t i = 0; i < fNOpDets; i++) {
               opDet_visDirect[i] /= n_samplings;
               opDet_visReflct[i] /= n_samplings;
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
 
