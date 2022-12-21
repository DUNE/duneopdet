// Dan Pershey
//
// OpHit clusterer, using y,z,t information to isolate OpHit's from a common
// origin.  Algorithm uses DBScan - but slightly modified, by calculating a
// cluster centroid, to allow for better clustering of mutiple-density
// clusters
//
// Using OpFlashFinder_module.cc as template
//


#ifndef OpSlicer_H
#define OpSlicer_H 1

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"


// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// ROOT includes

// C++ Includes
#include <vector>
#include <string>
#include <memory>
#include <limits>

namespace opdet {

  bool sortOpHitByTime(const art::Ptr<recob::OpHit> &left,
                       const art::Ptr<recob::OpHit> &right)
  {
    return left->PeakTime() < right->PeakTime();
  }

  class OpSlicer : public art::EDProducer{
  public:

    // Standard constructor and destructor for an ART module.
    explicit OpSlicer(const fhicl::ParameterSet&);
    virtual ~OpSlicer();

    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& pset);

    // The producer routine, called once per event.
    void produce(art::Event&);

  private:

    double YZDist(art::Ptr<recob::OpHit>, art::Ptr<recob::OpHit>) const;
    double Dist(art::Ptr<recob::OpHit>, art::Ptr<recob::OpHit>) const;
    int YZCentroid(std::vector< art::Ptr<recob::OpHit> >,
                   std::vector<int>) const;

    void GetHitYZ(std::vector< art::Ptr<recob::OpHit> >,
                  std::vector<int>,
                  std::vector<double>&,
                  std::vector<double>&) const;
    void ClusterHits(std::vector<art::Ptr<recob::OpHit>>,
                     std::vector<recob::OpFlash>&,
                     std::vector< std::vector<int> >&,
                     detinfo::DetectorClocksData const&)const;


    // The parameters we'll read from the .fcl file.
    std::string fOpHitModuleLabel; // Input tag for OpHit collection

    double fTScale;
    double fRScale;
    double fR0;
    double fBreakTime;
    int fMinN;

    double fTrigCoinc;

    art::ServiceHandle<geo::Geometry> geo;
  };

}

namespace opdet {
  DEFINE_ART_MODULE(OpSlicer)
}

#endif

namespace opdet {

  //--------------------------------------------------------------------------
  // Constructor
  OpSlicer::OpSlicer(const fhicl::ParameterSet & pset) : EDProducer{pset}
  {

    reconfigure(pset);

    produces< std::vector< recob::OpFlash > >();
    produces< art::Assns< recob::OpFlash, recob::OpHit > >();

  }

  //--------------------------------------------------------------------------
  void OpSlicer::reconfigure(fhicl::ParameterSet const& pset)
  {

    // Indicate that the Input Module comes from .fcl
    fOpHitModuleLabel = pset.get< std::string >("OpHitModuleLabel");

    fTScale    = pset.get<double>("TScale");
    fRScale    = pset.get<double>("RScale");
    fR0        = pset.get<double>("R0");
    fBreakTime = pset.get<double>("BreakTime");
    fMinN      = pset.get<int>   ("MinN");

    fTrigCoinc      = pset.get< double >("TrigCoinc");

  }

  //--------------------------------------------------------------------------
  // Destructor
  OpSlicer::~OpSlicer()
  {
  }
  //--------------------------------------------------------------------------
  void OpSlicer::beginJob()
  {
  }

  //--------------------------------------------------------------------------
  void OpSlicer::endJob()
  {
  }

  //--------------------------------------------------------------------------
  void OpSlicer::produce(art::Event& evt)
  {

    // These are the storage pointers we will put in the event
    std::unique_ptr< std::vector< recob::OpFlash > >
                      flashPtr(new std::vector< recob::OpFlash >);
    std::unique_ptr< art::Assns< recob::OpFlash, recob::OpHit > >
                      assnPtr(new art::Assns< recob::OpFlash, recob::OpHit >);

    // This will keep track of what flashes will assoc to what ophits
    // at the end of processing
    std::vector< std::vector< int > > assocList;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    // Get OpHits from the event
    auto opHitHandle = evt.getHandle< std::vector< recob::OpHit > >(fOpHitModuleLabel);

    std::vector< art::Ptr<recob::OpHit> > ohits;
    for (int i = 0; i < int(opHitHandle->size()); i++){
      art::Ptr<recob::OpHit> opHitPtr(opHitHandle,i);
      ohits.push_back(opHitPtr);
    }


    // Run the clustering
    ClusterHits(ohits, *flashPtr, assocList, clockData);

    // Make the associations which we noted we need
    for (size_t i = 0; i != assocList.size(); ++i)
    {
      art::PtrVector< recob::OpHit > opHitPtrVector;
      for (int const& hitIndex : assocList.at(i))
      {
        art::Ptr< recob::OpHit > opHitPtr(opHitHandle, hitIndex);
        opHitPtrVector.push_back(opHitPtr);
      }

      util::CreateAssn(*this, evt, *flashPtr, opHitPtrVector,
                       *(assnPtr.get()), i);
    }

    // Store results into the event
    evt.put(std::move(flashPtr));
    evt.put(std::move(assnPtr));
  }

  double OpSlicer::YZDist(art::Ptr<recob::OpHit> a,
                          art::Ptr<recob::OpHit> b) const
  {
    // First need to ask the geometry for the y-z location of both OpHits
    int channela = a->OpChannel();
    auto const xyza = geo->OpDetGeoFromOpChannel(channela).GetCenter();
    double ay = xyza.Y();
    double az = xyza.Z();

    int channelb = b->OpChannel();
    auto const xyzb = geo->OpDetGeoFromOpChannel(channelb).GetCenter();
    double by = xyzb.Y();
    double bz = xyzb.Z();

    double r2 = pow((ay-by),2) + pow((az-bz),2);

    double r = sqrt(r2);
    return r;
  }
  double OpSlicer::Dist(art::Ptr<recob::OpHit> a,
                        art::Ptr<recob::OpHit> b) const
  {
    double r = YZDist(a,b);
    return sqrt(pow((a->PeakTime()-b->PeakTime())/fTScale,2)+pow(r/fRScale,2));
  }

  int OpSlicer::YZCentroid(std::vector< art::Ptr<recob::OpHit> > ohits,
                           std::vector<int> curN) const
  {
    double maxDens = 0;
    int maxIdx = -1;
    for (int i = 0; i < int(curN.size()); i++){
      double dens = 0;
      for (int j = 0; j < int(curN.size()); j++){
        double r = YZDist(ohits[i],ohits[j]);
        r = sqrt(r)/100;
        dens += ohits[curN[j]]->PE() * exp(-r*r);
      }
      if (dens > maxDens){
        maxDens = dens;  maxIdx = curN[i];
      }
    }
    return maxIdx;
  }



  void OpSlicer::GetHitYZ(std::vector< art::Ptr<recob::OpHit> > ohits,
                          std::vector<int> curN,
                          std::vector<double> &ys,
                          std::vector<double> &zs) const
  {
    for (int cur : curN){
      art::Ptr<recob::OpHit> oh = ohits[cur];
      int channel = oh->OpChannel();
      auto const xyz = geo->OpDetGeoFromOpChannel(channel).GetCenter();
      ys.push_back(xyz.Y());
      zs.push_back(xyz.Z());
    }
  }


  void OpSlicer::ClusterHits(std::vector< art::Ptr<recob::OpHit> > ohits,
                             std::vector< recob::OpFlash>& oflashes,
                             std::vector< std::vector<int> >& assoc,
                             detinfo::DetectorClocksData const &ts)const
  {
    std::sort(ohits.begin(),ohits.end(),sortOpHitByTime);

    std::vector<bool> isClust;
    for (int i = 0; i < int(ohits.size()); i++) isClust.push_back(false);

    std::vector<int> neigh;

    for (int i = 0; i < int(ohits.size()); i++){

      if (isClust[i]) continue; // Don't base clusts off of clustered hits!
      neigh.erase(neigh.begin(),neigh.end()); // Start from scratch every time

      // check nearby hits in time for coincidence
      for (int j = i-1; j > 0; j--){
        if ( Dist(ohits[i],ohits[j]) < fR0 && !isClust[j] ){
          neigh.push_back(j);
        }
        if (abs(ohits[i]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2) break;
      }
      for (int j = i+1; j < int(ohits.size()); j++){
        if ( Dist(ohits[i],ohits[j]) < fR0 && !isClust[j] ){
          neigh.push_back(j);
        }
        if (abs(ohits[i]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2) break;
      }
      if (int(neigh.size())<fMinN) continue;
      neigh.erase(neigh.begin(),neigh.end());


      std::vector<int> cands;
      for (int j = i; j < int(ohits.size()); j++){
        if (isClust[j]) continue;
        if (ohits[j]->PeakTimeAbs()-ohits[i]->PeakTimeAbs()>fBreakTime) break;
        cands.push_back(j);
      }
      int centroidIdx = YZCentroid(ohits,cands);
      if (centroidIdx < 0) // No centroid found
        continue;
      
      art::Ptr<recob::OpHit> centroid = ohits[centroidIdx];
      neigh.push_back(centroidIdx);

      std::vector<int> curN;
      curN.push_back(centroidIdx);

      // check nearby hits in time for coincidence
      for (int j = centroidIdx-1; j > 0; j--){
        if ( Dist(ohits[centroidIdx],ohits[j]) < fR0 && !isClust[j] ){
          neigh.push_back(j);
          curN.push_back(j);
        }
        if (abs(ohits[centroidIdx]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2) break;
      }
      for (int j = centroidIdx+1; j < int(ohits.size()); j++){
        if ( Dist(ohits[centroidIdx],ohits[j]) < fR0 && !isClust[j] ){
          neigh.push_back(j);
          curN.push_back(j);
        }
        if (abs(ohits[centroidIdx]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2) break;
      }
      double totPE = 0;
      for (int idx : curN) totPE += ohits[idx]->PE();
      if (int(curN.size())<fMinN) continue;

      // Loop through neighboring hits, chck if it's a core hit
      while (neigh.size() > 0){
        std::vector<int> curNeigh;  curNeigh.push_back(neigh[0]);
        for (int j = neigh[0]-1; j > 0; j--){
          if ( Dist(ohits[neigh[0]],ohits[j]) < fR0 && !isClust[j] ){
            curNeigh.push_back(j);
          }
          if (abs(ohits[centroidIdx]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2)
            break;
        }
        for (int j = neigh[0]+1; j < int(ohits.size()); j++){
          if ( Dist(ohits[neigh[0]],ohits[j]) < fR0 && !isClust[j] ){
            curNeigh.push_back(j);
          }
          if (abs(ohits[centroidIdx]->PeakTimeAbs()-ohits[j]->PeakTimeAbs())>2)
            break;
        }
        // If this is a core point, add in all reachable hits to neighborhood
        if (int(curNeigh.size())>=fMinN){
          for (int cur : curNeigh){
            if (std::find(curN.begin(),curN.end(),cur)==curN.end()){
              curN.push_back(cur);
              if (Dist(ohits[centroidIdx],ohits[cur]) < fR0)
                neigh.push_back(cur);
            }
          }
        }
        neigh.erase(neigh.begin());
      }

      if (int(curN.size())<fMinN) continue;

      // Time to make the OpFlash;
      // Y-Z coordinates come from the centroid
      int channelcentroid = centroid->OpChannel();
      auto const xyzcentroid = geo->OpDetGeoFromOpChannel(channelcentroid).GetCenter();
      double yCenter = xyzcentroid.Y();
      double zCenter = xyzcentroid.Z();
      double tCenter = centroid->PeakTimeAbs();


      // Now that we have centroid coordinates, include ana delayed light
      for (int j = i; j < int(ohits.size()); j++){
        if (std::find(curN.begin(),curN.end(),j)!=curN.end()) continue;
        double r = YZDist(ohits[j],centroid);
        if ( r < fRScale*fR0){
          curN.push_back(j);
        }
        if (abs(ohits[j]->PeakTimeAbs()-tCenter)>fBreakTime) break;
      }

      double finE = 0;
      for (int idx : curN) finE += ohits[idx]->PE();

      // Grab the y-z information from the geometry
      std::vector<double> ys;
      std::vector<double> zs;
      GetHitYZ(ohits,curN,ys,zs);

      double minT = std::numeric_limits<double>::max(); double maxT = -std::numeric_limits<double>::max();
      double minY = 1e6; double maxY = -1e6;
      double minZ = 1e6; double maxZ = -1e6;

      std::vector<double> PEs (geo->MaxOpChannel() + 1,0.0);
      std::vector<double> PE2s (geo->MaxOpChannel() + 1,0.0);
      double fastToTotal = 0;
      for (int hIdx = 0; hIdx < int(ys.size()); hIdx++){
        int cIdx = curN[hIdx];

        minT = std::min(minT,ohits[cIdx]->PeakTimeAbs());
        maxT = std::max(maxT,ohits[cIdx]->PeakTimeAbs());
        minY = std::min(minY,ys[hIdx]);
        maxY = std::min(maxY,ys[hIdx]);
        minZ = std::min(minZ,zs[hIdx]);
        maxZ = std::min(maxZ,zs[hIdx]);
        PEs[ohits[cIdx]->OpChannel()] += ohits[cIdx]->PE();
        PE2s[ohits[cIdx]->OpChannel()] += ohits[cIdx]->PE();
        fastToTotal += ohits[hIdx]->FastToTotal();
      }
      double yWidth = maxY-minY;
      double zWidth = maxZ-minZ;

      double tot1 = 0;
      double tot2 = 0;
      for (double PE : PEs) tot1 += PE;
      for (double PE : PE2s) tot2 += PE;

      // From OpFlashAlg
      int Frame = ts.OpticalClock().Frame(tCenter - 18.1);
      if (Frame == 0) Frame = 1;

      int BeamFrame = ts.OpticalClock().Frame(ts.TriggerTime());
      bool InBeamFrame = false;
      if (!(ts.TriggerTime() < 0)) InBeamFrame = (Frame == BeamFrame);

      double tWidth = (maxT-minT)/2;

      int OnBeamTime = 0;
      if (InBeamFrame && (std::abs(tCenter) < fTrigCoinc)) OnBeamTime = 1;

      oflashes.emplace_back(tCenter,tWidth,tCenter,Frame,
                            PEs,InBeamFrame,OnBeamTime,fastToTotal,
                            yCenter,yWidth,zCenter,zWidth);
      assoc.emplace_back(curN);


      // And finally, indicate that current hits have been clustered
      for (int cur : curN) isClust[cur] = true;

    }
  }

} // namespace opdet
