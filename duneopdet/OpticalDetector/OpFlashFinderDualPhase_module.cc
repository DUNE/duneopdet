// -*- mode: c++; c-basic-offset: 2; -*-
//
// This module groups OpHits into OpFlashes,
// based in the time and space distribution 
// of the OpHits.
// by José Soto, CIEMAT (2019)
//


#ifndef OpFlashFinderDualPhase_H
#define OpFlashFinderDualPhase_H 1

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

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

namespace opdet {
 
  class OpFlashFinderDualPhase : public art::EDProducer{
  public:
 
    // Standard constructor and destructor for an ART module.
    explicit OpFlashFinderDualPhase(const fhicl::ParameterSet&);
    virtual ~OpFlashFinderDualPhase();

    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& pset);

    // The producer routine, called once per event. 
    void produce(art::Event&); 
    

    void RunFlashFinder(std::vector< recob::OpHit > const& HitVector,
                      std::vector< recob::OpFlash >&     FlashVector,
                      std::vector< std::vector< int > >& AssocList,
                      geo::GeometryCore const&           geom,
                      detinfo::DetectorClocksData const& ts,
                      float const&                       TrigCoinc);
    void AddHitContribution(recob::OpHit const&    currentHit,
                          double&                MaxTime,
                          double&                MinTime,
                          double&                AveTime,
                          double&                FastToTotal,
                          double&                AveAbsTime,
                          double&                TotalPE,
                          std::vector< double >& PEs);
    void GetHitGeometryInfo(recob::OpHit const&      currentHit,
                          geo::GeometryCore const& geom,
                          std::vector< double >&   sumw,
                          std::vector< double >&   sumw2,
                          double&                  sumy, 
                          double&                  sumy2,
                          double&                  sumz, 
                          double&                  sumz2);
    double CalculateWidth(double const& sum, 
                        double const& sum_squared, 
                        double const& weights_sum);
    std::vector<int> getNeighbors( std::vector< recob::OpHit > const&       HitVector,
                          int hitnumber,
                          std::vector<bool> &processed, float initimecluster, std::vector<int> &sorted,geo::GeometryCore const& geom);
    void AssignHitsToFlash( std::vector< recob::OpHit > const&       HitVector,
                         std::vector< std::vector< int > >&       HitsPerFlash,geo::GeometryCore const& geom);
  void ConstructFlash(std::vector< int > const&          HitsPerFlashVec,
                      std::vector< recob::OpHit > const& HitVector,
                      std::vector< recob::OpFlash >&     FlashVector,
                      geo::GeometryCore const&           geom,
                      detinfo::DetectorClocksData const& ts,
                      float const&                       TrigCoinc);
  private:

    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpHit collection
    

    Double_t fMaximumDistance;
    Double_t fMaximumTimeDistance;
    Double_t fMaximumTimeWindow;
    Double_t fTrigCoinc;

  };

} 

namespace opdet {
  DEFINE_ART_MODULE(OpFlashFinderDualPhase)
}

#endif 

namespace opdet {

  //----------------------------------------------------------------------------
  // Constructor
  OpFlashFinderDualPhase::OpFlashFinderDualPhase(const fhicl::ParameterSet & pset)
    : EDProducer{pset}
  {

    reconfigure(pset);

    produces< std::vector< recob::OpFlash > >();
    produces< art::Assns< recob::OpFlash, recob::OpHit > >();

  }

  //----------------------------------------------------------------------------
  void OpFlashFinderDualPhase::reconfigure(fhicl::ParameterSet const& pset)
  {

    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get< std::string >("InputModule");
      
    fTrigCoinc      = pset.get< double >("TrigCoinc");
    fMaximumDistance     = pset.get< double >("MaximumDistance");
    fMaximumTimeDistance = pset.get< double >("MaximumTimeDistance");
    fMaximumTimeWindow   = pset.get< double >("MaximumTimeWindow");

  }

  //----------------------------------------------------------------------------
  // Destructor
  OpFlashFinderDualPhase::~OpFlashFinderDualPhase() 
  {
  }
   
  //----------------------------------------------------------------------------
  void OpFlashFinderDualPhase::beginJob()
  {
  }

  //----------------------------------------------------------------------------
  void OpFlashFinderDualPhase::endJob()
  { 
  }

  //----------------------------------------------------------------------------
  void OpFlashFinderDualPhase::produce(art::Event& evt) 
  {

    // These are the storage pointers we will put in the event
    std::unique_ptr< std::vector< recob::OpFlash > > 
                      flashPtr(new std::vector< recob::OpFlash >);
    std::unique_ptr< art::Assns< recob::OpFlash, recob::OpHit > >  
                      assnPtr(new art::Assns< recob::OpFlash, recob::OpHit >);

    // This will keep track of what flashes will assoc to what ophits
    // at the end of processing
    std::vector< std::vector< int > > assocList;

    auto const& geometry(*lar::providerFrom< geo::Geometry >());

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    // Get OpHits from the event
    auto opHitHandle = evt.getHandle< std::vector< recob::OpHit > >(fInputModule);

    RunFlashFinder(*opHitHandle,
                   *flashPtr,
                   assocList,
                   geometry,
                   clockData,
                   fTrigCoinc);


   

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

  void OpFlashFinderDualPhase::RunFlashFinder(std::vector< recob::OpHit > const& HitVector,
                      std::vector< recob::OpFlash >&     FlashVector,
                      std::vector< std::vector< int > >& AssocList,
                      geo::GeometryCore const&           geom,
                      detinfo::DetectorClocksData const& ts,
                      float const&                       TrigCoinc) {

    // Now start to create flashes.
    // First, need vector to keep track of which hits belong to which flashes
    std::vector< std::vector< int > > HitsPerFlash;
    
    AssignHitsToFlash(HitVector,
                      HitsPerFlash,
                      geom);
    
    // Now we have all our hits assigned to a flash. 
    // Make the recob::OpFlash objects
    for (auto const& HitsPerFlashVec : HitsPerFlash)
      ConstructFlash(HitsPerFlashVec,
                     HitVector,
                     FlashVector,
                     geom,
                     ts,
                     TrigCoinc);

    // Finally, write the association list.
    // back_inserter tacks the result onto the end of AssocList
    for (auto& HitIndicesThisFlash : HitsPerFlash)
      AssocList.push_back(HitIndicesThisFlash);
    
  } // End RunFlashFinder

  //----------------------------------------------------------------------------
  void OpFlashFinderDualPhase::ConstructFlash(std::vector< int > const&          HitsPerFlashVec,
                      std::vector< recob::OpHit > const& HitVector,
                      std::vector< recob::OpFlash >&     FlashVector,
                      geo::GeometryCore const&           geom,
                      detinfo::DetectorClocksData const& ts,
                      float const&                       TrigCoinc) {

    double MaxTime = -1e9;
    double MinTime =  1e9;
    
    std::vector< double > PEs(geom.MaxOpChannel() + 1, 0.0);
    unsigned int Nplanes = geom.Nplanes();
    std::vector< double > sumw (Nplanes, 0.0);
    std::vector< double > sumw2(Nplanes, 0.0);
    
    double TotalPE     = 0;
    double AveTime     = 0;
    double AveAbsTime  = 0;
    double FastToTotal = 0;
    double sumy        = 0;
    double sumz        = 0;
    double sumy2       = 0;
    double sumz2       = 0;

    for (auto const& HitID : HitsPerFlashVec) {
      AddHitContribution(HitVector.at(HitID),
                         MaxTime,
                         MinTime,
                         AveTime,
                         FastToTotal,
                         AveAbsTime,
                         TotalPE,
                         PEs);
      GetHitGeometryInfo(HitVector.at(HitID),
                         geom,
                         sumw,
                         sumw2,
                         sumy, 
                         sumy2,
                         sumz, 
                         sumz2);
    }

    AveTime     /= TotalPE;
    AveAbsTime  /= TotalPE;
    FastToTotal /= TotalPE;
    
    double meany = sumy/TotalPE;
    double meanz = sumz/TotalPE;
    
    double widthy = CalculateWidth(sumy, sumy2, TotalPE);
    double widthz = CalculateWidth(sumz, sumz2, TotalPE);
    
    std::vector< double > WireCenters(Nplanes, 0.0);
    std::vector< double > WireWidths(Nplanes, 0.0);
    
    for (size_t p = 0; p != Nplanes; ++p) {
      WireCenters.at(p) = sumw.at(p)/TotalPE;
      WireWidths.at(p)  = CalculateWidth(sumw.at(p), sumw2.at(p), TotalPE);
    }

    // Emprical corrections to get the Frame right.
    // Eventual solution - remove frames
    int Frame = ts.OpticalClock().Frame(AveAbsTime - 18.1);
    if (Frame == 0) Frame = 1;
    
    int BeamFrame = ts.OpticalClock().Frame(ts.TriggerTime());
    bool InBeamFrame = false;
    if (!(ts.TriggerTime() < 0)) InBeamFrame = (Frame == BeamFrame);

    double TimeWidth = (MaxTime - MinTime)/2.0;

    int OnBeamTime = 0; 
    if (InBeamFrame && (std::abs(AveTime) < TrigCoinc)) OnBeamTime = 1;
    
    FlashVector.emplace_back(AveTime,
                             TimeWidth,
                             AveAbsTime,
                             Frame,
                             PEs, 
                             InBeamFrame,
                             OnBeamTime,
                             FastToTotal,
                             meany, 
                             widthy, 
                             meanz, 
                             widthz, 
                             WireCenters, 
                             WireWidths);

  }
 void OpFlashFinderDualPhase::AddHitContribution(recob::OpHit const&    currentHit,
                          double&                MaxTime,
                          double&                MinTime,
                          double&                AveTime,
                          double&                FastToTotal,
                          double&                AveAbsTime,
                          double&                TotalPE,
                          std::vector< double >& PEs) {

    double PEThisHit   = currentHit.PE();
    double TimeThisHit = currentHit.PeakTime();
    if (TimeThisHit > MaxTime) MaxTime = TimeThisHit;
    if (TimeThisHit < MinTime) MinTime = TimeThisHit;
    
    // These quantities for the flash are defined 
    // as the weighted averages over the hits
    AveTime     += PEThisHit*TimeThisHit;
    FastToTotal += PEThisHit*currentHit.FastToTotal();
    AveAbsTime  += PEThisHit*currentHit.PeakTimeAbs();
    
    // These are totals
    TotalPE     += PEThisHit;
    PEs.at(static_cast< unsigned int >(currentHit.OpChannel())) += PEThisHit;

  }

  //----------------------------------------------------------------------------
  void OpFlashFinderDualPhase::GetHitGeometryInfo(recob::OpHit const&      currentHit,
                          geo::GeometryCore const& geom,
                          std::vector< double >&   sumw,
                          std::vector< double >&   sumw2,
                          double&                  sumy, 
                          double&                  sumy2,
                          double&                  sumz, 
                          double&                  sumz2) {

    auto const xyz = geom.OpDetGeoFromOpChannel(currentHit.OpChannel()).GetCenter();
    double PEThisHit = currentHit.PE();
    
    geo::TPCID tpc = geom.FindTPCAtPosition(xyz);
    // if the point does not fall into any TPC,
    // it does not contribute to the average wire position
    if (tpc.isValid) {
      for (size_t p = 0; p != geom.Nplanes(); ++p) {
        geo::PlaneID const planeID(tpc, p);
        unsigned int w = geom.NearestWireID(xyz, planeID).Wire;
        sumw.at(p)  += PEThisHit*w;
        sumw2.at(p) += PEThisHit*w*w;
      }
    } // if we found the TPC
    sumy  += PEThisHit*xyz.Y();
    sumy2 += PEThisHit*xyz.Y()*xyz.Y();
    sumz  += PEThisHit*xyz.Z();
    sumz2 += PEThisHit*xyz.Z()*xyz.Z();
    
  }

  //----------------------------------------------------------------------------
  double OpFlashFinderDualPhase::CalculateWidth(double const& sum, 
                        double const& sum_squared, 
                        double const& weights_sum) {

    if (sum_squared*weights_sum - sum*sum < 0) return 0;
    else return std::sqrt(sum_squared*weights_sum - sum*sum)/weights_sum;

  }


 std::vector<int> OpFlashFinderDualPhase::getNeighbors( std::vector< recob::OpHit > const&       HitVector,int hitnumber, std::vector<bool> &processed, float initimecluster, std::vector<int> &sorted,geo::GeometryCore const& geom)
 {

  int itTmin=0;
  int itTmax=HitVector.size();
  for (int h=hitnumber;h>0;h--) if(HitVector[sorted[h]].PeakTime() < HitVector[sorted[hitnumber]].PeakTime()-fMaximumTimeDistance)
  {
    itTmin=h; break;
  }
  for (unsigned int h=hitnumber;h<HitVector.size();h++) if(HitVector[sorted[h]].PeakTime() > HitVector[sorted[hitnumber]].PeakTime()+fMaximumTimeDistance)
  {
    itTmax=h; break;
  }

   std::vector<int> neighbors;
   float dt, dtmax;
   auto const xyz_hitnumber = geom.OpDetGeoFromOpChannel(HitVector[sorted[hitnumber]].OpChannel()).GetCenter();

   for (int h=itTmin;h<itTmax;h++)
   {
     if(!processed[h])
     {
       dt  = TMath::Abs(  HitVector[sorted[hitnumber]].PeakTime() - HitVector[sorted[h]].PeakTime() );
       dtmax  = TMath::Abs(initimecluster - HitVector[sorted[h]].PeakTime() );
       auto const xyz_h = geom.OpDetGeoFromOpChannel(HitVector[sorted[h]].OpChannel()).GetCenter();
       float distance = (xyz_h - xyz_hitnumber).R();
//       std::cout << distance << " " << dt << " " << dtmax << std::endl; lets_pause();
       if(distance<fMaximumDistance*fMaximumDistance && dt<fMaximumTimeDistance && dtmax<fMaximumTimeWindow){ neighbors.push_back(h);processed[h]=true;}
     }
   }
   return neighbors;
 }
// AssignHitsToFlash(HitVector, HitsPerFlash, MaximumDistance, MaximumTimeDistance, MaximumTimeWindow);
 void OpFlashFinderDualPhase::AssignHitsToFlash( std::vector< recob::OpHit > const&       HitVector,
                         std::vector< std::vector< int > >&       HitsPerFlash,
                         geo::GeometryCore const& geom)
 {

   int ntothits=HitVector.size();
   std::size_t n(0);
   std::vector<int> sorted;
   sorted.resize(ntothits);
   std::generate(std::begin(sorted), std::end(sorted), [&]{ return n++; });
   std::sort(  std::begin(sorted), std::end(sorted), [&](int i1, int i2) { return HitVector[i1].PeakTime() < HitVector[i2].PeakTime(); }); 

   HitsPerFlash.clear();

   std::vector<bool> processed;
   processed.resize(HitVector.size(),false);

   for (unsigned int h=0;h<HitVector.size();h++)
   { 
     if(!processed[h])
     {
       processed[h]=true;
       std::vector<int> N; N.clear();
       std::vector<int> nb =getNeighbors(HitVector,h, processed, HitVector[sorted[h]].PeakTime(),sorted,geom); 
       N.insert( N.end(), nb.begin(), nb.end() );

       for (unsigned int i=0;i<N.size();i++)
       {
           std::vector<int> nb2 =getNeighbors(HitVector,N[i], processed, HitVector[sorted[h]].PeakTime(),sorted,geom); 
           N.insert( N.end(), nb2.begin(), nb2.end() );
       }
           N.push_back(sorted[h]);
           HitsPerFlash.push_back(N);
     } 
   }
    
 } // End AssignHitsToFlash

  //----------------------------------------------------------------------------

} // namespace opdet
