// Vertical drift flash finder based on the relative 
// time and tiles location distances between the OpHits used
// to compose flashes considering the particular 3 detecting 
// XArapuca planes geometry of the Vertical Drift.
//
// Author: F. Marinho (2022)
// 
// Based on modification on José Soto's code for DP 
// *Upper volume only/adapt if lower volume is implemented


#ifndef OpFlashFinderVerticalDrift_H
#define OpFlashFinderVerticalDrift_H 1

#include "larcore/CoreUtils/ServiceUtil.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
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
 
  class OpFlashFinderVerticalDrift : public art::EDProducer{
  public:
 
    // Standard constructor and destructor for an ART module.
    explicit OpFlashFinderVerticalDrift(const fhicl::ParameterSet&);
    virtual ~OpFlashFinderVerticalDrift();

    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& pset);

    // The producer routine, called once per event. 
    void produce(art::Event&); 
    

    void RunFlashFinder(std::vector< recob::OpHit > const& HitVector,
                        std::vector< recob::OpFlash >&     FlashVector,
                        std::vector< std::vector< int > >& AssocList,
                        detinfo::DetectorClocksData const& ts,
                        float const&                       TrigCoinc);
    void AddHitContribution(recob::OpHit const&    currentHit,
                            double&                MaxTime,
                            double&                MinTime,
                            double&                AveTime,
                            double&                FastToTotal,
                            double&                AveAbsTime,
                            double&                TotalPE,
                            std::vector<double> & PEs);
    void GetHitGeometryInfo(recob::OpHit const&      currentHit,
                            std::vector< double >&   sumw,
                            std::vector< double >&   sumw2,
                            double&                  sumy, 
                            double&                  sumy2,
                            double&                  sumz, 
                            double&                  sumz2);
    double CalculateWidth(double const& sum, 
                          double const& sum_squared, 
                          double const& weights_sum);
    std::vector<int> getNeighbors(std::vector< recob::OpHit > const&       HitVector,
                                  int hitnumber,
                                  std::vector<bool> &processed,
                                  float initimecluster,
                                  std::vector<int> &sorted);
    void AssignHitsToFlash(std::vector< recob::OpHit > const&       HitVector,
                           std::vector< std::vector< int > >&       HitsPerFlash);
    void ConstructFlash(std::vector< int > const&          HitsPerFlashVec,
                        std::vector< recob::OpHit > const& HitVector,
                        std::vector< recob::OpFlash >&     FlashVector,
                        detinfo::DetectorClocksData const& ts,
                        float const&                       TrigCoinc);
  private:

    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpHit collection    
    const geo::GeometryCore * geometry;
    geo::WireReadoutGeom const * wireReadout;


    Double_t fMaximumTimeDistance;        // time parameter for selecting time neighbouring hits
    Double_t fMaximumTimeWindow;          // time parameter for flash composition wrt attributed initial flash time
    Double_t fTrigCoinc;                  // Beam frame coincidence criteria
    Double_t fLMy;                        // position variables for identifying in which detecting plan a sensor is located 
    Double_t fRMy;
    Double_t fCx;
    Double_t fCdyz;                       // distance criteria for position neighbouring hits according to the different combinations of sensor planes with hits
    Double_t fCMdz;
    Double_t fMMdz;


  };

} 

namespace opdet {
  DEFINE_ART_MODULE(OpFlashFinderVerticalDrift)
}

#endif 

namespace opdet {

  //----------------------------------------------------------------------------
  // Constructor
  OpFlashFinderVerticalDrift::OpFlashFinderVerticalDrift(const fhicl::ParameterSet & pset)
    : EDProducer{pset}
  {

    reconfigure(pset);

    produces< std::vector< recob::OpFlash > >();
    produces< art::Assns< recob::OpFlash, recob::OpHit > >();

    geometry = &*(art::ServiceHandle<geo::Geometry>());  
    wireReadout = &art::ServiceHandle<geo::WireReadout>()->Get();
    
//    geometry = *lar::providerFrom< geo::Geometry >();
//    wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

  }

  //----------------------------------------------------------------------------
  void OpFlashFinderVerticalDrift::reconfigure(fhicl::ParameterSet const& pset)
  {

    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get< std::string >("InputModule");
    //These values need to be stablished   
    fTrigCoinc      = pset.get< double >("TrigCoinc",1.0);
    fMaximumTimeDistance = pset.get< double >("MaximumTimeDistance",1.0);
    fMaximumTimeWindow   = pset.get< double >("MaximumTimeWindow",1.0);
    fLMy  = pset.get< double >("LMy",-700.0);
    fRMy = pset.get< double >("RMy",700.0);
    fCx = pset.get< double >("Cx",-320.0);
    fCdyz = pset.get< double >("Cdyz",400.0);
    fCMdz = pset.get< double >("CMdz",450.0);
    fMMdz = pset.get< double >("MMdz",450.0);

  }

  //----------------------------------------------------------------------------
  // Destructor
  OpFlashFinderVerticalDrift::~OpFlashFinderVerticalDrift() 
  {
  }
   
  //----------------------------------------------------------------------------
  void OpFlashFinderVerticalDrift::beginJob()
  {
  }

  //----------------------------------------------------------------------------
  void OpFlashFinderVerticalDrift::endJob()
  { 
  }

  //----------------------------------------------------------------------------
  void OpFlashFinderVerticalDrift::produce(art::Event& evt) 
  {

    // These are the storage pointers we will put in the event
    std::unique_ptr< std::vector< recob::OpFlash > > 
                      flashPtr(new std::vector< recob::OpFlash >);
    std::unique_ptr< art::Assns< recob::OpFlash, recob::OpHit > >  
                      assnPtr(new art::Assns< recob::OpFlash, recob::OpHit >);

    // This will keep track of what flashes will assoc to what ophits
    // at the end of processing
    std::vector< std::vector< int > > assocList;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    // Get OpHits from the event
    auto opHitHandle = evt.getHandle< std::vector< recob::OpHit > >(fInputModule);

    RunFlashFinder(*opHitHandle,
                   *flashPtr,
                   assocList,
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

  void OpFlashFinderVerticalDrift::RunFlashFinder(std::vector< recob::OpHit > const& HitVector,
                                                  std::vector< recob::OpFlash >&     FlashVector,
                                                  std::vector< std::vector< int > >& AssocList,
                                                  detinfo::DetectorClocksData const& ts,
                                                  float const&                       TrigCoinc)
  {
    // Now start to create flashes.
    // First, need vector to keep track of which hits belong to which flashes
    std::vector< std::vector< int > > HitsPerFlash;
    
    AssignHitsToFlash(HitVector,
                      HitsPerFlash);
    
    // Now we have all our hits assigned to a flash. 
    // Make the recob::OpFlash objects
    for (auto const& HitsPerFlashVec : HitsPerFlash)
      ConstructFlash(HitsPerFlashVec,
                     HitVector,
                     FlashVector,
                     ts,
                     TrigCoinc);

    // Finally, write the association list.
    // back_inserter tacks the result onto the end of AssocList
    for (auto& HitIndicesThisFlash : HitsPerFlash)
      AssocList.push_back(HitIndicesThisFlash);
    
  } // End RunFlashFinder

  //----------------------------------------------------------------------------
  void OpFlashFinderVerticalDrift::ConstructFlash(std::vector< int > const&          HitsPerFlashVec,
                                                  std::vector< recob::OpHit > const& HitVector,
                                                  std::vector< recob::OpFlash >&     FlashVector,
                                                  detinfo::DetectorClocksData const& ts,
                                                  float const&                       TrigCoinc)
  {
    double MaxTime = -1e9;
    double MinTime =  1e9;
    
    std::vector<double > PEs(geometry->NOpDets(),0.0);
    unsigned int Nplanes = wireReadout->Nplanes();
    //std::cout << "N planes: " << Nplanes << std::endl;
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

    // Empirical corrections to get the Frame right.
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
 void OpFlashFinderVerticalDrift::AddHitContribution(recob::OpHit const&    currentHit,
                          double&                MaxTime,
                          double&                MinTime,
                          double&                AveTime,
                          double&                FastToTotal,
                          double&                AveAbsTime,
                          double&                TotalPE,
                          std::vector<double >& PEs) {

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
    unsigned int thisOpDet = wireReadout->OpDetFromOpChannel(currentHit.OpChannel());
    PEs.at(thisOpDet)+=PEThisHit;

  }

  //----------------------------------------------------------------------------
  void OpFlashFinderVerticalDrift::GetHitGeometryInfo(recob::OpHit const&      currentHit,
                                                      std::vector< double >&   sumw,
                                                      std::vector< double >&   sumw2,
                                                      double&                  sumy, 
                                                      double&                  sumy2,
                                                      double&                  sumz, 
                                                      double&                  sumz2) {

    auto const xyz = wireReadout->OpDetGeoFromOpChannel(currentHit.OpChannel()).GetCenter();

    double PEThisHit = currentHit.PE();
    
    geo::TPCID tpc = geometry->FindTPCAtPosition(xyz);
    // if the point does not fall into any TPC,
    // it does not contribute to the average wire position
    if (tpc.isValid) {
      for (size_t p = 0; p != wireReadout->Nplanes(); ++p) {
        geo::PlaneID const planeID(tpc, p);
        unsigned int w = wireReadout->NearestWireID(xyz, planeID).Wire;
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
  double OpFlashFinderVerticalDrift::CalculateWidth(double const& sum, 
                        double const& sum_squared, 
                        double const& weights_sum) {

    if (sum_squared*weights_sum - sum*sum < 0) return 0;
    else return std::sqrt(sum_squared*weights_sum - sum*sum)/weights_sum;

  }


 std::vector<int> OpFlashFinderVerticalDrift::getNeighbors(std::vector< recob::OpHit > const&       HitVector,
                                                           int hitnumber, 
                                                           std::vector<bool> &processed, 
                                                           float initimecluster, 
                                                           std::vector<int> &sorted)
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
   float dx, dy, dz;
   float dt, dtmax;
   auto const xyz_hitnumber = wireReadout->OpDetGeoFromOpChannel(HitVector[sorted[hitnumber]].OpChannel()).GetCenter();

   //std::cout <<"Earliest: "<< HitVector[sorted[hitnumber]].PeakTime() <<" " << xyz_hitnumber.X() << " " << xyz_hitnumber.Y() << " " << xyz_hitnumber.Z() << std::endl;

   //Hits clustering for flash formation considering selection requirements according to which planes hits belong to (cathode & membrane) 

   for (int h=itTmin;h<itTmax;h++)
   {
     if(!processed[h])
     {
       auto const xyz_h = wireReadout->OpDetGeoFromOpChannel(HitVector[sorted[h]].OpChannel()).GetCenter();
      // std::cout << "Next:    " << HitVector[sorted[h]].PeakTime() << " "<< xyz_h.X() << " " << xyz_h.Y() << " " << xyz_h.Z() << std::endl; 
       dx  = (xyz_h - xyz_hitnumber).X();
       dy  = (xyz_h - xyz_hitnumber).Y();
       dz  = (xyz_h - xyz_hitnumber).Z();
       dt  = TMath::Abs(  HitVector[sorted[hitnumber]].PeakTime() - HitVector[sorted[h]].PeakTime() );
       dtmax  = TMath::Abs(initimecluster - HitVector[sorted[h]].PeakTime() );       
       double dist = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
       
       if(xyz_h.X() < fCx && xyz_hitnumber.X() < fCx){
         if(dt<fMaximumTimeDistance && dtmax<fMaximumTimeWindow && dist < fCdyz){ neighbors.push_back(h);processed[h]=true;}  	
         //std::cout <<" C C "<< fCdyz <<" "<< dist << std::endl; 

       }
       if(xyz_h.Y() > fRMy && xyz_hitnumber.Y() > fRMy){
         if(dt<fMaximumTimeDistance && dtmax<fMaximumTimeWindow && dist < fMMdz){ neighbors.push_back(h);processed[h]=true;}  
         //std::cout <<" R R "<< fMMdz <<" "<< dist << std::endl;
       }
       if(xyz_h.Y() < fLMy && xyz_hitnumber.Y() < fLMy){
         if(dt<fMaximumTimeDistance && dtmax<fMaximumTimeWindow && dist < fMMdz){ neighbors.push_back(h);processed[h]=true;}
         //std::cout <<" L L "<< fMMdz <<" "<< dist << std::endl;
       }
       if(TMath::Abs(xyz_h.Y()) > fRMy && xyz_hitnumber.X() < fCx){
         dist = TMath::Abs(dz);
         if(dt<fMaximumTimeDistance && dtmax<fMaximumTimeWindow && dist < fCMdz){ neighbors.push_back(h);processed[h]=true;} 
         //std::cout <<" C RL "<< fCMdz <<" "<< dist << std::endl;
       }
       if(xyz_h.X() < fCx && TMath::Abs(xyz_hitnumber.Y()) > fRMy){
         dist = TMath::Abs(dz);
         if(dt<fMaximumTimeDistance && dtmax<fMaximumTimeWindow && dist < fCMdz){ neighbors.push_back(h);processed[h]=true;}
         //std::cout <<" RL C "<< fCMdz <<" "<< dist << std::endl;
       }
       if(xyz_h.Y() < fLMy && xyz_hitnumber.Y() > fRMy){
         dist = TMath::Abs(dz);
         if(dt<fMaximumTimeDistance && dtmax<fMaximumTimeWindow && dist < fMMdz){ neighbors.push_back(h);processed[h]=true;}
         //std::cout <<" R L "<< fMMdz <<" "<< dist << std::endl;
       }
       if(xyz_h.Y() > fRMy && xyz_hitnumber.Y() < fLMy){
         dist = TMath::Abs(dz);
         if(dt<fMaximumTimeDistance && dtmax<fMaximumTimeWindow && dist < fMMdz){ neighbors.push_back(h);processed[h]=true;}
         //std::cout <<" L R "<< fMMdz <<" "<< dist << std::endl;
       }       

     }
   }
   return neighbors;
 }
// AssignHitsToFlash(HitVector, HitsPerFlash, MaximumDistance, MaximumTimeDistance, MaximumTimeWindow);
 void OpFlashFinderVerticalDrift::AssignHitsToFlash(std::vector< recob::OpHit > const&       HitVector,
                                                    std::vector< std::vector< int > >&       HitsPerFlash)
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
       //std::cout << "Start h :" << h << std::endl;  
       processed[h]=true;
       std::vector<int> N; N.clear();
       std::vector<int> nb =getNeighbors(HitVector,h, processed, HitVector[sorted[h]].PeakTime(),sorted); 
       N.insert( N.end(), nb.begin(), nb.end() );
       //std::cout << "Main go: " << nb.size() << std::endl;
       for (unsigned int i=0;i<N.size();i++)
       {
           std::vector<int> nb2 =getNeighbors(HitVector,N[i], processed, HitVector[sorted[h]].PeakTime(),sorted); 
           //std::cout << "Go: " << i << nb2.size() << std::endl;
           N.insert( N.end(), nb2.begin(), nb2.end() );
       }
           N.push_back(sorted[h]);
           //std::cout << "N size:  " << N.size() << std::endl;         
           HitsPerFlash.push_back(N);
     } 
   }
    
 } // End AssignHitsToFlash

  //----------------------------------------------------------------------------

} // namespace opdet
