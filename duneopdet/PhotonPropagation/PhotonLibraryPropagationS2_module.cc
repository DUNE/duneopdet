////////////////////////////////////////////////////////////////////////
// Class:       PhotonLibraryPropagationS2
// Plugin Type: producer (art v2_05_00)
// File:        PhotonLibraryPropagationS2_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "nutools/RandomUtils/NuRandomService.h"


#include <memory>
#include <iostream>

#include "dune/PhotonPropagation/PhotonVisibilityServiceS2.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/LArG4/OpDetPhotonTable.h"
#include "larsim/LArG4/OpDetSensitiveDetector.h"
#include "larsim/LArG4/OpDetReadoutGeometry.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larsim/Simulation/PhotonVoxels.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimDriftedElectronCluster.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
// ROOT Includes
#include "TGeoManager.h"

namespace phot {
  class PhotonLibraryPropagationS2;
}


class phot::PhotonLibraryPropagationS2 : public art::EDProducer {
public:
  explicit PhotonLibraryPropagationS2(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonLibraryPropagationS2(PhotonLibraryPropagationS2 const &) = delete;
  PhotonLibraryPropagationS2(PhotonLibraryPropagationS2 &&) = delete;
  PhotonLibraryPropagationS2 & operator = (PhotonLibraryPropagationS2 const &) = delete;
  PhotonLibraryPropagationS2 & operator = (PhotonLibraryPropagationS2 &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob() override;
  void Print(std::map<int, std::map<int, int>>* StepPhotonTable);

private:

  bool fUseLitePhotons;

  std::string fDriftEModuleLabel;
  double fGain;//Number of photons created per drifted electron.

};


phot::PhotonLibraryPropagationS2::PhotonLibraryPropagationS2(fhicl::ParameterSet const & p)
{

  art::ServiceHandle<sim::LArG4Parameters> lgp;

  fUseLitePhotons=lgp->UseLitePhotons();
  if(fUseLitePhotons)
  {
     produces< std::vector<sim::OpDetBacktrackerRecord> >();
     produces< std::vector<sim::SimPhotonsLite> >();
  }
  else     produces< std::vector<sim::SimPhotons> >();  

  art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "photon",    p, "SeedPhoton");
  art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "scinttime", p, "SeedScintTime");  
  this->reconfigure(p);
}


void phot::PhotonLibraryPropagationS2::produce(art::Event & e)
{

  mf::LogInfo("PhotonLibraryPropagationS2") << "Producing S2 photons."<<std::endl;


  art::ServiceHandle<PhotonVisibilityServiceS2> pvs;


  if (!pvs->IncludeParPropTime())  mf::LogInfo("PhotonLibraryPropagationS2") << "Parametrized Propagation Time is not set up for S2. The propagation of S2 light will be instantaneous."<<std::endl;

  art::ServiceHandle<geo::Geometry> geom;
  //  const detinfo::LArProperties* larp = lar::providerFrom<detinfo::LArPropertiesService>();
  
  art::ServiceHandle<art::RandomNumberGenerator> rng;  
  CLHEP::HepRandomEngine &engine_photon = rng->getEngine("photon");
  CLHEP::RandPoissonQ randpoisphot(engine_photon);
  CLHEP::HepRandomEngine &engine_scinttime = rng->getEngine("scinttime");
  CLHEP::RandFlat randflatscinttime(engine_scinttime);


    larg4::OpDetPhotonTable::Instance()->ClearTable(geom->NOpDets());

  // Get the pointer to the fast scintillation table
  // larg4::OpDetPhotonTable * fst = larg4::OpDetPhotonTable::Instance();
  larg4::OpDetPhotonTable* litefst = larg4::OpDetPhotonTable::Instance();
  
  const size_t NOpChannels = pvs->NOpChannels();
  double nphot;

//  std::unique_ptr< std::vector<sim::SimPhotons>  >               PhotonCol                  (new std::vector<sim::SimPhotons>);
  std::unique_ptr< std::vector<sim::SimPhotonsLite>  >           LitePhotonCol              (new std::vector<sim::SimPhotonsLite>);
  std::unique_ptr< std::vector< sim::OpDetBacktrackerRecord > >  cOpDetBacktrackerRecordCol (new std::vector<sim::OpDetBacktrackerRecord>);

  art::ValidHandle< std::vector<sim::SimDriftedElectronCluster> > ElectronClusters_handle = e.getValidHandle<std::vector<sim::SimDriftedElectronCluster>>(fDriftEModuleLabel);

  //std::cout << "... GETTING ELECTRON CLUSTER HANDLE"<<std::endl;


  // int counter=0;

  if(fUseLitePhotons)  mf::LogDebug("PhotonLibraryPropagationS2") << "Creating S2 photons as a SimPhotonsLite data product."<<std::endl;
  else  mf::LogError("PhotonLibraryPropagationS2") << "Error creating S2 light, SimPhotons data product is not supported."<<std::endl;

  int counter=0;

  for (sim::SimDriftedElectronCluster const& ElectronCluster: *ElectronClusters_handle)
  {  
  //std::cout << "reading a cluster " << counter <<std::endl; //counter ++;
  int detectedcounter=0;

  double const xyz[3] = { ElectronCluster.FinalPositionX(), ElectronCluster.FinalPositionY(), ElectronCluster.FinalPositionZ() };
  float const* Visibilities = pvs->GetAllVisibilities(xyz);
  if(!Visibilities)
  continue;

  TF1 *ParPropTimeTF1 = nullptr;

  if(pvs->IncludeParPropTime())
  {
    ParPropTimeTF1 = pvs->GetTimingTF1(xyz);
  }

  //std::cout <<"\tPosition> " <<ElectronCluster.FinalPositionX()<<" " << ElectronCluster.FinalPositionY()<<" " <<ElectronCluster.FinalPositionZ() <<std::endl;  
  //std::cout <<"\tTiming> " <<ElectronCluster.getTime()<<std::endl;
  //std::cout <<"\tEnergy> " <<ElectronCluster.getEnergy()<<std::endl;
  //std::cout <<"\telectrons> " <<ElectronCluster.getNumberOfElectrons()<<std::endl;
  //std::cout <<"\tTrack> " <<ElectronCluster.getTrackID()<<std::endl;  
      
  nphot =ElectronCluster.NumberOfElectrons()*fGain; // # photons generated in the Gas Ar phase
  //mf::LogInfo("PhotonLibraryPropagationS2") << ElectronCluster.getNumberOfElectrons() <<" electron arrives."<< std::endl;
  //mf::LogInfo("PhotonLibraryPropagationS2") << nphot <<" S2 photons have been created in " << ElectronCluster.FinalPositionX() << " " << ElectronCluster.FinalPositionY() << " "<< ElectronCluster.FinalPositionZ()<< " at time " << ElectronCluster.getTime() << std::endl;

  //std::cout << ElectronCluster.getNumberOfElectrons() << " electrons." << std::endl;

  if(!Visibilities)
  {
  }
  else
  {  
    //std::cout <<"\t\tVisibilities loaded. Let's calculate the number of detected photons. OpChannels: " <<NOpChannels<<std::endl;  

    std::map<int, int> DetectedNum;

    for(size_t OpDet=0; OpDet!=NOpChannels; OpDet++)
    {
      G4int DetThisPMT = G4int(randpoisphot.fire(Visibilities[OpDet] * nphot));
      //mf::LogDebug("PhotonLibraryPropagationS2") << DetThisPMT <<" photons arrive to OpChannel " << OpDet<<std::endl;
      //std::cout <<"\t\tVisibilities loaded. Let's calculate the number of detected photons. OpChannels: " <<Visibilities[OpDet] << " "<< nphot << " " << Visibilities[OpDet] * nphot << std::endl;
      if(DetThisPMT>0) 
      {
        DetectedNum[OpDet]=DetThisPMT;
        //std::cout <<"\t\t\t#photons per pmt"<<OpDet<<": "<<DetThisPMT<<std::endl;
      }
    }
      // Now we run through each PMT figuring out num of detected photons
    if(fUseLitePhotons)
    {
      //std::cout << "\t\tLet's create our StepPhotonTable Map." << std::endl;  

      std::map<int, std::map<int, int>> StepPhotonTable;
      // And then add these to the total collection for the event  
      Print(&StepPhotonTable);

      for(std::map<int,int>::const_iterator itdetphot = DetectedNum.begin();
        itdetphot!=DetectedNum.end(); ++itdetphot)
      {

        //std::cout << "... iterating! over " <<itdetphot->second <<" photons. "<<  std::endl;
        std::map<int, int>  StepPhotons;

/*
        float maxarrivaltimerange=0.0;

        if(pvs->IncludeParPropTime()) maxarrivaltimerange+=static_cast<int>(10*PropParameters[itdetphot->first][0]);

        if(pvs->IncludeParPropTime()&&!pvs->IncludeMCParPropTime()&&itdetphot->second>50)
        {

          G4double deltaTime = ElectronCluster.getTime();
          //just normalizing the TF1 by the number of photons
          double integral=0;
          int ccounter=0;

          if(pvs->GetfParPropTime_nonanalyticalfunction())
          {
            for(size_t i=1; i<pvs->ParPropTimeNpar();i++)
            {
              ParPropTimeFunctionTF1->SetParameter(i, PropParameters[itdetphot->first][i]);
            }
             ParPropTimeFunctionTF1->SetRange(PropParameters[itdetphot->first][0],maxarrivaltimerange);

            integral=ParPropTimeFunctionTF1->Integral(PropParameters[itdetphot->first][0],PropParameters[itdetphot->first][0]+maxarrivaltimerange);
          }
          else
          {
            integral=pvs->GetParPropTimeFunctionIntegral(PropParameters[itdetphot->first],PropParameters[itdetphot->first][0]+maxarrivaltimerange);
          }

          float leftovers=0;
          for (int binnumber=0; binnumber < maxarrivaltimerange; binnumber++)
          {
           int ticks=static_cast<int> (deltaTime + binnumber); 
            if(pvs->GetfParPropTime_nonanalyticalfunction())
            {
              StepPhotons[ticks] = static_cast<int>(std::round(itdetphot->second*ParPropTimeFunctionTF1->Eval(binnumber+0.5+PropParameters[itdetphot->first][0])/integral));
              leftovers+=itdetphot->second*ParPropTimeFunctionTF1->Eval(binnumber+0.5+PropParameters[itdetphot->first][0])/integral-StepPhotons[ticks];
            }
            else
            {
              StepPhotons[ticks] = static_cast<int>(std::round(itdetphot->second*ParPropTimeFunction(PropParameters[itdetphot->first],binnumber+0.5+PropParameters[itdetphot->first][0])/integral));
              leftovers+=itdetphot->second*ParPropTimeFunction(PropParameters[itdetphot->first],binnumber+0.5+PropParameters[itdetphot->first][0])/integral-StepPhotons[ticks];
            }
            if (leftovers>1.0){leftovers--;StepPhotons[ticks]++;}
            else if (leftovers<-1.0){leftovers++;StepPhotons[ticks]--;}
            ccounter+=StepPhotons[ticks];

          }
          //std::cout << "CHECK COUNTER!!! " << ccounter << " =? " << itdetphot->second << " if they are not equal we have problems!"<<std::endl; 
          //Print(&StepPhotonTable);

        }
        else
        {
*/         
          for (G4int i = 0; i < itdetphot->second; ++i)
          {
            G4double deltaTime = ElectronCluster.Time();
            //standard propagation time, montecarlo simulation photon per photon
            if(pvs->IncludeParPropTime())
            {
              deltaTime += ParPropTimeTF1[itdetphot->first].GetRandom();
            }
            G4double aSecondaryTime = deltaTime;
            float Time = aSecondaryTime;
            int ticks = static_cast<int>(Time);
            StepPhotons[ticks]++;detectedcounter++;
          }
        //std::cout <<"itdetphot->first "  << itdetphot->first << std::endl;
        //std::cout <<"ElectronCluster.getTrackID() "  << ElectronCluster.getTrackID() <<",  ElectronCluster.getEnergy()/CLHEP::MeV "<< ElectronCluster.getEnergy()/CLHEP::MeV<< std::endl;
         StepPhotonTable[itdetphot->first] = StepPhotons;
         //Iterate over Step Photon Table to add photons to OpDetBacktrackerRecords.

         sim::OpDetBacktrackerRecord tmpOpDetBTRecord(itdetphot->first);
         //int thisG4TrackID = (aStep.GetTrack())->GetTrackID();
         int thisG4TrackID = ElectronCluster.TrackID();
         double xO = ( ElectronCluster.FinalPositionX() / CLHEP::cm );
         double yO = ( ElectronCluster.FinalPositionY() / CLHEP::cm );
         double zO = ( ElectronCluster.FinalPositionZ() / CLHEP::cm );
         double const xyzPos[3] = {xO,yO,zO};
        //         double energy  = ( aStep.GetTotalEnergyDeposit() / CLHEP::MeV );
         double energy  = ElectronCluster.Energy()/CLHEP::MeV;
         //Loop over StepPhotons to get number of photons detected at each time for this channel and G4Step.

         for(std::map<int,int>::iterator stepPhotonsIt = StepPhotons.begin(); stepPhotonsIt != StepPhotons.end(); ++stepPhotonsIt)
         {
             int photonTime = stepPhotonsIt->first;
             int numPhotons = stepPhotonsIt->second;
             tmpOpDetBTRecord.AddScintillationPhotons(thisG4TrackID, photonTime, numPhotons, xyzPos, energy);
         }
         litefst->AddOpDetBacktrackerRecord(tmpOpDetBTRecord);
      }//endloop per optdet
      //Print(&StepPhotonTable);
      litefst->AddPhoton(&StepPhotonTable);
    }//endif litephotons
    else
    {
      mf::LogError("PhotonLibraryPropagationS2") << "Error creating S2 light, SimPhotons data product is not supported."<<std::endl;
    }
  }//endif visibilities

  mf::LogInfo("PhotonLibraryPropagationS2") << counter <<" electron clusters processed. "<< ElectronCluster.NumberOfElectrons() <<" electron arrives. "<< nphot <<" S2 photons have been created in " << ElectronCluster.FinalPositionX() << " " << ElectronCluster.FinalPositionY() << " "<< ElectronCluster.FinalPositionZ()<< " at time " << ElectronCluster.Time() << std::endl;counter++;

//  OpDetSensitiveDetector *theOpDetDet = dynamic_cast<OpDetSensitiveDetector*>(sdManager->FindSensitiveDetector("OpDetSensitiveDetector"));
//  if (OpDetSensitiveDetector)
//  {

  }// endif electroncluster loop

  if(!fUseLitePhotons)
  {
     mf::LogError("PhotonLibraryPropagationS2") << "Error creating S2 light, SimPhotons data product is not supported."<<std::endl;
     //std::cout << "SIMPHOTONS NOT SUPPORTED!!! "<< std::endl;
      /*
    LOG_DEBUG("Optical") << "Storing OpDet Hit Collection in Event";
    std::vector<sim::SimPhotons>& ThePhotons = larg4::OpDetPhotonTable::Instance()->GetPhotons();
    PhotonCol->reserve(ThePhotons.size());
    for(auto& it : ThePhotons)
    PhotonCol->push_back(std::move(it));*/
  }
  else
  {
    mf::LogDebug("PhotonLibraryPropagationS2") << "Converting Photon Map in SimPhotonLite data product";
    std::map<int, std::map<int, int> > ThePhotons = litefst->GetLitePhotons();
    if(ThePhotons.size() > 0)
    {
      LitePhotonCol->reserve(ThePhotons.size());

      for(auto const& it : ThePhotons)
      {
              
        sim::SimPhotonsLite ph;
        ph.OpChannel = it.first;
        ph.DetectedPhotons = it.second;
        LitePhotonCol->push_back(ph);
      //std::cout << "... \t" << it.first << std::endl;
      }
  }
  //*cOpDetBacktrackerRecordCol = larg4::OpDetPhotonTable::Instance()->YieldOpDetBacktrackerRecords();
  *cOpDetBacktrackerRecordCol = litefst->YieldOpDetBacktrackerRecords();
  }

  if(!fUseLitePhotons)
  {
    mf::LogError("PhotonLibraryPropagationS2") << "Error creating S2 light, SimPhotons data product is not supported."<<std::endl;
    //e.put(std::move(PhotonCol));
  }
  else
  {
    mf::LogDebug("PhotonLibraryPropagationS2") << "Storing S2 Photon Collection in event " ;
    e.put(std::move(LitePhotonCol));
    e.put(std::move(cOpDetBacktrackerRecordCol));
  }
  
}

void phot::PhotonLibraryPropagationS2::Print(std::map<int, std::map<int, int>>* StepPhotonTable)
{
  for(auto it = StepPhotonTable->begin(); it!=StepPhotonTable->end(); it++)
  {
    for(auto in_it = it->second.begin(); in_it!=it->second.end(); in_it++)
    {
      std::cout << in_it->second << " ";
    }
    std::cout << std::endl;
  }
}


void phot::PhotonLibraryPropagationS2::reconfigure(fhicl::ParameterSet const & p)
{
  fGain = p.get<double>("Gain",500);
  fDriftEModuleLabel= p.get< std::string >("DriftEModuleLabel");
}

void phot::PhotonLibraryPropagationS2::beginJob()
{

}

DEFINE_ART_MODULE(phot::PhotonLibraryPropagationS2)
