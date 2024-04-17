////////////////////////////////////////////////////////////////////////
// Class:       MyAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        MyAnalysis_module.cc
//
// Generated at Tue Mar 21 10:26:24 2023 by Laura Pérez-Molina
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"

// Utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h" // simb::MCParticle
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h" //top level information form Pandora inLArSoft in recob::PFParticle
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h" // DetSim Analysis
#include "larcorealg/Geometry/GeometryCore.h" // DetSim Analysis
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"// DetSim Analysis
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "lardata/ArtDataHelper/HitCreator.h" // RawDigit
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"

// ROOT includes.
#include "TH1.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// C++ includes
#include <TTree.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <vector>
#include <cmath>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

// #include "duneana/DAQSimAna/TriggerPrimitiveFinderTool.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"


namespace ana { class MyAnalysis; }


class ana::MyAnalysis : public art::EDAnalyzer 
{
  public:
    explicit MyAnalysis(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    MyAnalysis(MyAnalysis const &) = delete;
    MyAnalysis(MyAnalysis &&) = delete;
    MyAnalysis& operator = (MyAnalysis const&) = delete;
    MyAnalysis& operator = (MyAnalysis&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override; //bulk of action
    // void reconfigure(fhicl::ParameterSet const & p);

    // Selected optional functions.
    void beginJob() override;
    // void endJob()   override;

    void reset(bool deepClean=false);

  private:

    // --- Fucntions --- //
    // long unsigned int WhichParType(int TrID);
    // void FillMyMaps(std::map<int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle<std::vector<simb::MCTruth>> Hand);
    // void CalcAdjHits(std::vector<recob::Hit> MyVec, std::vector<std::vector<recob::Hit>> &Clusters, TH1I *MyHist, TH1F *MyADCIntHist, bool HeavDebug);
    bool InMyMap(int TrID, std::map<int, simb::MCParticle> ParMap);
    void PrintInColor(std::string MyString, int MyColor, std::string Type = "Info");
    int GetColor(std::string MyString);
    std::string str(int MyInt);
    std::string str(unsigned int MyUnsignedInt);
    std::string str(float MyFloat);
    std::string str(double MyDouble);
    std::string str(std::vector<int> MyVec);
    std::string str(std::vector<float> MyVec);
    std::string str(std::vector<double> MyVec);
    int supress_stdout();
    void resume_stdout(int fd);

    // --- Declare member data here --- //

    // --- Our fcl parameter labels for the modules that made the data products
    std::vector<std::vector<std::string>> fProductsToDump;
    std::vector<std::vector<std::string>> fTreesToWrite;
    // --- Declare our services
    // const art::ServiceHandle<cheat::BackTrackerService> btServ;
    // art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    // art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // --- Input settings imported from the fcl
    const geo::Geometry* fGeom;
    bool fRollUpUnsavedIDs;
    bool debug = true;
    // --- Our TTrees, and its associated variables.
    TTree *fTruth;
    TTree *fDetSim;
    TTree *fReco;

    // Event Information //
    unsigned int fFlag, fEvent, fRun, fSubRun, fNMCParticles, fNPFParticles;
    float GenPDG;

    // TRUTH //
    std::vector<int>    fGenParticlePdgCode,fGenParticleTrackID,fMCParticlePdgCode,fMCParticleTrackID,fMCParticleParentTrackID,fSimEdepTrackID,fSimEdepPDGCode;
    std::vector<double> fGenParticleStartPositionX,fGenParticleStartPositionY,fGenParticleStartPositionZ,fGenParticleStartPositionT;
    std::vector<double> fGenParticleEndPositionX,fGenParticleEndPositionY,fGenParticleEndPositionZ,fGenParticleEndPositionT;
    std::vector<double> fGenParticleTrueEnergy,fGenParticleTrueMomentum,fGenParticleStartMomentumX,fGenParticleStartMomentumY,fGenParticleStartMomentumZ,fGenParticleStartMomentumE;
    std::vector<double> fGenParticleEndMomentumX,fGenParticleEndMomentumY,fGenParticleEndMomentumZ,fGenParticleEndMomentumE;
    std::vector<double> fMCParticleStartPositionX,fMCParticleStartPositionY,fMCParticleStartPositionZ,fMCParticleStartPositionT;
    std::vector<double> fMCParticleEndPositionX,fMCParticleEndPositionY,fMCParticleEndPositionZ,fMCParticleEndPositionT;
    std::vector<double> fMCParticleTrueEnergy,fMCParticleTrueMomentum;
    std::vector<double> fSimEdepE,fSimEdepX,fSimEdepY,fSimEdepZ;

    std::vector<bool>  fMCIsPrimary;

    // std::vector<std::map<int, simb::MCParticle>> Parts = {};

    // RECO //
    std::vector<int> fMCParticleNHits,fPFPID,fPFPPdgCode,fPFPNChildren,fPFPParentID,fPFPNClusters,fPFPTrackID;
    std::vector<double> fPFPTrackLength,fPFPTrackStartX,fPFPTrackStartY,fPFPTrackStartZ,fPFPTrackVertexX,fPFPTrackVertexY,fPFPTrackVertexZ;
    std::vector<double> fPFPTrackEndX,fPFPTrackEndY,fPFPTrackEndZ,fPFPTrackTheta,fPFPTrackPhi,fPFPTrackZenithAngle,fPFPTrackAzimuthAngle;
    std::vector<double> fPFPTrackStartDirectionX,fPFPTrackStartDirectionY,fPFPTrackStartDirectionZ,fPFPTrackVertexDirectionX,fPFPTrackVertexDirectionY,fPFPTrackVertexDirectionZ;
    std::vector<double> fPFPTrackEndDirectionX,fPFPTrackEndDirectionY,fPFPTrackEndDirectionZ,fPFPTrackChi2,fPFPTrackNdof;
    std::vector<int> fPFPTrueParticleMatchedID,fPFPTrueParticleMatchedPosition,fPFPNSharedTrueParticleHits;
    std::vector<bool> fPFPIsPrimary,fPFPIsTrack,fPFPIsShower;
    std::vector<int> fPFPNHits,fPFPShowerID,fPFPShowerBestPlane;
    std::vector<double> fPFPShowerDirectionX,fPFPShowerDirectionY,fPFPShowerDirectionZ,fPFPShowerDirectionErrX,fPFPShowerDirectionErrY,fPFPShowerDirectionErrZ;
    std::vector<double> fPFPShowerStartX,fPFPShowerStartY,fPFPShowerStartZ,fPFPShowerStartErrX,fPFPShowerStartErrY,fPFPShowerStartErrZ;
    std::vector<double> fPFPShowerLength,fPFPShowerOpenAngle,fPFPCompleteness,fPFPPurity;
    std::vector<std::vector<double>> fPFPShowerdEdx;
    
    std::vector<std::vector<double>> fPFPCompletenessView    {{}, {}, {}};
    std::vector<std::vector<double>> fPFPPurityView {{}, {}, {}};
    std::vector<std::vector<double>> fPFPTrueParticleMatchedPositionView  {{}, {}, {}};
    std::vector<std::vector<double>> fPFPTrueParticleMatchedIDView {{}, {}, {}};
    std::vector<std::vector<double>> fPFPNSharedTrueParticleHitsView {{}, {}, {}};
    std::vector<std::vector<double>> fPFPNHitsView   {{}, {}, {}};
    std::vector<std::vector<double>> fMCParticleNHitsView         {{}, {}, {}};
    std::vector<std::vector<double>> fPFPCluPlane    ;
    std::vector<std::vector<double>> fPFPCluView     ;
    std::vector<std::vector<double>> fPFPCluNHits    ;
    std::vector<std::vector<double>> fPFPCluIntegral ;

    static const int kNMaxMCParticles = 1000000;
    static const int kNMaxPFParticles = 2000;
    static const int kNMaxPFPClusters = 100;
    static const int kNViews          = 3;

    // // MC Particles //
    // unsigned int fNMCParticles;

    // int    fGenParticlePdgCode          [kNMaxMCParticles];
    // int    fGenParticleTrackID          [kNMaxMCParticles]; 
    // double fGenParticleStartPositionX   [kNMaxMCParticles];
    // double fGenParticleStartPositionY   [kNMaxMCParticles];
    // double fGenParticleStartPositionZ   [kNMaxMCParticles];
    // double fGenParticleStartPositionT   [kNMaxMCParticles];
    // double fGenParticleEndPositionX     [kNMaxMCParticles];
    // double fGenParticleEndPositionY     [kNMaxMCParticles];
    // double fGenParticleEndPositionZ     [kNMaxMCParticles];
    // double fGenParticleEndPositionT     [kNMaxMCParticles];
    // double fGenParticleTrueEnergy       [kNMaxMCParticles]; 
    // double fGenParticleTrueMomentum     [kNMaxMCParticles]; 
    // double fGenParticleStartMomentumX   [kNMaxMCParticles];
    // double fGenParticleStartMomentumY   [kNMaxMCParticles];
    // double fGenParticleStartMomentumZ   [kNMaxMCParticles];
    // double fGenParticleStartMomentumE   [kNMaxMCParticles];
    // double fGenParticleEndMomentumX     [kNMaxMCParticles];
    // double fGenParticleEndMomentumY     [kNMaxMCParticles];
    // double fGenParticleEndMomentumZ     [kNMaxMCParticles];
    // double fGenParticleEndMomentumE     [kNMaxMCParticles];

    // bool   fMCIsPrimary                [kNMaxMCParticles];
    // int    fMCParticlePdgCode          [kNMaxMCParticles];
    // int    fMCParticleTrackID          [kNMaxMCParticles]; 
    // int    fMCParticleParentTrackID    [kNMaxMCParticles]; 
    // double fMCParticleTrueEnergy       [kNMaxMCParticles]; 
    // double fMCParticleTrueMomentum     [kNMaxMCParticles]; 
    // int    fMCParticleNTrajectoryPoint [kNMaxMCParticles]; 
    // //Position
    // double fMCParticleStartPositionX   [kNMaxMCParticles];
    // double fMCParticleStartPositionY   [kNMaxMCParticles];
    // double fMCParticleStartPositionZ   [kNMaxMCParticles];
    // double fMCParticleStartPositionT   [kNMaxMCParticles];
    // double fMCParticleEndPositionX     [kNMaxMCParticles];
    // double fMCParticleEndPositionY     [kNMaxMCParticles];
    // double fMCParticleEndPositionZ     [kNMaxMCParticles];
    // double fMCParticleEndPositionT     [kNMaxMCParticles];
    // //Momentum
    // double fMCParticleStartMomentumX   [kNMaxMCParticles];
    // double fMCParticleStartMomentumY   [kNMaxMCParticles];
    // double fMCParticleStartMomentumZ   [kNMaxMCParticles];
    // double fMCParticleStartMomentumE   [kNMaxMCParticles];
    // double fMCParticleEndMomentumX     [kNMaxMCParticles];
    // double fMCParticleEndMomentumY     [kNMaxMCParticles];
    // double fMCParticleEndMomentumZ     [kNMaxMCParticles];
    // double fMCParticleEndMomentumE     [kNMaxMCParticles];

    // double fMCParticleVertexTime       [kNMaxMCParticles];
    // double fMCParticleEndTime          [kNMaxMCParticles];
    // int    fMCParticleNHits            [kNMaxMCParticles];
    // int    fMCParticleNHitsView        [kNMaxMCParticles][kNViews];


    // //SimEnergyDeposited //
    // int fSimEdepTrackID[kNMaxMCParticles];
    // int fSimEdepPDGCode[kNMaxMCParticles];
    // double fSimEdepE   [kNMaxMCParticles];
    // double fSimEdepX   [kNMaxMCParticles];
    // double fSimEdepY   [kNMaxMCParticles];
    // double fSimEdepZ   [kNMaxMCParticles];


    // "Particle Flow" Particle (PFP) information 
    // unsigned int fNPFParticles;
    // int    fPFPID                             [kNMaxPFParticles];
    // bool   fPFPIsPrimary                      [kNMaxPFParticles];
    // int    fPFPTrueParticleMatchedID          [kNMaxPFParticles];
    // int    fPFPTrueParticleMatchedPosition    [kNMaxPFParticles];
    // int    fPFPParentID                       [kNMaxPFParticles];
    // int    fPFPPdgCode                        [kNMaxPFParticles];
    // int    fPFPNChildren                      [kNMaxPFParticles];
    // int    fPFPNHits                          [kNMaxPFParticles];
    // int    fPFPNHitsView                      [kNMaxPFParticles][kNViews];
    // int    fPFPNSharedTrueParticleHits        [kNMaxPFParticles];
    // int    fPFPNSharedTrueParticleHitsView    [kNMaxPFParticles][kNViews];
    // int    fPFPTrueParticleMatchedIDView      [kNMaxPFParticles][kNViews];
    // int    fPFPTrueParticleMatchedPositionView[kNMaxPFParticles][kNViews];
    // double fPFPCompleteness                   [kNMaxMCParticles];
    // double fPFPCompletenessView               [kNMaxMCParticles][kNViews];
    // double fPFPPurity                         [kNMaxMCParticles];
    // double fPFPPurityView                     [kNMaxMCParticles][kNViews];
    // //////////
    // //Clusters
    // int    fPFPNClusters                      [kNMaxPFParticles];
    // int    fPFPCluPlane                       [kNMaxPFParticles][kNMaxPFPClusters];
    // int    fPFPCluView                        [kNMaxPFParticles][kNMaxPFPClusters];
    // int    fPFPCluNHits                       [kNMaxPFParticles][kNMaxPFPClusters];
    // double fPFPCluIntegral                    [kNMaxPFParticles][kNMaxPFPClusters];
    // //////////
    // //Tracks
    // bool   fPFPIsTrack                        [kNMaxPFParticles];
    // int    fPFPTrackID                        [kNMaxPFParticles];
    // double fPFPTrackLength                    [kNMaxPFParticles];
    // double fPFPTrackStartX                    [kNMaxPFParticles];
    // double fPFPTrackStartY                    [kNMaxPFParticles];
    // double fPFPTrackStartZ                    [kNMaxPFParticles];
    // double fPFPTrackVertexX                   [kNMaxPFParticles];
    // double fPFPTrackVertexY                   [kNMaxPFParticles];
    // double fPFPTrackVertexZ                   [kNMaxPFParticles];
    // double fPFPTrackEndX                      [kNMaxPFParticles];
    // double fPFPTrackEndY                      [kNMaxPFParticles];
    // double fPFPTrackEndZ                      [kNMaxPFParticles];
    // double fPFPTrackTheta                     [kNMaxPFParticles];
    // double fPFPTrackPhi                       [kNMaxPFParticles];
    // double fPFPTrackZenithAngle               [kNMaxPFParticles];
    // double fPFPTrackAzimuthAngle              [kNMaxPFParticles];
    // double fPFPTrackStartDirectionX           [kNMaxPFParticles];
    // double fPFPTrackStartDirectionY           [kNMaxPFParticles];
    // double fPFPTrackStartDirectionZ           [kNMaxPFParticles];
    // double fPFPTrackVertexDirectionX          [kNMaxPFParticles];
    // double fPFPTrackVertexDirectionY          [kNMaxPFParticles];
    // double fPFPTrackVertexDirectionZ          [kNMaxPFParticles];
    // double fPFPTrackEndDirectionX             [kNMaxPFParticles];
    // double fPFPTrackEndDirectionY             [kNMaxPFParticles];
    // double fPFPTrackEndDirectionZ             [kNMaxPFParticles];
    // float  fPFPTrackChi2                      [kNMaxPFParticles];
    // int    fPFPTrackNdof                      [kNMaxPFParticles];
    // //////////
    // //Showers
    // bool   fPFPIsShower                       [kNMaxPFParticles];
    // int    fPFPShowerID                       [kNMaxPFParticles];
    // int    fPFPShowerBestPlane                [kNMaxPFParticles];
    // double fPFPShowerDirectionX               [kNMaxPFParticles];
    // double fPFPShowerDirectionY               [kNMaxPFParticles];
    // double fPFPShowerDirectionZ               [kNMaxPFParticles];
    // double fPFPShowerDirectionErrX            [kNMaxPFParticles];
    // double fPFPShowerDirectionErrY            [kNMaxPFParticles];
    // double fPFPShowerDirectionErrZ            [kNMaxPFParticles];
    // double fPFPShowerStartX                   [kNMaxPFParticles];
    // double fPFPShowerStartY                   [kNMaxPFParticles];
    // double fPFPShowerStartZ                   [kNMaxPFParticles];
    // double fPFPShowerStartErrX                [kNMaxPFParticles];
    // double fPFPShowerStartErrY                [kNMaxPFParticles];
    // double fPFPShowerStartErrZ                [kNMaxPFParticles];
    // double fPFPShowerLength                   [kNMaxPFParticles];
    // double fPFPShowerOpenAngle                [kNMaxPFParticles];

    // std::vector< double > fPFPShowerdEdx      [kNMaxPFParticles];

    //Detector Hits Info
    int channel;
    int tdc;
    int adc;
    int view;
    int track_id;
    double energy;
    std::vector<int> track_ids;
    std::vector<double> energies;
}; //End class definition


ana::MyAnalysis::MyAnalysis(fhicl::ParameterSet const & p) : EDAnalyzer{p} //, 
  // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;

  std::vector<std::vector<std::string>> TreesToWrite   = p.get<std::vector<std::vector<std::string>>>("TreesToWrite");
  std::vector<std::vector<std::string>> ProductsToDump = p.get<std::vector<std::vector<std::string>>>("ProductsToDump");
  fTreesToWrite   = TreesToWrite;
  fProductsToDump = ProductsToDump;

  std::vector<std::string> trees; 
  int nTrees = fTreesToWrite.size();
  for (auto t = 0; t < nTrees; t++) { auto tree = fTreesToWrite[t][0]; trees.push_back(tree); }  


  if (std::find(trees.begin(), trees.end(), std::string("Truth")) != trees.end()) 
  {
    fTruth  = tfs->make<TTree>("Truth", "Truth");

    for (auto i : fProductsToDump)
    {
      auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
      TString out_label = o_label;

      // Set branches: Ifs over the different types of data
      // if (i_type == "simb::MCParticle")
      if (i_type == "simb::MCParticle")
      {
        //Event branches//
        fTruth->Branch("Event",      &fEvent,     "Event/I");
        fTruth->Branch("Run",        &fRun,       "Run/I");
        fTruth->Branch("Subrun",     &fSubRun,    "subRun/I");
        fTruth->Branch("Flag",       &fFlag,      "Flag/I");

        //MC truth branches//
        fTruth->Branch("nMCParticles",                &fNMCParticles,               "nMCParticles/I");
        // fTruth->Branch("genPdgCode",                  &fGenParticlePdgCode,         "genPdgCode[nMCParticles]/i");
        fTruth->Branch("genPdgCode",                  &fGenParticlePdgCode         );
        fTruth->Branch("genTrackID",                  &fGenParticleTrackID         );
        fTruth->Branch("genStartPositionX",           &fGenParticleStartPositionX  );
        fTruth->Branch("genStartPositionY",           &fGenParticleStartPositionY  );
        fTruth->Branch("genStartPositionZ",           &fGenParticleStartPositionZ  );
        fTruth->Branch("genStartPositionT",           &fGenParticleStartPositionT  );
        fTruth->Branch("genEndPositionX",             &fGenParticleEndPositionX    );
        fTruth->Branch("genEndPositionY",             &fGenParticleEndPositionY    );
        fTruth->Branch("genEndPositionZ",             &fGenParticleEndPositionZ    );
        fTruth->Branch("genEndPositionT",             &fGenParticleEndPositionT    );
        fTruth->Branch("genTrueEnergy",               &fGenParticleTrueEnergy      );
        fTruth->Branch("genTrueMomentum",             &fGenParticleTrueMomentum    );
        fTruth->Branch("genStartMomentumX",           &fGenParticleStartMomentumX  );
        fTruth->Branch("genStartMomentumY",           &fGenParticleStartMomentumY  );
        fTruth->Branch("genStartMomentumZ",           &fGenParticleStartMomentumZ  );
        fTruth->Branch("genEndMomentumX",             &fGenParticleEndMomentumX    );
        fTruth->Branch("genEndMomentumY",             &fGenParticleEndMomentumY    );
        fTruth->Branch("genEndMomentumZ",             &fGenParticleEndMomentumZ    );
        // fTruth->Branch("mcParticlePdgCode",           &fMCParticlePdgCode,          "MCParticlePdgCode[nMCParticles]/I");
        fTruth->Branch("mcParticlePdgCode",           &fMCParticlePdgCode          );
        fTruth->Branch("mcParticleTrackID",           &fMCParticleTrackID          );
        fTruth->Branch("mcParticleMotherID",          &fMCParticleParentTrackID    );
        fTruth->Branch("mcParticleTrueEnergy",        &fMCParticleTrueEnergy       );
        fTruth->Branch("mcParticleTrueMomentum",      &fMCParticleTrueMomentum     );
        fTruth->Branch("mcParticleStartPositionX",    &fMCParticleStartPositionX   );
        fTruth->Branch("mcParticleStartPositionY",    &fMCParticleStartPositionY   );
        fTruth->Branch("mcParticleStartPositionZ",    &fMCParticleStartPositionZ   );
        fTruth->Branch("mcParticleEndPositionX",      &fMCParticleEndPositionX     );
        fTruth->Branch("mcParticleEndPositionY",      &fMCParticleEndPositionY     );
        fTruth->Branch("mcParticleEndPositionZ",      &fMCParticleEndPositionZ     );

        // fTruth->Branch("mcIsMCPrimary",               &fMCIsPrimary,                "MCIsPrimary[nMCParticles]/O");
        // fTruth->Branch("mcParticleNTrajectoryPoints", &fMCParticleNTrajectoryPoint, "MCParticleNTrajectoryPoint[nMCParticles]/I");

        // fTruth->Branch("mcParticleStartPositionX",    &fMCParticleStartPositionX,    "MCParticleStartPositionX[nMCParticles]/D");
        // fTruth->Branch("mcParticleStartPositionY",    &fMCParticleStartPositionY,    "MCParticleStartPositionY[nMCParticles]/D");
        // fTruth->Branch("mcParticleStartPositionZ",    &fMCParticleStartPositionZ,    "MCParticleStartPositionZ[nMCParticles]/D");
        // fTruth->Branch("mcParticleStartPositionT",    &fMCParticleStartPositionT,    "MCParticleStartPositionT[nMCParticles]/D");
        // fTruth->Branch("mcParticleStartMomentumX",    &fMCParticleStartMomentumX,    "MCParticleStartMomentumX[nMCParticles]/D");
        // fTruth->Branch("mcParticleStartMomentumY",    &fMCParticleStartMomentumY,    "MCParticleStartMomentumY[nMCParticles]/D");
        // fTruth->Branch("mcParticleStartMomentumZ",    &fMCParticleStartMomentumZ,    "MCParticleStartMomentumZ[nMCParticles]/D");
        // fTruth->Branch("mcParticleStartMomentumE",    &fMCParticleStartMomentumE,    "MCParticleStartMomentumE[nMCParticles]/D");
        // fTruth->Branch("mcParticleEndPositionX",      &fMCParticleEndPositionX,      "MCParticleEndPositionX[nMCParticles]/D");
        // fTruth->Branch("mcParticleEndPositionY",      &fMCParticleEndPositionY,      "MCParticleEndPositionY[nMCParticles]/D");
        // fTruth->Branch("mcParticleEndPositionZ",      &fMCParticleEndPositionZ,      "MCParticleEndPositionZ[nMCParticles]/D");
        // fTruth->Branch("mcParticleEndPositionT",      &fMCParticleEndPositionT,      "MCParticleEndPositionT[nMCParticles]/D");
        // fTruth->Branch("mcParticleEndMomentumX",      &fMCParticleEndMomentumX,      "MCParticleEndMomentumX[nMCParticles]/D");
        // fTruth->Branch("mcParticleEndMomentumY",      &fMCParticleEndMomentumY,      "MCParticleEndMomentumY[nMCParticles]/D");
        // fTruth->Branch("mcParticleEndMomentumZ",      &fMCParticleEndMomentumZ,      "MCParticleEndMomentumZ[nMCParticles]/D");
        // fTruth->Branch("mcParticleEndMomentumE",      &fMCParticleEndMomentumE,      "MCParticleEndMomentumE[nMCParticles]/D");
        
        // fTruth->Branch("mcParticleVertexTime",        &fMCParticleVertexTime,        "MCParticleVertexTime[nMCParticles]/D");
        // fTruth->Branch("mcParticleEndTime",           &fMCParticleEndTime,           "MCParticleEndTime[nMCParticles]/D");
        // fTruth->Branch("mcParticleNHits",             &fMCParticleNHits,             "MCParticleNHits[nMCParticles]/I");
        // fTruth->Branch("mcParticleNHitsView",         &fMCParticleNHitsView,         "MCParticleNHitsView[nMCParticles][3]/I");

      }
      else if (i_type == "sim::SimEnergyDeposit")
      {
        // ENERGY DEPOSITED //
        fTruth->Branch("SimEdepTrackID",         &fSimEdepTrackID );
        fTruth->Branch("SimEdepPDGCode",         &fSimEdepPDGCode );
        fTruth->Branch("SimEdepEnergy",          &fSimEdepE );
        fTruth->Branch("SimEdepMiddlePositionX", &fSimEdepX );
        fTruth->Branch("SimEdepMiddlePositionY", &fSimEdepY );
        fTruth->Branch("SimEdepMiddlePositionZ", &fSimEdepZ );
      }
    } //for loop over products to dump
  } //if truth tree

  if (std::find(trees.begin(), trees.end(), std::string("Reco")) != trees.end()) 
  {
    fReco   = tfs->make<TTree>("Reco",  "Reco");

    //Event branches//
    fReco->Branch("Event",      &fEvent,     "Event/i");
    fReco->Branch("Run",        &fRun,       "Run/i");
    fReco->Branch("subRun",     &fSubRun,    "subRun/i");
    fReco->Branch("Flag",       &fFlag,      "Flag/I");

    //PFP branches
    fReco->Branch("nPFParticles",                       &fNPFParticles                  );                       
    fReco->Branch("pfpTrueParticleMatchedID",           &fPFPTrueParticleMatchedID      );
    fReco->Branch("pfpTrueParticleMatchedPosition",     &fPFPTrueParticleMatchedPosition);
    fReco->Branch("pfpIsPrimary",                       &fPFPIsPrimary                  );
    fReco->Branch("pfpID",                              &fPFPID                         );
    fReco->Branch("pfpParentID",                        &fPFPParentID                   );
    fReco->Branch("pfpPdgCode",                         &fPFPPdgCode                    );
    fReco->Branch("pfpNChildren",                       &fPFPNChildren                 );
    // fReco->Branch("pfpNHits",                           &fPFPNHits,                           "PFPNHits[nPFParticles]/I");
    // fReco->Branch("pfpNHitsView",                       &fPFPNHitsView,                       "PFPNHitsView[nPFParticles][3]/I");
    // fReco->Branch("pfpNSharedTrueParticleHits",         &fPFPNSharedTrueParticleHits,         "PFPNSharedTrueParticleHits[nPFParticles]/I");
    // fReco->Branch("pfpNSharedTrueParticleHitsView",     &fPFPNSharedTrueParticleHitsView,     "PFPNSharedTrueParticleHitsView[nPFParticles][3]/I");
    // fReco->Branch("pfpTrueParticleMatchedIDView",       &fPFPTrueParticleMatchedIDView,       "PFPTrueParticleMatchedIDView[nPFParticles][3]/I");
    // fReco->Branch("pfpTrueParticleMatchedPositionView", &fPFPTrueParticleMatchedPositionView, "PFPTrueParticleMatchedPositionView[nPFParticles][3]/I");
    
    fReco->Branch("pfpIsTrack",               &fPFPIsTrack               );
    fReco->Branch("pfpTrackID",               &fPFPTrackID               );
    fReco->Branch("pfpTrackLength",           &fPFPTrackLength           );
    fReco->Branch("pfpTrackStartX",           &fPFPTrackStartX           );
    fReco->Branch("pfpTrackStartY",           &fPFPTrackStartY           );
    fReco->Branch("pfpTrackStartZ",           &fPFPTrackStartZ           );
    fReco->Branch("pfpTrackVertexX",          &fPFPTrackVertexX          );
    fReco->Branch("pfpTrackVertexY",          &fPFPTrackVertexY          );
    fReco->Branch("pfpTrackVertexZ",          &fPFPTrackVertexZ          );
    fReco->Branch("pfpTrackEndX",             &fPFPTrackEndX             );
    fReco->Branch("pfpTrackEndY",             &fPFPTrackEndY             );
    fReco->Branch("pfpTrackEndZ",             &fPFPTrackEndZ             );
    fReco->Branch("pfpTrackTheta",            &fPFPTrackTheta            );
    fReco->Branch("pfpTrackPhi",              &fPFPTrackPhi              );
    fReco->Branch("pfpTrackZenithAngle",      &fPFPTrackZenithAngle      );
    fReco->Branch("pfpTrackAzimuthAngle",     &fPFPTrackAzimuthAngle     );
    fReco->Branch("pfpTrackStartDirectionX",  &fPFPTrackStartDirectionX  );
    fReco->Branch("pfpTrackStartDirectionY",  &fPFPTrackStartDirectionY  );
    fReco->Branch("pfpTrackStartDirectionZ",  &fPFPTrackStartDirectionZ  );
    fReco->Branch("pfpTrackVertexDirectionX", &fPFPTrackVertexDirectionX );
    fReco->Branch("pfpTrackVertexDirectionY", &fPFPTrackVertexDirectionY );
    fReco->Branch("pfpTrackVertexDirectionZ", &fPFPTrackVertexDirectionZ );
    fReco->Branch("pfpTrackEndDirectionX",    &fPFPTrackEndDirectionX    );
    fReco->Branch("pfpTrackEndDirectionY",    &fPFPTrackEndDirectionY    );
    fReco->Branch("pfpTrackEndDirectionZ",    &fPFPTrackEndDirectionZ    );
    fReco->Branch("pfpTrackChi2",             &fPFPTrackChi2             );
    fReco->Branch("pfpTrackStartNdof",        &fPFPTrackNdof             );

    fReco->Branch("pfpNClusters",             &fPFPNClusters             );
    fReco->Branch("pfpCluNHits",              &fPFPCluNHits              );
    // fReco->Branch("pfpCluPlane",    &fPFPCluPlane,    "PFPCluPlane[nPFParticles][100]/I");
    // fReco->Branch("pfpCluView",     &fPFPCluView,     "PFPCluView[nPFParticles][100]/I");
    // fReco->Branch("pfpCluIntegral", &fPFPCluIntegral, "PFPCluIntegral[nPFParticles][100]/D");

    fReco->Branch("pfpIsShower",            &fPFPIsShower            );
    fReco->Branch("pfpShowerID",            &fPFPShowerID            );
    fReco->Branch("pfpShowerBestPlane",     &fPFPShowerBestPlane     );
    fReco->Branch("pfpShowerDirectionX",    &fPFPShowerDirectionX    );
    fReco->Branch("pfpShowerDirectionY",    &fPFPShowerDirectionY    );
    fReco->Branch("pfpShowerDirectionZ",    &fPFPShowerDirectionZ    );
    // fReco->Branch("pfpShowerDirectionErrX", &fPFPShowerDirectionErrX, "PFPShowerDirectionErrX[nPFParticles]/D");
    // fReco->Branch("pfpShowerDirectionErrY", &fPFPShowerDirectionErrY, "PFPShowerDirectionErrY[nPFParticles]/D");
    // fReco->Branch("pfpShowerDirectionErrZ", &fPFPShowerDirectionErrZ, "PFPShowerDirectionErrZ[nPFParticles]/D");
    fReco->Branch("pfpShowerStartX",        &fPFPShowerStartX        );
    fReco->Branch("pfpShowerStartY",        &fPFPShowerStartY        );
    fReco->Branch("pfpShowerStartZ",        &fPFPShowerStartZ        );
    // fReco->Branch("pfpShowerStartErrX",     &fPFPShowerStartErrX,     "PFPShowerStartErrX[nPFParticles]/D");
    // fReco->Branch("pfpShowerStartErrY",     &fPFPShowerStartErrY,     "PFPShowerStartErrY[nPFParticles]/D");
    // fReco->Branch("pfpShowerStartErrZ",     &fPFPShowerStartErrZ,     "PFPShowerStartErrZ[nPFParticles]/D");
    fReco->Branch("pfpShowerLength",        &fPFPShowerLength        );
    fReco->Branch("pfpShowerOpeningAngle",  &fPFPShowerOpenAngle     );
    
    fReco->Branch("pfpShowerdEdx",          &fPFPShowerdEdx          );

    fReco->Branch("pfpCompleteness",      &fPFPCompleteness          );
    // fReco->Branch("pfpCompletenessView",  &fPFPCompletenessView, "PFPCompletenessView[nPFParticles][3]/D");
    fReco->Branch("pfpPurity",            &fPFPPurity                );
    // fReco->Branch("pfpPurityView",        &fPFPPurityView,       "PFPPurityView[nPFParticles][3]/D");
  } //if reco tree

  if (std::find(trees.begin(), trees.end(), std::string("DetSim")) != trees.end()) 
  {
    fDetSim = tfs->make<TTree>("DetSim","DetSim");
    fDetSim->Branch("Event",    &fEvent,     "Event/i");
    fDetSim->Branch("Channels", &channel);
    fDetSim->Branch("TDCs", &tdc);
    fDetSim->Branch("ADCs", &adc);
    fDetSim->Branch("View", &view);
    fDetSim->Branch("TrackID", &track_id);
    fDetSim->Branch("TrackIDs", &track_ids);
    fDetSim->Branch("Energy", &energy);
    fDetSim->Branch("Energies", &energies);
  } //if detsim tree

  fRollUpUnsavedIDs = p.get<bool>("RollUpUnsavedIDs"); 
  fGeom             = &*art::ServiceHandle<geo::Geometry>();
} 


void ana::MyAnalysis::beginJob()
{
  // Implementation of optional member function here.
  // reset(true); //deep clean the variables

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}//beginJob



void ana::MyAnalysis::analyze(const art::Event & evt)
{
  reset(true); //deep clean the variables
  fEvent  = evt.id().event();
  fRun    = evt.id().run();
  fSubRun = evt.id().subRun();
  std::map<int, int> fTrackIDMap; // store out of the loop for the truth_tree a map of trackID to MCParticle index

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt); //Contains all timing reference information for the detector. 
  fFlag = rand() % 10000000000;
  std::string sHead = "####################################################";
  sHead = sHead + "\nEVENT: " + str(fEvent) + " RUN: " + str(fRun) + " SUBRUN: " + str(fSubRun) + " FLAG: " + str(fFlag);
  sHead = sHead + "\nTPC Frequency in [MHz]: " + str(clockData.TPCClock().Frequency());
  sHead = sHead + "\nTPC Tick in [us]: " + str(clockData.TPCClock().TickPeriod());
  sHead = sHead + "\nSuccesfull reset of variables";
  sHead = sHead + "\n####################################################";
  PrintInColor(sHead, GetColor("magenta"));

  const sim::ParticleList &PartList = pi_serv->ParticleList();
  std::string sMcTruth = "";
  sMcTruth = sMcTruth + "\nThere are a total of " + str(int(PartList.size())) + " Particles in the event\n";
  PrintInColor(sMcTruth, GetColor("yellow"));
  std::vector<std::string> trees; 
  int nTrees = fTreesToWrite.size();
  for (auto t = 0; t < nTrees; t++) 
  {
    auto tree = fTreesToWrite[t][0];
    trees.push_back(tree);
  }  

  art::Handle<std::vector<simb::MCTruth>> mctruths;
  //=================== =============================== Tree Truth ==================================================//
  if (std::find(trees.begin(), trees.end(), std::string("Truth")) != trees.end()) 
  {
    std::string sMC_ini = "... Filling MC Truth Tree ... ";
    PrintInColor(sMC_ini, GetColor("magenta"));

    for (auto i : fProductsToDump)
    {
      std::string sInfoTruth = "--------------------------- TRUTH ----------------------------------";
      auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
    
      if (i_type == "simb::MCTruth") //"simb::MCParticle"
      {
        if (debug){sInfoTruth = sInfoTruth + "\n---> Module_Label: " + i_label + ";\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t Output_Name: " + o_label;}
        if(!evt.isRealData())
        {
          art::ValidHandle<std::vector<simb::MCTruth>> MCTrues = evt.getValidHandle<std::vector<simb::MCTruth>>(i_label);
          
          // --- Loop over the MCTruth vector (have your primary particles stored in the generator) --- //
          int i_gen(0);
          for (auto const &truth : *MCTrues)
          {
            int N = truth.NParticles(); // Number of descendant particles in this MCTruth (fPartList=n particles in the event)
            const simb::MCParticle gen_particle = truth.GetParticle(i_gen);
            GenPDG = gen_particle.PdgCode();
            float GenE = gen_particle.E();
            float GenX = gen_particle.Vx();
            float GenY = gen_particle.Vy();
            float GenZ = gen_particle.Vz();

            std::cout << "PDG: " << GenPDG << std::endl;
            // if (abs(GenPDG) == 13)
            // { 
            fGenParticlePdgCode       .push_back(gen_particle.PdgCode());
            fGenParticleTrackID       .push_back(gen_particle.TrackId());
            fGenParticleStartPositionX.push_back(gen_particle.Position().X());
            fGenParticleStartPositionY.push_back(gen_particle.Position().Y());
            fGenParticleStartPositionZ.push_back(gen_particle.Position().Z());
            fGenParticleStartPositionT.push_back(gen_particle.T());
            fGenParticleEndPositionX  .push_back(gen_particle.EndX());
            fGenParticleEndPositionY  .push_back(gen_particle.EndY());
            fGenParticleEndPositionZ  .push_back(gen_particle.EndZ());
            fGenParticleEndPositionT  .push_back(gen_particle.EndT());
            fGenParticleTrueEnergy    .push_back(gen_particle.E());
            fGenParticleTrueMomentum  .push_back(gen_particle.P()); //std::sqrt(std::pow(Momentum(i).E(),2.) - std::pow(fmass,2.)
            fGenParticleStartMomentumX.push_back(gen_particle.Px());
            fGenParticleStartMomentumY.push_back(gen_particle.Py());
            fGenParticleStartMomentumZ.push_back(gen_particle.Pz());
            fGenParticleEndMomentumX  .push_back(gen_particle.EndPx());
            fGenParticleEndMomentumY  .push_back(gen_particle.EndPy());
            fGenParticleEndMomentumZ  .push_back(gen_particle.EndPz());
            
            sInfoTruth = sInfoTruth + "\nPrimary particle (" + str(GenPDG) + "," + str(fGenParticleTrackID.back()) + ") energy: " + str(GenE) + " GeV";
            sInfoTruth = sInfoTruth + "\nVertex position (" + str(GenX) + ", " + str(GenY) + ", " + str(GenZ) + ") cm";
            fNMCParticles = N; // generalizar
            art::FindManyP<simb::MCParticle> mcParticles(MCTrues, evt, "largeant");
            auto trueParticles = mcParticles.at(i_gen); // Get the list of MCParticle daughters ("largeant" label) of i_gen primary particle
            int N_all = trueParticles.size();           // Number of descendant particles in this MCParticle
            sInfoTruth = sInfoTruth + "\nPrimary MCParticle " + str(i_gen) + " has " + str(N) + " daughters. (Total particles: " + str(N_all) +")";
            sInfoTruth = sInfoTruth + "\nGen.\tPdgCode\t\tEnergy\t\tEndPosition";
            sInfoTruth = sInfoTruth + "\n--------------------------------------------------------------------";
            // --- Loop over the MCParticle ALL descendants stored in trueParticle --- //
            for (int iDau = 1; iDau<N_all; iDau++) // all largeant objects (possible daughters)
            {
              auto mcParticlePtr = trueParticles.at(iDau);
              auto trueParticle  = *mcParticlePtr;
              if (trueParticle.Mother() == 0) //fGenParticleTrackID[i_gen]
              // if (trueParticle.Mother() == fGenParticleTrackID[i_gen])
              {
                // trueParticle.Process() == "primary"? isMCPrimary = true:isMCPrimary = false;
                fMCParticlePdgCode      .push_back(trueParticle.PdgCode());
                fMCParticleTrackID      .push_back(trueParticle.TrackId());
                fMCParticleParentTrackID.push_back(trueParticle.Mother());
                fMCParticleTrueEnergy   .push_back(trueParticle.E());
                fMCParticleTrueMomentum .push_back(trueParticle.P()); //std::sqrt(std::pow(Momentum(i).E(),2.) - std::pow(fmass,2.)
                fMCParticleStartPositionX .push_back(trueParticle.Position().X());
                fMCParticleStartPositionY .push_back(trueParticle.Position().Y());
                fMCParticleStartPositionZ .push_back(trueParticle.Position().Z());
                fMCParticleEndPositionX .push_back(trueParticle.EndX());
                fMCParticleEndPositionY .push_back(trueParticle.EndY());
                fMCParticleEndPositionZ .push_back(trueParticle.EndZ());
                // OLD WAY // int    fMCParticlePdgCode[kNMaxMCParticles] //
                // fMCParticlePdgCode      [iMc] = trueParticle.PdgCode();
                
                // if (trueParticle.PdgCode() < 1000000)
                // {
                //   sInfoTruth = sInfoTruth + "\n" + str(iMc) + " with PDG " + gen_particle.PdgCode() + "\t" + str(trueParticle.PdgCode()) + "\t\t" + str(trueParticle.E()) + "\t(" + str(trueParticle.EndX()) + ", " + str(trueParticle.EndY()) + ", " + str(trueParticle.EndZ()) + ") cm";
                // }
                // else
                // {
                //   sInfoTruth = sInfoTruth + "\n" + str(iMc) + " with PDG " + gen_particle.PdgCode() + "\t" + str(trueParticle.PdgCode()) + "\t" + str(trueParticle.E()) + "\t\t(" + str(trueParticle.EndX()) + ", " + str(trueParticle.EndY()) + ", " + str(trueParticle.EndZ()) + ") cm";
                // }
                // Position().X() = Vx()
                // std:: cout << "CHECK X (" <<  gen_particle.EndX() - trueParticle.Vx() << ")"<< " & (" <<  gen_particle.EndX() - trueParticle.Position().X() << ")" <<  std::endl;
                // std:: cout << "CHECK Y (" <<  gen_particle.EndY() - trueParticle.Vy() << ")"<< " & (" <<  gen_particle.EndY() - trueParticle.Position().Y() << ")" << std::endl;
                // std:: cout << "CHECK Z (" <<  gen_particle.EndZ() - trueParticle.Vz() << ")"<< " & (" <<  gen_particle.EndZ() - trueParticle.Position().Z() << ")" << std::endl;
              } //if (trueParticle.Mother() == 0)
            } //for loop filling MCInfo  
            std::cout << fMCParticlePdgCode.size() << std::endl;
            PrintInColor(sInfoTruth, GetColor("blue"));
            // } //if GenPDG == 13 or GenPDG != -13
            i_gen++;
          } //for loop over MCTruth vector
        } //if !evt.isRealData()
      } //if i_type == "simb::MCTruth"

      // Access Deposited Energy 
      if (i_type == "sim::SimEnergyDeposit")
      {
        if (debug)
        { 
          std::string sInfoEdep = "\n---> Module_Label: " + i_label + ";\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t Output_Name: " + o_label;
          PrintInColor(sInfoEdep, GetColor("magenta"));
        }

        art::Handle<std::vector<sim::SimEnergyDeposit>> EdepHandle;
        evt.getByLabel(i_label, i_instance, EdepHandle);
        if(!EdepHandle.isValid()) 
        { 
          std::string sWarningEdep = "\nUnable to find std::vector<sim::SimEnergyDeposit> with module label: " + i_label;
          PrintInColor(sWarningEdep, GetColor("red"));
        } // !EdepHandle
        if(EdepHandle.isValid())
        {
          for(auto edep : *EdepHandle) 
          {
            // int particle_id = fTrackIDMap[edep.TrackID()];
            fSimEdepTrackID.push_back(edep.TrackID());
            fSimEdepPDGCode.push_back(edep.PdgCode());
            fSimEdepX      .push_back(edep.MidPointX());
            fSimEdepY      .push_back(edep.MidPointY());
            fSimEdepZ      .push_back(edep.MidPointZ());
            fSimEdepE      .push_back(edep.E());
          } // for(auto edep : *EdepHandle)
        } // EdepHandle.isValid()  
      } //if i_type == "sim::SimEnergyDeposit"
    } //for loop over products to dump

    fTruth->Fill(); // FILL THE TREE

  } //if Truth is in the vector


  //================================================== Tree DetSim ==================================================//
  if (std::find(trees.begin(), trees.end(), std::string("DetSim")) != trees.end()) 
  {
    std::string sDS_ini = "... Filling DetSim Tree ... ";
    struct detector_output
          {
            int channel   = {  };
            int tdc       = {  };
            int adc       = {  };
            int view      = {  };
            int track_id  = {-1};
            double energy = { 0};

            std::vector<int> track_ids   = {};
            std::vector<double> energies = {};

            detector_output(int _channel, int _tdc, int _adc, int _view, std::vector<sim::IDE> ides)
            { 
              channel = _channel;
              tdc     = _tdc;
              adc     = _adc;
              view    = _view;
              for (auto ide : ides)
              {
                track_ids.emplace_back(ide.trackID);
                energies.emplace_back(ide.energy);
                if (ide.energy > energy) 
                {
                  energy   = ide.energy;
                  track_id = ide.trackID;
                }
              }
            }
          };

    std::vector<detector_output> output_array; // vector of detector_output objects

    art::Handle<std::vector<raw::RawDigit>> RawDigHandle;
    art::Handle<std::vector<sim::SimChannel>> SimChaHandle;
    for (auto i : fProductsToDump)
    {
      auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
    
      if (i_type == "raw::RawDigit")
      {
        if (debug){sDS_ini = sDS_ini + "\n---> Module_Label: " + i_label + ";\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t Output_Name: " + o_label;}

        evt.getByLabel(i_label, i_instance, RawDigHandle);

        if(!RawDigHandle.isValid())
        {
          std::string sWarningEdep = "\nUnable to find std::vector<raw::RawDigit> with module label: " + i_label + "; Instance: " + i_instance;
          PrintInColor(sWarningEdep, GetColor("red"));
          return;
        } // !RawDigHandle
      }

      if (i_type == "sim::SimChannel")
      {
        if (debug)
        {
          sDS_ini = sDS_ini + "\n---> Module_Label: " + i_label + ";\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t Output_Name: " + o_label;
        }
        evt.getByLabel(i_label, i_instance, SimChaHandle);

        if(!SimChaHandle.isValid())
        {
          std::string sWarningEdep = "\nUnable to find std::vector<sim::SimChannel> with module label: " + i_label + "; Instance: " + i_instance;
          PrintInColor(sWarningEdep, GetColor("red"));
          return;
        } // !DetSimHandle
      }
    } //for loop over products to dump

    PrintInColor(sDS_ini, GetColor("magenta"));

    std::map<std::tuple<int, int, int, int>, sim::SimChannel> channels; //map with keys as a tuple of three integers and values of type sim::SimChannel
    if(SimChaHandle.isValid() && RawDigHandle.isValid())
    {
      std::vector<sim::SimChannel> sim_channels; 
      std::vector<raw::RawDigit> digits;

      for(auto sim_channel : *SimChaHandle) { sim_channels.push_back(sim_channel); }
      for(auto digit : *RawDigHandle)       { digits.push_back(digit); }  // loop over all raw digits (i.e ADC counts for each channel for each time tick)
      // There are 15360 channels in the ProtoDUNE-SP detector.
      // RawDigit is the raw data from the detector. SimChannel is the simulated data from the detector.
      // RawDigit is the output of the electronics chain. SimChannel is the output of the simulation chain.
      // RawDigit is the ADC counts for each channel for each time tick. SimChannel is the number of electrons
      // that were deposited in each channel for each time tick.
      
      int DigsSize = digits.size(); // number of digits [DUNE-FD: 30720]
      std::cout << "Digit has " << DigsSize << " size" << std::endl;
      std::cout << "Digit has " << digits[0].Samples() << " samples" << std::endl;
      std::cout << "Digit has " << (int)digits[0].GetPedestal() << " pedestal" << std::endl;

      for (int dig=0; dig<DigsSize; dig++)
      {
        // Processing each channel and tick.
        int channel  = digits[dig].Channel();          // channel number
        int NSamples = digits[dig].Samples();          // number of ADC samples (TDC ticks) [DUNE-FD: 6000]
        int pedestal = (int)digits[dig].GetPedestal(); // pedestal value [DUNE-FD: 2350]

        // uncompress the digits and remove the pedestal
        std::vector<short> uncompressed(NSamples); // uncompressed ADC values
        raw::Uncompress( digits[dig].ADCs(), uncompressed, pedestal, digits[dig].Compression()); 

        for (int tick=0; tick<NSamples; tick++)
        {
          geo::GeometryCore const* geometry_core = lar::providerFrom<geo::Geometry>(); // geometry provider --> which wire plane
          if(std::abs(uncompressed[tick]) > 20.0) // ADC threshold > 20
          {
            auto sim_channel = sim_channels[channel]; // sim channel
            std::vector<sim::IDE> ides = sim_channel.TrackIDsAndEnergies(tick, tick);
            for (auto const& ide : ides) 
            {
              int particle_id = fTrackIDMap[ide.trackID];
              view = geometry_core->View(channel);
              std::tuple<int, int, int, int> key(channel, tick, particle_id, view);
              channels[key] = sim_channel;
              // Print the values being added to the channels map
              std::cout << "Added to channels map: Channel=" << channel << " Tick=" << tick << " ParticleID=" << particle_id << std::endl;
            } // for (auto const& ide : ides)
          } // if(std::abs(uncompressed[tick]) > 20.0)
        } // for (int tick=0; tick<NSamples; tick++)
      } // for (int dig=0; dig<DigsSize; dig++)

      if (channels.empty())  { std::cout << "WARNING channels empty"   << std::endl; }
      if (!channels.empty()) { std::cout << "GREAT channels not empty" << std::endl; }

      for (const auto& entry : channels) //iterates over each entry in the channels map
      {
        int channel, tick, particle_id, view;
        std::tie(channel, tick, particle_id, view) = entry.first;  //unpacks the tuple into the three variables
        const sim::SimChannel& sim_channel         = entry.second; //unpacks the sim channel into the variable

        output_array.emplace_back(channel, tick, particle_id, view, sim_channel.TrackIDsAndEnergies(tick, tick)); //creates a new entry in the output array
      }
    } // if(SimChaHandle.isValid() && RawDigHandle.isValid())
    
    for (const auto& output : output_array) 
    {
      // Access the components of each detector_output struct
      channel   = output.channel;
      tdc       = output.tdc;
      adc       = output.adc;
      view      = output.view;
      track_id  = output.track_id;
      energy    = output.energy;
      track_ids = output.track_ids;
      energies  = output.energies;
      
      fDetSim->Fill();
      
      // Print or use the components as needed
      // std::cout << "Channel: " << channel << std::endl;
      // std::cout << "TDC: " << tdc << std::endl;
      // std::cout << "ADC: " << adc << std::endl;
      // std::cout << "Track ID: " << track_id << std::endl;
      // std::cout << "Energy: " << energy << std::endl;
      // // Print the track IDs
      // std::cout << "Track IDs: ";
      // for (int track : track_ids) { std::cout << track << " "; }
      // std::cout << std::endl;
      // // Print the energies
      // std::cout << "Energies: ";
      // for (double e : energies) { std::cout << e << " "; }
      // std::cout << std::endl;
    } // for (const auto& output : output_array)

  } // if (tree == "DetSim")

  // ================================================== Tree Reco ==================================================//
  if (std::find(trees.begin(), trees.end(), std::string("Reco")) != trees.end()) 
  {
    // if (abs(GenPDG) == 13)
    // { 
    std::string sReco_ini = "... Filling Reco Tree ... ";
    PrintInColor(sReco_ini, GetColor("magenta"));
    std::string sInfoReco  = "------------------------ RECO ---------------------------------------";
    std::string label_pfp;      std::string instance_pfp;     std::string type_pfp; 
    std::string label_track;    std::string instance_track;   std::string type_track;
    std::string label_cluster;  std::string instance_cluster; std::string type_cluster;
    std::string label_shower;   std::string instance_shower;  std::string type_shower;

    for (auto i : fProductsToDump)
    {
      auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
      if (i_type == "recob::PFParticle")
      {
        if (debug){sInfoReco = sInfoReco + "\n---> Module_Label: " + i_label + ";\t\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + "; Output_Name: " + o_label;}
        label_pfp = i_label; instance_pfp = i_instance; type_pfp = i_type; 
      }

      if (i_type == "recob::Track")
      {
        if (debug){sInfoReco = sInfoReco + "\n---> Module_Label: " + i_label + ";\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t Output_Name: " + o_label;}
        label_track = i_label; instance_track = i_instance; type_track = i_type;
      }

      if (i_type == "recob::Cluster")
      {
        if (debug){sInfoReco = sInfoReco +  "\n---> Module_Label: " + i_label + ";\t\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t Output_Name: " + o_label;}
        label_cluster = i_label; instance_cluster = i_instance; type_cluster = i_type;
      }

      if (i_type == "recob::Shower")
      {
        if (debug){sInfoReco = sInfoReco + "\n---> Module_Label: " + i_label + ";\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t Output_Name: " + o_label;}
        label_shower = i_label; instance_shower = i_instance; type_shower = i_type;
      }
    } //for loop over products to dump

    //Access the PFParticles from Pandora
    const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, label_pfp);
    fNPFParticles = pfparticleVect.size();
    if(!fNPFParticles){std::string sWarningReco =  "No PFParticles found!"; PrintInColor(sWarningReco, GetColor("yellow"));}
    
    sInfoReco = sInfoReco + "\nNumber of PFParticles: " + str(fNPFParticles);
    sInfoReco = sInfoReco + "\n---------------------------------------------------------------------";
    PrintInColor(sInfoReco, GetColor("blue"));

    //Access the Clusters
    std::vector<art::Ptr<recob::Cluster>> clusterVect;
    auto clusterHandle = evt.getHandle<std::vector<recob::Cluster> >(label_cluster);
    if (clusterHandle){ art::fill_ptr_vector(clusterVect,clusterHandle); } 
    art::FindManyP<recob::Cluster> clusterParticleAssoc(pfparticleVect, evt, label_cluster);

    //Access the Tracks
    auto trackHandle = evt.getHandle<std::vector<recob::Track> >(label_track);
    if (!trackHandle){std::string sWarningReco = "\nUnable to find std::vector<recob::Track> with module label: " + label_track; PrintInColor(sWarningReco, GetColor("yellow"));}

    std::vector<art::Ptr<recob::Track> > trackList;
    art::fill_ptr_vector(trackList, trackHandle);

    int iPfp(0);
    for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect) // loop over all PFParticles
    {
      fPFPID       .push_back(iPfp);
      fPFPIsPrimary.push_back(pfp->IsPrimary());
      fPFPPdgCode  .push_back(pfp->PdgCode());
      fPFPNChildren.push_back(pfp->NumDaughters());
      if (pfp->IsPrimary()){fPFPParentID.push_back(-1)           ;}
      else                 {fPFPParentID.push_back(pfp->Parent());} // if primary, parent ID = -1

      std::vector<art::Ptr<recob::Cluster>> pfpClusters = clusterParticleAssoc.at(pfp.key()); // get clusters associated to the PFParticle
      fPFPNClusters.push_back(pfpClusters.size()); // number of clusters associated to the PFParticle
      
      if(!pfpClusters.empty())
      {
        int iClu(0);
        for(const art::Ptr<recob::Cluster> &clu:pfpClusters)
        {
          std::vector<double> fPFPCluPlane_in,fPFPCluView_in,fPFPCluNHits_in,fPFPCluIntegral_in;
          fPFPCluPlane_in   .push_back(clu->Plane().asPlaneID().Plane);
          fPFPCluView_in    .push_back(clu->View());
          fPFPCluNHits_in   .push_back(clu->NHits());
          fPFPCluIntegral_in.push_back(clu->Integral());
          
          fPFPCluPlane.push_back(fPFPCluPlane_in); fPFPCluView.push_back(fPFPCluView_in); fPFPCluNHits.push_back(fPFPCluNHits_in); fPFPCluIntegral.push_back(fPFPCluIntegral_in);

          iClu++;
          if (iClu == kNMaxPFPClusters){ break; }
        }

      }

      std::vector<art::Ptr<recob::Hit>> pfpHits; 
      if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt, label_pfp, label_track)) // if PFParticle is a track
      {
        art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, evt, label_pfp, label_track); 

        fPFPIsTrack              .push_back(true);
        fPFPTrackID              .push_back(track->ID());
        fPFPTrackLength          .push_back(track->Length());
        fPFPTrackStartX          .push_back(track->Start().X());
        fPFPTrackStartY          .push_back(track->Start().Y());
        fPFPTrackStartZ          .push_back(track->Start().Z());
        fPFPTrackVertexX         .push_back(track->Vertex().X());
        fPFPTrackVertexY         .push_back(track->Vertex().Y());
        fPFPTrackVertexZ         .push_back(track->Vertex().Z());
        fPFPTrackEndX            .push_back(track->End().X());
        fPFPTrackEndY            .push_back(track->End().Y());
        fPFPTrackEndZ            .push_back(track->End().Z());
        fPFPTrackTheta           .push_back(track->Theta());
        fPFPTrackPhi             .push_back(track->Phi());
        fPFPTrackZenithAngle     .push_back(track->ZenithAngle());
        fPFPTrackAzimuthAngle    .push_back(track->AzimuthAngle());
        fPFPTrackStartDirectionX .push_back(track->StartDirection().X());
        fPFPTrackStartDirectionY .push_back(track->StartDirection().Y());
        fPFPTrackStartDirectionZ .push_back(track->StartDirection().Z());
        fPFPTrackVertexDirectionX.push_back(track->VertexDirection().X());
        fPFPTrackVertexDirectionY.push_back(track->VertexDirection().Y());
        fPFPTrackVertexDirectionZ.push_back(track->VertexDirection().Z());
        fPFPTrackEndDirectionX   .push_back(track->EndDirection().X());
        fPFPTrackEndDirectionY   .push_back(track->EndDirection().Y());
        fPFPTrackEndDirectionZ   .push_back(track->EndDirection().Z());
        fPFPTrackChi2            .push_back(track->Chi2());
        fPFPTrackNdof            .push_back(track->Ndof());

        pfpHits = dune_ana::DUNEAnaTrackUtils::GetHits(track, evt, label_track); // get hits associated to the track
        
        std::vector<art::Ptr<recob::Hit>> pfpHitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, label_pfp, 0);
        std::vector<art::Ptr<recob::Hit>> pfpHitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, label_pfp, 1);
        std::vector<art::Ptr<recob::Hit>> pfpHitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, label_pfp, 2);
        
        fPFPNHits       .push_back(pfpHits.size());
        fPFPNHitsView[0].push_back(pfpHitsView0.size());
        fPFPNHitsView[1].push_back(pfpHitsView1.size());
        fPFPNHitsView[2].push_back(pfpHitsView2.size());

        if(!evt.isRealData()) // if MC
        {
          TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits, fRollUpUnsavedIDs));
          if (TruthMatchUtils::Valid(g4ID))
          {
            fPFPTrueParticleMatchedID.push_back(g4ID); int pos(999999); 
            for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
              if(fMCParticleTrackID.at(ipos)==g4ID) {pos = ipos;} // get the position of the MCParticle in the MCParticle tree
            }
            fPFPTrueParticleMatchedPosition.push_back(pos);
          }

          TruthMatchUtils::G4ID g4IDView0(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView0,fRollUpUnsavedIDs));
          if (TruthMatchUtils::Valid(g4ID))
          {
            fPFPTrueParticleMatchedIDView[0].push_back(g4IDView0);
            for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
              if (fMCParticleTrackID.at(ipos) != g4IDView0){ continue; }
              fPFPTrueParticleMatchedPositionView[0].push_back(ipos);
              break;
            }
          }

          TruthMatchUtils::G4ID g4IDView1(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHitsView1, fRollUpUnsavedIDs));
          if (TruthMatchUtils::Valid(g4ID))
          {
            fPFPTrueParticleMatchedIDView[1].push_back(g4IDView1);
            for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
              if (fMCParticleTrackID.at(ipos) != g4IDView1){ continue; }
              fPFPTrueParticleMatchedPositionView[1].push_back(ipos);
              break;
            }
          }

          TruthMatchUtils::G4ID g4IDView2(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView2,fRollUpUnsavedIDs));
          if (TruthMatchUtils::Valid(g4ID))
          {
            fPFPTrueParticleMatchedIDView[2].push_back(g4IDView2);
            for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
              if (fMCParticleTrackID.at(ipos) != g4IDView2){ continue; }
              fPFPTrueParticleMatchedPositionView[2].push_back(ipos);
              break;
            }
          }
        } // if !evt.isRealData()
      } // if IsTrack
      
      if(dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp, evt, label_pfp, label_shower)) // if PFParticle is a shower
      {
        art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp, evt, label_pfp, label_shower); 
        fPFPIsShower           .push_back(true);
        fPFPShowerID           .push_back(shower->ID());
        fPFPShowerBestPlane    .push_back(shower->best_plane());
        fPFPShowerDirectionX   .push_back(shower->Direction().X());
        fPFPShowerDirectionY   .push_back(shower->Direction().Y());
        fPFPShowerDirectionZ   .push_back(shower->Direction().Z());
        fPFPShowerDirectionErrX.push_back(shower->DirectionErr().X());
        fPFPShowerDirectionErrY.push_back(shower->DirectionErr().Y());
        fPFPShowerDirectionErrZ.push_back(shower->DirectionErr().Z());
        fPFPShowerStartX       .push_back(shower->ShowerStart().X());
        fPFPShowerStartY       .push_back(shower->ShowerStart().Y());
        fPFPShowerStartZ       .push_back(shower->ShowerStart().Z());
        fPFPShowerStartErrX    .push_back(shower->ShowerStartErr().X());
        fPFPShowerStartErrY    .push_back(shower->ShowerStartErr().Y());
        fPFPShowerStartErrZ    .push_back(shower->ShowerStartErr().Z());
        fPFPShowerLength       .push_back(shower->Length());
        fPFPShowerOpenAngle    .push_back(shower->OpenAngle());
        fPFPShowerdEdx         .push_back(shower->dEdx());
        

        pfpHits = dune_ana::DUNEAnaShowerUtils::GetHits(shower, evt, label_shower);
        std::vector<art::Ptr<recob::Hit> > pfpHitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, label_pfp, 0);
        std::vector<art::Ptr<recob::Hit> > pfpHitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, label_pfp, 1);
        std::vector<art::Ptr<recob::Hit> > pfpHitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, label_pfp, 2);
        fPFPNHits       .push_back(pfpHits.size());
        fPFPNHitsView[0].push_back(pfpHitsView0.size());
        fPFPNHitsView[1].push_back(pfpHitsView1.size());
        fPFPNHitsView[2].push_back(pfpHitsView2.size());


        if(!evt.isRealData()) // if MC
        {
          TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits, fRollUpUnsavedIDs));
          if (TruthMatchUtils::Valid(g4ID))
          {
            fPFPTrueParticleMatchedID.push_back(g4ID); int pos(999999); 
            for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            { 
              if(fMCParticleTrackID.at(ipos)==g4ID){ pos=ipos; }
            }
            fPFPTrueParticleMatchedPosition.push_back(pos);
          }

          TruthMatchUtils::G4ID g4IDView0(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHitsView0, fRollUpUnsavedIDs));
          if (TruthMatchUtils::Valid(g4ID))
          {
            fPFPTrueParticleMatchedIDView[0].push_back(g4IDView0);
            for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
              if (fMCParticleTrackID.at(ipos) != g4IDView0) { continue; }
              fPFPTrueParticleMatchedPositionView[0].push_back(ipos);
              break;
            }
          }

          TruthMatchUtils::G4ID g4IDView1(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHitsView1, fRollUpUnsavedIDs));
          if (TruthMatchUtils::Valid(g4ID))
          {
            fPFPTrueParticleMatchedIDView[1].push_back(g4IDView1);
            for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
              if (fMCParticleTrackID.at(ipos) != g4IDView1){ continue; }
              fPFPTrueParticleMatchedPositionView[1].push_back(ipos);
              break;
            }
          }

          TruthMatchUtils::G4ID g4IDView2(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHitsView2, fRollUpUnsavedIDs));
          if (TruthMatchUtils::Valid(g4ID))
          {
            fPFPTrueParticleMatchedIDView[2].push_back(g4IDView2);
            for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
            {
              if (fMCParticleTrackID.at(ipos) != g4IDView2){ continue; }
              fPFPTrueParticleMatchedPositionView[2].push_back(ipos);
              break;
            }
          }
        } // if(!evt.isRealData())
      } // if(dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp, evt, "pandora", "pandoraShower"))

      if(!evt.isRealData()) // if MC
      {
        //Fill shared MC particle to hits map for this specific PFP

        for (const auto& hit: pfpHits)
        {
          TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleID(clockData, hit, fRollUpUnsavedIDs));
          if (TruthMatchUtils::Valid(g4ID))
          {
            if(g4ID==fPFPTrueParticleMatchedID.back())                  
            { 
              if (!fPFPNSharedTrueParticleHits.empty()) {++fPFPNSharedTrueParticleHits.back();} 
            }
            if(g4ID==fPFPTrueParticleMatchedID.back() && hit->View()==0)
            { 
              if (!fPFPNSharedTrueParticleHitsView[0].empty()){++fPFPNSharedTrueParticleHitsView[0].back();}
            }
            if(g4ID==fPFPTrueParticleMatchedID.back() && hit->View()==1)
            { 
              if (!fPFPNSharedTrueParticleHitsView[1].empty()){++fPFPNSharedTrueParticleHitsView[1].back();}
            }
            if(g4ID==fPFPTrueParticleMatchedID.back() && hit->View()==2)
            { 
              if (!fPFPNSharedTrueParticleHitsView[2].empty()){++fPFPNSharedTrueParticleHitsView[2].back();}
            }
          }
        }

        // PURITY
        if (!fPFPNHits.empty() && fPFPNHits.back()<999999 && !fPFPNSharedTrueParticleHits.empty())
        {
          fPFPPurity.push_back(fPFPNSharedTrueParticleHits.back() / fPFPNHits.back());
        }
        if (!fPFPNHitsView[0].empty() && fPFPNHitsView[0].back()<999999 && !fPFPNSharedTrueParticleHitsView[0].empty())
        {
          fPFPPurityView[0].push_back(fPFPNSharedTrueParticleHitsView[0].back() / fPFPNHitsView[0].back());
        }
        if (!fPFPNHitsView[1].empty() && fPFPNHitsView[1].back()<999999 && !fPFPNSharedTrueParticleHitsView[1].empty())
        {
          fPFPPurityView[1].push_back(fPFPNSharedTrueParticleHitsView[1].back() / fPFPNHitsView[1].back());
        }
        if (!fPFPNHitsView[2].empty() && fPFPNHitsView[2].back()<999999 && !fPFPNSharedTrueParticleHitsView[2].empty())
        {
          fPFPPurityView[2].push_back(fPFPNSharedTrueParticleHitsView[2].back() / fPFPNHitsView[2].back());
        }

        // COMPLETENESS
        if (!fPFPTrueParticleMatchedPosition.empty())
        {
          std::vector<int>::size_type last_pos = fPFPTrueParticleMatchedPosition.back();
          if (last_pos != 999999)
          {
            if(last_pos < fMCParticleNHits       .size() && fMCParticleNHits       .at(last_pos) > 0 && fMCParticleNHits       .at(last_pos) < 999999) fPFPCompleteness       .push_back(fPFPNSharedTrueParticleHits       .back() / fMCParticleNHits       .at(last_pos));
            if(last_pos < fMCParticleNHitsView[0].size() && fMCParticleNHitsView[0].at(last_pos) > 0 && fMCParticleNHitsView[0].at(last_pos) < 999999) fPFPCompletenessView[0].push_back(fPFPNSharedTrueParticleHitsView[0].back() / fMCParticleNHitsView[0].at(last_pos));
            if(last_pos < fMCParticleNHitsView[1].size() && fMCParticleNHitsView[1].at(last_pos) > 0 && fMCParticleNHitsView[1].at(last_pos) < 999999) fPFPCompletenessView[1].push_back(fPFPNSharedTrueParticleHitsView[1].back() / fMCParticleNHitsView[1].at(last_pos));
            if(last_pos < fMCParticleNHitsView[2].size() && fMCParticleNHitsView[2].at(last_pos) > 0 && fMCParticleNHitsView[2].at(last_pos) < 999999) fPFPCompletenessView[2].push_back(fPFPNSharedTrueParticleHitsView[2].back() / fMCParticleNHitsView[2].at(last_pos));
          } // if(last_pos != 999999)
        }  // if (!fPFPTrueParticleMatchedPosition.empty())
      } // if(!evt.isRealData())
      iPfp++;
    } // for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect)
    // } // if (GenPDG == 13 or GenPDG != -13)
    // else {std::cout << "WARNING: GenPDG is not muon! (not computing Reco) \n" << std::endl;}
    fReco->Fill();

  } // if (tree == "Reco")


}// analyze()

// ########################################################################################################################################//
//_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_FUNCTION_TIME!_//
// ########################################################################################################################################//

void ana::MyAnalysis::reset(bool deepClean)
{
    fNMCParticles = 0;
    fGenParticlePdgCode       .clear();
    fGenParticleTrackID       .clear();
    fGenParticleStartPositionX.clear();
    fGenParticleStartPositionY.clear();
    fGenParticleStartPositionZ.clear();
    fGenParticleStartPositionT.clear();
    fGenParticleEndPositionX  .clear();
    fGenParticleEndPositionY  .clear();
    fGenParticleEndPositionZ  .clear();
    fGenParticleEndPositionT  .clear();
    fGenParticleTrueEnergy    .clear();
    fGenParticleTrueMomentum  .clear();
    fGenParticleStartMomentumX.clear();
    fGenParticleStartMomentumY.clear();
    fGenParticleStartMomentumZ.clear();
    fGenParticleEndMomentumX  .clear();
    fGenParticleEndMomentumY  .clear();
    fGenParticleEndMomentumZ  .clear();
    fMCParticlePdgCode      .clear();
    fMCParticleTrackID      .clear();
    fMCParticleParentTrackID.clear();
    fMCParticleTrueEnergy   .clear();
    fMCParticleTrueMomentum .clear();
    fMCParticleEndPositionX .clear();
    fMCParticleEndPositionY .clear();
    fMCParticleEndPositionZ .clear();


    // for(unsigned int iPfp=0; iPfp<(deepClean ? kNMaxPFParticles : fNPFParticles); iPfp++){

    //  fPFPID[iPfp]=999999;
    //  fPFPTrueParticleMatchedID[iPfp]=999999;
    //  fPFPTrueParticleMatchedPosition[iPfp]=999999;

    //  fPFPNHits[iPfp]=999999;
    //  fPFPNSharedTrueParticleHits[iPfp]=0;

    //  fPFPNClusters[iPfp]=999999;
    //  fPFPIsTrack[iPfp]=0;
    //  fPFPIsShower[iPfp]=0;

    //  fPFPTrackID[iPfp]=999999;
    //  fPFPTrackLength[iPfp]=999999;
    //  fPFPTrackStartX[iPfp]=999999;
    //  fPFPTrackStartY[iPfp]=999999;
    //  fPFPTrackStartZ[iPfp]=999999;
    //  fPFPTrackVertexX[iPfp]=999999;
    //  fPFPTrackVertexY[iPfp]=999999;
    //  fPFPTrackVertexZ[iPfp]=999999;
    //  fPFPTrackEndX[iPfp]=999999;
    //  fPFPTrackEndY[iPfp]=999999;
    //  fPFPTrackEndZ[iPfp]=999999;
    //  fPFPTrackTheta[iPfp]=999999;
    //  fPFPTrackPhi[iPfp]=999999;
    //  fPFPTrackZenithAngle[iPfp]=999999;
    //  fPFPTrackAzimuthAngle[iPfp]=999999;
    //  fPFPTrackStartDirectionX[iPfp]=999999;
    //  fPFPTrackStartDirectionY[iPfp]=999999;
    //  fPFPTrackStartDirectionZ[iPfp]=999999;
    //  fPFPTrackVertexDirectionX[iPfp]=999999;
    //  fPFPTrackVertexDirectionY[iPfp]=999999;
    //  fPFPTrackVertexDirectionZ[iPfp]=999999;
    //  fPFPTrackEndDirectionX[iPfp]=999999;
    //  fPFPTrackEndDirectionY[iPfp]=999999;
    //  fPFPTrackEndDirectionZ[iPfp]=999999;
    //  fPFPTrackChi2[iPfp]=999999;
    //  fPFPTrackNdof[iPfp]=999999;
     
    //  fPFPShowerID[iPfp]=999999;
    //  fPFPShowerBestPlane[iPfp]=999999;
    //  fPFPShowerDirectionX[iPfp]=999999;
    //  fPFPShowerDirectionY[iPfp]=999999;
    //  fPFPShowerDirectionZ[iPfp]=999999;
    //  fPFPShowerDirectionErrX[iPfp]=999999;
    //  fPFPShowerDirectionErrY[iPfp]=999999;
    //  fPFPShowerDirectionErrZ[iPfp]=999999;
    //  fPFPShowerStartX[iPfp]=999999;
    //  fPFPShowerStartY[iPfp]=999999;
    //  fPFPShowerStartZ[iPfp]=999999;
    //  fPFPShowerStartErrX[iPfp]=999999;
    //  fPFPShowerStartErrY[iPfp]=999999;
    //  fPFPShowerStartErrZ[iPfp]=999999;
    //  fPFPShowerLength[iPfp]=999999;
    //  fPFPShowerOpenAngle[iPfp]=999999;

    //  for(int iClu=0; iClu<kNMaxPFPClusters; iClu++){
    //    fPFPCluPlane[iPfp][iClu]=999999; 
    //    fPFPCluView[iPfp][iClu]=999999; 
    //    fPFPCluNHits[iPfp][iClu]=999999; 
    //    fPFPCluIntegral[iPfp][iClu]=999999; 
    //  }

    //  fPFPCompleteness[iPfp]=999999;
    //  fPFPPurity[iPfp]=999999;

    //  for (unsigned int iView = 0; iView < kNViews; iView++)
    //  {
    //       fPFPNHitsView[iPfp][iView] = 999999;
    //       fPFPNSharedTrueParticleHitsView[iPfp][iView]=0;
    //       fPFPTrueParticleMatchedIDView[iPfp][iView] = 999999;
    //       fPFPTrueParticleMatchedPositionView[iPfp][iView] = 999999;
    //       fPFPPurityView[iPfp][iView]=999999;
    //       fPFPCompletenessView[iPfp][iView]=999999;
    //  }
    // }
    // fNPFParticles = 0;

    return;
}

//......................................................
// This function checks if a given TrackID is in a given map
bool ana::MyAnalysis::InMyMap(int TrID, std::map<int, simb::MCParticle> ParMap)
{
  std::map<int, simb::MCParticle>::iterator ParIt;
  ParIt = ParMap.find(TrID);
  if (ParIt != ParMap.end())
  {
    return true;
  }
  else
    return false;
}

//......................................................
// This function creates a terminal color printout
void ana::MyAnalysis::PrintInColor(std::string MyString, int Color, std::string Type)
{
  if (Type == "Info")
  {
    mf::LogInfo("ana") << "\033[" << Color << "m" << MyString << "\033[0m";
  }
  if (Type == "Degub")
  {
    mf::LogDebug("ana") << "\033[" << Color << "m" << MyString << "\033[0m";
  }
  if (Type == "Error")
  {
    mf::LogError("ana") << "\033[" << Color << "m" << MyString << "\033[0m";
  }
  return;
}

// ......................................................
// This function returns an integer that corresponds to a given color name
int ana::MyAnalysis::GetColor(std::string ColorName)
{
  if (ColorName == "black")
    return 30;
  else if (ColorName == "red")
    return 31;
  else if (ColorName == "green")
    return 32;
  else if (ColorName == "yellow")
    return 33;
  else if (ColorName == "blue")
    return 34;
  else if (ColorName == "magenta")
    return 35;
  else if (ColorName == "cyan")
    return 36;
  else if (ColorName == "white")
    return 37;
  else if (ColorName == "bright_black")
    return 90;
  else if (ColorName == "bright_red")
    return 91;
  else if (ColorName == "bright_green")
    return 92;
  else if (ColorName == "bright_yellow")
    return 93;
  else if (ColorName == "bright_blue")
    return 94;
  else if (ColorName == "bright_magenta")
    return 95;
  else if (ColorName == "bright_cyan")
    return 96;
  else if (ColorName == "bright_white")
    return 97;
  else
  {
    mf::LogError("ana::MyAnalysis") << "Color " << ColorName << " not recognized. Returning white.";
    return 37;
  }
  return 0;
}

std::string ana::MyAnalysis::str(int i)
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}
std::string ana::MyAnalysis::str(unsigned int i)
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}
std::string ana::MyAnalysis::str(double i)
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}
std::string ana::MyAnalysis::str(float i)
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}
std::string ana::MyAnalysis::str(std::vector<int> i)
{
  std::stringstream ss;
  for (int j = 0; j < int(i.size()); j++)
  {
    ss << i[j] << " ";
  }
  return ss.str();
}
std::string ana::MyAnalysis::str(std::vector<double> i)
{
  std::stringstream ss;
  for (int j = 0; j < int(i.size()); j++)
  {
    ss << i[j] << " ";
  }
  return ss.str();
}
std::string ana::MyAnalysis::str(std::vector<float> i)
{
  std::stringstream ss;
  for (int j = 0; j < int(i.size()); j++)
  {
    ss << i[j] << " ";
  }
  return ss.str();
}

void ana::MyAnalysis::resume_stdout(int fd)
{
  std::fflush(stdout);
  dup2(fd, 1);
  close(fd);
}


// void ana::MyAnalysis::endJob()
// {
//   // Implementation of optional member function here.
// }

DEFINE_ART_MODULE(ana::MyAnalysis)
