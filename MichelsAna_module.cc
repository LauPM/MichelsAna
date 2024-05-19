////////////////////////////////////////////////////////////////////////////////
// Class:MichelsAna
// Plugin Type: analyzer (art v3_05_01)
// File:MichelsAna_module.cc
//
// Use: lar -c run_MichelsAna.fcl -s /pnfs/dune/.../file.root (-n1)
// Output: root file with the tree as configured in MichelsAna.fcl
// 
// Written by Laura Pérez-Molina (Tue Mar 21 10:26:24 2023) with guidence 
// of Sergio Manthey Corchado and from Aleena Rafique MichelsAnalysis for PD-SP
////////////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/RootIOPolicy.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCTruth.h"    // simb::MCTruth
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h" // simb::MCParticle
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RecoBase/PFParticle.h" //top level information form Pandora inLArSoft in recob::PFParticle
#include "lardataobj/AnalysisBase/Calorimetry.h"
// #include "lardataobj/RecoBase/SpacePoint.h"
// #include "lardataobj/RecoBase/PointCharge.h"
// #include "lardataobj/RecoBase/Wire.h"
// #include "lardataobj/RecoBase/EndPoint2D.h"
// #include "lardataobj/RecoBase/Vertex.h"
// #include "lardataobj/RawData/BeamInfo.h"
// #include "lardataobj/AnalysisBase/CosmicTag.h"
// #include "lardataobj/AnalysisBase/FlashMatch.h"
// #include "lardataobj/AnalysisBase/T0.h"
// #include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
// Includes with particular tools and functions for this module
#include "MichelsTools.h"

#define MVA_LENGTH 4

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TTimeStamp.h"
#include "TProfile.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TMinuit.h"

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
#include <cstdlib>


namespace michels 
{ 

	class MichelsAna : public art::EDAnalyzer 
	{
		public:
			// --- Standard constructor and destructor for an ART module --- //
			explicit MichelsAna(fhicl::ParameterSet const& p);
			virtual ~MichelsAna();
			MichelsAna(MichelsAna const &) = delete;
			MichelsAna(MichelsAna &&) = delete;
			MichelsAna& operator = (MichelsAna const&) = delete;
			MichelsAna& operator = (MichelsAna&&) = delete;

			// Required functions.
			void beginJob() override;
			void endJob() override;
			void analyze(art::Event const& e) override; //bulk of action

			double activeBounds_eff[6] = {-359.4155, 359.4155, 0., 607.49875, -0.49375, 695.28625};
			double fThicknessIniVolume = 30;
			bool IsPointInBounds(double *v, TVector3 const & p)  
			{
				return ((p.Y()>=(v[3]-fThicknessIniVolume) && p.Y()<=v[3]) || (p.X()>=v[0] && p.X()<=(v[0]+fThicknessIniVolume)) || (p.X()<=v[1] && p.X()>=(v[1]-fThicknessIniVolume)) || (p.Z()>=v[4] && p.Z()<=(v[4]+fThicknessIniVolume)) || (p.Z()<=v[5] && p.Z()>=(v[5]-fThicknessIniVolume)));
			}
			int fiducialBounds[6]={int(activeBounds_eff[0])+50,int(activeBounds_eff[1])-50,int(activeBounds_eff[2])+80,int(activeBounds_eff[3])-50,int(activeBounds_eff[4])+80,int(activeBounds_eff[5])-85};
			bool IsPointInFV(int *v, TVector3 const & t)
			{
				return (t.X() >= v[0] && t.X() <= v[1] && t.Y() >= v[2] && t.Y() <= v[3] && t.Z() >= v[4] && t.Z() <= v[5]);
			}
			int fiducialBounds1[6]={int(activeBounds_eff[0])+50,int(activeBounds_eff[1])-50,int(activeBounds_eff[2])+50,int(activeBounds_eff[3])-50,int(activeBounds_eff[4])+50,int(activeBounds_eff[5])-50};
			
		private:
			// --- Our functions stored in MichelsBasicUtils --- //
			std::unique_ptr<michels::MichelsBasicUtils> utils;

			// --- Our fcl parameter labels for the modules that made the data products
			std::vector<std::vector<std::string>> fProductsToDump;
			std::vector<std::vector<std::string>> fTreesToWrite;
      double fRadiusThreshold; // Radius for hits in event selection
      double fCNNThreshold;    // CNN output threshold for selection
      int fCloseHitsThreshold; // Selection threshold for number of close hits 
      bool debug;              // Debugging flag
			// --- Declare our services
			art::ServiceHandle<cheat::BackTrackerService> bt_serv;
			art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
			art::ServiceHandle<cheat::PhotonBackTrackerService> pbt_serv;
			// --- Input settings imported from the fcl
			const geo::Geometry* fGeom;
			bool fRollUpUnsavedIDs;
			// --- Our TTrees, and its associated variables.
			TTree *Gen, *Cry, *AllReco, *CutReco;

      //------------ VARIABLES DEFINED ------------//


      unsigned int fFlag, fEvent, fRun, fSubRun, genNDaughters, genOrigin, genPdgCode, genTrackID, cryNDaughters, cryOrigin, cryPdgCode, cryTrackID, fNPFParticles;
      unsigned int fcutFlag, fcutEvent, fcutRun, fcutSubRun, fcut_endhitkey,fcut_endwire,fcut_endchno,fcut_endtpcno;
      int    fNFlashes, fall_trksint, fReadOutWindowSize, fNumberTimeSamples;
      int    fPFP_trks, fT0_trks, fstartinbound_trks, fendinFV_trks, fccrosser_trks;
      int    fnoelecdiv_bound, funbroken_trks, flongcm_trks, fminhitpt_trks, fmaxhitpt_trks;
      int    fnearhits_trks, fcut_mu, funcont_trks, fcrossZ_trks, fbackward_trks, fEndinTPC_trks, NTrueMichels;
      int    ft0ana, fstartinboundana, fendinFVana, fccrossana, fstopzana, fdistana, fbrokcoana, ftrklenana, fminhitptana, fmaxhitptana;
      int    ftrkdistcolana, fhitctana, fshwrdisana, NTracks, ftrkcolhits, fshwrcolhits, fshwrallhits;
      int    fntrkhits,fhitsV,fhitsU, ftrueparhitallcount;
      int    fhitsY, fnearhitct, fmhitcount, ftrueparhitcolcount;
      int    fyear_mon_day, fhour_min_sec;
      int    fcut_endhitmult, fcut_ccrosse, fcut_trkI, fcut_ncolhit, fcut_nearhitcoun, fMichelcountcutan, fcutshwr_ke;
      int    fcutshwr_ID, fcutshwr_bestplane, fmccut_trkid, fmccut_endprocess, fmccut_pdg, fmccut_status_code, fmccut_ND;
      int    fmccut_mother, fmccut_origin, fmccut_process, fmccut_rescatter, fmchits_trkid, fmchits_endprocess, fmchits_pdg; 
      int    fmchits_status_code, fmchits_mother, fmchits_origin, fmchits_process, fmchits_rescatter, fmchits_ND;
      int    fmcconehits_trkid, fmcconehits_endprocess, fmcconehits_pdg, fmcconehits_status_code, fmcconehits_ND;
      int    fmcconehits_mother, fmcconehits_origin, fmcconehits_process, fmcconehits_rescatter, fmcd_trkid, fmcd_endprocess;
      int    fmcd_pdg, fmcd_status_code, fmcd_ND, fmcd_mother, fmcd_origin, fmcd_process, fmcd_rescatter, fmcshwr_trkid; 
      int    fmcshwr_endprocess, fmcshwr_pdg, fmcshwr_status_code, fmcshwr_ND, fmcshwr_mother, fmcshwr_origin, fmcshwr_process;
      int    fmcshwr_rescatter, fmichel_conesize, fmichel_wcount, totflash, ftotintwf;
      float  fcut_endptime, fcut_endhitchrg, fcut_endhitx, fcut_endhity, fcut_endhitz, fcut_endhitsigptime, fcut_endhitsigchrg;
      float  fcut_endhitsigpamp, fcut_endhitdof, fcut_endhitgof, fcut_endhitptplusRMS, fcut_endhitptminusRMS, fcut_dist_hit_end, fcut_dist_times;
      float  fcut_trkstartx, fcut_trkstarty, fcut_trkstartz, fcut_trkendx, fcut_trkendy, fcut_trkendz, fcut_trkstartcosx, fcut_trkstartcosy;
      float  fcut_trkstartcosz, fcut_trkendcosx, fcut_trkendcosy, fcut_trkendcosz, fcut_trackthetaxz, fcut_trackthetayz, fcut_trklen, fcut_trktheta;
      float  fcut_trkphi, fcut_minhitptime, fcut_maxhitptime, fcut_CorrWirePtime, fcut_trkrecotime, fcutshwr_length, fcutshwr_startx;
      float  fcutshwr_starty, fcutshwr_startz, fcutshwr_startdcosx, fcutshwr_startdcosy, fcutshwr_startdcosz, fcutshwr_openangle, fcutshwr_dist;
      float  fmccut_vx, fmccut_vy, fmccut_vz, fmccut_t, fmccut_endx, fmccut_endy, fmccut_endz;
      float  fmccut_endt, fmccut_px, fmccut_py, fmccut_pz, fmccut_momentum, fmccut_energy, fmccut_endpx;
      float  fmccut_endpy, fmccut_endpz, fmccut_endenergy, fmccut_length, fmccut_pathlen, fmccut_vxdrifted, fmccut_vydrifted;
      float  fmccut_vzdrifted, fmccut_tdrifted, fmccut_endxdrifted, fmccut_endydrifted, fmccut_endzdrifted, fmccut_endtdrifted, fmccut_pxdrifted;
      float  fmccut_pydrifted, fmccut_pzdrifted, fmccut_momentumdrifted, fmccut_energydrifted, fmccut_endpxdrifted, fmccut_endpydrifted, fmccut_endpzdrifted;
      float  fmccut_endenergydrifted, fmccut_pathlendrifted, fmccut_theta, fmccut_phi, fmccut_mass, fmchits_vx, fmchits_vy;
      float  fmchits_vz, fmchits_t, fmchits_endx, fmchits_endy, fmchits_endz, fmchits_endt, fmchits_px;
      float  fmchits_py, fmchits_pz, fmchits_momentum, fmchits_energy, fmchits_endpx, fmchits_endpy, fmchits_endpz;
      float  fmchits_endenergy, fmchits_pathlen, fmchits_vxdrifted, fmchits_vydrifted, fmchits_vzdrifted, fmchits_tdrifted, fmchits_endxdrifted;
      float  fmchits_endydrifted, fmchits_endzdrifted, fmchits_endtdrifted, fmchits_pxdrifted, fmchits_pydrifted, fmchits_pzdrifted, fmchits_momentumdrifted;
      float  fmchits_energydrifted, fmchits_endpxdrifted, fmchits_endpydrifted, fmchits_endpzdrifted, fmchits_endenergydrifted, fmchits_pathlendrifted, fmchits_theta;
      float  fmchits_phi, fmchits_mass, fmcconehits_vx, fmcconehits_vy, fmcconehits_vz, fmcconehits_t, fmcconehits_endx;
      float  fmcconehits_endy, fmcconehits_endz, fmcconehits_endt, fmcconehits_px, fmcconehits_py, fmcconehits_pz, fmcconehits_momentum;
      float  fmcconehits_energy, fmcconehits_endpx, fmcconehits_endpy, fmcconehits_endpz, fmcconehits_endenergy, fmcconehits_pathlen, fmcconehits_vxdrifted;
      float  fmcconehits_vydrifted, fmcconehits_vzdrifted, fmcconehits_tdrifted, fmcconehits_endxdrifted, fmcconehits_endydrifted, fmcconehits_endzdrifted, fmcconehits_endtdrifted, fmcconehits_pxdrifted;
      float  fmcconehits_pydrifted, fmcconehits_pzdrifted, fmcconehits_momentumdrifted, fmcconehits_energydrifted, fmcconehits_endpxdrifted, fmcconehits_endpydrifted, fmcconehits_endpzdrifted;
      float  fmcconehits_endenergydrifted, fmcconehits_pathlendrifted, fmcconehits_theta, fmcconehits_phi, fmcconehits_mass, fmcd_vx, fmcd_vy;
      float  fmcd_vz, fmcd_t, fmcd_endx, fmcd_endy, fmcd_endz, fmcd_endt, fmcd_px;
      float  fmcd_py, fmcd_pz, fmcd_momentum, fmcd_energy, fmcd_truecuthitsE, ftrueEdepo, fmcd_endpx;
      float  fmcd_endpy, fmcd_endpz, fmcd_endenergy, fmcd_pathlen, fmcd_vxdrifted, fmcd_vydrifted, fmcd_vzdrifted;
      float  fmcd_tdrifted, fmcd_endxdrifted, fmcd_endydrifted, fmcd_endzdrifted, fmcd_endtdrifted, fmcd_pxdrifted, fmcd_pydrifted;
      float  fmcd_pzdrifted, fmcd_momentumdrifted, fmcd_energydrifted, fmcd_endpxdrifted, fmcd_endpydrifted, fmcd_endpzdrifted, fmcd_endenergydrifted;
      float  fmcd_pathlendrifted, fmcd_theta, fmcd_phi, fmcd_mass, fmcshwr_vx, fmcshwr_vy, fmcshwr_vz;
      float  fmcshwr_t, fmcshwr_endx, fmcshwr_endy, fmcshwr_endz, fmcshwr_endt, fmcshwr_px, fmcshwr_py;
      float  fmcshwr_pz, fmcshwr_momentum, fmcshwr_energy, fmcshwr_endpx, fmcshwr_endpy, fmcshwr_endpz, fmcshwr_endenergy;
      float  fmcshwr_pathlen, fmcshwr_vxdrifted, fmcshwr_vydrifted, fmcshwr_vzdrifted, fmcshwr_tdrifted, fmcshwr_endxdrifted, fmcshwr_endydrifted;
      float  fmcshwr_endzdrifted, fmcshwr_endtdrifted, fmcshwr_pxdrifted, fmcshwr_pydrifted, fmcshwr_pzdrifted, fmcshwr_momentumdrifted, fmcshwr_energydrifted;
      float  fmcshwr_endpxdrifted, fmcshwr_endpydrifted, fmcshwr_endpzdrifted, fmcshwr_endenergydrifted, fmcshwr_pathlendrifted, fmcshwr_theta, fmcshwr_phi, fmcshwr_mass;
      float  fflashcut_track_dist_int, fflashcut_time_int, fflashcut_pe_int, fflashcut_ycenter_int, fflashcut_zcenter_int, fflashcut_ywidth_int, fflashcut_zwidth_int, fflashcut_timewidth_int, fflashcut_abstime_int, fflashcut_frame_int;
      float  favg_missing_energy, favg_missing_numelec;
      
      double genIniPosX,genIniPosY,genIniPosZ,genIniPosT;
      double genEndPosX,genEndPosY,genEndPosZ,genEndPosT;
      double cryIniPosX, cryIniPosY,cryIniPosZ,cryIniPosT;
      double cryEndPosX, cryEndPosY,cryEndPosZ,cryEndPosT;
      double genTrueEne,genTrueMom,genIniMomX,genIniMomY,genIniMomZ,genIniMomE;
      double genEndMomX,genEndMomY,genEndMomZ,genEndMomE;
      double cryTrueEne,cryTrueMom,cryIniMomX,cryIniMomY,cryIniMomZ,cryIniMomE;
      double cryEndMomX,cryEndMomY,cryEndMomZ,cryEndMomE;

      bool   fMCIsPrimary;

      double fEvtTime, fcutEvtTime, fDriftV, fTrigTime, fEField, fTemper;
      double _fPFP_trks=0, _ftruePFP_trks=0, _fT0_trks=0, _ftrueT0_trks=0, _ftrueshwr_trks=0;
      double _fstartinbound_trks=0, _ftruestartinbound_trks=0, _fendinFV_trks=0, _ftrueMichelcount=0;
      double _ftrueendinFV_trks=0, _fccrosser_trks=0, _ftrueccrosser_trks=0, _fnoelecdiv_bound=0;
      double _ftruenoelecdiv_bound=0, _funbroken_trks=0, _ftrueunbroken_trks=0, _flongcm_trks=0;
      double _ftruelongcm_trks=0, _fminhitpt_trks=0, _ftrueminhitpt_trks=0, _fmaxhitpt_trks=0;
      double _ftruemaxhitpt_trks=0, _fnearhits_trks=0, _ftruenearhits_trks=0, _fshwr_trks=0;

      static constexpr int kNMaxMCParticles = 1000000;
      static constexpr int kNMaxPFParticles = 2000;
      static constexpr int kNMaxPFPClusters = 100;
      static constexpr int kNViews= 3; //constexpr int kMaxPlanes     = 3;
      static constexpr int kMaxCry= 5000;
      static constexpr int kMaxGen= 5000;
      static constexpr int kMaxTracks= 5000;
      static constexpr int kMaxPlanes= 3;
      static constexpr int kMaxHits= 10000;
      static constexpr int kMaxFlashes      = 10000;   //maximum number of flashes
      static constexpr int kMaxCh= 30000; //maximum number of channels;   
      static constexpr int kMaxTrackHits    = 3000;
      static constexpr int kMaxDaugh= 10000;
      static constexpr int kMaxWf= 90000;
      static constexpr int kMaxct= 1000;

      int fNCloseHitsIni[3];               // Number of hits within selection radius
      int fNCloseHitsEnd[3];               // Number of hits within selection radius

      // Purity and Efficiency Numbers
      int fNMichel_gen;     // Number of Michel electrons in the event
      int fNMichel_cry;     // Number of Michel electrons in the event
      int fYesSelected = 0; // Number of correctly selected events for each set of reco params
      int fNoSelected = 0;  // Number of wrongly selected events for each set of reco params

      //kMaxTracks
      std::vector<int>    NpfpsFromTrack, ft0sana, fstartinboundsana, fendinFVsana, fccrosserana, fbrokencountana, fnearhitcountana;
      std::vector<int>    fMichelcountt0ana, NMichelpfpFromTrack, fMichelcountendinFVana, fMichelcountstartinboundana, fMichelcountelecdivstopzana, fMichelcountccrosserana;
      std::vector<int>    fMichelcountbrokencountana, fMichelcountlenana, fMichelcountminhitptimeana, fMichelcountmaxhitptimeana, fMichelcountnearhitana, fMichelcountshwrdistana;
      std::vector<float>  felecdivstopzana, fdistanceana, ftrklengthana, fminhitptimeana, fmaxhitptimeana, fnshwrdistana, fshwrdistana;
      std::vector<float>  mcFromTrack_michelE, ftrueEt0ana, ftrueEstartinboundana, ftrueEendinFVana, ftrueEccrosserana, ftrueEelecdivstopzana, ftrueEbrokencountana;  
      std::vector<float>  ftrueElenana, ftrueEminhitptimeana, ftrueEmaxhitptimeana, ftrueEnearhitana, ftrueEshwrdistana;
      // kMaxFlashes
      std::vector<float>  flash_time, flash_abstime, flash_totalpe, flash_ycenter, flash_zcenter, flash_ywidth, flash_zwidth;
      std::vector<float>  flash_twidth, flash_frame, flash_distana, flash_peana;
      // kMaxHits
      std::vector<int>    fhits_TPC, fhits_mult, fshwrhits_chno, fshwrhits_wire, fshwrhits_TPC, fshwrhits_plane, fshwrhits_mult;
      std::vector<int>    fshwrallhits_chno, fshwrallhits_wire, fshwrallhits_TPC, fshwrallhits_plane, fshwrallhits_mult, fnearhits_key, fnearhits_chno;
      std::vector<int>    fnearhits_wire, fnearhits_TPC, fnearhits_plane, fnearhits_mult, fmhits_key, fmhits_chno, fmhits_wire;
      std::vector<int>    fmhits_TPC, fmhits_plane, fmhits_mult, fmhits_longtrk, fmhits_sametrk, fmhits_corrhit, ftrueparhitsall_key;
      std::vector<int>    ftrueparhitsall_chno, ftrueparhitsall_wire, ftrueparhitsall_TPC, ftrueparhitsall_plane, ftrueparhitsall_mult, ftrueparhitscol_key, ftrueparhitscol_chno;
      std::vector<int>    ftrueparhitscol_wire, ftrueparhitscol_TPC, ftrueparhitscol_plane, ftrueparhitscol_mult;
      std::vector<float>  fhits_key, fhits_charge, fhits_wire, fhits_peakT, fhits_chno, fhits_xpos, fhits_ypos, fhits_zpos;
      std::vector<float>  fhits_sigptime, fhits_sigchrg, fhits_sigpamp, fhits_dof, fhits_gof, fhits_ptplusRMS, fhits_ptminusRMS;
      std::vector<float>  fhits_cnnMichel, fhits_cnnEM, fhits_cnnTrack, fshwrhits_peakT, fshwrhits_charge, fshwrhits_xpos, fshwrhits_ypos;
      std::vector<float>  fshwrhits_zpos,fshwrhits_sigptime,fshwrhits_sigchrg,fshwrhits_sigpamp,fshwrhits_dof,fshwrhits_gof,fshwrhits_ptplusRMS;
      std::vector<float>  fshwrhits_ptminusRMS, fshwrallhits_peakT, fshwrallhits_charge, fshwrallhits_xpos, fshwrallhits_ypos, fshwrallhits_zpos, fshwrallhits_sigptime;
      std::vector<float>  fshwrallhits_sigchrg, fshwrallhits_sigpamp, fshwrallhits_dof, fshwrallhits_gof, fshwrallhits_ptplusRMS, fshwrallhits_ptminusRMS, ftrkdqdxU;
      std::vector<float>  ftrkdedxU, ftrkresrangeU, ftrkhitxU, ftrkhityU, ftrkhitzU, ftrkpitchU, ftrkdqdxV;
      std::vector<float>  ftrkdedxV, ftrkresrangeV, ftrkhitxV, ftrkhityV, ftrkhitzV, ftrkpitchV, ftrkdqdxY;
      std::vector<float>  ftrkdedxY, ftrkresrangeY, ftrkhitxY, ftrkhityY, ftrkhitzY, ftrkpitchY, fnearhits_peakT;
      std::vector<float>  fnearhits_charge, fnearhits_xpos, fnearhits_ypos, fnearhits_zpos, fnearhits_energy, fnearhits_sigptime, fnearhits_sigchrg;
      std::vector<float>  fnearhits_sigpamp, fnearhits_dof, fnearhits_gof, fnearhits_ptplusRMS, fnearhits_ptminusRMS, fnearhits_cnnMichel, fnearhits_cnnEM;
      std::vector<float>  fnearhits_cnnTrack, fmhits_peakT, fmhits_charge, fmhits_xpos, fmhits_ypos, fmhits_zpos, fmhits_sigptime;
      std::vector<float>  fmhits_sigchrg, fmhits_sigpamp, fmhits_dof, fmhits_gof, fmhits_angledeg, fmhits_maghitveccostheta, fmhits_distance;
      std::vector<float>  fmhits_ptplusRMS, fmhits_ptminusRMS, fmhits_cnnMichel, fmhits_cnnEM, fmhits_cnnTrack, ftrueparhitsall_peakT, ftrueparhitsall_charge;
      std::vector<float>  ftrueparhitsall_xpos, ftrueparhitsall_ypos, ftrueparhitsall_zpos, ftrueparhitsall_sigptime, ftrueparhitsall_sigchrg, ftrueparhitsall_sigpamp, ftrueparhitsall_dof;
      std::vector<float>  ftrueparhitsall_gof, ftrueparhitsall_ptplusRMS, ftrueparhitsall_ptminusRMS, ftrueMiEFrac, ftrueparhitscol_peakT, ftrueparhitscol_charge, ftrueparhitscol_xpos;
      std::vector<float>  ftrueparhitscol_ypos, ftrueparhitscol_zpos, ftrueparhitscol_sigptime, ftrueparhitscol_sigchrg, ftrueparhitscol_sigpamp, ftrueparhitscol_dof, ftrueparhitscol_gof;
      std::vector<float>  ftrueparhitscol_angledeg, ftrueparhitscol_maghitveccostheta, ftrueparhitscol_distance, ftrueparhitscol_ptplusRMS, ftrueparhitscol_ptminusRMS;
      // kMaxct
      std::vector<float>  fmichel_zpos, fmichel_ypos, fmichel_xpos;
      std::vector<float>  fmichel_chrg, fmichel_chno, fmichel_key, fmichel_wire;
      std::vector<float>  fmichel_chargehit, fmichel_tpc, fmichel_ptime, fmichel_angledeg;
      std::vector<float>  fmichel_maghitveccostheta, fmichel_distance, fmichel_mult, fmichel_sigptime;
      std::vector<float>  fmichel_sigchrg, fmichel_sigpamp, fmichel_dof, fmichel_gof;
      std::vector<float>  fmichel_ptminusRMS, fmichel_ptplusRMS, fmichel_status, fmichel_cnnMichel;
      std::vector<float>  fmichel_cnnEM, fmichel_cnnTrack;
      // kMaxWf
      std::vector<int>    fwfchan, fwfcut_endhitkey, fwfcut_endwire;
      std::vector<int>    fwfcut_endchno, fwfcut_endtpcno, fwfcut_endhitmult;
      std::vector<float>  fwfcut_endhitx, fwfcut_endhity, fwfcut_endhitz;
      std::vector<float>  fwfcut_endhitsigptime, fwfcut_endhitsigchrg, fwfcut_endhitsigpamp;
      std::vector<float>  fwfcut_endhitdof, fwfcut_endhitgof, fwfcut_endhitchrg, fwfcut_endptime;
      std::vector<double> fwftimeint, fwftimeext, ftrkrecotime, fwftime, fwftracktimediff;

    ///////////////////////////

      std::vector<int>    gen_g4PdgCode,gen_g4TrackID,gen_g4MumsTID,fSimEdepTID,fSimEdepPDG;
      std::vector<double> gen_g4IniPosX,gen_g4IniPosY,gen_g4IniPosZ,gen_g4IniPosT;
      std::vector<double> gen_g4EndPosX,gen_g4EndPosY,gen_g4EndPosZ,gen_g4EndPosT;
      std::vector<double> gen_g4TrueEne,gen_g4TrueMom;
      std::vector<double> fSimEdepE,fSimEdepX,fSimEdepY,fSimEdepZ;

      std::vector<int>    cry_g4PdgCode, cry_g4TrackID, cry_g4MumsTID;
      std::vector<double> cry_g4TrueEne, cry_g4TrueMom, cry_g4IniPosX, cry_g4IniPosY, cry_g4IniPosZ, cry_g4EndPosX, cry_g4EndPosY, cry_g4EndPosZ; 

      std::vector<int>    gen_g4NHits,fPFPID,fPFPPdgCode,fPFPNChildren,fPFPMumsTID,fPFPNClusters,fPFPTrackID;
      std::vector<int>    fPFPTrueParticleMatchedID,fPFPTrueParticleMatchedPosition,fPFPNSharedTrueParticleHits;
      std::vector<int>    fPFPNHits,fPFPShowerID,fPFPShowerBestPlane;
      std::vector<bool>   fPFPIsPrimary,fPFPIsTrack,fPFPIsShower;
      std::vector<double> fPFPTrackLength,fPFPTrackIniX,fPFPTrackIniY,fPFPTrackIniZ,fPFPTrackVertexX,fPFPTrackVertexY,fPFPTrackVertexZ;
      std::vector<double> fPFPTrackEndX,fPFPTrackEndY,fPFPTrackEndZ,fPFPTrackTheta,fPFPTrackPhi,fPFPTrackZenithAngle,fPFPTrackAzimuthAngle;
      std::vector<double> fPFPTrackIniDirectionX,fPFPTrackIniDirectionY,fPFPTrackIniDirectionZ,fPFPTrackVertexDirectionX,fPFPTrackVertexDirectionY,fPFPTrackVertexDirectionZ;
      std::vector<double> fPFPTrackEndDirectionX,fPFPTrackEndDirectionY,fPFPTrackEndDirectionZ,fPFPTrackChi2,fPFPTrackNdof;
      std::vector<double> fPFPShowerDirectionX,fPFPShowerDirectionY,fPFPShowerDirectionZ,fPFPShowerDirectionErrX,fPFPShowerDirectionErrY,fPFPShowerDirectionErrZ;
      std::vector<double> fPFPShowerIniX,fPFPShowerIniY,fPFPShowerIniZ,fPFPShowerIniErrX,fPFPShowerIniErrY,fPFPShowerIniErrZ;
      std::vector<double> fPFPShowerLength,fPFPShowerOpenAngle,fPFPCompleteness,fPFPPurity;
      std::vector<std::vector<double>> fPFPShowerdEdx;

      TCanvas * c;

      // Declare analysis utils
      protoana::ProtoDUNETruthUtils truthUtil;
      protoana::ProtoDUNETrackUtils trackUtil;
      protoana::ProtoDUNEPFParticleUtils pfpUtil;
	}; //End class definition


	MichelsAna::MichelsAna(fhicl::ParameterSet const & p) : EDAnalyzer{p}, utils(new michels::MichelsBasicUtils(p)) //, 
	// More initializers here.
	{
		art::ServiceHandle<art::TFileService> tfs;
		std::vector<std::vector<std::string>> TreesToWrite   = p.get<std::vector<std::vector<std::string>>>("TreesToWrite");
		std::vector<std::vector<std::string>> ProductsToDump = p.get<std::vector<std::vector<std::string>>>("ProductsToDump");
		fTreesToWrite   = TreesToWrite;
		fProductsToDump = ProductsToDump;
    debug               = p.get<bool>("debug");
    fRadiusThreshold    = p.get<double>("RadiusThreshold");
    fCNNThreshold       = p.get<double>("CNNThreshold");  
    fCloseHitsThreshold = p.get<int>("CloseHitsThreshold"); 


    c = new TCanvas( "canv", "canv" );


		std::vector<std::string> trees; 
		int nTrees = fTreesToWrite.size();
		for (auto t = 0; t < nTrees; t++) { auto tree = fTreesToWrite[t][0]; trees.push_back(tree); }  

		if (std::find(trees.begin(), trees.end(), std::string("Gen")) != trees.end()) 
		{
			Gen  = tfs->make<TTree>("Gen", "Gen");
			for (auto i : fProductsToDump)
			{
				auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
				TString out_label = o_label;

				// Set branches: Ifs over the different types of data
				if (o_label == "Gen_MCTruth")
				{
					//Event branches//
					Gen->Branch("Event",  &fEvent,  "Event/I");
					Gen->Branch("Run",    &fRun,    "Run/I");
					Gen->Branch("Subrun", &fSubRun, "subRun/I");
					Gen->Branch("Flag",   &fFlag,   "Flag/I");
					Gen->Branch("TriggerTime",   &fTrigTime,      "TrigTime/D");
					Gen->Branch("DriftVelocity", &fDriftV,        "DriftVel/D");
					Gen->Branch("EventTime",     &fEvtTime,       "fEvtTime/D");
					Gen->Branch("YearMonDay",    &fyear_mon_day,  "fyear_mon_day/I");
					Gen->Branch("HourMinSec",    &fhour_min_sec,  "fhour_min_sec/I");

					//Generator truth info branches//
					Gen->Branch("genNMichel",    &fNMichel_gen,  "fNMichel_gen/I");
					Gen->Branch("genOrigin",     &genOrigin,     "genOrigin/I");
					Gen->Branch("genNDaughters", &genNDaughters, "genNDaughters/I");
					Gen->Branch("genPdgCode",    &genPdgCode);
					Gen->Branch("genTrackID",    &genTrackID);
          // The following branches will be filled when the primary=muon and have Michels candidates
					Gen->Branch("genIniPosX",   &genIniPosX);
					Gen->Branch("genIniPosY",   &genIniPosY);
					Gen->Branch("genIniPosZ",   &genIniPosZ);
					Gen->Branch("genIniPosT",   &genIniPosT);
					Gen->Branch("genEndPosX",   &genEndPosX);
					Gen->Branch("genEndPosY",   &genEndPosY);
					Gen->Branch("genEndPosZ",   &genEndPosZ);
					Gen->Branch("genEndPosT",   &genEndPosT);
					Gen->Branch("genTrueEne",   &genTrueEne);
					Gen->Branch("genTrueMom",   &genTrueMom);
					Gen->Branch("genIniMomX",   &genIniMomX);
					Gen->Branch("genIniMomY",   &genIniMomY);
					Gen->Branch("genIniMomZ",   &genIniMomZ);
					Gen->Branch("genEndMomX",   &genEndMomX);
					Gen->Branch("genEndMomY",   &genEndMomY);
					Gen->Branch("genEndMomZ",   &genEndMomZ);
					Gen->Branch("gen_g4PdgCode", &gen_g4PdgCode);
					Gen->Branch("gen_g4TrackID", &gen_g4TrackID);
					Gen->Branch("gen_g4MumsTID", &gen_g4MumsTID);
					Gen->Branch("gen_g4TrueEne", &gen_g4TrueEne);
					Gen->Branch("gen_g4TrueMom", &gen_g4TrueMom);
					Gen->Branch("gen_g4IniPosX", &gen_g4IniPosX);
					Gen->Branch("gen_g4IniPosY", &gen_g4IniPosY);
					Gen->Branch("gen_g4IniPosZ", &gen_g4IniPosZ);
					Gen->Branch("gen_g4EndPosX", &gen_g4EndPosX);
					Gen->Branch("gen_g4EndPosY", &gen_g4EndPosY);
					Gen->Branch("gen_g4EndPosZ", &gen_g4EndPosZ);
				}//if Gen_MCTruth
				else if (i_type == "sim::SimEnergyDeposit")
				{
					// ENERGY DEPOSITED //
					Gen->Branch("SimEdepTID",     &fSimEdepTID);
					Gen->Branch("SimEdepPDG",     &fSimEdepPDG);
					Gen->Branch("SimEdepEnergy",  &fSimEdepE );
					Gen->Branch("SimEdepMiddleX", &fSimEdepX );
					Gen->Branch("SimEdepMiddleY", &fSimEdepY );
					Gen->Branch("SimEdepMiddleZ", &fSimEdepZ );
				}// if sim::SimEnergyDeposit
			} //for loop over products to dump
		} //if truth tree

		if (std::find(trees.begin(), trees.end(), std::string("Cry")) != trees.end()) 
		{
			Cry  = tfs->make<TTree>("Cry", "Cry");
			for (auto i : fProductsToDump)
			{
				auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
				TString out_label = o_label;

				// Set branches: Ifs over the different types of data
				if (o_label == "Cry_MCTruth")
				{
					//Event branches//
					Cry->Branch("Event",  &fEvent,  "Event/I");
					Cry->Branch("Run",    &fRun,    "Run/I");
					Cry->Branch("Subrun", &fSubRun, "subRun/I");
					Cry->Branch("TruthFlag",   &fFlag,   "Flag/I");
					Cry->Branch("TriggerTime",    &fTrigTime,     "TrigTime/D");
					Cry->Branch("DriftVelocity",  &fDriftV,       "DriftVel/D");
					Cry->Branch("EventTime",      &fEvtTime,      "fEvtTime/D");
					Cry->Branch("YearMonDay",     &fyear_mon_day, "fyear_mon_day/I");
					Cry->Branch("HourMinSec",     &fhour_min_sec, "fhour_min_sec/I");

					//MC truth branches//
					Cry->Branch("cryNMichel",    &fNMichel_cry,   "fNMichel_cry/I");
					Cry->Branch("cryNDaughters", &cryNDaughters, "cryNDaughters/I");
					Cry->Branch("cryOrigin",     &cryOrigin,     "cryOrigin/I");
					Cry->Branch("cryPdgCode",    &cryPdgCode);
					Cry->Branch("cryTrackID",    &cryTrackID);
          // The following branches will be filled when the primary=muon and have Michels candidates
					Cry->Branch("cryIniPosX",   &cryIniPosX);
					Cry->Branch("cryIniPosY",   &cryIniPosY);
					Cry->Branch("cryIniPosZ",   &cryIniPosZ);
					Cry->Branch("cryIniPosT",   &cryIniPosT);
					Cry->Branch("cryEndPosX",   &cryEndPosX);
					Cry->Branch("cryEndPosY",   &cryEndPosY);
					Cry->Branch("cryEndPosZ",   &cryEndPosZ);
					Cry->Branch("cryEndPosT",   &cryEndPosT);
					Cry->Branch("cryTrueEne",   &cryTrueEne);
					Cry->Branch("cryTrueMom",   &cryTrueMom);
					Cry->Branch("cryIniMomX",   &cryIniMomX);
					Cry->Branch("cryIniMomY",   &cryIniMomY);
					Cry->Branch("cryIniMomZ",   &cryIniMomZ);
					Cry->Branch("cryEndMomX",   &cryEndMomX);
					Cry->Branch("cryEndMomY",   &cryEndMomY);
					Cry->Branch("cryEndMomZ",   &cryEndMomZ);
					Cry->Branch("cry_g4PdgCode", &cry_g4PdgCode);
					Cry->Branch("cry_g4TrackID", &cry_g4TrackID);
					Cry->Branch("cry_g4MumsTID", &cry_g4MumsTID);
					Cry->Branch("cry_g4TrueEne", &cry_g4TrueEne);
					Cry->Branch("cry_g4TrueMom", &cry_g4TrueMom);
					Cry->Branch("cry_g4IniPosX", &cry_g4IniPosX);
					Cry->Branch("cry_g4IniPosY", &cry_g4IniPosY);
					Cry->Branch("cry_g4IniPosZ", &cry_g4IniPosZ);
					Cry->Branch("cry_g4EndPosX", &cry_g4EndPosX);
					Cry->Branch("cry_g4EndPosY", &cry_g4EndPosY);
					Cry->Branch("cry_g4EndPosZ", &cry_g4EndPosZ);
				} //if Cry_MCTruth
			} //for loop over products to dump
		} //if cry tree

		if (std::find(trees.begin(), trees.end(), std::string("AllReco")) != trees.end()) 
		{
			AllReco   = tfs->make<TTree>("AllReco",  "AllReco");

			//Event branches//
			AllReco->Branch("Event",  &fEvent,  "Event/I");
			AllReco->Branch("Run",    &fRun,    "Run/I");
			AllReco->Branch("Subrun", &fSubRun, "subRun/I");
			AllReco->Branch("Flag",   &fFlag,   "Flag/I");
			AllReco->Branch("TriggerTime",   &fTrigTime, "TrigTime/D");
			AllReco->Branch("DriftVelocity", &fDriftV,   "DriftVel/D");
			AllReco->Branch("EventTime",     &fEvtTime,  "fEvtTime/D");
			AllReco->Branch("YearMonDay",    &fyear_mon_day,  "fyear_month_day/I");
			AllReco->Branch("HourMinSec",    &fhour_min_sec,  "fhour_min_sec/I");

			//PFP branches
			// AllReco->Branch("nPFParticles",&fNPFParticles   );
			// AllReco->Branch("pfpTrueParticleMatchedID",&fPFPTrueParticleMatchedID      );
			// AllReco->Branch("pfpTrueParticleMatchedPosition",     &fPFPTrueParticleMatchedPosition);
			// AllReco->Branch("pfpIsPrimary",&fPFPIsPrimary   );
			// AllReco->Branch("pfpID",&fPFPID);
			// AllReco->Branch("pfpMumsTID",&fPFPMumsTID    );
			// AllReco->Branch("pfpPdgCode",&fPFPPdgCode     );
			// AllReco->Branch("pfpNChildren",&fPFPNChildren   );
			
			// AllReco->Branch("pfpIsTrack",&fPFPIsTrack);
			// AllReco->Branch("pfpTrackID",&fPFPTrackID);
			// AllReco->Branch("pfpTrackLength",&fPFPTrackLength);
			// AllReco->Branch("pfpTrackIniX",&fPFPTrackIniX);
			// AllReco->Branch("pfpTrackIniY",&fPFPTrackIniY);
			// AllReco->Branch("pfpTrackIniZ",&fPFPTrackIniZ);
			// AllReco->Branch("pfpTrackVertexX",&fPFPTrackVertexX);
			// AllReco->Branch("pfpTrackVertexY",&fPFPTrackVertexY);
			// AllReco->Branch("pfpTrackVertexZ",&fPFPTrackVertexZ);
			// AllReco->Branch("pfpTrackEndX",&fPFPTrackEndX);
			// AllReco->Branch("pfpTrackEndY",&fPFPTrackEndY);
			// AllReco->Branch("pfpTrackEndZ",&fPFPTrackEndZ);
			// AllReco->Branch("pfpTrackTheta",&fPFPTrackTheta);
			// AllReco->Branch("pfpTrackPhi",&fPFPTrackPhi);
			// AllReco->Branch("pfpTrackZenithAngle",      &fPFPTrackZenithAngle      );
			// AllReco->Branch("pfpTrackAzimuthAngle",     &fPFPTrackAzimuthAngle     );
			// AllReco->Branch("pfpTrackIniDirectionX",  &fPFPTrackIniDirectionX  );
			// AllReco->Branch("pfpTrackIniDirectionY",  &fPFPTrackIniDirectionY  );
			// AllReco->Branch("pfpTrackIniDirectionZ",  &fPFPTrackIniDirectionZ  );
			// AllReco->Branch("pfpTrackVertexDirectionX", &fPFPTrackVertexDirectionX );
			// AllReco->Branch("pfpTrackVertexDirectionY", &fPFPTrackVertexDirectionY );
			// AllReco->Branch("pfpTrackVertexDirectionZ", &fPFPTrackVertexDirectionZ );
			// AllReco->Branch("pfpTrackEndDirectionX",    &fPFPTrackEndDirectionX    );
			// AllReco->Branch("pfpTrackEndDirectionY",    &fPFPTrackEndDirectionY    );
			// AllReco->Branch("pfpTrackEndDirectionZ",    &fPFPTrackEndDirectionZ    );
			// AllReco->Branch("pfpTrackChi2",&fPFPTrackChi2);
			// AllReco->Branch("pfpTrackIniNdof",&fPFPTrackNdof);

			// AllReco->Branch("pfpNClusters",&fPFPNClusters);
			// AllReco->Branch("pfpCluNHits",&fPFPCluNHits);

			// AllReco->Branch("pfpIsShower",&fPFPIsShower);
			// AllReco->Branch("pfpShowerID",&fPFPShowerID);
			// AllReco->Branch("pfpShowerBestPlane",     &fPFPShowerBestPlane     );
			// AllReco->Branch("pfpShowerDirectionX",    &fPFPShowerDirectionX    );
			// AllReco->Branch("pfpShowerDirectionY",    &fPFPShowerDirectionY    );
			// AllReco->Branch("pfpShowerDirectionZ",    &fPFPShowerDirectionZ    );
			// AllReco->Branch("pfpShowerIniX",&fPFPShowerIniX);
			// AllReco->Branch("pfpShowerIniY",&fPFPShowerIniY);
			// AllReco->Branch("pfpShowerIniZ",&fPFPShowerIniZ);
			// AllReco->Branch("pfpShowerLength",&fPFPShowerLength);
			// AllReco->Branch("pfpShowerOpeningAngle",  &fPFPShowerOpenAngle     );
			
			// AllReco->Branch("pfpShowerdEdx",&fPFPShowerdEdx);
			// AllReco->Branch("pfpCompleteness",      &fPFPCompleteness);
			// AllReco->Branch("pfpPurity",&fPFPPurity );

			AllReco->Branch("fNFlashes",       &fNFlashes    );
			AllReco->Branch("flash_time",      &flash_time   );
			AllReco->Branch("flash_totalpe",   &flash_totalpe);
			AllReco->Branch("flash_ycenter",   &flash_ycenter);
			AllReco->Branch("flash_zcenter",   &flash_zcenter);
			AllReco->Branch("flash_ywidth",    &flash_ywidth );
			AllReco->Branch("flash_zwidth",    &flash_zwidth );
			AllReco->Branch("flash_twidth",    &flash_twidth );
			AllReco->Branch("flash_abstime",   &flash_abstime);
			AllReco->Branch("flash_frame",     &flash_frame  );

			AllReco->Branch("NTracks"       , &NTracks       , "NTracks/I");
			AllReco->Branch("NpfpsFromTrack", &NpfpsFromTrack);
			// AllReco->Branch("fPFP_trks",&fPFP_trks,"fPFP_trks/I");
			// AllReco->Branch("fT0_trks",&fT0_trks,"fT0_trks/I");
			// AllReco->Branch("fstartinbound_trks", &fstartinbound_trks, "fstartinbound_trks/I");
			// AllReco->Branch("fendinFV_trks",      &fendinFV_trks,      "fendinFV_trks/I");
			// AllReco->Branch("fccrosser_trks",     &fccrosser_trks,     "fccrosser_trks/I");
			// AllReco->Branch("fnoelecdiv_bound",   &fnoelecdiv_bound,   "fnoelecdiv_bound/I");
			// AllReco->Branch("funbroken_trks",     &funbroken_trks,     "funbroken_trks/I");
			// AllReco->Branch("flongcm_trks",&flongcm_trks,"flongcm_trks/I");
			// AllReco->Branch("fminhitpt_trks",     &fminhitpt_trks,     "fminhitpt_trks/I");
			// AllReco->Branch("fmaxhitpt_trks",     &fmaxhitpt_trks,     "fmaxhitpt_trks/I");
			// //    fEventTree->Branch("fmuendy_trks",&fmuendy_trks,"fmuendy_trks/I");
			// //    fEventTree->Branch("fmuendz_trks",&fmuendz_trks,"fmuendz_trks/I");
			// //    fEventTree->Branch("fdistmorehits_trks",&fdistmorehits_trks,"fdistmorehits_trks/I");
			// //    fEventTree->Branch("fPHtest_trks",&fPHtest_trks,"fPHtest_trks/I");
			// AllReco->Branch("fnearhits_trks", &fnearhits_trks, "fnearhits_trks/I");
			// AllReco->Branch("fcut_mu",&fcut_mu,"fcut_mu/I");
			// AllReco->Branch("funcont_trks",   &funcont_trks,   "funcont_trks/I");
			// AllReco->Branch("fcrossZ_trks",   &fcrossZ_trks,   "fcrossZ_trks/I");
			// AllReco->Branch("fbackward_trks", &fbackward_trks, "fbackward_trks/I");
			// AllReco->Branch("fEndinTPC_trks", &fEndinTPC_trks, "fEndinTPC_trks/I");
			// AllReco->Branch("NTrueMichels",   &NTrueMichels,   "NTrueMichels/I");

			// AllReco->Branch("fpfpana",&fpfpana,"fpfpana/I");  
			// AllReco->Branch("ft0ana",&ft0ana,"ft0ana/I");  
			// AllReco->Branch("fstartinboundana",&fstartinboundana,"fstartinboundana/I");  
			// AllReco->Branch("fendinFVana",     &fendinFVana,     "fendinFVana/I");    
			// AllReco->Branch("fccrossana",      &fccrossana,      "fccrossana/I");
			// AllReco->Branch("fstopzana",&fstopzana,"fstopzana/I");
			// AllReco->Branch("fdistana",&fdistana,"fdistana/I");
			// AllReco->Branch("fbrokcoana",      &fbrokcoana,      "fbrokcoana/I");
			// AllReco->Branch("fminhitptana",    &fminhitptana,    "fminhitptana/I");
			// AllReco->Branch("fmaxhitptana",    &fmaxhitptana,    "fmaxhitptana/I");
			// //    fEventTree->Branch("fmuendyana",&fmuendyana,"fmuendyana/I");
			// //    fEventTree->Branch("fmuendzana",&fmuendzana,"fmuendzana/I");
			// //    fEventTree->Branch("fPHana",&fPHana,"fPHana/I");
			// AllReco->Branch("ftrklenana",     &ftrklenana,     "ftrklenana/I");
			// AllReco->Branch("ftrkdistcolana", &ftrkdistcolana, "ftrkdistcolana/I");
			// AllReco->Branch("fhitctana",      &fhitctana,   "fhitctana/I");
			// AllReco->Branch("fshwrdisana",    &fshwrdisana, "fshwrdisana/I");
			
			// AllReco->Branch("fpfpsana",fpfpsana,"fpfpsana[fpfpana]/I");  
			// AllReco->Branch("ft0sana",ft0sana,"ft0sana[ft0ana]/I");  
			// AllReco->Branch("fstartinboundsana",fstartinboundsana,"fstartinboundsana[fstartinboundana]/I");  
			// AllReco->Branch("fendinFVsana",     fendinFVsana,     "fendinFVsana[fendinFVana]/I");    
			// AllReco->Branch("fccrosserana",     fccrosserana,     "fccrosserana[fccrossana]/I");
			// AllReco->Branch("felecdivstopzana", felecdivstopzana, "felecdivstopzana[fstopzana]/F");
			// AllReco->Branch("fdistanceana",     fdistanceana,     "fdistanceana[fdistana]/F");
			// AllReco->Branch("fbrokencountana",  fbrokencountana,  "fbrokencountana[fbrokcoana]/I");
			// AllReco->Branch("ftrklengthana",    ftrklengthana,    "ftrklengthana[ftrklenana]/F");
			// AllReco->Branch("fminhitptimeana",  fminhitptimeana,  "fminhitptimeana[fminhitptana]/F");
			// AllReco->Branch("fmaxhitptimeana",  fmaxhitptimeana,  "fmaxhitptimeana[fmaxhitptana]/F");
			// //    fEventTree->Branch("fmuonendyana",fmuonendyana,"fmuonendyana[fmuendyana]/F");
			// //    fEventTree->Branch("fmuonendzana",fmuonendzana,"fmuonendzana[fmuendzana]/F");
			// //    fEventTree->Branch("ftrkdistcollhitsana",ftrkdistcollhitsana,"ftrkdistcollhitsana[ftrkdistcolana]/F");
			// //    fEventTree->Branch("fPHtestana",fPHtestana,"fPHtestana[fPHana]/I");
			// AllReco->Branch("fnearhitcountana",fnearhitcountana, "fnearhitcountana[fhitctana]/I");
			// AllReco->Branch("fnshwrdistana",   fnshwrdistana,    "fnshwrdistana[fhitctana]/F");
			// AllReco->Branch("fshwrdistana",    fshwrdistana,     "fshwrdistana[fshwrdisana]/F");
			
			// AllReco->Branch("NMichelpfpFromTrack",NMichelpfpFromTrack,"NMichelpfpFromTrack[fpfpana]/I");
			// AllReco->Branch("fMichelcountt0ana",fMichelcountt0ana,"fMichelcountt0ana[ft0ana]/I");
			// AllReco->Branch("fMichelcountstartinboundana", fMichelcountstartinboundana,"fMichelcountstartinboundana[fstartinboundana]/I");
			// AllReco->Branch("fMichelcountendinFVana",      fMichelcountendinFVana,     "fMichelcountendinFVana[fendinFVana]/I");
			// AllReco->Branch("fMichelcountccrosserana",     fMichelcountccrosserana,    "fMichelcountccrosserana[fccrossana]/I");
			// AllReco->Branch("fMichelcountelecdivstopzana", fMichelcountelecdivstopzana,"fMichelcountelecdivstopzana[fstopzana]/I");
			// AllReco->Branch("fMichelcountbrokencountana", fMichelcountbrokencountana, "fMichelcountbrokencountana[fbrokcoana]/I");
			// AllReco->Branch("fMichelcountlenana",fMichelcountlenana,"fMichelcountlenana[ftrklenana]/I");
			// AllReco->Branch("fMichelcountminhitptimeana", fMichelcountminhitptimeana, "fMichelcountminhitptimeana[fminhitptana]/I");
			// AllReco->Branch("fMichelcountmaxhitptimeana", fMichelcountmaxhitptimeana, "fMichelcountmaxhitptimeana[fmaxhitptana]/I");
			// //    fEventTree->Branch("fMichelcountdistana",fMichelcountdistana,"fMichelcountdistana[fdistana]/I");
			// //    fEventTree->Branch("fMichelcountmuonendyana",fMichelcountmuonendyana,"fMichelcountmuonendyana[fmuendyana]/I");
			// //    fEventTree->Branch("fMichelcountmuonendzana",fMichelcountmuonendzana,"fMichelcountmuonendzana[fmuendzana]/I");
			// //    fEventTree->Branch("fMichelcountdistcollana",fMichelcountdistcollana,"fMichelcountdistcollana[ftrkdistcolana]/I");
			// //    fEventTree->Branch("fMichelcountPHtestana",fMichelcountPHtestana,"fMichelcountPHtestana[fPHana]/I");
			// AllReco->Branch("fMichelcountnearhitana",  fMichelcountnearhitana, "fMichelcountnearhitana[fhitctana]/I");
			// AllReco->Branch("fMichelcountshwrdistana", fMichelcountshwrdistana,"fMichelcountshwrdistana[fshwrdisana]/I");
			//     // CutReco->Branch("favg_missing_energy",&favg_missing_energy,"favg_missing_energy/F");
			//     // CutReco->Branch("favg_missing_numelec",&favg_missing_numelec,"favg_missing_numelec/F");
			
			// AllReco->Branch("mcFromTrack_michelE",mcFromTrack_michelE,"mcFromTrack_michelE[fpfpana]/F");
			// AllReco->Branch("ftrueEt0ana",ftrueEt0ana,"ftrueEt0ana[ft0ana]/F");
			// AllReco->Branch("ftrueEstartinboundana", ftrueEstartinboundana, "ftrueEstartinboundana[fstartinboundana]/F");
			// AllReco->Branch("ftrueEendinFVana",      ftrueEendinFVana,      "ftrueEendinFVana[fendinFVana]/F");
			// AllReco->Branch("ftrueEccrosserana",     ftrueEccrosserana,     "ftrueEccrosserana[fccrossana]/F");
			// AllReco->Branch("ftrueEelecdivstopzana", ftrueEelecdivstopzana, "ftrueEtelecdivstopzana[fstopzana]/F");
			// //    fEventTree->Branch("ftrueEdistana",ftrueEdistana,"ftrueEdistana[fdistana]/F");
			// AllReco->Branch("ftrueEbrokencountana", ftrueEbrokencountana, "ftrueEbrokencountana[fbrokcoana]/F");
			// AllReco->Branch("ftrueElenana",ftrueElenana,"ftrueElenana[ftrklenana]/F");
			// AllReco->Branch("ftrueEminhitptimeana", ftrueEminhitptimeana, "ftrueEminhitptimeana[fminhitptana]/F");
			// AllReco->Branch("ftrueEmaxhitptimeana", ftrueEmaxhitptimeana, "ftrueEmaxhitptimeana[fmaxhitptana]/F");
			// //    fEventTree->Branch("fMichelcountmuonendyana",fMichelcountmuonendyana,"fMichelcountmuonendyana[fmuendyana]/I");
			// //    fEventTree->Branch("fMichelcountmuonendzana",fMichelcountmuonendzana,"fMichelcountmuonendzana[fmuendzana]/I");
			// //    fEventTree->Branch("fMichelcountdistcollana",fMichelcountdistcollana,"fMichelcountdistcollana[ftrkdistcolana]/I");
			// //    fEventTree->Branch("fMichelcountPHtestana",fMichelcountPHtestana,"fMichelcountPHtestana[fPHana]/I");
			// AllReco->Branch("ftrueEnearhitana", ftrueEnearhitana, "ftrueEnearhitana[fhitctana]/F");
			// AllReco->Branch("ftrueEshwrdistana",ftrueEshwrdistana,"ftrueEshwrdistana[fshwrdisana]/F");
		} //if reco tree

		if (std::find(trees.begin(), trees.end(), std::string("CutReco")) != trees.end()) 
		{
			CutReco   = tfs->make<TTree>("CutReco",  "CutReco");

			//Event branches//
			CutReco->Branch("Event",         &fcutEvent,     "fcutEvent/I");
			CutReco->Branch("Run",           &fcutRun,       "fcutRun/I");
			CutReco->Branch("Subrun",        &fcutSubRun,    "fcutsubRun/I");
			CutReco->Branch("Flag",          &fcutFlag,      "fcutFlag/I");
			CutReco->Branch("EventTime",     &fcutEvtTime,   "fcutEvtTime/D");
			CutReco->Branch("TriggerTime",   &fTrigTime,     "TrigTime/D");
			CutReco->Branch("DriftVelocity", &fDriftV,       "DriftVel/D");
			CutReco->Branch("EvtTime",       &fEvtTime,      "EvtTime/D");
			CutReco->Branch("YearMonDay",    &fyear_mon_day, "fyear_month_day/I");
			CutReco->Branch("HourMinSec",    &fhour_min_sec, "fhour_min_sec/I");

			//     CutReco->Branch("favg_missing_energy",&favg_missing_energy,"favg_missing_energy/F");
			//     CutReco->Branch("favg_missing_numelec",&favg_missing_numelec,"favg_missing_numelec/F");
			//     CutReco->Branch("fcut_endhitkey",&fcut_endhitkey,"fcut_endhitkey/I");
			//     CutReco->Branch("fcut_endwire",&fcut_endwire,"fcut_endwire/I");
			//     CutReco->Branch("fcut_endchno",&fcut_endchno,"fcut_endchno/I");
			//     CutReco->Branch("fcut_endtpcno",&fcut_endtpcno,"fcut_endtpcno/I");
			//     CutReco->Branch("fcut_endptime",&fcut_endptime,"fcut_endptime/F");
			//     CutReco->Branch("fcut_endhitchrg",&fcut_endhitchrg,"fcut_endhitchrg/F");
			//     CutReco->Branch("fcut_endhitx",&fcut_endhitx,"fcut_endhitx/F");
			//     CutReco->Branch("fcut_endhity",&fcut_endhity,"fcut_endhity/F");
			//     CutReco->Branch("fcut_endhitz",&fcut_endhitz,"fcut_endhitz/F");
			//     CutReco->Branch("fcut_ccrosser",&fcut_ccrosser,"fcut_ccrosser/I");
			//     CutReco->Branch("fcut_endhitmult",&fcut_endhitmult,"fcut_endhitmult/I");
			//     CutReco->Branch("fcut_endhitsigptime",&fcut_endhitsigptime,"fcut_endhitsigptime/F");
			//     CutReco->Branch("fcut_endhitsigchrg",&fcut_endhitsigchrg,"fcut_endhitsigchrg/F");
			//     CutReco->Branch("fcut_endhitsigpamp",&fcut_endhitsigpamp,"fcut_endhitsigpamp/F");
			//     CutReco->Branch("fcut_endhitdof",&fcut_endhitdof,"fcut_endhitdof/F");
			//     CutReco->Branch("fcut_endhitgof",&fcut_endhitgof,"fcut_endhitgof/F");
			//     CutReco->Branch("fcut_endhitptminusRMS",&fcut_endhitptminusRMS,"fcut_endhitptminusRMS/F");
			//     CutReco->Branch("fcut_endhitptplusRMS",&fcut_endhitptplusRMS,"fcut_endhitptplusRMS/F");

			//     CutReco->Branch("fMichelcountcutana",&fMichelcountcutana,"fMichelcountcutana/I");
			//     CutReco->Branch("fcut_dist_hit_end",&fcut_dist_hit_end,"fcut_dist_hit_end/F");
			//     CutReco->Branch("fcut_dist_times",&fcut_dist_times,"fcut_dist_times/F");
			//     CutReco->Branch("fcut_trackthetaxz",&fcut_trackthetaxz,"fcut_trackthetaxz/F");
			//     CutReco->Branch("fcut_trackthetayz",&fcut_trackthetayz,"fcut_trackthetayz/F");
			//     CutReco->Branch("fcut_trkstartx",&fcut_trkstartx,"fcut_trkstartx/F");
			//     CutReco->Branch("fcut_trkstarty",&fcut_trkstarty,"fcut_trkstarty/F");
			//     CutReco->Branch("fcut_trkstartz",&fcut_trkstartz,"fcut_trkstartz/F");
			//     CutReco->Branch("fcut_trkendx",&fcut_trkendx,"fcut_trkendx/F");
			//     CutReco->Branch("fcut_trkendy",&fcut_trkendy,"fcut_trkendy/F");
			//     CutReco->Branch("fcut_trkendz",&fcut_trkendz,"fcut_trkendz/F");
			//     CutReco->Branch("fcut_trkstartcosx",&fcut_trkstartcosx,"fcut_trkstartcosx/F");
			//     CutReco->Branch("fcut_trkstartcosy",&fcut_trkstartcosy,"fcut_trkstartcosy/F");
			//     CutReco->Branch("fcut_trkstartcosz",&fcut_trkstartcosz,"fcut_trkstartcosz/F");
			//     CutReco->Branch("fcut_trkendcosx",&fcut_trkendcosx,"fcut_trkendcosx/F");
			//     CutReco->Branch("fcut_trkendcosy",&fcut_trkendcosy,"fcut_trkendcosy/F");
			//     CutReco->Branch("fcut_trkendcosz",&fcut_trkendcosz,"fcut_trkendcosxz/F");
			//     CutReco->Branch("fcut_trklen",&fcut_trklen,"fcut_trklen/F");
			//     CutReco->Branch("fcut_trktheta",&fcut_trktheta,"fcut_trktheta/F");
			//     CutReco->Branch("fcut_trkphi",&fcut_trkphi,"fcut_trkphi/F");
			//     CutReco->Branch("fcut_trkID",&fcut_trkID,"fcut_trkID/I");
			// //    CutReco->Branch("fcut_PHratio",&fcut_PHratio,"fcut_PHratio/F");
			// //    CutReco->Branch("fcut_PH",&fcut_PH,"fcut_PH/I");
			//     CutReco->Branch("fcut_minhitptime",&fcut_minhitptime,"fcut_minhitptime/F");
			//     CutReco->Branch("fcut_maxhitptime",&fcut_maxhitptime,"fcut_maxhitptime/F");
			//     CutReco->Branch("fcut_ncolhits",&fcut_ncolhits,"fcut_ncolhits/I");
			//     CutReco->Branch("fcut_nearhitcount",&fcut_nearhitcount,"fcut_nearhitcount/I");
			//     CutReco->Branch("fcut_CorrWirePtime",&fcut_CorrWirePtime,"fcut_CorrWirePtime/F");
			//     CutReco->Branch("fcut_trkrecotime",&fcut_trkrecotime,"fcut_trkrecotime/F");
			// //    CutReco->Branch("fmorechargeend",&fmorechargeend,"fmorechargeend/I");
				
			//     CutReco->Branch("fcutshwr_key",&fcutshwr_key,"fcutshwr_key/I");
			//     CutReco->Branch("fcutshwr_ID",&fcutshwr_ID,"fcutshwr_ID/I");
			//     CutReco->Branch("fcutshwr_length",&fcutshwr_length,"fcutshwr_length/F");
			//     CutReco->Branch("fcutshwr_startx",&fcutshwr_startx,"fcutshwr_startx/F");
			//     CutReco->Branch("fcutshwr_starty",&fcutshwr_starty,"fcutshwr_starty/F");
			//     CutReco->Branch("fcutshwr_startz",&fcutshwr_startz,"fcutshwr_startz/F");
			//     CutReco->Branch("fcutshwr_bestplane",&fcutshwr_bestplane,"fcutshwr_bestplane/I");
			//     CutReco->Branch("fcutshwr_startdcosx",&fcutshwr_startdcosx,"fcutshwr_startdcosx/F");
			//     CutReco->Branch("fcutshwr_startdcosy",&fcutshwr_startdcosy,"fcutshwr_startdcosy/F");
			//     CutReco->Branch("fcutshwr_startdcosz",&fcutshwr_startdcosz,"fcutshwr_startdcosz/F");
			//     CutReco->Branch("fcutshwr_openangle",&fcutshwr_openangle,"fcutshwr_openangle/F");
			// //    CutReco->Branch("fcutall_shwrEnergy",&fcutall_shwrEnergy,"fcutall_shwrEnergy/F");
			// //    CutReco->Branch("fcutcol_shwrEnergy",&fcutcol_shwrEnergy,"fcutcol_shwrEnergy/F");
			//     CutReco->Branch("fcutshwr_dist",&fcutshwr_dist,"fcutshwr_dist/F");
			// //    CutReco->Branch("fcutshwr_dEdx",&fcutshwr_dEdx,"fcutshwr_dEdx/F");
			// //    CutReco->Branch("fcutshwr_energy",&fcutshwr_energy,"fcutshwr_energy/F");
			// //    CutReco->Branch("fcutshwr_mipenergy",&fcutshwr_mipenergy,"fcutshwr_energy/F");

			//     CutReco->Branch("ftrkcolhits",&ftrkcolhits,"ftrkcolhits/I");
			//     CutReco->Branch("fhits_key",fhits_key,"fhits_key[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_charge",fhits_charge,"fhits_charge[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_wire",fhits_wire,"fhits_wire[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_peakT",fhits_peakT,"fhits_peakT[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_TPC",fhits_TPC,"fhits_TPC[ftrkcolhits]/I");
			//     CutReco->Branch("fhits_chno",fhits_chno,"fhits_chno[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_xpos",fhits_xpos,"fhits_xpos[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_ypos",fhits_ypos,"fhits_ypos[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_zpos",fhits_zpos,"fhits_zpos[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_mult",fhits_mult,"fhits_mult[ftrkcolhits]/I");
			//     CutReco->Branch("fhits_sigptime",fhits_sigptime,"fhits_sigptime[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_sigchrg",fhits_sigchrg,"fhits_sigchrg[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_sigpamp",fhits_sigpamp,"fhits_sigpamp[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_dof",fhits_dof,"fhits_dof[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_gof",fhits_gof,"fhits_gof[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_ptminusRMS",fhits_ptminusRMS,"fhits_ptminusRMS[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_ptplusRMS",fhits_ptplusRMS,"fhits_ptplusRMS[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_cnnMichel",fhits_cnnMichel,"fhits_cnnMichel[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_cnnEM",fhits_cnnEM,"fhits_cnnEM[ftrkcolhits]/F");
			//     CutReco->Branch("fhits_cnnTrack",fhits_cnnTrack,"fhits_cnnTrack[ftrkcolhits]/F");
			// //    fEventTree->Branch("fhits_chrg",&fhits_chrg);

			//     CutReco->Branch("fshwrcolhits",&fshwrcolhits,"fshwrcolhits/I");
			//     CutReco->Branch("fshwrhits_chno",fshwrhits_chno,"fshwrhits_chno[fshwrcolhits]/I");
			//     CutReco->Branch("fshwrhits_peakT",fshwrhits_peakT,"fshwrhits_peakT[fshwrcolhits]/F");
			//     CutReco->Branch("fshwrhits_charge",fshwrhits_charge,"fshwrhits_charge[fshwrcolhits]/F");
			//     CutReco->Branch("fshwrhits_wire",fshwrhits_wire,"fshwrhits_wire[fshwrcolhits]/I");
			//     CutReco->Branch("fshwrhits_plane",fshwrhits_plane,"fshwrhits_plane[fshwrcolhits]/I");
			//     CutReco->Branch("fshwrhits_TPC",fshwrhits_TPC,"shwrhits_TPC[fshwrcolhits]/I");
			//     CutReco->Branch("fshwrhits_xpos",fshwrhits_xpos,"fshwrhits_xpos[fshwrcolhits]/F");
			//     CutReco->Branch("fshwrhits_ypos",fshwrhits_ypos,"fshwrhits_ypos[fshwrcolhits]/F");
			//     CutReco->Branch("fshwrhits_zpos",fshwrhits_zpos,"fshwrhits_zpos[fshwrcolhits]/F");
			//     CutReco->Branch("fshwrhits_mult",fshwrhits_mult,"fshwrhits_mult[fshwrcolhits]/I");
			//     CutReco->Branch("fshwrhits_sigptime",fshwrhits_sigptime,"fshwrhits_sigptime[fshwrcolhits]/F");
			//     CutReco->Branch("fshwrhits_sigchrg",fshwrhits_sigchrg,"fshwrhits_sigchrg[fshwrcolhits]/F");
			//     CutReco->Branch("fshwrhits_sigpamp",fshwrhits_sigpamp,"fshwrhits_sigpamp[fshwrcolhits]/F");
			//     CutReco->Branch("fshwrhits_dof",fshwrhits_dof,"fshwrhits_dof[fshwrcolhits]/F");
			//     CutReco->Branch("fshwrhits_gof",fshwrhits_gof,"fshwrhits_gof[fshwrcolhits]/F");
			//     CutReco->Branch("fshwrhits_ptminusRMS",fshwrhits_ptminusRMS,"fshwrhits_ptminusRMS[fshwrcolhits]/F");
			//     CutReco->Branch("fshwrhits_ptplusRMS",fshwrhits_ptplusRMS,"fshwrhits_ptplusRMS[fshwrcolhits]/F");

			//     CutReco->Branch("fshwrallhits",&fshwrallhits,"fshwrallhits/I");
			//     CutReco->Branch("fshwrallhits_chno",fshwrallhits_chno,"fshwrallhits_chno[fshwrallhits]/I");
			//     CutReco->Branch("fshwrallhits_peakT",fshwrallhits_peakT,"fshwrallhits_peakT[fshwrallhits]/F");
			//     CutReco->Branch("fshwrallhits_charge",fshwrallhits_charge,"fshwrallhits_charge[fshwrallhits]/F");
			//     CutReco->Branch("fshwrallhits_wire",fshwrallhits_wire,"fshwrallhits_wire[fshwrallhits]/I");
			//     CutReco->Branch("fshwrallhits_plane",fshwrallhits_plane,"fshwrallhits_plane[fshwrallhits]/I");
			//     CutReco->Branch("fshwrallhits_TPC",fshwrallhits_TPC,"fshwrallhits_TPC[fshwrallhits]/I");
			//     CutReco->Branch("fshwrallhits_xpos",fshwrallhits_xpos,"fshwrallhits_xpos[fshwrallhits]/F");
			//     CutReco->Branch("fshwrallhits_ypos",fshwrallhits_ypos,"fshwrallhits_ypos[fshwrallhits]/F");
			//     CutReco->Branch("fshwrallhits_zpos",fshwrallhits_zpos,"fshwrallhits_zpos[fshwrallhits]/F");
			//     CutReco->Branch("fshwrallhits_mult",fshwrallhits_mult,"fshwrallhits_mult[fshwrallhits]/I");
			//     CutReco->Branch("fshwrallhits_sigptime",fshwrallhits_sigptime,"fshwrallhits_sigptime[fshwrallhits]/F");
			//     CutReco->Branch("fshwrallhits_sigchrg",fshwrallhits_sigchrg,"fshwrallhits_sigchrg[fshwrallhits]/F");
			//     CutReco->Branch("fshwrallhits_sigpamp",fshwrallhits_sigpamp,"fshwrallhits_sigpamp[fshwrallhits]/F");
			//     CutReco->Branch("fshwrallhits_dof",fshwrallhits_dof,"fshwrallhits_dof[fshwrallhits]/F");
			//     CutReco->Branch("fshwrallhits_gof",fshwrallhits_gof,"fshwrallhits_gof[fshwrallhits]/F");
			//     CutReco->Branch("fshwrallhits_ptminusRMS",fshwrallhits_ptminusRMS,"fshwrallhits_ptminusRMS[fshwrallhits]/F");
			//     CutReco->Branch("fshwrallhits_ptplusRMS",fshwrallhits_ptplusRMS,"fshwrallhits_ptplusRMS[fshwrallhits]/F");

			//     CutReco->Branch("fntrkhits",&fntrkhits,"fntrkhits/I");
			//     CutReco->Branch("fhitsU",&fhitsU,"fhitsU/I");
			//     CutReco->Branch("ftrkdqdxU",ftrkdqdxU,"ftrkdqdxU[fhitsU]/F");
			//     CutReco->Branch("ftrkdedxU",ftrkdedxU,"ftrkdedxU[fhitsU]/F");
			//     CutReco->Branch("ftrkresrangeU",ftrkresrangeU,"ftrkresrangeU[fhitsU]/F");
			//     CutReco->Branch("ftrkhitxU",ftrkhitxU,"ftrkhitxU[fhitsU]/F");
			//     CutReco->Branch("ftrkhityU",ftrkhityU,"ftrkhityU[fhitsU]/F");
			//     CutReco->Branch("ftrkhitzU",ftrkhitzU,"ftrkhitzU[fhitsU]/F");
			//     CutReco->Branch("ftrkpitchU",ftrkpitchU,"ftrkpitchU[fhitsU]/F");
			//     CutReco->Branch("fhitsV",&fhitsV,"fhitsV/I");
			//     CutReco->Branch("ftrkdqdxV",ftrkdqdxV,"ftrkdqdxV[fhitsV]/F");
			//     CutReco->Branch("ftrkdedxV",ftrkdedxV,"ftrkdedxV[fhitsV]/F");
			//     CutReco->Branch("ftrkresrangeV",ftrkresrangeV,"ftrkresrangeV[fhitsV]/F");
			//     CutReco->Branch("ftrkhitxV",ftrkhitxV,"ftrkhitxV[fhitsV]/F");
			//     CutReco->Branch("ftrkhityV",ftrkhityV,"ftrkhityV[fhitsV]/F");
			//     CutReco->Branch("ftrkhitzV",ftrkhitzV,"ftrkhitzV[fhitsV]/F");
			//     CutReco->Branch("ftrkpitchV",ftrkpitchV,"ftrkpitchV[fhitsV]/F");
			//     CutReco->Branch("fhitsY",&fhitsY,"fhitsY/I");
			//     CutReco->Branch("ftrkdqdxY",ftrkdqdxY,"ftrkdqdxY[fhitsY]/F");
			//     CutReco->Branch("ftrkdedxY",ftrkdedxY,"ftrkdedxY[fhitsY]/F");
			//     CutReco->Branch("ftrkresrangeY",ftrkresrangeY,"ftrkresrangeY[fhitsY]/F");
			//     CutReco->Branch("ftrkhitxY",ftrkhitxY,"ftrkhitxY[fhitsY]/F");
			//     CutReco->Branch("ftrkhityY",ftrkhityY,"ftrkhityY[fhitsY]/F");
			//     CutReco->Branch("ftrkhitzY",ftrkhitzY,"ftrkhitzY[fhitsY]/F");
			//     CutReco->Branch("ftrkpitchY",ftrkpitchY,"ftrkpitchY[fhitsY]/F");

			//     CutReco->Branch("fnearhitct",&fnearhitct,"fnearhitct/I");
			//     CutReco->Branch("fnearhits_key",fnearhits_key,"fnearhits_key[fnearhitct]/I");
			//     CutReco->Branch("fnearhits_chno",fnearhits_chno,"fnearhits_chno[fnearhitct]/I");
			//     CutReco->Branch("fnearhits_peakT",fnearhits_peakT,"fnearhits_peakT[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_charge",fnearhits_charge,"fnearhits_charge[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_wire",fnearhits_wire,"fnearhits_wire[fnearhitct]/I");
			//     CutReco->Branch("fnearhits_plane",fnearhits_plane,"fnearhits_plane[fnearhitct]/I");
			//     CutReco->Branch("fnearhits_TPC",fnearhits_TPC,"fnearhits_TPC[fnearhitct]/I");
			//     CutReco->Branch("fnearhits_xpos",fnearhits_xpos,"fnearhits_xpos[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_ypos",fnearhits_ypos,"fnearhits_ypos[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_zpos",fnearhits_zpos,"fnearhits_zpos[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_mult",fnearhits_mult,"fnearhits_mult[fnearhitct]/I");
			//     CutReco->Branch("fnearhits_sigptime",fnearhits_sigptime,"fnearhits_sigptime[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_sigchrg",fnearhits_sigchrg,"fnearhits_sigchrg[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_sigpamp",fnearhits_sigpamp,"fnearhits_sigpamp[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_dof",fnearhits_dof,"fnearhits_dof[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_gof",fnearhits_gof,"fnearhits_gof[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_ptplusRMS",fnearhits_ptplusRMS,"fnearhits_ptplusRMS[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_ptminusRMS",fnearhits_ptminusRMS,"fnearhits_ptminusRMS[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_cnnMichel",fnearhits_cnnMichel,"fnearhits_cnnMichel[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_cnnEM",fnearhits_cnnEM,"fnearhits_cnnEM[fnearhitct]/F");
			//     CutReco->Branch("fnearhits_cnnTrack",fnearhits_cnnTrack,"fnearhits_cnnTrack[fnearhitct]/F");

			//     CutReco->Branch("fmhitcount",&fmhitcount,"fmhitcount/I");
			//     CutReco->Branch("fmhits_key",fmhits_key,"fmhits_key[fmhitcount]/I");
			//     CutReco->Branch("fmhits_chno",fmhits_chno,"fmhits_chno[fmhitcount]/I");
			//     CutReco->Branch("fmhits_peakT",fmhits_peakT,"fmhits_peakT[fmhitcount]/F");
			//     CutReco->Branch("fmhits_charge",fmhits_charge,"fmhits_charge[fmhitcount]/F");
			//     CutReco->Branch("fmhits_wire",fmhits_wire,"fmhits_wire[fmhitcount]/I");
			//     CutReco->Branch("fmhits_plane",fmhits_plane,"fmhits_plane[fmhitcount]/I");
			//     CutReco->Branch("fmhits_TPC",fmhits_TPC,"fmhits_TPC[fmhitcount]/I");
			//     CutReco->Branch("fmhits_xpos",fmhits_xpos,"fmhits_xpos[fmhitcount]/F");
			//     CutReco->Branch("fmhits_ypos",fmhits_ypos,"fmhits_ypos[fmhitcount]/F");
			//     CutReco->Branch("fmhits_zpos",fmhits_zpos,"fmhits_zpos[fmhitcount]/F");
			//     CutReco->Branch("fmhits_mult",fmhits_mult,"fmhits_mult[fmhitcount]/I");
			//     CutReco->Branch("fmhits_sigptime",fmhits_sigptime,"fmhits_sigptime[fmhitcount]/F");
			//     CutReco->Branch("fmhits_sigchrg",fmhits_sigchrg,"fmhits_sigchrg[fmhitcount]/F");
			//     CutReco->Branch("fmhits_sigpamp",fmhits_sigpamp,"fmhits_sigpamp[fmhitcount]/F");
			//     CutReco->Branch("fmhits_dof",fmhits_dof,"fmhits_dof[fmhitcount]/F");
			//     CutReco->Branch("fmhits_gof",fmhits_gof,"fmhits_gof[fmhitcount]/F");
			//     CutReco->Branch("fmhits_angledeg",fmhits_angledeg,"fmhits_angledeg[fmhitcount]/F");
			//     CutReco->Branch("fmhits_maghitveccostheta",fmhits_maghitveccostheta,"fmhits_maghitveccostheta[fmhitcount]/F");
			//     CutReco->Branch("fmhits_distance",fmhits_distance,"fmhits_distance[fmhitcount]/F");
			//     CutReco->Branch("fmhits_longtrk",fmhits_longtrk,"fmhits_longtrk[fmhitcount]/I");
			//     CutReco->Branch("fmhits_sametrk",fmhits_sametrk,"fmhits_sametrk[fmhitcount]/I");
			//     CutReco->Branch("fmhits_corrhit",fmhits_corrhit,"fmhits_corrhit[fmhitcount]/I");
			//     CutReco->Branch("fmhits_ptminusRMS",fmhits_ptminusRMS,"fmhits_ptminusRMS[fmhitcount]/F");
			//     CutReco->Branch("fmhits_ptplusRMS",fmhits_ptplusRMS,"fmhits_ptplusRMS[fmhitcount]/F");
			//     CutReco->Branch("fmhits_cnnMichel",fmhits_cnnMichel,"fmhits_cnnMichel[fmhitcount]/F");
			//     CutReco->Branch("fmhits_cnnEM",fmhits_cnnEM,"fmhits_cnnEM[fmhitcount]/F");
			//     CutReco->Branch("fmhits_cnnTrack",fmhits_cnnTrack,"fmhits_cnnTrack[fmhitcount]/F");
				
			//     CutReco->Branch("ftrueparhitallcount",&ftrueparhitallcount,"ftrueparhitallcount/I");
			//     CutReco->Branch("ftrueparhitsall_key",ftrueparhitsall_key,"ftrueparhitsall_key[ftrueparhitallcount]/I");
			//     CutReco->Branch("ftrueparhitsall_chno",ftrueparhitsall_chno,"ftrueparhitsall_chno[ftrueparhitallcount]/I");
			//     CutReco->Branch("ftrueparhitsall_peakT",ftrueparhitsall_peakT,"ftrueparhitsall_peakT[ftrueparhitallcount]/F");
			//     CutReco->Branch("ftrueparhitsall_charge",ftrueparhitsall_charge,"ftrueparhitsall_charge[ftrueparhitallcount]/F");
			//     CutReco->Branch("ftrueparhitsall_wire",ftrueparhitsall_wire,"ftrueparhitsall_wire[ftrueparhitallcount]/I");
			//     CutReco->Branch("ftrueparhitsall_plane",ftrueparhitsall_plane,"ftrueparhitsall_plane[ftrueparhitallcount]/I");
			//     CutReco->Branch("ftrueparhitsall_TPC",ftrueparhitsall_TPC,"ftrueparhitsall_TPC[ftrueparhitallcount]/I");
			//     CutReco->Branch("ftrueparhitsall_xpos",ftrueparhitsall_xpos,"ftrueparhitsall_xpos[ftrueparhitallcount]/F");
			//     CutReco->Branch("ftrueparhitsall_ypos",ftrueparhitsall_ypos,"ftrueparhitsall_ypos[ftrueparhitallcount]/F");
			//     CutReco->Branch("ftrueparhitsall_zpos",ftrueparhitsall_zpos,"ftrueparhitsall_zpos[ftrueparhitallcount]/F");
			//     CutReco->Branch("ftrueparhitsall_mult",ftrueparhitsall_mult,"ftrueparhitsall_mult[ftrueparhitallcount]/I");
			//     CutReco->Branch("ftrueparhitsall_sigptime",ftrueparhitsall_sigptime,"ftrueparhitsall_sigptime[ftrueparhitallcount]/F");
			//     CutReco->Branch("ftrueparhitsall_sigchrg",ftrueparhitsall_sigchrg,"ftrueparhitsall_sigchrg[ftrueparhitallcount]/F");
			//     CutReco->Branch("ftrueparhitsall_sigpamp",ftrueparhitsall_sigpamp,"ftrueparhitsall_sigpamp[ftrueparhitallcount]/F");
			//     CutReco->Branch("ftrueparhitsall_dof",ftrueparhitsall_dof,"ftrueparhitsall_dof[ftrueparhitallcount]/F");
			//     CutReco->Branch("ftrueparhitsall_gof",ftrueparhitsall_gof,"ftrueparhitsall_gof[ftrueparhitallcount]/F");
			//     CutReco->Branch("ftrueparhitsall_ptminusRMS",ftrueparhitsall_ptminusRMS,"ftrueparhitsall_ptminusRMS[ftrueparhitallcount]/F");
			//     CutReco->Branch("ftrueparhitsall_ptplusRMS",ftrueparhitsall_ptplusRMS,"ftrueparhitsall_ptplusRMS[ftrueparhitallcount]/F");

			//     CutReco->Branch("ftrueparhitcolcount",&ftrueparhitcolcount,"ftrueparhitcolcount/I");
			//     CutReco->Branch("ftrueMiEFrac",ftrueMiEFrac,"ftrueMiEFrac[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_key",ftrueparhitscol_key,"ftrueparhitscol_key[ftrueparhitcolcount]/I");
			//     CutReco->Branch("ftrueparhitscol_chno",ftrueparhitscol_chno,"ftrueparhitscol_chno[ftrueparhitcolcount]/I");
			//     CutReco->Branch("ftrueparhitscol_peakT",ftrueparhitscol_peakT,"ftrueparhitscol_peakT[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_charge",ftrueparhitscol_charge,"ftrueparhitscol_charge[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_wire",ftrueparhitscol_wire,"ftrueparhitscol_wire[ftrueparhitcolcount]/I");
			//     CutReco->Branch("ftrueparhitscol_plane",ftrueparhitscol_plane,"ftrueparhitscol_plane[ftrueparhitcolcount]/I");
			//     CutReco->Branch("ftrueparhitscol_TPC",ftrueparhitscol_TPC,"ftrueparhitscol_TPC[ftrueparhitcolcount]/I");
			//     CutReco->Branch("ftrueparhitscol_xpos",ftrueparhitscol_xpos,"ftrueparhitscol_xpos[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_ypos",ftrueparhitscol_ypos,"ftrueparhitscol_ypos[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_zpos",ftrueparhitscol_zpos,"ftrueparhitscol_zpos[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_angledeg",ftrueparhitscol_angledeg,"ftrueparhitscol_angledeg[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_mult",ftrueparhitscol_mult,"ftrueparhitscol_mult[ftrueparhitcolcount]/I");
			//     CutReco->Branch("ftrueparhitscol_sigptime",ftrueparhitscol_sigptime,"ftrueparhitscol_sigptime[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_sigchrg",ftrueparhitscol_sigchrg,"ftrueparhitscol_sigchrg[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_sigpamp",ftrueparhitscol_sigpamp,"ftrueparhitscol_sigpamp[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_dof",ftrueparhitscol_dof,"ftrueparhitscol_dof[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_gof",ftrueparhitscol_gof,"ftrueparhitscol_gof[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_ptminusRMS",ftrueparhitscol_ptminusRMS,"ftrueparhitscol_ptminusRMS[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_ptplusRMS",ftrueparhitscol_ptplusRMS,"ftrueparhitscol_ptplusRMS[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_maghitveccostheta",ftrueparhitscol_maghitveccostheta,"ftrueparhitscol_maghitveccostheta[ftrueparhitcolcount]/F");
			//     CutReco->Branch("ftrueparhitscol_distance",ftrueparhitscol_distance,"ftrueparhitscol_distance[ftrueparhitcolcount]/F");

			// //    CutReco->Branch("fshwrhits_charge",fshwrhits_charge,"fshwrhits_charge[fcut_mu][3][]/I");

			//     CutReco->Branch("fmccut_trkid",&fmccut_trkid,"fmccut_trkid/I");
			//     CutReco->Branch("fmccut_vx",&fmccut_vx,"fmccut_vx/F");
			//     CutReco->Branch("fmccut_vy",&fmccut_vy,"fmccut_vy/F");
			//     CutReco->Branch("fmccut_vz",&fmccut_vz,"fmccut_vz/F");
			//     CutReco->Branch("fmccut_t",&fmccut_t,"fmccut_t/F");
			//     CutReco->Branch("fmccut_endx",&fmccut_endx,"fmccut_endx/F");
			//     CutReco->Branch("fmccut_endy",&fmccut_endy,"fmccut_endy/F");
			//     CutReco->Branch("fmccut_endz",&fmccut_endz,"fmccut_endz/F");
			//     CutReco->Branch("fmccut_endt",&fmccut_endt,"fmccut_endt/F");
			//     CutReco->Branch("fmccut_px",&fmccut_px,"fmccut_px/F");
			//     CutReco->Branch("fmccut_py",&fmccut_py,"fmccut_py/F");
			//     CutReco->Branch("fmccut_pz",&fmccut_pz,"fmccut_pz/F");
			//     CutReco->Branch("fmccut_momentum",&fmccut_momentum,"fmccut_momentum/F");
			//     CutReco->Branch("fmccut_energy",&fmccut_energy,"fmccut_energy/F");
			//     CutReco->Branch("fmccut_endpx",&fmccut_endpx,"fmccut_endpx/F");
			//     CutReco->Branch("fmccut_endpy",&fmccut_endpy,"fmccut_endpy/F");
			//     CutReco->Branch("fmccut_endpz",&fmccut_endpz,"fmccut_endpz/F");
			//     CutReco->Branch("fmccut_endenergy",&fmccut_endenergy,"fmccut_endenergy/F");
			//     CutReco->Branch("fmccut_pathlen",&fmccut_pathlen,"fmccut_pathlen/F");
			//     CutReco->Branch("fmccut_length",&fmccut_length,"fmccut_length/F");
			//     CutReco->Branch("fmccut_vxdrifted",&fmccut_vxdrifted,"fmccut_vxdrifted/F");
			//     CutReco->Branch("fmccut_vydrifted",&fmccut_vydrifted,"fmccut_vydrifted/F");
			//     CutReco->Branch("fmccut_vzdrifted",&fmccut_vzdrifted,"fmccut_vzdrifted/F");
			//     CutReco->Branch("fmccut_tdrifted",&fmccut_tdrifted,"fmccut_tdrifted/F");
			//     CutReco->Branch("fmccut_endxdrifted",&fmccut_endxdrifted,"fmccut_endxdrifted/F");
			//     CutReco->Branch("fmccut_endydrifted",&fmccut_endydrifted,"fmccut_endydrifted/F");
			//     CutReco->Branch("fmccut_endzdrifted",&fmccut_endzdrifted,"fmccut_endzdrifted/F");
			//     CutReco->Branch("fmccut_endtdrifted",&fmccut_endtdrifted,"fmccut_endtdrifted/F");
			//     CutReco->Branch("fmccut_pxdrifted",&fmccut_pxdrifted,"fmccut_pxdrifted/F");
			//     CutReco->Branch("fmccut_pydrifted",&fmccut_pydrifted,"fmccut_pydrifted/F");
			//     CutReco->Branch("fmccut_pzdrifted",&fmccut_pzdrifted,"fmccut_pzdrifted/F");
			//     CutReco->Branch("fmccut_momentumdrifted",&fmccut_momentumdrifted,"fmccut_momentumdrifted/F");
			//     CutReco->Branch("fmccut_energydrifted",&fmccut_energydrifted,"fmccut_energydrifted/F");
			//     CutReco->Branch("fmccut_endpxdrifted",&fmccut_endpxdrifted,"fmccut_endpxdrifted/F");
			//     CutReco->Branch("fmccut_endpydrifted",&fmccut_endpydrifted,"fmccut_endpydrifted/F");
			//     CutReco->Branch("fmccut_endpzdrifted",&fmccut_endpzdrifted,"fmccut_endpzdrifted/F");
			//     CutReco->Branch("fmccut_endenergydrifted",&fmccut_endenergydrifted,"fmccut_endenergydrifted/F");
			//     CutReco->Branch("fmccut_pathlendrifted",&fmccut_pathlendrifted,"fmccut_pathlendrifted/F");
			//     CutReco->Branch("fmccut_endprocess",&fmccut_endprocess,"fmccut_endprocess/I");
			//     CutReco->Branch("fmccut_theta",&fmccut_theta,"fmccut_theta/F");
			//     CutReco->Branch("fmccut_phi",&fmccut_phi,"fmccut_phi/F");
			//     CutReco->Branch("fmccut_pdg",&fmccut_pdg,"fmccut_pdg/I");
			//     CutReco->Branch("fmccut_status_code",&fmccut_status_code,"fmccut_status_code/I");
			//     CutReco->Branch("fmccut_mass",&fmccut_mass,"fmccut_mass/F");
			//     CutReco->Branch("fmccut_ND",&fmccut_ND,"fmccut_ND/I");
			//     CutReco->Branch("fmccut_mother",&fmccut_mother,"fmccut_mother/I");
			//     CutReco->Branch("fmccut_origin",&fmccut_origin,"fmccut_origin/I");
			//     CutReco->Branch("fmccut_process",&fmccut_process,"fmccut_process/I");
			//     CutReco->Branch("fmccut_rescatter",&fmccut_rescatter,"fmccut_rescatter/I");

			// //    CutReco->Branch("fhasElect",&fhasElect,"fhasElect/I"); 
			// //    CutReco->Branch("fcutdau_mu",&fcutdau_mu,"fcutdau_mu/I");
			//     CutReco->Branch("fmcd_trkid",&fmcd_trkid,"fmcd_trkid/I");
			//     CutReco->Branch("fmcd_vx",&fmcd_vx,"fmcd_vx/F");
			//     CutReco->Branch("fmcd_vy",&fmcd_vy,"fmcd_vy/F");
			//     CutReco->Branch("fmcd_vz",&fmcd_vz,"fmcd_vz/F");
			//     CutReco->Branch("fmcd_t",&fmcd_t,"fmcd_t/F");
			//     CutReco->Branch("fmcd_endx",&fmcd_endx,"fmcd_endx/F");
			//     CutReco->Branch("fmcd_endy",&fmcd_endy,"fmcd_endy/F");
			//     CutReco->Branch("fmcd_endz",&fmcd_endz,"fmcd_endz/F");
			//     CutReco->Branch("fmcd_endt",&fmcd_endt,"fmcd_endt/F");
			//     CutReco->Branch("fmcd_px",&fmcd_px,"fmcd_px/F");
			//     CutReco->Branch("fmcd_py",&fmcd_py,"fmcd_py/F");
			//     CutReco->Branch("fmcd_pz",&fmcd_pz,"fmcd_pz/F");
			//     CutReco->Branch("fmcd_momentum",&fmcd_momentum,"fmcd_momentum/F");
			//     CutReco->Branch("fmcd_energy",&fmcd_energy,"fmcd_energy/F");
			//     CutReco->Branch("fmcd_truecuthitsE",&fmcd_truecuthitsE,"fmcd_truecuthitsE/F");
			//     CutReco->Branch("ftrueEdepo",&ftrueEdepo,"ftrueEdepo/F");
			//     CutReco->Branch("fmcd_endpx",&fmcd_endpx,"fmcd_endpx/F");
			//     CutReco->Branch("fmcd_endpy",&fmcd_endpy,"fmcd_endpy/F");
			//     CutReco->Branch("fmcd_endpz",&fmcd_endpz,"fmcd_endpz/F");
			//     CutReco->Branch("fmcd_endenergy",&fmcd_endenergy,"fmcd_endenergy/F");
			//     CutReco->Branch("fmcd_pathlen",&fmcd_pathlen,"fmcd_pathlen/F");
			//     CutReco->Branch("fmcd_vxdrifted",&fmcd_vxdrifted,"fmcd_vxdrifted/F");
			//     CutReco->Branch("fmcd_vydrifted",&fmcd_vydrifted,"fmcd_vydrifted/F");
			//     CutReco->Branch("fmcd_vzdrifted",&fmcd_vzdrifted,"fmcd_vzdrifted/F");
			//     CutReco->Branch("fmcd_tdrifted",&fmcd_tdrifted,"fmcd_tdrifted/F");
			//     CutReco->Branch("fmcd_endxdrifted",&fmcd_endxdrifted,"fmcd_endxdrifted/F");
			//     CutReco->Branch("fmcd_endydrifted",&fmcd_endydrifted,"fmcd_endydrifted/F");
			//     CutReco->Branch("fmcd_endzdrifted",&fmcd_endzdrifted,"fmcd_endzdrifted/F");
			//     CutReco->Branch("fmcd_endtdrifted",&fmcd_endtdrifted,"fmcd_endtdrifted/F");
			//     CutReco->Branch("fmcd_pxdrifted",&fmcd_pxdrifted,"fmcd_pxdrifted/F");
			//     CutReco->Branch("fmcd_pydrifted",&fmcd_pydrifted,"fmcd_pydrifted/F");
			//     CutReco->Branch("fmcd_pzdrifted",&fmcd_pzdrifted,"fmcd_pzdrifted/F");
			//     CutReco->Branch("fmcd_momentumdrifted",&fmcd_momentumdrifted,"fmcd_momentumdrifted/F");
			//     CutReco->Branch("fmcd_energydrifted",&fmcd_energydrifted,"fmcd_energydrifted/F");
			//     CutReco->Branch("fmcd_endpxdrifted",&fmcd_endpxdrifted,"fmcd_endpxdrifted/F");
			//     CutReco->Branch("fmcd_endpydrifted",&fmcd_endpydrifted,"fmcd_endpydrifted/F");
			//     CutReco->Branch("fmcd_endpzdrifted",&fmcd_endpzdrifted,"fmcd_endpzdrifted/F");
			//     CutReco->Branch("fmcd_endenergydrifted",&fmcd_endenergydrifted,"fmcd_endenergydrifted/F");
			//     CutReco->Branch("fmcd_pathlendrifted",&fmcd_pathlendrifted,"fmcd_pathlendrifted/F");
			//     CutReco->Branch("fmcd_endprocess",&fmcd_endprocess,"fmcd_endprocess/I");
			//     CutReco->Branch("fmcd_theta",&fmcd_theta,"fmcd_theta/F");
			//     CutReco->Branch("fmcd_phi",&fmcd_phi,"fmcd_phi/F");
			//     CutReco->Branch("fmcd_pdg",&fmcd_pdg,"fmcd_pdg/I");
			//     CutReco->Branch("fmcd_status_code",&fmcd_status_code,"fmcd_status_code/I");
			//     CutReco->Branch("fmcd_mass",&fmcd_mass,"fmcd_mass/F");
			//     CutReco->Branch("fmcd_ND",&fmcd_ND,"fmcd_ND/I");
			//     CutReco->Branch("fmcd_mother",&fmcd_mother,"fmcd_mother/I");
			//     CutReco->Branch("fmcd_origin",&fmcd_origin,"fmcd_origin/I");
			//     CutReco->Branch("fmcd_process",&fmcd_process,"fmcd_process/I");
			//     CutReco->Branch("fmcd_rescatter",&fmcd_rescatter,"fmcd_rescatter/I");

			//     CutReco->Branch("fmchits_trkid",&fmchits_trkid,"fmchits_trkid/I");
			//     CutReco->Branch("fmchits_vx",&fmchits_vx,"fmchits_vx/F");
			//     CutReco->Branch("fmchits_vy",&fmchits_vy,"fmchits_vy/F");
			//     CutReco->Branch("fmchits_vz",&fmchits_vz,"fmchits_vz/F");
			//     CutReco->Branch("fmchits_t",&fmchits_t,"fmchits_t/F");
			//     CutReco->Branch("fmchits_endx",&fmchits_endx,"fmchits_endx/F");
			//     CutReco->Branch("fmchits_endy",&fmchits_endy,"fmchits_endy/F");
			//     CutReco->Branch("fmchits_endz",&fmchits_endz,"fmchits_endz/F");
			//     CutReco->Branch("fmchits_endt",&fmchits_endt,"fmchits_endt/F");
			//     CutReco->Branch("fmchits_px",&fmchits_px,"fmchits_px/F");
			//     CutReco->Branch("fmchits_py",&fmchits_py,"fmchits_py/F");
			//     CutReco->Branch("fmchits_pz",&fmchits_pz,"fmchits_pz/F");
			//     CutReco->Branch("fmchits_momentum",&fmchits_momentum,"fmchits_momentum/F");
			//     CutReco->Branch("fmchits_energy",&fmchits_energy,"fmchits_energy/F");
			//     CutReco->Branch("fmchits_endpx",&fmchits_endpx,"fmchits_endpx/F");
			//     CutReco->Branch("fmchits_endpy",&fmchits_endpy,"fmchits_endpy/F");
			//     CutReco->Branch("fmchits_endpz",&fmchits_endpz,"fmchits_endpz/F");
			//     CutReco->Branch("fmchits_endenergy",&fmchits_endenergy,"fmchits_endenergy/F");
			//     CutReco->Branch("fmchits_pathlen",&fmchits_pathlen,"fmchits_pathlen/F");
			//     CutReco->Branch("fmchits_vxdrifted",&fmchits_vxdrifted,"fmchits_vxdrifted/F");
			//     CutReco->Branch("fmchits_vydrifted",&fmchits_vydrifted,"fmchits_vydrifted/F");
			//     CutReco->Branch("fmchits_vzdrifted",&fmchits_vzdrifted,"fmchits_vzdrifted/F");
			//     CutReco->Branch("fmchits_tdrifted",&fmchits_tdrifted,"fmchits_tdrifted/F");
			//     CutReco->Branch("fmchits_endxdrifted",&fmchits_endxdrifted,"fmchits_endxdrifted/F");
			//     CutReco->Branch("fmchits_endydrifted",&fmchits_endydrifted,"fmchits_endydrifted/F");
			//     CutReco->Branch("fmchits_endzdrifted",&fmchits_endzdrifted,"fmchits_endzdrifted/F");
			//     CutReco->Branch("fmchits_endtdrifted",&fmchits_endtdrifted,"fmchits_endtdrifted/F");
			//     CutReco->Branch("fmchits_pxdrifted",&fmchits_pxdrifted,"fmchits_pxdrifted/F");
			//     CutReco->Branch("fmchits_pydrifted",&fmchits_pydrifted,"fmchits_pydrifted/F");
			//     CutReco->Branch("fmchits_pzdrifted",&fmchits_pzdrifted,"fmchits_pzdrifted/F");
			//     CutReco->Branch("fmchits_momentumdrifted",&fmchits_momentumdrifted,"fmchits_momentumdrifted/F");
			//     CutReco->Branch("fmchits_energydrifted",&fmchits_energydrifted,"fmchits_energydrifted/F");
			//     CutReco->Branch("fmchits_endpxdrifted",&fmchits_endpxdrifted,"fmchits_endpxdrifted/F");
			//     CutReco->Branch("fmchits_endpydrifted",&fmchits_endpydrifted,"fmchits_endpydrifted/F");
			//     CutReco->Branch("fmchits_endpzdrifted",&fmchits_endpzdrifted,"fmchits_endpzdrifted/F");
			//     CutReco->Branch("fmchits_endenergydrifted",&fmchits_endenergydrifted,"fmchits_endenergydrifted/F");
			//     CutReco->Branch("fmchits_pathlendrifted",&fmchits_pathlendrifted,"fmchits_pathlendrifted/F");
			//     CutReco->Branch("fmchits_endprocess",&fmchits_endprocess,"fmchits_endprocess/I");
			//     CutReco->Branch("fmchits_theta",&fmchits_theta,"fmchits_theta/F");
			//     CutReco->Branch("fmchits_phi",&fmchits_phi,"fmchits_phi/F");
			//     CutReco->Branch("fmchits_pdg",&fmchits_pdg,"fmchits_pdg/I");
			//     CutReco->Branch("fmchits_status_code",&fmchits_status_code,"fmchits_status_code/I");
			//     CutReco->Branch("fmchits_mass",&fmchits_mass,"fmchits_mass/F");
			//     CutReco->Branch("fmchits_ND",&fmchits_ND,"fmchits_ND/I");
			//     CutReco->Branch("fmchits_mother",&fmchits_mother,"fmchits_mother/I");
			//     CutReco->Branch("fmchits_origin",&fmchits_origin,"fmchits_origin/I");
			//     CutReco->Branch("fmchits_process",&fmchits_process,"fmchits_process/I");
			//     CutReco->Branch("fmchits_rescatter",&fmchits_rescatter,"fmchits_rescatter/I");
				
			//     CutReco->Branch("fmcconehits_trkid",&fmcconehits_trkid,"fmcconehits_trkid/I");
			//     CutReco->Branch("fmcconehits_vx",&fmcconehits_vx,"fmcconehits_vx/F");
			//     CutReco->Branch("fmcconehits_vy",&fmcconehits_vy,"fmcconehits_vy/F");
			//     CutReco->Branch("fmcconehits_vz",&fmcconehits_vz,"fmcconehits_vz/F");
			//     CutReco->Branch("fmcconehits_t",&fmcconehits_t,"fmcconehits_t/F");
			//     CutReco->Branch("fmcconehits_endx",&fmcconehits_endx,"fmcconehits_endx/F");
			//     CutReco->Branch("fmcconehits_endy",&fmcconehits_endy,"fmcconehits_endy/F");
			//     CutReco->Branch("fmcconehits_endz",&fmcconehits_endz,"fmcconehits_endz/F");
			//     CutReco->Branch("fmcconehits_endt",&fmcconehits_endt,"fmcconehits_endt/F");
			//     CutReco->Branch("fmcconehits_px",&fmcconehits_px,"fmcconehits_px/F");
			//     CutReco->Branch("fmcconehits_py",&fmcconehits_py,"fmcconehits_py/F");
			//     CutReco->Branch("fmcconehits_pz",&fmcconehits_pz,"fmcconehits_pz/F");
			//     CutReco->Branch("fmcconehits_momentum",&fmcconehits_momentum,"fmcconehits_momentum/F");
			//     CutReco->Branch("fmcconehits_energy",&fmcconehits_energy,"fmcconehits_energy/F");
			//     CutReco->Branch("fmcconehits_endpx",&fmcconehits_endpx,"fmcconehits_endpx/F");
			//     CutReco->Branch("fmcconehits_endpy",&fmcconehits_endpy,"fmcconehits_endpy/F");
			//     CutReco->Branch("fmcconehits_endpz",&fmcconehits_endpz,"fmcconehits_endpz/F");
			//     CutReco->Branch("fmcconehits_endenergy",&fmcconehits_endenergy,"fmcconehits_endenergy/F");
			//     CutReco->Branch("fmcconehits_pathlen",&fmcconehits_pathlen,"fmcconehits_pathlen/F");
			//     CutReco->Branch("fmcconehits_vxdrifted",&fmcconehits_vxdrifted,"fmcconehits_vxdrifted/F");
			//     CutReco->Branch("fmcconehits_vydrifted",&fmcconehits_vydrifted,"fmcconehits_vydrifted/F");
			//     CutReco->Branch("fmcconehits_vzdrifted",&fmcconehits_vzdrifted,"fmcconehits_vzdrifted/F");
			//     CutReco->Branch("fmcconehits_tdrifted",&fmcconehits_tdrifted,"fmcconehits_tdrifted/F");
			//     CutReco->Branch("fmcconehits_endxdrifted",&fmcconehits_endxdrifted,"fmcconehits_endxdrifted/F");
			//     CutReco->Branch("fmcconehits_endydrifted",&fmcconehits_endydrifted,"fmcconehits_endydrifted/F");
			//     CutReco->Branch("fmcconehits_endzdrifted",&fmcconehits_endzdrifted,"fmcconehits_endzdrifted/F");
			//     CutReco->Branch("fmcconehits_endtdrifted",&fmcconehits_endtdrifted,"fmcconehits_endtdrifted/F");
			//     CutReco->Branch("fmcconehits_pxdrifted",&fmcconehits_pxdrifted,"fmcconehits_pxdrifted/F");
			//     CutReco->Branch("fmcconehits_pydrifted",&fmcconehits_pydrifted,"fmcconehits_pydrifted/F");
			//     CutReco->Branch("fmcconehits_pzdrifted",&fmcconehits_pzdrifted,"fmcconehits_pzdrifted/F");
			//     CutReco->Branch("fmcconehits_momentumdrifted",&fmcconehits_momentumdrifted,"fmcconehits_momentumdrifted/F");
			//     CutReco->Branch("fmcconehits_energydrifted",&fmcconehits_energydrifted,"fmcconehits_energydrifted/F");
			//     CutReco->Branch("fmcconehits_endpxdrifted",&fmcconehits_endpxdrifted,"fmcconehits_endpxdrifted/F");
			//     CutReco->Branch("fmcconehits_endpydrifted",&fmcconehits_endpydrifted,"fmcconehits_endpydrifted/F");
			//     CutReco->Branch("fmcconehits_endpzdrifted",&fmcconehits_endpzdrifted,"fmcconehits_endpzdrifted/F");
			//     CutReco->Branch("fmcconehits_endenergydrifted",&fmcconehits_endenergydrifted,"fmcconehits_endenergydrifted/F");
			//     CutReco->Branch("fmcconehits_pathlendrifted",&fmcconehits_pathlendrifted,"fmcconehits_pathlendrifted/F");
			//     CutReco->Branch("fmcconehits_endprocess",&fmcconehits_endprocess,"fmcconehits_endprocess/I");
			//     CutReco->Branch("fmcconehits_theta",&fmcconehits_theta,"fmcconehits_theta/F");
			//     CutReco->Branch("fmcconehits_phi",&fmcconehits_phi,"fmcconehits_phi/F");
			//     CutReco->Branch("fmcconehits_pdg",&fmcconehits_pdg,"fmcconehits_pdg/I");
			//     CutReco->Branch("fmcconehits_status_code",&fmcconehits_status_code,"fmcconehits_status_code/I");
			//     CutReco->Branch("fmcconehits_mass",&fmcconehits_mass,"fmcconehits_mass/F");
			//     CutReco->Branch("fmcconehits_ND",&fmcconehits_ND,"fmcconehits_ND/I");
			//     CutReco->Branch("fmcconehits_mother",&fmcconehits_mother,"fmcconehits_mother/I");
			//     CutReco->Branch("fmcconehits_origin",&fmcconehits_origin,"fmcconehits_origin/I");
			//     CutReco->Branch("fmcconehits_process",&fmcconehits_process,"fmcconehits_process/I");
			//     CutReco->Branch("fmcconehits_rescatter",&fmcconehits_rescatter,"fmcconehits_rescatter/I");

			//     CutReco->Branch("fmcshwr_trkid",&fmcshwr_trkid,"fmcshwr_trkid/I");
			//     CutReco->Branch("fmcshwr_vx",&fmcshwr_vx,"fmcshwr_vx/F");
			//     CutReco->Branch("fmcshwr_vy",&fmcshwr_vy,"fmcshwr_vy/F");
			//     CutReco->Branch("fmcshwr_vz",&fmcshwr_vz,"fmcshwr_vz/F");
			//     CutReco->Branch("fmcshwr_t",&fmcshwr_t,"fmcshwr_t/F");
			//     CutReco->Branch("fmcshwr_endx",&fmcshwr_endx,"fmcshwr_endx/F");
			//     CutReco->Branch("fmcshwr_endy",&fmcshwr_endy,"fmcshwr_endy/F");
			//     CutReco->Branch("fmcshwr_endz",&fmcshwr_endz,"fmcshwr_endz/F");
			//     CutReco->Branch("fmcshwr_endt",&fmcshwr_endt,"fmcshwr_endt/F");
			//     CutReco->Branch("fmcshwr_px",&fmcshwr_px,"fmcshwr_px/F");
			//     CutReco->Branch("fmcshwr_py",&fmcshwr_py,"fmcshwr_py/F");
			//     CutReco->Branch("fmcshwr_pz",&fmcshwr_pz,"fmcshwr_pz/F");
			//     CutReco->Branch("fmcshwr_momentum",&fmcshwr_momentum,"fmcshwr_momentum/F");
			//     CutReco->Branch("fmcshwr_energy",&fmcshwr_energy,"fmcshwr_energy/F");
			//     CutReco->Branch("fmcshwr_endpx",&fmcshwr_endpx,"fmcshwr_endpx/F");
			//     CutReco->Branch("fmcshwr_endpy",&fmcshwr_endpy,"fmcshwr_endpy/F");
			//     CutReco->Branch("fmcshwr_endpz",&fmcshwr_endpz,"fmcshwr_endpz/F");
			//     CutReco->Branch("fmcshwr_endenergy",&fmcshwr_endenergy,"fmcshwr_endenergy/F");
			//     CutReco->Branch("fmcshwr_pathlen",&fmcshwr_pathlen,"fmcshwr_pathlen/F");
			//     CutReco->Branch("fmcshwr_vxdrifted",&fmcshwr_vxdrifted,"fmcshwr_vxdrifted/F");
			//     CutReco->Branch("fmcshwr_vydrifted",&fmcshwr_vydrifted,"fmcshwr_vydrifted/F");
			//     CutReco->Branch("fmcshwr_vzdrifted",&fmcshwr_vzdrifted,"fmcshwr_vzdrifted/F");
			//     CutReco->Branch("fmcshwr_tdrifted",&fmcshwr_tdrifted,"fmcshwr_tdrifted/F");
			//     CutReco->Branch("fmcshwr_endxdrifted",&fmcshwr_endxdrifted,"fmcshwr_endxdrifted/F");
			//     CutReco->Branch("fmcshwr_endydrifted",&fmcshwr_endydrifted,"fmcshwr_endydrifted/F");
			//     CutReco->Branch("fmcshwr_endzdrifted",&fmcshwr_endzdrifted,"fmcshwr_endzdrifted/F");
			//     CutReco->Branch("fmcshwr_endtdrifted",&fmcshwr_endtdrifted,"fmcshwr_endtdrifted/F");
			//     CutReco->Branch("fmcshwr_pxdrifted",&fmcshwr_pxdrifted,"fmcshwr_pxdrifted/F");
			//     CutReco->Branch("fmcshwr_pydrifted",&fmcshwr_pydrifted,"fmcshwr_pydrifted/F");
			//     CutReco->Branch("fmcshwr_pzdrifted",&fmcshwr_pzdrifted,"fmcshwr_pzdrifted/F");
			//     CutReco->Branch("fmcshwr_momentumdrifted",&fmcshwr_momentumdrifted,"fmcshwr_momentumdrifted/F");
			//     CutReco->Branch("fmcshwr_energydrifted",&fmcshwr_energydrifted,"fmcshwr_energydrifted/F");
			//     CutReco->Branch("fmcshwr_endpxdrifted",&fmcshwr_endpxdrifted,"fmcshwr_endpxdrifted/F");
			//     CutReco->Branch("fmcshwr_endpydrifted",&fmcshwr_endpydrifted,"fmcshwr_endpydrifted/F");
			//     CutReco->Branch("fmcshwr_endpzdrifted",&fmcshwr_endpzdrifted,"fmcshwr_endpzdrifted/F");
			//     CutReco->Branch("fmcshwr_endenergydrifted",&fmcshwr_endenergydrifted,"fmcshwr_endenergydrifted/F");
			//     CutReco->Branch("fmcshwr_pathlendrifted",&fmcshwr_pathlendrifted,"fmcshwr_pathlendrifted/F");
			//     CutReco->Branch("fmcshwr_endprocess",&fmcshwr_endprocess,"fmcshwr_endprocess/I");
			//     CutReco->Branch("fmcshwr_theta",&fmcshwr_theta,"fmcshwr_theta/F");
			//     CutReco->Branch("fmcshwr_phi",&fmcshwr_phi,"fmcshwr_phi/F");
			//     CutReco->Branch("fmcshwr_pdg",&fmcshwr_pdg,"fmcshwr_pdg/I");
			//     CutReco->Branch("fmcshwr_status_code",&fmcshwr_status_code,"fmcshwr_status_code/I");
			//     CutReco->Branch("fmcshwr_mass",&fmcshwr_mass,"fmcshwr_mass/F");
			//     CutReco->Branch("fmcshwr_ND",&fmcshwr_ND,"fmcshwr_ND/I");
			//     CutReco->Branch("fmcshwr_mother",&fmcshwr_mother,"fmcshwr_mother/I");
			//     CutReco->Branch("fmcshwr_origin",&fmcshwr_origin,"fmcshwr_origin/I");
			//     CutReco->Branch("fmcshwr_process",&fmcshwr_process,"fmcshwr_process/I");
			//     CutReco->Branch("fmcshwr_rescatter",&fmcshwr_rescatter,"fmcshwr_rescatter/I");

			//     CutReco->Branch("totflash",&totflash,"totflash/I");
			//     CutReco->Branch("flash_distana",flash_distana,"flash_distana[totflash]/F");
			//     CutReco->Branch("flash_peana",flash_peana,"flash_peana[totflash]/F");
				
			//     CutReco->Branch("fmichel_conesize",&fmichel_conesize,"fmichel_conesize/I");
			//     CutReco->Branch("fmichel_wcount",&fmichel_wcount,"fmichel_wcount/I");
			//     CutReco->Branch("fmichel_zpos",fmichel_zpos,"fmichel_zpos[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_ypos",fmichel_ypos,"fmichel_ypos[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_xpos",fmichel_xpos,"fmichel_xpos[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_chrg",fmichel_chrg,"fmichel_chrg[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_chno",fmichel_chno,"fmichel_chno[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_key",fmichel_key,"fmichel_keyo[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_wire",fmichel_wire,"fmichel_wire[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_chargehit",fmichel_chargehit,"fmichel_chargehit[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_tpc",fmichel_tpc,"fmichel_tpc[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_ptime",fmichel_ptime,"fmichel_ptime[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_angledeg",fmichel_angledeg,"fmichel_angledeg[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_maghitveccostheta",fmichel_maghitveccostheta,"fmichel_maghitveccostheta[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_distance",fmichel_distance,"fmichel_distance[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_mult",fmichel_mult,"fmichel_mult[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_sigptime",fmichel_sigptime,"fmichel_sigptime[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_sigchrg",fmichel_sigchrg,"fmichel_sigchrg[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_sigpamp",fmichel_sigpamp,"fmichel_sigpamp[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_dof",fmichel_dof,"fmichel_dof[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_gof",fmichel_gof,"fmichel_gof[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_ptminusRMS",fmichel_ptminusRMS,"fmichel_ptminusRMS[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_ptplusRMS",fmichel_ptplusRMS,"fmichel_ptplusRMS[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_status",fmichel_status,"fmichel_status[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_cnnMichel",fmichel_cnnMichel,"fmichel_cnnMichel[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_cnnEM",fmichel_cnnEM,"fmichel_cnnEM[fmichel_wcount]/F");
			//     CutReco->Branch("fmichel_cnnTrack",fmichel_cnnTrack,"fmichel_cnnTrack[fmichel_wcount]/F");
		}//if CutReco

		fRollUpUnsavedIDs = p.get<bool>("RollUpUnsavedIDs"); 
		fGeom= &*art::ServiceHandle<geo::Geometry>();
	} 

	// --- Our Histograms --- //
	// hDriftTime      = tfs->make<TH2F>("hDriftTime",     "hDriftTime", 100, -400., 400., 100, 0., 10000.);

	MichelsAna::~MichelsAna()
	{}

	void MichelsAna::beginJob()
	{
		utils->reset(true); //deep clean the variables
	}//beginJob



	void MichelsAna::analyze(const art::Event & evt)
	{
		fEvent  = evt.id().event();
		fRun    = evt.id().run();
		fSubRun = evt.id().subRun();
		art::Timestamp ts = evt.time();
		TTimeStamp tts(ts.timeHigh(), ts.timeLow());
		fEvtTime = tts.AsDouble();

		UInt_t year=0, mon=0, day=0, hour=0, min=0, sec=0;
		fyear_mon_day = tts.GetDate(kTRUE,0,&year,&mon,&day);
		fhour_min_sec = tts.GetTime(kTRUE,0,&hour,&min,&sec);
		
		auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt); //All timing information
		auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData); //Detector properties
		fEField = detProp.Efield(); // Electric field in the detector
		fTemper = detProp.Temperature(); // Temperature in the detector
		fDriftV = detProp.DriftVelocity(fEField,fTemper); // Drift velocity in the detector cm/us

		std::map<int, int> fTrackIDMap; // store out of the loop for the truth_tree a map of trackID to MCParticle index

		fFlag = rand() % 10000000000;
		std::string sHead = "################################################################";
		sHead = sHead + "\n  YYYY/MM/DD: " + utils->str(year) + "/" + utils->str(mon) + "/" + utils->str(day) + " HH:MM:SS: " + utils->str(hour) + ":" + utils->str(min) + ":" + utils->str(sec);
		sHead = sHead + "\n  EVENT: " + utils->str(fEvent) + " RUN: " + utils->str(fRun) + " SUBRUN: " + utils->str(fSubRun) + " FLAG: " + utils->str(fFlag);
		sHead = sHead + "\n  TPC Frequency in [MHz]: " + utils->str(clockData.TPCClock().Frequency());
		sHead = sHead + "\n  TPC Tick in [us]: " + utils->str(clockData.TPCClock().TickPeriod());
		sHead = sHead + "\n  Electric Field in [kV/cm]: " + utils->str(fEField) + " | Temperature in [K]: " + utils->str(fTemper);
		sHead = sHead + "\n  Drift Velocity in [cm/us]: " + utils->str(fDriftV);

		sHead = sHead + "\n  Succesfull reset of variables";
		sHead = sHead + "\n################################################################";
		utils->PrintInColor(sHead, utils->GetColor("magenta"));

		const sim::ParticleList &PartList = pi_serv->ParticleList();
		std::string sMcTruth = "";
		sMcTruth = sMcTruth + "\nThere are a total of " + utils->str(int(PartList.size())) + " Particles in the event\n";
		utils->PrintInColor(sMcTruth, utils->GetColor("yellow"));
		// Load the configured trees to store
		std::vector<std::string> trees; 
		int nTrees = fTreesToWrite.size();
		for (auto t = 0; t < nTrees; t++) 
		{
			auto tree = fTreesToWrite[t][0];
			trees.push_back(tree);
		}  

		//=================== =============================== Tree Gen ==================================================//
		if (std::find(trees.begin(), trees.end(), std::string("Gen")) != trees.end()) 
		{
			std::string sMC_ini = "... Filling Generator Truth Info [GEN TREE] ... ";
			utils->PrintInColor(sMC_ini, utils->GetColor("magenta"));

			for (auto i : fProductsToDump)
			{
				std::string sInfoGen = "--------------------------------------------- GENIE (MC TRUTH) ---------------------------------------------------";
				auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
				
				if (o_label == "Gen_MCTruth") 
				{
					if (debug){sInfoGen = sInfoGen + "\n---> Module_Label: " + i_label + ";\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t Output_Name: " + o_label;}
					if(!evt.isRealData())
					{
						// Vector with all the generator/GENIE particles of the event
						art::ValidHandle<std::vector<simb::MCTruth>> MCTrues = evt.getValidHandle<std::vector<simb::MCTruth>>(i_label);
						// --- Loop over the MCTruth vector (have your primary particles stored in the generator) --- //
						int i_gen(0);
						for (auto const &truth : *MCTrues)
						{
              auto gen_result = utils->MuonCandidateAnalysis(evt, fNMichel_gen, i_gen, MCTrues, truth);
              fNMichel_gen    = std::get<0>(gen_result); // Saving the number of Michel electrons
              genPdgCode      = std::get<1>(gen_result); // Saving the PDG code of the generator particle
              genTrackID      = std::get<2>(gen_result); // Saving the TrackID of the generator particle
              genOrigin       = std::get<3>(gen_result); // Saving the Origin of the generator particle
              genNDaughters   = std::get<4>(gen_result); // Saving the number of daughters of the generator particle
							
              ParticleData mumData = std::get<5>(gen_result);
              genIniPosX = mumData.iniPosX;
              genIniPosY = mumData.iniPosY;
              genIniPosZ = mumData.iniPosZ;
              genIniPosT = mumData.iniPosT;
              genEndPosX = mumData.endPosX;
              genEndPosY = mumData.endPosY;
              genEndPosZ = mumData.endPosZ;
              genEndPosT = mumData.endPosT;
              genTrueEne = mumData.Energy;
              genTrueMom = mumData.Moment;
              genIniMomX = mumData.iniMomX;
              genIniMomY = mumData.iniMomY;
              genIniMomZ = mumData.iniMomZ;
              genEndMomX = mumData.endMomX;
              genEndMomY = mumData.endMomY;
              genEndMomZ = mumData.endMomZ;

              ParticleData dauData = std::get<6>(gen_result);
              gen_g4PdgCode.push_back(dauData.PdgCode);
              gen_g4TrackID.push_back(dauData.TrackId);
              gen_g4MumsTID.push_back(dauData.MumTID);
              gen_g4IniPosX.push_back(dauData.iniPosX);
              gen_g4IniPosY.push_back(dauData.iniPosY);
              gen_g4IniPosZ.push_back(dauData.iniPosZ);
              gen_g4IniPosT.push_back(dauData.iniPosT);
              gen_g4EndPosX.push_back(dauData.endPosX);
              gen_g4EndPosY.push_back(dauData.endPosY);
              gen_g4EndPosZ.push_back(dauData.endPosZ);
              gen_g4EndPosT.push_back(dauData.endPosT);
              gen_g4TrueEne.push_back(dauData.Energy);
              gen_g4TrueMom.push_back(dauData.Moment);
              // gen_g4IniMomX.push_back(iniMomX);
              // gen_g4IniMomY.push_back(iniMomY);
              // gen_g4IniMomZ.push_back(iniMomZ);
              // gen_g4EndMomX.push_back(endMomX);
              // gen_g4EndMomY.push_back(endMomY);
              // gen_g4EndMomZ.push_back(endMomZ);

              std::string auxInfoGen = std::get<7>(gen_result);
              sInfoGen = sInfoGen + auxInfoGen;
              sInfoGen = sInfoGen + "\n------------------------------------------------------------------------------------------------------------------";
							utils->PrintInColor(sInfoGen, utils->GetColor("blue"));
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
						utils->PrintInColor(sInfoEdep, utils->GetColor("magenta"));
					}

					art::Handle<std::vector<sim::SimEnergyDeposit>> EdepHandle;
					evt.getByLabel(i_label, i_instance, EdepHandle);
					if(!EdepHandle.isValid()) 
					{ 
						std::string sWarningEdep = "\nUnable to find std::vector<sim::SimEnergyDeposit> with module label: " + i_label;
						utils->PrintInColor(sWarningEdep, utils->GetColor("red"));
					} // !EdepHandle
					if(EdepHandle.isValid())
					{
						for(auto edep : *EdepHandle) 
						{
							fSimEdepTID.push_back(edep.TrackID());
							fSimEdepPDG.push_back(edep.PdgCode());
							fSimEdepX  .push_back(edep.MidPointX());
							fSimEdepY  .push_back(edep.MidPointY());
							fSimEdepZ  .push_back(edep.MidPointZ());
							fSimEdepE  .push_back(edep.E());
						} // for(auto edep : *EdepHandle)
					} // EdepHandle.isValid()  
					// Sim channels
					// art::Handle< std::vector<sim::SimChannel> > simchannelListHandle;
					// std::vector<art::Ptr<sim::SimChannel> > simchannellist;
					// if (evt.getByLabel(fSimChannelModuleLabel,simchannelListHandle)) art::fill_ptr_vector(simchannellist, simchannelListHandle);
				} //if i_type == "sim::SimEnergyDeposit"
			} //for loop over products to dump

			Gen->Fill(); // FILL THE TREE

		} //if Gen is in the vector

		// ================================================== Tree Cry ==================================================//
		if (std::find(trees.begin(), trees.end(), std::string("Cry")) != trees.end()) 
		{
			std::string sCry_ini = "... Filling Cosmic Rays generator Info [CRY TREE] ... ";
			utils->PrintInColor(sCry_ini, utils->GetColor("magenta"));
			for (auto i : fProductsToDump)
			{
				std::string sInfoCry = "---------------------------------------- COSMIC RAY GENERATION (MC) ----------------------------------------------";
				auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
				if (o_label == "Cry_MCTruth") 
				{
					if (debug){sInfoCry = sInfoCry + "\n---> Module_Label: " + i_label + ";\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t Output_Name: " + o_label;}
					if(!evt.isRealData())
					{
						// Vector with all the Cosmic generator particles of the event
						art::ValidHandle<std::vector<simb::MCTruth>> CryTrues = evt.getValidHandle<std::vector<simb::MCTruth>>(i_label);
						// --- Loop over the MCTruth vector (have your primary particles stored in the generator) --- //
						int i_cry(0);
						for (auto const &cry : *CryTrues)
						{

              auto cry_result = utils->MuonCandidateAnalysis(evt, fNMichel_cry, i_cry, CryTrues, cry);
              fNMichel_cry    = std::get<0>(cry_result); // Saving the number of Michel electrons
              cryPdgCode      = std::get<1>(cry_result); // Saving the PDG code of the generator particle
              cryTrackID      = std::get<2>(cry_result); // Saving the TrackID of the generator particle
              cryOrigin       = std::get<3>(cry_result); // Saving the Origin of the generator particle
              cryNDaughters   = std::get<4>(cry_result); // Saving the number of daughters of the generator particle
							
              ParticleData mumData = std::get<5>(cry_result);
              genIniPosX = mumData.iniPosX;
              genIniPosY = mumData.iniPosY;
              genIniPosZ = mumData.iniPosZ;
              genIniPosT = mumData.iniPosT;
              genEndPosX = mumData.endPosX;
              genEndPosY = mumData.endPosY;
              genEndPosZ = mumData.endPosZ;
              genEndPosT = mumData.endPosT;
              genTrueEne = mumData.Energy;
              genTrueMom = mumData.Moment;
              genIniMomX = mumData.iniMomX;
              genIniMomY = mumData.iniMomY;
              genIniMomZ = mumData.iniMomZ;
              genEndMomX = mumData.endMomX;
              genEndMomY = mumData.endMomY;
              genEndMomZ = mumData.endMomZ;

              ParticleData dauData = std::get<6>(cry_result);
              cry_g4PdgCode.push_back(dauData.PdgCode);
              cry_g4TrackID.push_back(dauData.TrackId);
              cry_g4MumsTID.push_back(dauData.MumTID);
              gen_g4IniPosX.push_back(dauData.iniPosX);
              gen_g4IniPosY.push_back(dauData.iniPosY);
              gen_g4IniPosZ.push_back(dauData.iniPosZ);
              gen_g4IniPosT.push_back(dauData.iniPosT);
              gen_g4EndPosX.push_back(dauData.endPosX);
              gen_g4EndPosY.push_back(dauData.endPosY);
              gen_g4EndPosZ.push_back(dauData.endPosZ);
              gen_g4EndPosT.push_back(dauData.endPosT);
              gen_g4TrueEne.push_back(dauData.Energy);
              gen_g4TrueMom.push_back(dauData.Moment);
              // gen_g4IniMomX.push_back(iniMomX);
              // gen_g4IniMomY.push_back(iniMomY);
              // gen_g4IniMomZ.push_back(iniMomZ);
              // gen_g4EndMomX.push_back(endMomX);
              // gen_g4EndMomY.push_back(endMomY);
              // gen_g4EndMomZ.push_back(endMomZ);

              std::string auxInfoCry = std::get<7>(cry_result);
              sInfoCry = sInfoCry + auxInfoCry;
              sInfoCry = sInfoCry + "\n------------------------------------------------------------------------------------------------------------------";
							utils->PrintInColor(sInfoCry, utils->GetColor("blue"));
							i_cry++;
						} //for loop over MCTruth vector
					} //if !evt.isRealData()
				} //if Cry_MCTruth
			} //for loop over products to dump

			Cry->Fill(); // FILL THE TREE

		} //if Cry is in the vector


		// ================================================== Tree AllReco ==================================================//
		if (std::find(trees.begin(), trees.end(), std::string("AllReco")) != trees.end()) 
		{
			std::string sReco_ini = "... Filling AllReco Info [RECO TREE] ... ";
			utils->PrintInColor(sReco_ini, utils->GetColor("magenta"));

			std::string label_pfp,     instance_pfp,     type_pfp; 
			std::string label_tracks,  instance_track,   type_track;
			std::string label_cluster, instance_cluster, type_cluster;
			std::string label_shower,  instance_shower,  type_shower;
			std::string label_hit,     instance_hit,     type_hit;
			std::string label_spacept, instance_spacept, type_spacept;
			std::string label_opflash, instance_opflash, type_opflash;
			
			std::string sInfoReco  = "----------------------------------------------------- RECO -----------------------------------------------------";
			std::string sErrsReco  = "";
			int load_reco_info = 0;
			for (auto i : fProductsToDump)
			{
				auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
				if (o_label == "Reco_PFP")
				{
					if (debug){sInfoReco = sInfoReco + "\n---> Module_Label: " + i_label + ";\t\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t Output_Name: " + o_label;}
					label_pfp = i_label; instance_pfp = i_instance; type_pfp = i_type; 
					sInfoReco = sInfoReco + "\nFounded PFP to analyze.";
					// Number of PFParticles: " + utils->str(fNPFParticles);
					++load_reco_info;
					// } catch (const std::exception& e) {sErrsReco = sErrsReco + "An error occurred: " + e.what();}
				}

				if (o_label == "Reco_Track")
				{
					if (debug){sInfoReco = sInfoReco + "\n---> Module_Label: " + i_label + ";\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t\t Output_Name: " + o_label;}
					label_tracks = i_label; instance_track = i_instance; type_track = i_type;
					sInfoReco = sInfoReco + "\nFounded tracks to analyze.";
					++load_reco_info;
					// } catch (const std::exception& e) {sErrsReco = sErrsReco + "An error occurred: " + e.what();}
				}

				if (o_label == "Reco_Shower")
				{
					if (debug){sInfoReco = sInfoReco + "\n---> Module_Label: " + i_label + ";\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t\t Output_Name: " + o_label;}
					label_shower = i_label; instance_shower = i_instance; type_shower = i_type;
					sInfoReco = sInfoReco + "\nFounded showers to analyze.";
					++load_reco_info;
					// } catch (const std::exception& e) {sErrsReco = sErrsReco + "An error occurred: " + e.what();}
				}

				if (o_label == "Reco_Hit")
				{
					if (debug){sInfoReco = sInfoReco + "\n---> Module_Label: " + i_label + ";\t\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t\t Output_Name: " + o_label;}
					label_hit = i_label; instance_hit = i_instance; type_hit = i_type;
					sInfoReco = sInfoReco + "\nFounded hits to analyze.";
					++load_reco_info;
					// } catch (const std::exception& e) {sErrsReco = sErrsReco + "An error occurred: " + e.what();}
				}

				if (o_label == "Reco_Cluster")
				{
					if (debug){sInfoReco = sInfoReco +  "\n---> Module_Label: " + i_label + ";\t\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t\t Output_Name: " + o_label;}
					label_cluster = i_label; instance_cluster = i_instance; type_cluster = i_type;
					sInfoReco = sInfoReco + "\nFounded clusters to analyze.";
					++load_reco_info;
					// } catch (const std::exception& e) {sErrsReco = sErrsReco + "An error occurred: " + e.what();}
				}

				if (o_label == "Reco_SpacePoint")
				{
					if (debug){sInfoReco = sInfoReco +  "\n---> Module_Label: " + i_label + ";\t\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t Output_Name: " + o_label;}
					label_spacept = i_label; instance_spacept = i_instance; type_spacept = i_type;
					sInfoReco = sInfoReco + "\nFounded spacepoints to analyze.";
					++load_reco_info;
					// } catch (const std::exception& e) {sErrsReco = sErrsReco + "An error occurred: " + e.what();}
				}

				if (o_label == "Reco_OpFlash")
				{
					if (debug){sInfoReco = sInfoReco +  "\n---> Module_Label: " + i_label + ";\t\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t\t Output_Name: " + o_label;}
					label_opflash = i_label; instance_opflash = i_instance; type_opflash = i_type;
					sInfoReco = sInfoReco + "\nFounded External OpFlashes to analyze.";
					++load_reco_info;
					// } catch (const std::exception& e) {sErrsReco = sErrsReco + "An error occurred: " + e.what();}
				}
				//TO DO: add waveforms
			} //for loop over products to dump

			// --- We confirm that the basic RecoObjects are present in the event we are about to analyse
			if (load_reco_info != 7) 
			{
				sErrsReco = sErrsReco + "An error occurred: Not all the necessary information was loaded [" + utils->str(load_reco_info) + "/7]\n";
				utils->PrintInColor(sErrsReco, utils->GetColor("red"));
				std::exit(EXIT_FAILURE);
			}// if load_reco_info != 7 not all the necessary information was loaded
			else
			{
				sInfoReco = sInfoReco + "\n\nAll the necessary information was loaded [" + utils->str(load_reco_info) + "/7]";
				
        // Once we have checked that all the needed information is present, we can proceed to the analysis
        // Now we extract the reco variables to be analysed.

				// --- PFParticles --- //
				std::vector<art::Ptr<recob::PFParticle>  > pfpartVect;
				auto pfpHandle = evt.getHandle<std::vector<recob::PFParticle> >(label_pfp);
				art::fill_ptr_vector(pfpartVect, pfpHandle);
				// pfpartVect  = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, label_pfp);
				fNPFParticles = pfpartVect.size();
				// --- Tracks --- //
				std::vector<art::Ptr<recob::Track>  > tracksVect;
				auto trackHandle = evt.getHandle<std::vector<recob::Track> >(label_tracks);
				art::fill_ptr_vector(tracksVect, trackHandle);
				// --- Showers --- //
				std::vector<art::Ptr<recob::Shower> > showerVect;
				auto showerHandle = evt.getHandle<std::vector<recob::Shower> >(label_shower);
				art::fill_ptr_vector(showerVect, showerHandle);
				// --- Hits --- //
				std::vector<art::Ptr<recob::Hit> > hitVect;
				auto hitHandle = evt.getHandle<std::vector<recob::Hit> >(label_hit);
				art::fill_ptr_vector(hitVect,hitHandle);
				// --- Clusters --- //
				std::vector<art::Ptr<recob::Cluster> > clusterVect;
				auto clusterHandle = evt.getHandle<std::vector<recob::Cluster> >(label_cluster);
				art::fill_ptr_vector(clusterVect,clusterHandle);
				// ---  OpFlashes --- //
				std::vector<art::Ptr<recob::OpFlash> > opflashVect;
				auto opflashHandle = evt.getHandle<std::vector<recob::OpFlash> >(label_opflash);
				art::fill_ptr_vector(opflashVect,opflashHandle);
				// --- SpacePoints --- //
				std::vector<art::Ptr<recob::SpacePoint> > spaceptVect;
				auto spaceptHandle = evt.getHandle<std::vector<recob::SpacePoint> >(label_spacept);
				art::fill_ptr_vector(spaceptVect,spaceptHandle);

				sInfoReco = sInfoReco + "\n-----------------------------------------------------------------------------------------------------------";
				utils->PrintInColor(sInfoReco, utils->GetColor("blue"));

				// --- Associations --- //
        auto const & tracks = * trackHandle;
        art::FindManyP<recob::Hit> hitsFromTracks( trackHandle, evt, label_tracks );

				art::FindManyP<recob::PFParticle> tracks_pfp_assn (tracksVect,   evt,label_tracks);     //to associate tracks with PFParticles
				// art::FindManyP<recob::PFParticle> shower_pfp_assn (showerVect,   evt,label_shower);     //to associate showers with PFParticles
				// art::FindManyP<recob::Cluster>    pfp_cluster     (pfpartVect,   evt, label_cluster);   //to associate PFParticles and clusters
				// art::FindMany<recob::Hit>         tracks_hits_assn(tracksVect,   evt, label_tracks);    //to associate tracks and hits
				// // art::FindMany<recob::Hit>shower_hits_assn(showerVect,   evt, label_shower);          //to associate showers and hits
				// art::FindManyP<recob::Hit>        shower_hits_assn(showerVect,   evt, label_shower);    //to associate showers and hits
				// art::FindManyP<recob::Hit>        cluster_hits    (clusterVect,  evt, label_cluster);   //to associate clusters and hits
				// art::FindManyP<recob::Hit>        tracks_hits_p   (tracksVect,   evt, label_tracks);    //to associate tracks and hits
				// art::FindManyP<recob::Track>      hits_tracks_p   (hitVect,      evt, label_tracks);    //to associate hits with tracks
				// art::FindManyP<recob::Cluster>    hits_clustr_p   (hitVect,      evt, label_cluster);   //to associate hits with clusters
				// art::FindManyP<anab::Calorimetry> tracks_caloSCE  (tracksVect,   evt, "pandoracaloSCE"); //to associate tracks with calorimetry
				// art::FindManyP<recob::Wire>       wire_hits       (hitVect,      evt, label_hit);
				// art::FindManyP<recob::SpacePoint> hits_spacept    (hitVect,      evt, label_spacept);
				// art::FindManyP<recob::OpHit>      in_op_flashs    (opflashVect,  evt, label_opflash);
				art::FindManyP<anab::T0>          t0_pfp_tracks   (pfpartVect,   evt, label_pfp); //to associate PFParticles with T0
				// art::FindManyP<anab::T0>          t0_pfp_shower   (pfpartVect,   evt, label_pfp); //to associate PFParticles with T0
				// art::FindManyP<anab::T0>          t0_tracks       (tracksVect,   evt, "pmtrack"); //to associate tracks with T0
				art::FindManyP<recob::Hit, recob::TrackHitMeta>    tracks_hits_p_meta(tracksVect, evt, label_tracks); // to associate tracks and hits
        anab::MVAReader<recob::Hit,MVA_LENGTH> hitResults( evt, "emtrkmichelid:emtrkmichel" ); // get the hit results helper


        // ------------------------------------------------- Loop over ALL RECO OPTICAL FLASHES ------------------------------------------------- //
				// We are going to store the information of the optical flashes in the event

        // std::cout<< "opflashVect.size() "<<opflashVect.size()<<std::endl;
				fNFlashes =  0; fTrigTime = 0; 
        double high_NPEs = 0;
				for (size_t i = 0; i < opflashVect.size() && i < kMaxFlashes; ++i)
				{
					art::Ptr<recob::OpFlash> pflash(opflashHandle, i);
					const recob::OpFlash& flash = *pflash;
					flash_time   .push_back(flash.Time());
					flash_abstime.push_back(flash.AbsTime()); // Time by PMT readout clock (?)
					flash_totalpe.push_back(flash.TotalPE());
					flash_ycenter.push_back(flash.YCenter());
					flash_zcenter.push_back(flash.ZCenter());
					flash_ywidth .push_back(flash.YWidth());
					flash_zwidth .push_back(flash.ZWidth());
					flash_twidth .push_back(flash.TimeWidth());
					flash_frame  .push_back(flash.Frame()); // Frame number 
					
					if(flash.TotalPE()>high_NPEs)
					{ 
					high_NPEs = flash.TotalPE();
					fTrigTime = flash.Time();
					}

					fNFlashes++;
				}//loop over flashes
        // std::cout<<"fNFlashes "<<fNFlashes<<std::endl;

        // ------------------------------------------------- Loop over ALL RECO HITS ------------------------------------------------- //
        // (We are going to classify the hits in the collection plane of the long tracks.
        // Three categories depending if the hit is contain in the same/another track or not.
        std::vector <double> HitsTracksKey;
        for(size_t tdx=0; tdx<tracksVect.size(); tdx++)
        {
          art::Ptr<recob::Track> aux_track(trackHandle, tdx);
          const recob::Track& track = *aux_track;
          auto hit_track  = tracks_hits_p_meta.at(tdx);
          auto meta_track = tracks_hits_p_meta.data(tdx);

          if(track.Length()> hits_track_length )
          {
            for (size_t hdx = 0; hdx<hit_track.size(); ++hdx)
            {  
              if(hit_track[hdx]->WireID().Plane == 2) HitsTracksKey.push_back(hit_track[hdx].key()); // only hits in collection plane
            }//looping over the hits of the selected long track	    
          }//if(track.Length()> hits_track_length)
        }//loop over tracks   
        std::cout<<"HitsTracksKey.size() "<<HitsTracksKey.size()<<std::endl;

        // ------------------------------------------------- Loop over ALL RECO TRAKS ------------------------------------------------- //
        // We are going to see 
        // 1) if the truth associated info has a Michel on it 
        // 2) by looping on the PFPs associated to the track and after making the selection cuts we are going to analyse the candidate Michel hits
        
        // std::cout<<"tracksVect.size() "<<tracksVect.size()<<std::endl;
        NTracks = 0;
        for(auto const & track : tracks)
        {
          // std::cout<<"NTracks "<<NTracks<<std::endl;
          std::vector<art::Ptr<recob::PFParticle>> pfpsFromTrack=tracks_pfp_assn.at(NTracks);
          // TRUTH INFO:  associated to the reco track
          int fNMichel_mcFromTrack = 0;
          double fmcFromTrackDau_E = 0;
          if(!evt.isRealData()) // Look on the MC truth associated to the track
          {
            // TRUTH MATCHING
            const simb::MCParticle* mcFromTrack = truthUtil.GetMCParticleFromRecoTrack(clockData,track,evt,label_tracks);
            if(!mcFromTrack) continue;

            const art::Ptr<simb::MCTruth> mcFromTrackTID = pi_serv->TrackIdToMCTruth_P(mcFromTrack->TrackId());
            if(!mcFromTrackTID) continue;

            bool muonDecay = false;
            if((abs(mcFromTrack->PdgCode())==13) && (mcFromTrack->NumberDaughters())>0)
            {		    
              muonDecay = utils->isMuonDecaying(mcFromTrack);
            }//if((abs(mcFromTrack->PdgCode())==13) && (mcFromTrack->NumberDaughters())>0)
            if(muonDecay)  
            { 
              std::cout<<"Reco MuonDecay!"<<std::endl;
              if((abs(mcFromTrack->PdgCode())==13) && (mcFromTrack->NumberDaughters())>0)
              {
                for (int ii=0; ii<(mcFromTrack->NumberDaughters());++ii) 
                {
                  const simb::MCParticle* mcFromTrackDau = pi_serv->TrackIdToParticle_P((mcFromTrack->Daughter(ii)));
                  if(!mcFromTrackDau) continue;
                  if(abs(mcFromTrackDau->PdgCode())==11 && mcFromTrackDau->E()>fmcFromTrackDau_E && mcFromTrackDau->Vx()>fiducialBounds1[0] && mcFromTrackDau->Vx()<fiducialBounds1[1] &&
                    mcFromTrackDau->Vy()>fiducialBounds1[2] && mcFromTrackDau->Vy()<fiducialBounds1[3] && 
                    mcFromTrackDau->Vz()>fiducialBounds1[4] && mcFromTrackDau->Vz()<fiducialBounds1[5]) 
                  {     
                    fNMichel_mcFromTrack = 1;    
                    fmcFromTrackDau_E = mcFromTrackDau->E();
                    mcFromTrack_michelE.push_back(fmcFromTrackDau_E);
                  } 
                }//for (int ii=0; ii<(mcFromTrack->NumberDaughters());++ii)	 		       
              }//if((abs(mcFromTrack->PdgCode())==13) && (mcFromTrack->NumberDaughters())>0)
            }//if(hasElectron && hasNuMu && hasNuE) //MCFromTrack // TRUTH INFO //
          }//if(!evt.isRealData())
          // std::cout<<"NTracks "<<NTracks<<std::endl;
          NpfpsFromTrack.push_back(pfpsFromTrack.size());
          NMichelpfpFromTrack.push_back(fNMichel_mcFromTrack);
          
          // RECO INFO
          // find best true track ID for this track
          // int bestTrackId = utils->trackMatching(clockData, &track - &tracks[0], hitsFromTracks );
          // std::cout<<"bestTrackId "<<bestTrackId<<std::endl;

          int fNpfpsFromTrack =  0;
          double fNMichel_pfpFromTrack = 0;

          //at least one PFP association
          if (pfpsFromTrack.size())
          {
            // if (fNMichel_mcFromTrack) fNMichelFromTrack++;
            fNpfpsFromTrack++;

            // Tagged tracks (T0)
            double fNT0 = 0;
            std::vector<art::Ptr<anab::T0>> T0_pfpsFromTrack_vect=t0_pfp_tracks.at(pfpsFromTrack[0].key());
            ft0sana.push_back(T0_pfpsFromTrack_vect.size());
            fMichelcountt0ana.push_back(fNMichel_mcFromTrack);
            ftrueEt0ana.push_back(fmcFromTrackDau_E);
            ft0ana++;

            // T0 tagged tracks
            double T0_value0 = 0;
            if(T0_pfpsFromTrack_vect.size())
            {
              fNT0++;
              fNMichel_pfpFromTrack++;

              T0_value0 = double (T0_pfpsFromTrack_vect.at(0)->Time()); //nano seconds
              std::cout<<"T0_value0 "<<T0_value0<<std::endl;
              TVector3 dir_ini = track.DirectionAtPoint<TVector3>(track.FirstValidPoint());
              TVector3 dir_end = track.DirectionAtPoint<TVector3>(track.LastValidPoint());
              TVector3 pos_end = track.End<TVector3>();
              TVector3 pos_ini(track.LocationAtPoint(track.FirstValidPoint()).X(), track.LocationAtPoint(track.FirstValidPoint()).Y(), track.LocationAtPoint(track.FirstValidPoint()).Z());
          
              fMichelcountstartinboundana[fstartinboundana] = fNMichel_mcFromTrack;
              ftrueEstartinboundana[fstartinboundana] = fmcFromTrackDau_E;
              fstartinboundana++;
              
              // Track ini_pos and end_pos inside the active volume
              if(utils->IsPointInBounds(activeBounds_eff,pos_ini) || utils->IsPointInBounds(activeBounds_eff,pos_end)) // Start or end point within the edges of Active volume bounds.
              {
                // Track start and end points
                auto const [ini_pt, ini] = std::make_tuple(track.Vertex(), track.Vertex<TVector3>());
                auto const [end_pt, end] = std::make_tuple(track.End(), track.End<TVector3>());
                // Check if the track starts or ends in the fiducial volume
                bool iniInFidVol = utils->insideFidVol(ini_pt);
                bool endInFidVol = utils->insideFidVol(end_pt);

                // Check if the track starts or ends in the fiducial volume
                if (iniInFidVol || endInFidVol) 
                {
                  // reset hit counts
                  for (int plane = 0; plane < 3; plane++) {fNCloseHitsIni[plane] = 0; fNCloseHitsEnd[plane] = 0;}
                  // loop over hits --> store the selected
                  std::vector< recob::Hit > taggedHitsIni;
                  std::vector< recob::Hit > taggedHitsEnd;
                  for (size_t hdx = 0; hdx < hitResults.size(); ++hdx) 
                  {
                    // cnn output vector: [em,trk,none,michel]
                    std::array<float,MVA_LENGTH> cnn_out = hitResults.getOutput(hdx);
                    // Looking at hits with large CNN values classified as michel
                    // std::cout << "Your selected fCNNThreshold is: " << fCNNThreshold << std::endl; //TO_DO: Find a way to save all the values given in the configuration file and print them just once
                    if (cnn_out[hitResults.getIndex("michel")] > fCNNThreshold) 
                    {
                      // need to check at both ends of the track
                      if (iniInFidVol) 
                      {
                        for (int plane = 0; plane < 3; plane++) // for each plane check if michel tagged hits are close to the track end
                        {
                          // project the track onto the plane
                          geo::GeometryCore const & geom = *art::ServiceHandle<geo::Geometry>();
                          geo::TPCID tpcid = geom.FindTPCAtPosition(ini_pt);
                          if (!tpcid.isValid) continue;
                          auto const & trackIni2D = pma::GetProjectionToPlane(ini,
                                                                              plane,
                                                                              geom.HasTPC(tpcid),
                                                                              geom.PositionToCryostatID(ini_pt).Cryostat);
                          double ini2D[2] = {trackIni2D.X(), trackIni2D.Y()};
                          const recob::Hit & hit = hitResults.item(hdx); // get the hit itself
                          if (hit.View() == plane) 
                          {
                              // std::cout << "hit.View() " << std::endl;
                              // std::cout << ini_pt.X() << " " << ini_pt.Y() << " " << ini_pt.Z() << std::endl;
                              // std::cout << ini2D[0] << " " << ini2D[1] << std::endl;

                            // check that the hit is close to the track endpoint
                            // and add to tagged hits if so
                            if ( utils->hitCloseToTrackEnd( detProp, fRadiusThreshold, ini_pt, ini2D, hit, geom ) ) 
                            {
                              fNCloseHitsIni[plane] += 1;
                              if (plane == 2) { taggedHitsIni.push_back( hit ); };
                            }//if ( utils->hitCloseToTrackEnd( detProp, fRadiusThreshold, ini_pt, ini2D, hit, geom ) )
                          }//if (hit.View() == plane)
                        }//for ( int plane = 0; plane < 3; plane++ )
                      }//if ( iniInFidVol )
                      if (endInFidVol) 
                      {
                        // In each plane check if the michel tagged hits are close to the end of the track
                        for (int plane = 0; plane < 3; plane++) 
                        {
                          // project the track onto the plane 
                          geo::GeometryCore const & geom = *art::ServiceHandle<geo::Geometry>();
                          geo::TPCID tpcid = geom.FindTPCAtPosition(end_pt);
                          if (!tpcid.isValid) continue;
                          auto const & trackEnd2D = pma::GetProjectionToPlane(end,
                                                                              plane,
                                                                              geom.HasTPC(tpcid),
                                                                              geom.PositionToCryostatID(end_pt).Cryostat);
                          double end2D[2] = {trackEnd2D.X(), trackEnd2D.Y()};
                          const recob::Hit & hit = hitResults.item(hdx); // get the hit itself
                          if (hit.View() == plane) {

                            // check that the hit is close to the track endpoint
                            // and add to tagged hits if so
                            if ( utils->hitCloseToTrackEnd( detProp, fRadiusThreshold, end_pt, end2D, hit, geom ) ) 
                            {
                              fNCloseHitsEnd[plane] += 1;
                              if (plane == 2) { taggedHitsEnd.push_back(hit); }
                            }//if ( utils->hitCloseToTrackEnd( detProp, fRadiusThreshold, end_pt, end2D, hit, geom ) )
                          }//if (hit.View() == plane)
                        }//for ( int plane = 0; plane < 3; plane++ )
                      }//if (endInFidVol)
                    }//if (cnn_out[hitResults.getIndex("michel")] > fCNNThreshold)
                  } // end of loop over hits

                  std::cout << "fNCloseHitsIni[0] " << fNCloseHitsIni[0] << std::endl;
                  std::cout << "fNCloseHitsIni[1] " << fNCloseHitsIni[1] << std::endl;
                  std::cout << "fNCloseHitsIni[2] " << fNCloseHitsIni[2] << std::endl;
                  std::cout << "fNCloseHitsEnd[0] " << fNCloseHitsEnd[0] << std::endl;
                  std::cout << "fNCloseHitsEnd[1] " << fNCloseHitsEnd[1] << std::endl;
                  std::cout << "fNCloseHitsEnd[2] " << fNCloseHitsEnd[2] << std::endl;
                  // event selection decision: clusters of hits in all planes at one end of the track
                  bool iniSelected = ( fNCloseHitsIni[0] > fCloseHitsThreshold && fNCloseHitsIni[1] > fCloseHitsThreshold && fNCloseHitsIni[2] > fCloseHitsThreshold );
                  bool endSelected = ( fNCloseHitsEnd[0] > fCloseHitsThreshold && fNCloseHitsEnd[1] > fCloseHitsThreshold && fNCloseHitsEnd[2] > fCloseHitsThreshold );
                  bool evtSelected = ( iniSelected || endSelected );
                  if ( iniSelected && endSelected ) { evtSelected = false; }
                  
                  // check if the event was correclty tagged
                  //  i.e. do the tagged hits correspond to a michel 
                  bool particleIsMichel(false);
                  // std::cout << "evtSelected " << evtSelected << std::endl;
                  // std::cout << "taggedHitsIni.size() " << taggedHitsIni.size() << std::endl;
                  // std::cout << "taggedHitsEnd.size() " << taggedHitsEnd.size() << std::endl;
                  if (evtSelected) 
                  {
                    if ( iniSelected ) { particleIsMichel = utils->areHitsMichel(clockData, taggedHitsIni); }
                    if ( endSelected ) { particleIsMichel = utils->areHitsMichel(clockData, taggedHitsEnd); }
                  }//if (evtSelected)
                  // std::cout << "particleIsMichel " << particleIsMichel << std::endl;
                  // update purity and efficiecny numbers and draw example events
                  if (  particleIsMichel && evtSelected ) { fYesSelected += 1; }
                  if ( !particleIsMichel && evtSelected ) { fNoSelected  += 1; }
                }//if (iniInFidVol || endInFidVol)
                NTracks++;
              }//if(T0_pfpsFromTrack_vect.size())
            }//if(IsPointInBounds(activeBounds_eff,pos) || IsPointInBounds(activeBounds_eff,end))
          }//if (pfpsFromTrack.size())
        }//for(auto const & track : tracks)

        HitsTracksKey.clear();
        AllReco->Fill();
        CutReco->Fill();
      } //if load_reco_info == 7 we can continue the analysis

      std::cout << "END but not going to print anything" << std::endl;

    } // if (tree == "AllReco")
  }// analyze()

	void MichelsAna::endJob()
	{
		// --- Print summary at the end of the job
		std::string sEndJob = "--------------------------- END JOB ----------------------------------";
    sEndJob = sEndJob + "\nNTracks      " + utils->str(NTracks);
    sEndJob = sEndJob + "\nfYesSelected " + utils->str(fYesSelected);
    sEndJob = sEndJob + "\nfNoSelected  " + utils->str(fNoSelected);
		utils->PrintInColor(sEndJob, utils->GetColor("magenta"));


	}// endJob()

} // namespace michels
DEFINE_ART_MODULE(michels::MichelsAna)
