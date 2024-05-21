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
      bool debug;              // Debugging flag
      double fFidVolCut, fCNNThreshold, fRadiusThreshold, fMuonTrackLengthCut, fMuonTrackLengthHitsCut, fMuonMinHitPeakCut, fMuonMaxHitPeakCut, fHitsDistance, fShowerDistance, fDistHitEnd, fAngleHitEnd, fConeLength, fConeAngle;
      int fCloseHitsThreshold, fMinHitCountMichel, fMaxHitCountMichel;

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
      unsigned int NTracks, fcutFlag, fcutEvent, fcutRun, fcutSubRun, fcut_endhitkey,fcut_endwire,fcut_endchno,fcut_endtpcno;
      int    fNFlashes, fall_trks, fall_trksint, fReadOutWindowSize, fNumberTimeSamples;
      int    fPFP_trks, fT0_trks, fstartinbound_trks, finiinFV_trks, fendinFV_trks, fccrosser_trks;
      int    fnoelecdiv_bound, funbroken_trks, flongcm_trks, fminhitpt_trks, fmaxhitpt_trks;
      int    fnearhits_trks, fcut_mu, funcont_trks, fcrossZ_trks, fbackward_trks, fEndinTPC_trks, NTrueMichels;
      int    fpfpana, ft0ana, fstartinboundana, fendinFVana, fccrossana, fstopzana, fdistana, fbrokcoana, ftrklenana, fminhitptana, fmaxhitptana;
      int    ftrkdistcolana, fhitctana, fshwrdisana, ftrkcolhits, fshwrcolhits, fshwrallhits;
      int    fntrkhits,fhitsV,fhitsU, ftrueparhitallcount;
      int    fhitsY, fnearhitct, fmhitcount, ftrueparhitcolcount;
      int    fyear_mon_day, fhour_min_sec;
      int    fcut_endhitmult, fcut_ccrosser, fcut_trkID, fcut_ncolhits, fcut_nearhitcount, fMichelcountcutana, fcutshwr_key;
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
      double _fPFP_trks=0,_fall_trks=0, _ftruePFP_trks=0, _fT0_trks=0, _ftrueT0_trks=0, _ftrueshwr_trks=0;
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
    debug                   = p.get<bool>("debug");
    fCloseHitsThreshold     = p.get<int>("CloseHitsThreshold"); 
    fMinHitCountMichel      = p.get<int>("MinHitCountMichel");    
    fMaxHitCountMichel      = p.get<int>("MaxHitCountMichel");
    fFidVolCut              = p.get<double>("FidVolCut");
    fCNNThreshold           = p.get<double>("CNNThreshold");   
    fRadiusThreshold        = p.get<double>("RadiusThreshold");
    fMuonTrackLengthCut     = p.get<double>("MuonTrackLengthCut"); 
    fMuonTrackLengthHitsCut = p.get<double>("MuonTrackLengthHitsCut"); 
    fMuonMinHitPeakCut      = p.get<double>("MuonMinHitPeakCut"); 
    fMuonMaxHitPeakCut      = p.get<double>("MuonMaxHitPeakCut");
    fHitsDistance           = p.get<double>("HitsDistance");      
    fShowerDistance         = p.get<double>("ShowerDistance");     
    fDistHitEnd             = p.get<double>("DistHitEnd");         
    fAngleHitEnd            = p.get<double>("AngleHitEnd");
    fConeLength             = p.get<double>("ConeLength");        
    fConeAngle              = p.get<double>("ConeAngle");          

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
              genTrueEne = mumData.iniEne;
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
              gen_g4TrueEne.push_back(dauData.iniEne);
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
              genTrueEne = mumData.iniEne;
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
              gen_g4TrueEne.push_back(dauData.iniEne);
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
      std::string label_wires,   instance_wires,   type_wires;
			
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

        if (o_label == "Reco_Wire")
				{
					if (debug){sInfoReco = sInfoReco +  "\n---> Module_Label: " + i_label + ";\t\t Instance_Name: " + i_instance + ";\t Product_Type: " + i_type + ";\t\t Output_Name: " + o_label;}
					label_wires = i_label; instance_wires = i_instance; type_wires = i_type;
					sInfoReco = sInfoReco + "\nFounded Wires to analyze.";
					++load_reco_info;
					// } catch (const std::exception& e) {sErrsReco = sErrsReco + "An error occurred: " + e.what();}
				}
				//TO DO: add waveforms
				//TO DO: add 
			} //for loop over products to dump

			// --- We confirm that the basic RecoObjects are present in the event we are about to analyse
			if (load_reco_info != 8) 
			{
				sErrsReco = sErrsReco + "An error occurred: Not all the necessary information was loaded [" + utils->str(load_reco_info) + "/7]\n";
				utils->PrintInColor(sErrsReco, utils->GetColor("red"));
				std::exit(EXIT_FAILURE);
			}// if load_reco_info != 7 not all the necessary information was loaded
			else
			{
				sInfoReco = sInfoReco + "\n\nAll the necessary information was loaded [" + utils->str(load_reco_info) + "/8]";
				
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
        // --- Wires --- //
        std::vector<art::Ptr<recob::Wire> > wireVect;
        auto wireHandle = evt.getHandle<std::vector<recob::Wire> >(label_wires);
        art::fill_ptr_vector(wireVect,wireHandle);

				sInfoReco = sInfoReco + "\n-----------------------------------------------------------------------------------------------------------";
				utils->PrintInColor(sInfoReco, utils->GetColor("blue"));

				// --- Associations --- //
        auto const & tracks = * trackHandle;
        art::FindManyP<recob::Hit> hitsFromTracks( trackHandle, evt, label_tracks );

				art::FindManyP<recob::PFParticle> tracks_pfp_assn (tracksVect,   evt,label_tracks);     //to associate tracks with PFParticles
				art::FindManyP<recob::PFParticle> shower_pfp_assn (showerVect,   evt,label_shower);     //to associate showers with PFParticles
				// art::FindManyP<recob::Cluster>    pfp_cluster     (pfpartVect,   evt, label_cluster);   //to associate PFParticles and clusters
				// art::FindMany<recob::Hit>         tracks_hits_assn(tracksVect,   evt, label_tracks);    //to associate tracks and hits
				// // art::FindMany<recob::Hit>shower_hits_assn(showerVect,   evt, label_shower);          //to associate showers and hits
				art::FindManyP<recob::Hit>        shower_hits_assn(showerVect,   evt, label_shower);    //to associate showers and hits
				// art::FindManyP<recob::Hit>        cluster_hits    (clusterVect,  evt, label_cluster);   //to associate clusters and hits
				// art::FindManyP<recob::Hit>        tracks_hits_p   (tracksVect,   evt, label_tracks);    //to associate tracks and hits
				// art::FindManyP<recob::Track>      hits_tracks_p   (hitVect,      evt, label_tracks);    //to associate hits with tracks
				// art::FindManyP<recob::Cluster>    hits_clustr_p   (hitVect,      evt, label_cluster);   //to associate hits with clusters
				art::FindManyP<anab::Calorimetry> tracks_caloSCE  (tracksVect,   evt, "pandoracaloSCE"); //to associate tracks with calorimetry
				art::FindManyP<recob::Wire>       wire_hits       (hitVect,      evt, label_hit);
				art::FindManyP<recob::SpacePoint> hits_spacept    (hitVect,      evt, label_spacept);
				// art::FindManyP<recob::OpHit>      in_op_flashs    (opflashVect,  evt, label_opflash);
				art::FindManyP<anab::T0>          t0_pfp_tracks   (pfpartVect,   evt, label_pfp); //to associate PFParticles with T0
				art::FindManyP<anab::T0>          t0_pfp_shower   (pfpartVect,   evt, label_pfp); //to associate PFParticles with T0
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
        // We now store the hits that satisfy this condition and use them later.
        double dt_min = 99999.; // us
        std::vector <double> HitsTracksKey;

        for(size_t tdx=0; tdx<tracksVect.size(); tdx++)
        {
          art::Ptr<recob::Track> aux_track(trackHandle, tdx);
          const recob::Track& track = *aux_track;
          auto hit_track  = tracks_hits_p_meta.at(tdx);
          auto meta_track = tracks_hits_p_meta.data(tdx);

          if(track.Length()> fMuonTrackLengthHitsCut )
          {
            for (size_t hdx = 0; hdx<hit_track.size(); ++hdx)
            {  
              if(hit_track[hdx]->WireID().Plane == 2) HitsTracksKey.push_back(hit_track[hdx].key()); // only hits in collection plane
            }//looping over the hits of the selected long track	    
          }//if(track.Length()> fMuonTrackLengthHitsCut)
        }//loop over tracks   
        std::cout<<"HitsTracksKey.size() "<<HitsTracksKey.size()<<std::endl;

        // ------------------------------------------------- Loop over ALL RECO TRAKS ------------------------------------------------- //
        // We are going to see 
        // 1) if the truth associated info has a Michel on it 
        // 2) by looping on the PFPs associated to the track and after making the selection cuts we are going to analyse the candidate Michel hits
        
        // Counters
        // initialize all in the same line of code:
        fall_trks = fPFP_trks = fT0_trks = fstartinbound_trks = fendinFV_trks = funbroken_trks = fccrosser_trks = fnoelecdiv_bound = flongcm_trks = fminhitpt_trks = fmaxhitpt_trks = fnearhits_trks = fcut_mu = funcont_trks = fcrossZ_trks = fbackward_trks = fEndinTPC_trks = 0;

        // std::cout<<"tracksVect.size() "<<tracksVect.size()<<std::endl;
        NTracks = 0;
        for(auto const & track : tracks)
        {
          // std::cout<<"NTracks "<<NTracks<<std::endl;
          std::vector<art::Ptr<recob::PFParticle>> pfpsFromTrack=tracks_pfp_assn.at(NTracks);
          fall_trks++; //Count how many tracks we have
          _fall_trks++;

          // TRUTH INFO:  associated to the reco track
          int fNMichel_mcFromTrack = 0;
          double fmcFromTrackDau_E = 0;
          if(!evt.isRealData()) // Look on the MC truth associated to the track
          {
            // TRUTH MATCHING
            const simb::MCParticle* mcFromTrack_raw = truthUtil.GetMCParticleFromRecoTrack(clockData,track,evt,label_tracks);
            if(!mcFromTrack_raw) continue;

            const art::Ptr<simb::MCTruth> mcFromTrackTID_raw = pi_serv->TrackIdToMCTruth_P(mcFromTrack_raw->TrackId());
            if(!mcFromTrackTID_raw) continue;

            bool muonDecay = false;
            if((abs(mcFromTrack_raw->PdgCode())==13) && (mcFromTrack_raw->NumberDaughters())>0)
            {		    
              muonDecay = utils->isMuonDecaying(mcFromTrack_raw);
            }//if((abs(mcFromTrack_raw->PdgCode())==13) && (mcFromTrack_raw->NumberDaughters())>0)
            if(muonDecay)  
            { 
              std::cout<<"Reco MuonDecay!"<<std::endl;
              if((abs(mcFromTrack_raw->PdgCode())==13) && (mcFromTrack_raw->NumberDaughters())>0)
              {
                for (int ii=0; ii<(mcFromTrack_raw->NumberDaughters());++ii) 
                {
                  const simb::MCParticle* mcFromTrackDau = pi_serv->TrackIdToParticle_P((mcFromTrack_raw->Daughter(ii)));
                  if(!mcFromTrackDau) continue;
                  if(abs(mcFromTrackDau->PdgCode())==11 && mcFromTrackDau->E()>fmcFromTrackDau_E && mcFromTrackDau->Vx()>fiducialBounds1[0] && mcFromTrackDau->Vx()<fiducialBounds1[1] &&
                    mcFromTrackDau->Vy()>fiducialBounds1[2] && mcFromTrackDau->Vy()<fiducialBounds1[3] && 
                    mcFromTrackDau->Vz()>fiducialBounds1[4] && mcFromTrackDau->Vz()<fiducialBounds1[5]) 
                  {     
                    fNMichel_mcFromTrack = 1;    
                    fmcFromTrackDau_E = mcFromTrackDau->E();
                    mcFromTrack_michelE.push_back(fmcFromTrackDau_E);
                  } 
                }//for (int ii=0; ii<(mcFromTrack_raw->NumberDaughters());++ii)	 		       
              }//if((abs(mcFromTrack_raw->PdgCode())==13) && (mcFromTrack_raw->NumberDaughters())>0)
            }//if(hasElectron && hasNuMu && hasNuE) //MCFromTrack // TRUTH INFO //
          }//if(!evt.isRealData())
          // std::cout<<"NTracks "<<NTracks<<std::endl;
          NpfpsFromTrack.push_back(pfpsFromTrack.size());
          NMichelpfpFromTrack.push_back(fNMichel_mcFromTrack);
          fpfpana++;
          
          // RECO INFO
          // find best true track ID for this track
          // int bestTrackId = utils->trackMatching(clockData, &track - &tracks[0], hitsFromTracks );
          // std::cout<<"bestTrackId "<<bestTrackId<<std::endl;

          int fNpfpsFromTrack =  0;
          double fNMichel_pfpFromTrack = 0;

          //at least one PFP association
          if (pfpsFromTrack.size())
          {
            fPFP_trks++; //Count how many tracks have at least one PFP associated
            _fPFP_trks++;
            if (fNMichel_mcFromTrack) _ftruePFP_trks++;
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
              fT0_trks++;  // T0 tagged tracks per event count
              _fT0_trks++; // total T0 tagged track count
              fNT0++;
              fNMichel_pfpFromTrack++;
              if(fNMichel_mcFromTrack==1) _ftrueT0_trks++;

              T0_value0 = double (T0_pfpsFromTrack_vect.at(0)->Time()); //nano seconds
              std::cout<<"T0_value0 "<<T0_value0<<std::endl;
              TVector3 dir_ini = track.DirectionAtPoint<TVector3>(track.FirstValidPoint());
              TVector3 dir_end = track.DirectionAtPoint<TVector3>(track.LastValidPoint());
              TVector3 pos_end = track.End<TVector3>();
              TVector3 pos_ini(track.LocationAtPoint(track.FirstValidPoint()).X(), track.LocationAtPoint(track.FirstValidPoint()).Y(), track.LocationAtPoint(track.FirstValidPoint()).Z());
          
              fMichelcountstartinboundana[fstartinboundana] = fNMichel_mcFromTrack;
              ftrueEstartinboundana[fstartinboundana] = fmcFromTrackDau_E;
              fstartinboundana++;

              //TODO:simplify or decide best way of checking limits
              // Track start and end points
              auto const [ini_pt, ini] = std::make_tuple(track.Vertex(), track.Vertex<TVector3>());
              auto const [end_pt, end] = std::make_tuple(track.End(), track.End<TVector3>());
              // Check if the track starts or ends in the fiducial volume
              bool iniInFidVol = utils->insideFidVol(ini_pt); //TODO:maybe IsPointInBounds?
              bool endInFidVol = utils->insideFidVol(end_pt);
              
              // Track ini_pos and end_pos inside the active volume
              // Start or end point within the edges of Active volume bounds.
              // if(utils->IsPointInBounds(activeBounds_eff,pos_ini) || utils->IsPointInBounds(activeBounds_eff,pos_end)) 
              // Check if the track starts or ends in the fiducial volume
              if (iniInFidVol || endInFidVol) 
              {
                fstartinbound_trks++;  // start in bounds tracks per event count
                _fstartinbound_trks++; // total start in bounds track count
                if (fNMichel_mcFromTrack) _ftruestartinbound_trks++;
                // We switch the start and end points if the 
                if(utils->IsPointInBounds(activeBounds_eff,pos_end)){utils->SwitchEndPoints(pos_ini, pos_end, dir_ini,dir_end);}  


                float trkinix = pos_ini.X();
                float trkiniy = pos_ini.Y();
                float trkiniz = pos_ini.Z();
                float trkendx = pos_end.X();
                float trkendy = pos_end.Y();
                float trkendz = pos_end.Z();	 
                //-----------Checking which tracks are cathode crossers-------///
                int ccrosser = 0;
                if((ini.X()*end.X()) < 0) ccrosser = 1; //TODO: confir this does what we want
                  
                fendinFVsana[fendinFVana] = IsPointInFV(fiducialBounds,end);
                fMichelcountendinFVana[fendinFVana] = fNMichel_mcFromTrack;
                ftrueEendinFVana[fendinFVana] = fmcFromTrackDau_E;
                fendinFVana++;  

                // End Point within the fiducial volume.
                if (endInFidVol)
                {
                  fendinFV_trks++;  // start in bounds tracks per event count
                  _fendinFV_trks++; // total start in bounds track count

                  //************************Getting truth info of the tracks to count the true michels**********************//
                  NTrueMichels = 0;
                  double daugx = -999, daugy = -999, daugz = -999, mcd_vx = -999, mcd_vy = -999, mcd_vz = -999, mcd_t = -999, mcd_endx = -999, mcd_endy = -999, mcd_endz = -999, mcd_endt = -999, mcd_px = -999, mcd_py = -999, mcd_pz = -999, mcd_momentum = -999, mcd_energy = -999, mcd_endpx = -999, mcd_endpy = -999, mcd_endpz = -999, mcd_endenergy = -999, mcd_pathlen = -999, mcd_vxdrifted = -999, mcd_vydrifted = -999, mcd_vzdrifted = -999, mcd_tdrifted = -999, mcd_endxdrifted = -999, mcd_endydrifted = -999, mcd_endzdrifted = -999, mcd_endtdrifted = -999, mcd_pxdrifted = -999, mcd_pydrifted = -999, mcd_pzdrifted = -999, mcd_momentumdrifted = -999, mcd_energydrifted = -999, mcd_endpxdrifted = -999, mcd_endpydrifted = -999, mcd_endpzdrifted = -999, mcd_endenergydrifted = -999, mcd_pathlendrifted = -999, mcd_theta = -999, mcd_phi = -999, mcd_mass = -999;
                  int mcd_trkid = -999, mcd_endprocess = -999, mcd_pdg = -999, mcd_status_code = -999, mcd_ND = -999, mcd_mother = -999, mcd_origin = -999, mcd_process = -999, mcd_rescatter = -999;
                  
                  std::cout<<"\n"<<std::endl;
                  // TODO: same structure as the part before filtering. Maybe write a function to do the common things.
                  if(!evt.isRealData())
                  {
                    // TRUTH MATCHING
                    const simb::MCParticle* mcFromTrack_cut1 = truthUtil.GetMCParticleFromRecoTrack(clockData,track,evt,label_tracks);
                    if(!mcFromTrack_cut1) continue;

                    const art::Ptr<simb::MCTruth> mcFromTrackTID = pi_serv->TrackIdToMCTruth_P(mcFromTrack_cut1->TrackId());
                    if(!mcFromTrackTID) continue;
                    double distance = 99999/*, distance0 = 99999*/;
                    bool muonDecay_cut1 = false;
                    if((abs(mcFromTrack_cut1->PdgCode())==13) && (mcFromTrack_cut1->NumberDaughters())>0)
                    {		    
                      muonDecay_cut1 = utils->isMuonDecaying(mcFromTrack_cut1);
                    }//if((abs(mcFromTrack_cut1->PdgCode())==13) && (mcFromTrack_cut1->NumberDaughters())>0)
                    if(muonDecay_cut1)  
                    { 
                      std::cout<<"Reco MuonDecay (cut1)!"<<std::endl;
                      if((abs(mcFromTrack_cut1->PdgCode())==13) && (mcFromTrack_cut1->NumberDaughters())>0)
                      {
                        double daugE = 0;
                        for (int ii=0; ii<(mcFromTrack_cut1->NumberDaughters());++ii) 
                        {
                          const simb::MCParticle* mcFromTrackDau_cut1 = pi_serv->TrackIdToParticle_P((mcFromTrack_cut1->Daughter(ii)));
                          if(!mcFromTrackDau_cut1) continue;
                          
                          TLorentzVector dmcstart, dmcend, dmcstartdrifted, dmcenddrifted;
                          unsigned int dpstarti, dpendi, dpstartdriftedi, dpenddriftedi; //mcparticle indices for starts and ends in tpc or drifted volumes
                          double dplen = utils->length(*mcFromTrackDau_cut1, dmcstart, dmcend, dpstarti, dpendi);
                          double dplendrifted = utils->driftedLength(detProp, *mcFromTrackDau_cut1, dmcstartdrifted, dmcenddrifted, dpstartdriftedi, dpenddriftedi);
                          // bool isActive = plen != 0;
                          bool disDrifted = dplendrifted!= 0;
                          distance=sqrt(pow(mcFromTrackDau_cut1->Vx()-trkendx,2)+pow(mcFromTrackDau_cut1->Vz()-trkendz,2)); //only in collection plane

                          std::cout<<"mother pdg "<<mcFromTrack_cut1->PdgCode()<<" mother energy "<<mcFromTrack_cut1->E()<<" true end x, y, z "<<mcFromTrack_cut1->EndX()<<" "<<mcFromTrack_cut1->EndY()<<" "<<mcFromTrack_cut1->EndZ()<<" trackid "<<mcFromTrack_cut1->TrackId()<<std::endl;
                          std::cout<<"daughter start distance "<<distance<<" x, y, z "<<mcFromTrackDau_cut1->Vx()<<" "<<mcFromTrackDau_cut1->Vy()<<" "<<mcFromTrackDau_cut1->Vz()<<" end x, y, z "<<mcFromTrackDau_cut1->EndX()<<" "<<mcFromTrackDau_cut1->EndY()<<" "<<mcFromTrackDau_cut1->EndZ()<<std::endl;
                          std::cout<<"daughter->PdgCode() "<<mcFromTrackDau_cut1->PdgCode()<<" daughter energy "<<mcFromTrackDau_cut1->E()<<" daughter process "<<mcFromTrackDau_cut1->Process()<<" daughter track id "<<mcFromTrackDau_cut1->TrackId()<<std::endl;

                          //TODO: fhicl configure this value
                          if(distance<=30 && abs(mcFromTrackDau_cut1->PdgCode())==11 && mcFromTrackDau_cut1->E()>daugE && mcFromTrackDau_cut1->Vx()>fiducialBounds1[0] && mcFromTrackDau_cut1->Vx()<fiducialBounds1[1] &&
                            mcFromTrackDau_cut1->Vy()>fiducialBounds1[2] && mcFromTrackDau_cut1->Vy()<fiducialBounds1[3] && 
                            mcFromTrackDau_cut1->Vz()>fiducialBounds1[4] && mcFromTrackDau_cut1->Vz()<fiducialBounds1[5]) // z positions are slightly different between a true and reco particles 
                          {     
                            fdistanceana[fdistana] = distance;
                            fdistana++;
                            // distance0 = distance;
                            NTrueMichels  = 1;
                            ParticleData mcFromTrackDau_cut1_struct = utils->extractParticleData(mcFromTrackDau_cut1);
                            mcd_pathlen = dplen;
                            if(disDrifted)
                            {
                              ParticleData mcFromTrackDau_cut1_drifted_struct = utils->extractParticleData(mcFromTrackDau_cut1);
                              mcd_pathlendrifted = dplendrifted;
                            }//if(disDrifted)
                            mcd_endprocess = int(mcFromTrackDau_cut1->Trajectory().ProcessToKey(mcFromTrackDau_cut1->EndProcess()));
                            mcd_theta      = mcFromTrackDau_cut1->Momentum().Theta();
                            mcd_phi        = mcFromTrackDau_cut1->Momentum().Phi();
                            mcd_pdg        = mcFromTrackDau_cut1->PdgCode();
                            mcd_status_code= mcFromTrackDau_cut1->StatusCode();
                            mcd_mass       = mcFromTrackDau_cut1->Mass();
                            mcd_ND         = mcFromTrackDau_cut1->NumberDaughters();
                            mcd_mother     = mcFromTrackDau_cut1->Mother();
                            // mcd_origin     = mcFromTrack_cut1->Origin();
                            mcd_process    = int(mcFromTrackDau_cut1->Trajectory().ProcessToKey(mcFromTrackDau_cut1->Process()));
                            mcd_rescatter  = mcFromTrackDau_cut1->Rescatter();
          
                            std::cout<<"daughter->PdgCode() "<<mcFromTrackDau_cut1->PdgCode()<<" daughter track id "<<mcFromTrackDau_cut1->TrackId()<<std::endl;
                            // cout<<" daughter.Mother() "<<daughter->Mother()<<" vx, vy, vz "<<fmcd_vx<<" "<<fmcd_vy<<" "<<fmcd_vz<<" endx, endy, endz "<<fmcd_endx<<" "<<fmcd_endy<<" "<<fmcd_endz<<std::endl;

                          }//if distance < 30
                        }//for (int ii=0; ii<(mcFromTrack_cut1->NumberDaughters());++ii)	 		       
                      }//if((abs(mcFromTrack_cut1->PdgCode())==13) && (mcFromTrack_cut1->NumberDaughters())>0)
                    }//if(hasElectron && hasNuMu && hasNuE) //MCFromTrack // TRUTH INFO //
                  }//if(!evt.isRealData())
      
                  if(fNMichel_mcFromTrack==1) _ftrueendinFV_trks++;
                  //***************Consider only cathode crossers******************//
                  fccrosserana[fccrossana] = ccrosser;
                  fMichelcountccrosserana[fccrossana] = NTrueMichels;
                  ftrueEccrosserana[fccrossana] = mcd_energy;
                  fccrossana++;
                  if(ccrosser == 1)
                  {
                    std::cout<<"Cathode crosser!"<<std::endl;
                    fccrosser_trks++;
                    _fccrosser_trks++;
                    if(fNMichel_mcFromTrack==1) _ftrueccrosser_trks++;
                    
                    //***************Removing electron diverters boundaries******************//
                    felecdivstopzana[fstopzana] = trkendz;
                    fMichelcountelecdivstopzana[fstopzana] = NTrueMichels;
                    ftrueEelecdivstopzana[fstopzana] = mcd_energy;
                    fstopzana++;

                    if((trkendz<APAnegbound1 || trkendz>APAposbound1) && (trkendz<APAnegbound2 || trkendz>APAposbound2))
                    {
                      fnoelecdiv_bound++;
                      _fnoelecdiv_bound++;
                      if(NTrueMichels==1) _ftruenoelecdiv_bound++;

                      /********************************** Broken tracks removal ****************************************/
                      // We loop over the tracks to see if there is some other track near the main loop one, with same direction
                      // and if so remove it from the selection
                      size_t bcount = 0;
                      for(unsigned int k=0;k<tracksVect.size();++k)
                      {
                        art::Ptr<recob::Track> ptrack_k(trackHandle, k);
                        const recob::Track& track_k = *ptrack_k;

                        TVector3 dir_ini_k, dir_end_k;
                        TVector3 pos_k(track_k.LocationAtPoint(track_k.FirstValidPoint()).X(), track_k.LocationAtPoint(track_k.FirstValidPoint()).Y(), track_k.LocationAtPoint(track_k.FirstValidPoint()).Z());
                        TVector3 end_k = track_k.End<TVector3>();
                        dir_ini_k = track_k.DirectionAtPoint<TVector3>(track_k.FirstValidPoint());
                        dir_end_k   = track_k.DirectionAtPoint<TVector3>(track_k.LastValidPoint());

                        if(k==NTracks) continue;
              
                        double bwdiststart          = sqrt(pow(pos_k.Y()-end.Y(),2)+pow(pos_k.Z()-end.Z(),2)); //since the non-T0 tagged tracks have wrong x position associated with them
                        double cosopeningangleStart = cos(dir_end.Theta())*cos(dir_ini_k.Theta())+sin(dir_end.Theta())*sin(dir_ini_k.Theta())*cos(dir_ini_k.Phi()-dir_end.Phi());
                        double bwdistend            = sqrt(pow(end_k.Y()-end.Y(),2)+pow(end_k.Z()-end.Z(),2));
                        double cosopeningangleEnd   = cos(dir_end.Theta())*cos(dir_end_k.Theta())+sin(dir_end.Theta())*sin(dir_end_k.Theta())*cos(dir_end_k.Phi()-dir_end.Phi());

                        if((abs(bwdiststart)<30 && abs(cosopeningangleStart) > 0.97) || (abs(bwdistend)<30 && abs(cosopeningangleEnd) > 0.97)
                          ||(abs(bwdiststart)<50 && abs(cosopeningangleStart) > 0.998) || (abs(bwdistend)<50 && abs(cosopeningangleEnd) > 0.998))
                        { bcount++; } //if broken track
                      }//for(size_t k=0;k<tracksVect.size();++k)
                      
                      // fbrokencountana[fbrokcoana] = tracksVect.size()-count;
                      fbrokencountana[fbrokcoana] = bcount;
                      fMichelcountbrokencountana[fbrokcoana] = NTrueMichels;
                      ftrueEbrokencountana[fbrokcoana] = mcd_energy;
                      fbrokcoana++;

                      // if(count==tracksVect.size()-1) //unbroken tracks

                      // SELECT UNBROKEN TRACKS [hardware problem in SP in pple in HD should be solved]
                      if(bcount==0) //unbroken tracks
                      {
                        funbroken_trks++;
                        _funbroken_trks++;
                        if(NTrueMichels==1)_ftrueunbroken_trks++;
        
                        //------------------track length cut------------------// 
                        ftrklengthana[ftrklenana] = track.Length();
                        fMichelcountlenana[ftrklenana] = NTrueMichels;
                        ftrueElenana[ftrklenana] = mcd_energy;
                        ftrklenana++;

                        // Select tracks with at least X cm long --> fhicl input
                        if(track.Length()> fMuonTrackLengthCut ) 
                        { 
                          flongcm_trks++;
                          _flongcm_trks++;
                          if(NTrueMichels==1)_ftruelongcm_trks++;

                          //----------------Get Hits associated with track-----------//
                          int ncolhits = 0;
                          const std::vector<const recob::Hit*> Hits = trackUtil.GetRecoTrackHits(track,evt,label_tracks);
                          // Fill vector with Hit Peak Times and store minimum
                          std::vector<double> HitPeakTimes;
                          for (unsigned int hitIndex = 0; hitIndex < Hits.size(); ++hitIndex)  
                          {
                            HitPeakTimes.push_back(Hits[hitIndex]->PeakTime());
                            if(Hits[hitIndex]->WireID().Plane==2) ncolhits++;      
                          }//loop over hits
                          float MinHitPeakTime = *(std::min_element(HitPeakTimes.begin(), HitPeakTimes.end()));
                          float MaxHitPeakTime = *(std::max_element(HitPeakTimes.begin(), HitPeakTimes.end()));
            
                          HitPeakTimes.clear();

                          //--------------Min hit peak time cut--------------------// 
                          fminhitptimeana[fminhitptana] = MinHitPeakTime;
                          fMichelcountminhitptimeana[fminhitptana] = NTrueMichels;
                          ftrueEminhitptimeana[fminhitptana] = mcd_energy;
                          fminhitptana++;
                          if(MinHitPeakTime > fMinHitCountMichel)
                          {
                            fminhitpt_trks++;
                            _fminhitpt_trks++;
                            if(NTrueMichels==1)_ftrueminhitpt_trks++;
                            //-----------Max hit peak time cut------------------//
                            fmaxhitptimeana[fmaxhitptana] = MaxHitPeakTime;
                            fMichelcountmaxhitptimeana[fmaxhitptana] = NTrueMichels;
                            ftrueEmaxhitptimeana[fmaxhitptana] = mcd_energy;
                            fmaxhitptana++;
          
                            if(MaxHitPeakTime < fMaxHitCountMichel)
                            {
                              fmaxhitpt_trks++;
                              _fmaxhitpt_trks++;
                              if(NTrueMichels==1)_ftruemaxhitpt_trks++;

                              std::cout<<"NTrueMichels "<<NTrueMichels<<" daugx "<<daugx<<" daugy "<<daugy<<" daugz "<<daugz<<std::endl;
                              //*****************Storing the candidate muon track collection plane hits****************************//
                              //Define variables in same line code:
                              float dist=-999,enddist=-999,endhitkey=-999,endpeaktime=-999,endwireno=-999,endtpcno=-999,endchno=-999,endhitchrg=-999,endhitx=-999,endhity=-999,endhitz=-999,endhitmult=-999,endhitsigptime=-999,endhitsigchrg=-999,endhitsigpamp=-999,endhitdof=-999,endhitgof=-999,endhitptminusRMS=-999,endhitptplusRMS=-999;

                              double hits_key[kMaxCh]       ={-999};
                              double hits_charge[kMaxCh]    ={-999};
                              double hits_chno[kMaxCh]      ={-999};
                              double hits_wire[kMaxCh]      ={-999};
                              double hits_peakT[kMaxCh]     ={-999};
                              double hits_TPC[kMaxCh]       ={-999};
                              float hits_xpos[kMaxCh]       ={-999};
                              float hits_ypos[kMaxCh]       ={-999};
                              float hits_zpos[kMaxCh]       ={-999};
                              float hits_mult[kMaxCh]       ={-999};
                              float hits_sigptime[kMaxCh]   ={-999};
                              float hits_sigchrg[kMaxCh]    ={-999};
                              float hits_sigpamp[kMaxCh]    ={-999};
                              float hits_dof[kMaxCh]        ={-999};
                              float hits_gof[kMaxCh]        ={-999};
                              float hits_ptminusRMS[kMaxCh] ={-999};
                              float hits_ptplusRMS[kMaxCh]  ={-999};
                              float hits_cnnMichel[kMaxCh]  ={-999};
                              float hits_cnnEM[kMaxCh]      ={-999};
                              float hits_cnnTrack[kMaxCh]   ={-999};
                              
                              //*************Determine the position of the last hit and storing the collection plane hits for each track****************************//
      
                              int trkcolhits=0;
                              std::vector <double> longtrk_hitkey;
              
                              if(tracks_hits_p_meta.isValid())
                              {
                                auto vhit=tracks_hits_p_meta.at(NTracks);
                                auto vmeta=tracks_hits_p_meta.data(NTracks);

                                //loop over all meta data hit
                                for (size_t ii = 0; ii<vhit.size(); ++ii)
                                {
                                  bool fBadhit = false;
                                  if (vmeta[ii]->Index() == std::numeric_limits<int>::max())
                                  { fBadhit = true; }
                                  if (vmeta[ii]->Index()>=tracksVect[NTracks]->NumberTrajectoryPoints())
                                  { fBadhit = true; }
                                  if (!tracksVect[NTracks]->HasValidPoint(vmeta[ii]->Index()))
                                  { fBadhit = true; }

                                  TVector3 loc;	  
                                  if (fBadhit) continue;
                                  else loc = tracksVect[NTracks]->LocationAtPoint<TVector3>(vmeta[ii]->Index());
                        
                                  if(vhit[ii]->WireID().Plane==2)
                                  {
                                    longtrk_hitkey.push_back(vhit[ii].key());
                                    if(loc.Z()!=0 && (loc.Z()<APAnegbound1 || loc.Z()>APAposbound1) && (loc.Z()<APAnegbound2 || loc.Z()>APAposbound2))
                                    {
                                      hits_key[trkcolhits]       =   vhit[ii].key();
                                      hits_charge[trkcolhits]    =   vhit[ii]->Integral();///(TMath::Exp(-(hits_peakT[trkcolhits1]-800)*500/tau)); //multiplied by 500nsec to convert time ticks to actual generation time
                                      hits_wire[trkcolhits]      =   vhit[ii]->WireID().Wire;
                                      hits_peakT[trkcolhits]     =   vhit[ii]->PeakTime();
                                      hits_TPC[trkcolhits]       =   vhit[ii]->WireID().TPC;
                                      hits_chno[trkcolhits]      =   vhit[ii]->Channel();
                                      hits_xpos[trkcolhits]      =   loc.X();
                                      hits_ypos[trkcolhits]      =   loc.Y();
                                      hits_zpos[trkcolhits]      =   loc.Z();
                                      hits_mult[trkcolhits]      =   vhit[ii]->Multiplicity();
                                      hits_sigptime[trkcolhits]  =   vhit[ii]->SigmaPeakTime();
                                      hits_sigchrg[trkcolhits]   =   vhit[ii]->SigmaIntegral();
                                      hits_sigpamp[trkcolhits]   =   vhit[ii]->SigmaPeakAmplitude();
                                      hits_dof[trkcolhits]       =   vhit[ii]->DegreesOfFreedom();
                                      hits_gof[trkcolhits]       =   vhit[ii]->GoodnessOfFit();      
                                      hits_ptminusRMS[trkcolhits]=   vhit[ii]->PeakTimeMinusRMS(5.0);
                                      hits_ptplusRMS[trkcolhits] =   vhit[ii]->PeakTimePlusRMS(5.0);

                                      std::array<float,4> cnn_out = hitResults.getOutput(vhit[ii].key());
                                      hits_cnnMichel[trkcolhits] =   cnn_out[ hitResults.getIndex("michel") ]; 
                                      // hits_cnnEM[trkcolhits]     =   cnn_out[ hitResults.getIndex("em") ]; 
                                      // hits_cnnTrack[trkcolhits]  =   cnn_out[ hitResults.getIndex("track") ]; 
                                          
                                      dist=sqrt(pow(hits_xpos[trkcolhits]-trkendx,2)+pow(hits_ypos[trkcolhits]-trkendy,2)+pow(hits_zpos[trkcolhits]-trkendz,2));
                  
                                      trkcolhits++;
                                      //---------storing the information for the end point of the track----------//
                                      if(abs(dist)<enddist)
                                      {
                                        enddist=dist;
                                        endhitkey=vhit[ii].key();
                                        endwireno=vhit[ii]->WireID().Wire;
                                        endpeaktime=vhit[ii]->PeakTime();
                                        endtpcno=vhit[ii]->WireID().TPC;
                                        endchno=vhit[ii]->Channel();
                                        endhitchrg=vhit[ii]->Integral();
                                        endhitx=loc.X();
                                        endhity=loc.Y();
                                        endhitz=loc.Z();
                                        endhitmult=vhit[ii]->Multiplicity();
                                        endhitsigptime=vhit[ii]->SigmaPeakTime();
                                        endhitsigchrg=vhit[ii]->SigmaIntegral();
                                        endhitsigpamp=vhit[ii]->SigmaPeakAmplitude();
                                        endhitdof=vhit[ii]->DegreesOfFreedom();
                                        endhitgof=vhit[ii]->GoodnessOfFit();
                                        endhitptminusRMS=vhit[ii]->PeakTimeMinusRMS(5.0);
                                        endhitptplusRMS=vhit[ii]->PeakTimePlusRMS(5.0);
                                      }//if(abs(dist)<enddist)
                                    }//if(loc.Z()!=0 && (loc.Z()<APAnegbound1 || loc.Z()>APAposbound1) && (loc.Z()<APAnegbound2 || loc.Z()>APAposbound2))
                                  }//if(vhit[ii]->WireID().Plane==2)
                                }//loop over all meta data hit
                              }//if(tracks_hits_p_meta.isValid())

                              // MinHitPeakTime & MaxHitPeakTime conditions satifistied:
                              // We compute now the first (and last) non-zero distances from the end wire location
                              // For that we need to rearrange the wires such that the track end point is at the end of the wire order
                              
                              double dist0= 999, dist1 = 999;
                              // first non-zero hit distance
                              for(int jj = 0; jj<trkcolhits; jj++)
                              {
                                dist0 = sqrt(pow(hits_xpos[jj]-trkendx,2)+pow(hits_ypos[jj]-trkendy,2)+pow(hits_zpos[jj]-trkendz,2));
                                break;  
                              }
                              // last non-zero hit distance
                              for(int jj = trkcolhits-1; jj>=0; jj--)
                              {
                                dist1 = sqrt(pow(hits_xpos[jj]-trkendx,2)+pow(hits_ypos[jj]-trkendy,2)+pow(hits_zpos[jj]-trkendz,2));
                                break;
                              }

                              //------Rearranging the wires-----------//
                              std::vector<Hit> hits(trkcolhits);
                              if(abs(dist1) > abs(dist0)) 
                              {
                                for(int jj=0; jj<trkcolhits-1; jj++) 
                                {
                                  for (int kk = jj + 1; kk < trkcolhits; kk++) 
                                  {
                                      std::swap(hits[jj], hits[kk]);
                                  }//for (int kk = jj + 1; kk < trkcolhits; kk++)
                                }//for(int jj=0; jj<trkcolhits-1; jj++)
                              }//if(abs(dist1) > abs(dist0))

                              //--------Finding the hits close to the stopping point of the tracks-------//
                              //bool foundMichel = false;
                              int nearhitct=0;
                              int nearhits_key[kMaxHits]       = {-999};
                              float nearhits_peakT[kMaxHits]   = {-999};
                              float nearhits_charge[kMaxHits]  = {-999};
                              float nearhits_wire[kMaxHits]    = {-999};
                              float nearhits_chno[kMaxHits]    = {-999};
                              float nearhits_TPC[kMaxHits]     = {-999};
                              float nearhits_plane[kMaxHits]   = {-999};
                              float nearhits_xpos[kMaxHits]    = {-999};
                              float nearhits_ypos[kMaxHits]    = {-999};
                              float nearhits_zpos[kMaxHits]    = {-999};
                              float nearhits_mult[kMaxHits]    = {-999};
                              float nearhits_sigptime[kMaxHits]= {-999};
                              float nearhits_sigchrg[kMaxHits] = {-999};
                              float nearhits_sigpamp[kMaxHits]  = {-999};
                              float nearhits_dof[kMaxHits]     = {-999};
                              float nearhits_gof[kMaxHits]     = {-999};
                              float nearhits_ptminusRMS[kMaxHits]    = {-999};
                              float nearhits_ptplusRMS[kMaxHits]     = {-999};
                              float nearhits_cnnMichel[kMaxHits]     = {-999};
                              float nearhits_cnnEM[kMaxHits]         = {-999};
                              float nearhits_cnnTrack[kMaxHits]      = {-999};

                              //cout<<"trkcolhits "<<trkcolhits<<" longtrk_hitkey.size() "<<longtrk_hitkey.size()<<std::endl;

                              // auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
                              // art::ServiceHandle<geo::Geometry const> geometry;

                              for(size_t ll=0; ll<hitVect.size();++ll) //loop over all hits
                              {
                                //auto tracks = thass.at(hitVect[ll].key());
                                if(hitVect[ll]->WireID().Plane==2)
                                {
                                  //----Make sure that this hit does not belong to the candidate muon track----//
                                  // if(hitVect[ll].key()==38025)cout<<"hitVect[ll].ptime "<<hitVect[ll]->PeakTime()<<" charge "<<hitVect[ll]->Integral()<<std::endl;
                                  int same_trk = 0;
                                  for(size_t mm=0;mm<longtrk_hitkey.size();mm++)
                                  {
                                    if(longtrk_hitkey.at(mm)==hitVect[ll].key()) { same_trk = 1; }
                                  }
                        
                                  //---Make sure that these hits don't belong to another long track-----//
                                  int long_trk = 0;
                                  for(size_t o1=0;o1<HitsTracksKey.size();o1++)
                                  {
                                    if(HitsTracksKey.at(o1)==hitVect[ll].key())
                                    {
                                      long_trk = 1;
                                    }
                                  }

                                  double diffpeaktt0     = hitVect[ll]->PeakTime() - (T0_value0/1000)*2; //calculate the peak time in ticks
                                  // double diffpeaktt0     = hitVect[ll]->PeakTime();
                                  double allhitX = detProp.ConvertTicksToX(diffpeaktt0,hitVect[ll]->WireID().Plane,hitVect[ll]->WireID().TPC,hitVect[ll]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
                                  double Wirestart[3], Wireend[3];
                                  fGeom->WireEndPoints(hitVect[ll]->WireID());
                                  // fGeom->WireEndPoints(hitVect[ll]->WireID().Cryostat, hitVect[ll]->WireID().TPC, hitVect[ll]->WireID().Plane, hitVect[ll]->WireID().Wire, Wirestart, Wireend);
                                  double allhitZ = Wirestart[2];
                                  double allhitY = 0;
                                  if(hits_spacept.at(ll).size()) 
                                  allhitY = hits_spacept.at(ll)[0]->XYZ()[1];
                                  if(!hits_spacept.at(ll).size()) 
                                  allhitY = 303.5;	

                                  if(same_trk==0 && long_trk==0 && hitVect[ll]->WireID().TPC==endtpcno &&(allhitZ<APAnegbound1 || allhitZ>APAposbound1) && (allhitZ<APAnegbound2 || allhitZ>APAposbound2))
                                  {
                                    double hitdist = sqrt(pow(abs(allhitX-trkendx),2)+pow(abs(allhitZ-trkendz),2)); 
                                    if(hitdist <= _hitdist)
                                    {
                                      nearhits_key[nearhitct]     = hitVect[ll].key();
                                      nearhits_peakT[nearhitct]   = hitVect[ll]->PeakTime();
                                      nearhits_charge[nearhitct]  = hitVect[ll]->Integral();
                                      nearhits_wire[nearhitct]    = hitVect[ll]->WireID().Wire;
                                      nearhits_chno[nearhitct]    = hitVect[ll]->Channel();
                                      nearhits_TPC[nearhitct]     = hitVect[ll]->WireID().TPC;
                                      nearhits_plane[nearhitct]   = hitVect[ll]->WireID().Plane;
                                      nearhits_xpos[nearhitct]    = allhitX;
                                      nearhits_ypos[nearhitct]    = allhitY;
                                      nearhits_zpos[nearhitct]    = allhitZ;
                                      nearhits_mult[nearhitct]      = hitVect[ll]->Multiplicity();
                                      nearhits_sigptime[nearhitct]  = hitVect[ll]->SigmaPeakTime();
                                      nearhits_sigchrg[nearhitct]   = hitVect[ll]->SigmaIntegral();
                                      nearhits_sigpamp[nearhitct]	= hitVect[ll]->SigmaPeakAmplitude();
                                      nearhits_dof[nearhitct]	= hitVect[ll]->DegreesOfFreedom();
                                      nearhits_gof[nearhitct]       = hitVect[ll]->GoodnessOfFit();
                                      nearhits_ptminusRMS[nearhitct]= hitVect[ll]->PeakTimeMinusRMS(5.0);
                                      nearhits_ptplusRMS[nearhitct] = hitVect[ll]->PeakTimePlusRMS(5.0);

                                      std::array<float,4> cnn_out = hitResults.getOutput( hitVect[ll].key() );
                                      nearhits_cnnMichel[nearhitct] = cnn_out[ hitResults.getIndex("michel") ]; 
                                      nearhits_cnnEM[nearhitct]     = cnn_out[ hitResults.getIndex("em") ]; 
                                      nearhits_cnnTrack[nearhitct]  = cnn_out[ hitResults.getIndex("track") ]; 

                                      std::cout<<"hitVect[ll].key() "<<hitVect[ll].key()<<" peakT "<<hitVect[ll]->PeakTime()<<" wireID "<<hitVect[ll]->WireID().Wire<<" x, y, z "<<nearhits_xpos[nearhitct]<<" "<<nearhits_ypos[nearhitct]<<" "<<nearhits_zpos[nearhitct]<<std::endl;
                                      std::cout<<"Last hit: track ID "<<track.ID()<<" wireno: "<<endwireno<<" peaktime:  "<<endpeaktime<<" tpcno: "<<endtpcno<<" hitchrg: "<<endhitchrg<<" This hit: wire: "<<hitVect[ll]->WireID().Wire<<" peaktime "<<hitVect[ll]->PeakTime()<<" TPC: "<<hitVect[ll]->WireID().TPC<<std::endl;

                                      nearhitct++;
                                    }//if(hitdist <= _hitdist)
                                  }// if(same_trk==0 && long_trk==0 && hitVect[ll]->WireID().TPC==endtpcno &&(allhitZ<APAnegbound1 || allhitZ>APAposbound1) && (allhitZ<APAnegbound2 || allhitZ>APAposbound2))
                                }//if(hitVect[ll]->WireID().Plane==2)
                              }//for(size_t ll=0; ll<hitVect.size();++ll) //loop over all hits
    
                              std::cout<<"\nCANDIDATE MUON COLLECTION PLANE HITS "<<trkcolhits<<std::endl;
                              
                              //----------Looping over showers-------------//
                              double shwr_dist        = 99999, shwr_dist0 = 99999; 		 
                              int    shwr_key         = -999;
                              int    shwr_ID          = -999;
                              double shwr_length      = -999;
                              double shwr_startx      = -999;
                              double shwr_starty      = -999;
                              double shwr_startz      = -999;
                              int    shwr_bestplane   = -999;
                              double shwr_startdcosx  = -999;
                              double shwr_startdcosy  = -999;
                              double shwr_startdcosz  = -999;
                              double shwr_openangle   = -999;

                              std::cout<<"showerVect.size() "<<showerVect.size()<<std::endl;
                              for(size_t jjj=0; jjj<showerVect.size();++jjj)
                              {
                                art::Ptr<recob::Shower> pshower(showerHandle, jjj);
                                const recob::Shower& shower = *pshower;

                                std::vector<art::Ptr<recob::PFParticle>> pfpsh=shower_pfp_assn.at(jjj);
                                if(pfpsh.size()) 
                                { 
                                  //---Only consider T0 tagged shower-----//
                                  std::vector<art::Ptr<anab::T0>> t0sh=t0_pfp_shower.at(pfpsh[0].key());
                                  if(t0sh.size())
                                  {  
                                    TVector3 const& shwr_start = shower.ShowerStart();
                                    float shwr_startx1     = shwr_start.X();
                                    float shwr_starty1     = shwr_start.Y();
                                    float shwr_startz1     = shwr_start.Z();

                                    shwr_dist = sqrt(pow(shwr_startx1 - trkendx,2) + pow(shwr_startz1 - trkendz,2));
                                    if( shwr_dist<shwr_dist0)
                                    {			   
                                      //cout<<"jjj "<<jjj<<" showerVect[jjj].key() "<<showerVect[jjj].key()<<" shwr_ID "<<shower.ID()<<std::endl;
                                      shwr_dist0       = shwr_dist;
                                      shwr_key         = showerVect[jjj].key();
                                      shwr_ID          = shower.ID();
                                      shwr_length      = shower.Length();
                                      shwr_startx      = shwr_startx1;
                                      shwr_starty      = shwr_starty1;
                                      shwr_startz      = shwr_startz1;
                                      shwr_bestplane   = shower.best_plane();
                                      TVector3 const& shwrdir_start = shower.Direction();
                                      shwr_startdcosx  = shwrdir_start.X();
                                      shwr_startdcosy  = shwrdir_start.Y();
                                      shwr_startdcosz  = shwrdir_start.Z();  
                                      shwr_openangle   = shower.OpenAngle();
                                    }// if( shwr_dist<shwr_dist0)
                                  }//if(t0sh.size())
                                }//if(pfpsh.size())
                              }// for(size_t jjj=0; jjj<showerVect.size();++jjj)
      
                              fnearhitcountana[fhitctana] = nearhitct;   
                              fnshwrdistana[fhitctana] = shwr_dist0;
                              fMichelcountnearhitana[fhitctana] = NTrueMichels;
                              ftrueEnearhitana[fhitctana] = mcd_energy;
                              fhitctana++;

                              ///*****************Nearby hit count cut******************//   
                              std::cout<<"nearhitct "<<nearhitct<<std::endl;
                              if(nearhitct>= _minhitcountmichel && nearhitct< _maxhitcountmichel)
                              {
                                fnearhits_trks++;
                          
                                _fnearhits_trks++;
                                if(NTrueMichels==1)_ftruenearhits_trks++;
                                std::cout<<"Nearby hit count "<<nearhitct<<std::endl;	

                                std::vector<double> HitPeakTimesshwrall;
                                std::vector<double> HitPeakTimesshwrcol;
                                std::vector<art::Ptr<recob::Hit >> shwrhits = shower_hits_assn.at(shwr_key);		
                              
                                for(size_t itr = 0; itr<shwrhits.size(); ++itr)
                                {
                                  // cout<<"inside shower hits plane: "<<shwrhits[itr]->WireID().Plane<<std::endl;
                                  HitPeakTimesshwrall.push_back(shwrhits[itr]->PeakTime());
                                  // looping over the collection plane  hits of the selected shower
                                  if(shwrhits[itr]->WireID().Plane == 2){ HitPeakTimesshwrcol.push_back(shwrhits[itr]->PeakTime()); }
                                }// for(itr = shwrhits.begin(); itr < shwrhits.end(); itr++)
                                float MinHitPeakTimeshwrall = 0;
                                float MinHitPeakTimeshwrcol = 0;
                                float MaxHitPeakTimeshwrall = 0;
                                float MaxHitPeakTimeshwrcol = 0;
                                if(HitPeakTimesshwrall.size() && HitPeakTimesshwrcol.size())
                                {
                                  MinHitPeakTimeshwrall = *(std::min_element(HitPeakTimesshwrall.begin(), HitPeakTimesshwrall.end()));
                                  MinHitPeakTimeshwrcol = *(std::min_element(HitPeakTimesshwrcol.begin(), HitPeakTimesshwrcol.end()));
                                  MaxHitPeakTimeshwrall = *(std::max_element(HitPeakTimesshwrall.begin(), HitPeakTimesshwrall.end()));
                                  MaxHitPeakTimeshwrcol = *(std::max_element(HitPeakTimesshwrcol.begin(), HitPeakTimesshwrcol.end()));
                                }
                                std::cout<<"MinHitPeakTimeshwrall "<<MinHitPeakTimeshwrall<<" MinHitPeakTimeshwrcol "<<MinHitPeakTimeshwrcol<<std::endl;
                                std::cout<<"MaxHitPeakTimeshwrall "<<MaxHitPeakTimeshwrall<<" MaxHitPeakTimeshwrcol "<<MaxHitPeakTimeshwrcol<<std::endl;
                          
                                HitPeakTimesshwrall.clear();
                                HitPeakTimesshwrcol.clear();

                                fshwrdistana[fshwrdisana] = shwr_dist0;
                                fMichelcountshwrdistana[fshwrdisana] = NTrueMichels;
                                ftrueEshwrdistana[fshwrdisana] = mcd_energy;
                                fshwrdisana++;

                                //=====================================================//
                                //----------- FINAL SELECTION CUT ---------------------//
                                //=====================================================//

                                if(abs(shwr_dist0)< _shwr_dist && MinHitPeakTimeshwrcol > _minhitpeakcut  && MaxHitPeakTimeshwrcol < _maxhitpeakcut)
                                {
                                  _fshwr_trks++;
                                  if(NTrueMichels==1)_ftrueshwr_trks++;
                            
                                  std::cout<<"\n"<<"\n"<<"Selected the michel event!!!! run, subrun, event numbers: "<<fRun<<" "<<fSubRun<<" "<<fEvent<<" TrkackID "<<track.ID()<<" start x, y, z "<<trkinix<<" "<<trkiniy<<"  "<<trkiniz<<" end x, y, z "<<trkendx<<" "<<trkendy<<"  "<<trkendz<<"  Track length "<<track.Length()<<std::endl;
                                  std::cout<<"shwr_dist0 "<<shwr_dist0<<" fshwr_key "<<shwr_key<<" fshwr_startx "<<shwr_startx<<" fshwr_starty "<<shwr_starty<<" fshwr_startz "<<shwr_startz<<std::endl;

                                  //----------Order near hits such that the distance from the muon end point increases--------------------//
                                  std::cout<<" end xpos "<<endhitx<<" end zpos "<<endhitz<<std::endl;
                                  float hits_dis1 = 0, hits_dis2 = 0;
                                  float aaa=0, bbb=0, ccc=0, ddd=0, eee=0, fff=0, ggg=0, hhh=0, iii=0, mmm=0, jjj=0, kkk=0, lll=0, nnn=0, ooo=0, ppp=0, qqq=0, rrr=0, sss=0, ttt=0, uuu=0;
                                  std::vector<Hit> hits(nearhitct);

                                  //TODO: replace long varibles with structures to avoid more loops
                                  // Populate the hits vector...
                                  // for(int i = 0; i < nearhitct; i++) {
                                  //     hits[NTracks].key = nearhits_key[NTracks];
                                  //     hits[NTracks].charge = nearhits_charge[NTracks];
                                  //     hits[NTracks].wire = nearhits_wire[NTracks];
                                  //     // ... set other properties ...
                                  //     hits[NTracks].cnnTrack = nearhits_cnnTrack[NTracks];
                                  //     hits[NTracks].distance = sqrt(pow((nearhits_xpos[NTracks]-endhitx),2)+pow((nearhits_zpos[NTracks]-endhitz),2));
                                  // }

                                  // std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
                                  //     return a.distance < b.distance;
                                  // });
                                  for(int qq =0; qq < nearhitct; qq++)
                                  {
                                    for(int rr =qq+1; rr < nearhitct; rr++)
                                    {
                                      hits_dis1 = sqrt(pow((nearhits_xpos[qq]-endhitx),2)+pow((nearhits_zpos[qq]-endhitz),2));
                                      hits_dis2 = sqrt(pow((nearhits_xpos[rr]-endhitx),2)+pow((nearhits_zpos[rr]-endhitz),2));
                                      if(hits_dis1 > hits_dis2)
                                      {
                                        aaa = nearhits_key[qq];
                                        nearhits_key[qq] = nearhits_key[rr];
                                        nearhits_key[rr] = aaa;
                                
                                        bbb = nearhits_charge[qq];
                                        nearhits_charge[qq] = nearhits_charge[rr];
                                        nearhits_charge[rr] = bbb;
                                        
                                        ccc = nearhits_wire[qq];
                                        nearhits_wire[qq] = nearhits_wire[rr];
                                        nearhits_wire[rr] = ccc;
                                        
                                        ddd = nearhits_chno[qq];
                                        nearhits_chno[qq] = nearhits_chno[rr];
                                        nearhits_chno[rr] = ddd;
                                        
                                        eee = nearhits_TPC[qq];
                                        nearhits_TPC[qq] = nearhits_TPC[rr];
                                        nearhits_TPC[rr] = eee;
                                        
                                        fff = nearhits_plane[qq];
                                        nearhits_plane[qq] = nearhits_plane[rr];
                                        nearhits_plane[rr] = fff;
                                        
                                        ggg = nearhits_xpos[qq];
                                        nearhits_xpos[qq] = nearhits_xpos[rr];
                                        nearhits_xpos[rr] = ggg;
                                        
                                        mmm = nearhits_ypos[qq];
                                        nearhits_ypos[qq] = nearhits_ypos[rr];
                                        nearhits_ypos[rr] = mmm;

                                        hhh = nearhits_zpos[qq];
                                        nearhits_zpos[qq] = nearhits_zpos[rr];
                                        nearhits_zpos[rr] = hhh;
                                        
                                        iii = nearhits_peakT[qq];
                                        nearhits_peakT[qq] = nearhits_peakT[rr];
                                        nearhits_peakT[rr] = iii;

                                        jjj = nearhits_mult[qq];
                                        nearhits_mult[qq] = nearhits_mult[rr];
                                        nearhits_mult[rr] = jjj;

                                        kkk = nearhits_sigptime[qq];
                                        nearhits_sigptime[qq] = nearhits_sigptime[rr];
                                        nearhits_sigptime[rr] = kkk;

                                        lll = nearhits_sigchrg[qq];
                                        nearhits_sigchrg[qq] = nearhits_sigchrg[rr];
                                        nearhits_sigchrg[rr] = lll;

                                        nnn = nearhits_sigpamp[qq];
                                        nearhits_sigpamp[qq] = nearhits_sigpamp[rr];
                                        nearhits_sigpamp[rr] = nnn;

                                        ooo = nearhits_dof[qq];
                                        nearhits_dof[qq] = nearhits_dof[rr];
                                        nearhits_dof[rr] = ooo;

                                        ppp = nearhits_gof[qq];
                                        nearhits_gof[qq] = nearhits_gof[rr];
                                        nearhits_gof[rr] = ppp;

                                        qqq = nearhits_ptminusRMS[qq];
                                        nearhits_ptminusRMS[qq] = nearhits_ptminusRMS[rr];
                                        nearhits_ptminusRMS[rr] = qqq;

                                        rrr = nearhits_ptplusRMS[qq];
                                        nearhits_ptplusRMS[qq] = nearhits_ptplusRMS[rr];
                                        nearhits_ptplusRMS[rr] = rrr;

                                        sss = nearhits_cnnMichel[qq];
                                        nearhits_cnnMichel[qq] = nearhits_cnnMichel[rr];
                                        nearhits_cnnMichel[rr] = sss;

                                        ttt = nearhits_cnnEM[qq];
                                        nearhits_cnnEM[qq] = nearhits_cnnEM[rr];
                                        nearhits_cnnEM[rr] = ttt;

                                        uuu = nearhits_cnnTrack[qq];
                                        nearhits_cnnTrack[qq] = nearhits_cnnTrack[rr];
                                        nearhits_cnnTrack[rr] = uuu;
                                      }// if(hits_dis1 > hits_dis2)		       
                                    }// for(int rr =qq+1; rr < nearhitct; rr++)			 
                                  }// for (int qq =0; qq < nearhitct; qq++) 

                                  //-------------checking the correlation of the michel hits----------------//	       	 
                                  float sumwire=0, sumptime=0, meanwire=0, meanptime=0, diffWirePtime=0, diffWire=0, diffPtime=0, stDWire=0, stDPtime=0, CorrWirePtime=0;
                                  for(int qq =0; qq < nearhitct; qq++)
                                  {
                                    std::cout<<"nearby hits_dis "<<sqrt(pow((nearhits_xpos[qq]-endhitx),2)+pow((nearhits_zpos[qq]-endhitz),2))<<" TPC "<<nearhits_TPC[qq]<<" wire "<<nearhits_wire[qq]<<" peakt "<<nearhits_peakT[qq]<<" x, y, z "<<nearhits_xpos[qq]<<" "<<nearhits_ypos[qq]<<" "<<nearhits_zpos[qq]<<std::endl;
                                    sumwire += nearhits_wire[qq];
                                    sumptime += nearhits_peakT[qq];
                                  }
                        
                                  meanwire = sumwire/nearhitct;
                                  meanptime = sumptime/nearhitct;
                            
                                  for(int qq =0; qq < nearhitct; qq++)
                                  {
                                    diffWirePtime += (nearhits_wire[qq]-meanwire)*(nearhits_peakT[qq]-meanptime);
                                    diffWire      += (nearhits_wire[qq]-meanwire)*(nearhits_wire[qq]-meanwire);
                                    diffPtime     += (nearhits_peakT[qq]-meanptime)*(nearhits_peakT[qq]-meanptime);
                                  }  
                      
                                  stDWire         = sqrt(diffWire/nearhitct);
                                  stDPtime        = sqrt(diffPtime/nearhitct);
                                  CorrWirePtime   = diffWirePtime/((nearhitct)*stDWire*stDPtime);
                                  //cout<<"CorrWirePtime "<<CorrWirePtime<<std::endl;
                                    

                                  // -------- Storing selected candidate nuon quantities -------------//
                                  fcutRun          = evt.run();
                                  fcutSubRun       = evt.subRun();
                                  fcutEvent        = evt.id().event();
                                  fcutEvtTime      = tts.AsDouble();
                                  fcut_endhitkey    = endhitkey;	
                                  fcut_endwire      = endwireno;
                                  fcut_endchno      = endchno;
                                  fcut_endtpcno     = endtpcno;
                                  fcut_endhitchrg   = endhitchrg;
                                  fcut_endptime     = endpeaktime;
                                  fcut_endhitx      = endhitx;
                                  fcut_endhity      = endhity;
                                  fcut_endhitz      = endhitz;
                                  fcut_endhitmult      = endhitmult;
                                  fcut_endhitsigptime  = endhitsigptime;
                                  fcut_endhitsigchrg   = endhitsigchrg;
                                  fcut_endhitsigpamp    = endhitsigpamp;
                                  fcut_endhitdof       = endhitdof;
                                  fcut_endhitgof       = endhitgof;
                                  fcut_endhitptplusRMS = endhitptplusRMS;
                                  fcut_endhitptminusRMS = endhitptminusRMS;
                                  fcut_ccrosser     = ccrosser;
                                  fcut_minhitptime  = MinHitPeakTime;
                                  fcut_maxhitptime  = MaxHitPeakTime;
                                  fcut_ncolhits     = trkcolhits;
                                  fcut_nearhitcount = nearhitct;
                                  fcut_CorrWirePtime= CorrWirePtime;
                                  fcut_trkrecotime  = double(T0_value0/1000); //track reco time T0 (convert from nsec to usec)
                          
                                  fcut_dist_hit_end = enddist;
                                  fcut_dist_times   = dt_min;
                                  fcut_trackthetaxz = std::atan2(dir_ini.X(), dir_ini.Z());
                                  fcut_trackthetayz = std::atan2(dir_ini.Y(), dir_ini.Z());
                                  fcut_trkstartx    = trkinix;  
                                  fcut_trkstarty    = trkiniy;  
                                  fcut_trkstartz    = trkiniz;  
                                  fcut_trkendx      = trkendx;
                                  fcut_trkendy      = trkendy;
                                  fcut_trkendz      = trkendz;
                                  fcut_trkstartcosx = dir_ini.X();
                                  fcut_trkstartcosy = dir_ini.Y();
                                  fcut_trkstartcosz = dir_ini.Z();
                                  fcut_trkendcosx   = dir_end.X();
                                  fcut_trkendcosy   = dir_end.Y();
                                  fcut_trkendcosz   = dir_end.Z();
                                  fcut_trklen       = track.Length();
                                  fcut_trktheta     = dir_ini.Theta();
                                  fcut_trkphi       = dir_ini.Phi();
                                  fcut_trkID        = track.ID();
                                  //fcut_PHratio      = avgdownchrg/avgupchrg;
                                  //fcut_PH           = PH;
                          
                                  fMichelcountcutana= NTrueMichels;
                            
                                  //---------storing selected shower reco quantities -------------//
                                  fcutshwr_key        =  shwr_key;
                                  fcutshwr_ID         =  shwr_ID;
                                  fcutshwr_length     =  shwr_length;
                                  fcutshwr_startx     =  shwr_startx;
                                  fcutshwr_starty     =  shwr_starty;
                                  fcutshwr_startz     =  shwr_startz;
                                  fcutshwr_bestplane  =  shwr_bestplane;
                                  fcutshwr_startdcosx =  shwr_startdcosx;
                                  fcutshwr_startdcosx =  shwr_startdcosy;
                                  fcutshwr_startdcosx =  shwr_startdcosz;  
                                  fcutshwr_openangle  =  shwr_openangle;
                                  fcutshwr_dist       =  shwr_dist0;
                                  // fcutshwr_dEdx       =  shwr_dEdx;
                                  // fcutshwr_energy     =  shwr_energy;
                                  // fcutshwr_mipenergy  =  shwr_mipenergy;  


                                  //-------------Saving selected candidate muon hits info------------------------//
                            
                                  ftrkcolhits = 0;
                                  for(int mm=0; mm<trkcolhits; mm++)
                                  {
                                    fhits_key[ftrkcolhits]      =  hits_key[mm];
                                    fhits_charge[ftrkcolhits]   =  hits_charge[mm];
                                    fhits_wire[ftrkcolhits]     =  hits_wire[mm];
                                    fhits_peakT[ftrkcolhits]    =  hits_peakT[mm];
                                    fhits_TPC[ftrkcolhits]      =  hits_TPC[mm];
                                    fhits_chno[ftrkcolhits]     =  hits_chno[mm];
                                    fhits_xpos[ftrkcolhits]     =  hits_xpos[mm];
                                    fhits_ypos[ftrkcolhits]     =  hits_ypos[mm];
                                    fhits_zpos[ftrkcolhits]     =  hits_zpos[mm];
                                    fhits_mult[ftrkcolhits]     =  hits_mult[mm];
                                    fhits_sigptime[ftrkcolhits] =  hits_sigptime[mm];
                                    fhits_sigchrg[ftrkcolhits]  =  hits_sigchrg[mm];
                                    fhits_sigpamp[ftrkcolhits]   =  hits_sigpamp[mm];
                                    fhits_dof[ftrkcolhits]      =  hits_dof[mm];
                                    fhits_gof[ftrkcolhits]      =  hits_gof[mm];
                                    fhits_ptminusRMS[ftrkcolhits]= hits_ptminusRMS[mm];
                                    fhits_ptplusRMS[ftrkcolhits]	 = hits_ptplusRMS[mm];
                                    fhits_cnnMichel[ftrkcolhits]	 = hits_cnnMichel[mm];
                                    fhits_cnnEM[ftrkcolhits]	 = hits_cnnEM[mm];
                                    fhits_cnnTrack[ftrkcolhits]	 = hits_cnnTrack[mm];
                              
                                    ftrkcolhits++;
                                  }//for(int mm=0; mm<trkcolhits; mm++)
          
                                  //-----------Filling selected shower hits info---------------------//
                                  fshwrcolhits = 0;
                                  fshwrallhits = 0;
                      
                                  // std::vector<const recob::Hit* > shwrhits = shower_hits_assn.at(shwr_key);	
                      
                                  // for(std::vector<const recob::Hit* >::iterator itr = shwrhits.begin(); itr < shwrhits.end(); itr++)
                                  for(size_t itr = 0; itr<shwrhits.size(); ++itr)
                                  {  
                                    
                                    //fshwrallhits_key[fshwrallhits]     =   shwrhits[itr]->ID();
                                    fshwrallhits_chno[fshwrallhits]      =   shwrhits[itr]->Channel();
                                    fshwrallhits_peakT[fshwrallhits]     =   shwrhits[itr]->PeakTime();
                                    fshwrallhits_charge[fshwrallhits]    =   shwrhits[itr]->Integral();///(TMath::Exp(-(hits_peakT[longcolhits]-800)*500/tau)); //multiplied by 500nsec to convert time ticks to actual generation time
                                    fshwrallhits_wire[fshwrallhits]      =   shwrhits[itr]->WireID().Wire;
                                    fshwrallhits_plane[fshwrallhits]     =   shwrhits[itr]->WireID().Plane;
                                    fshwrallhits_TPC[fshwrallhits]       =   shwrhits[itr]->WireID().TPC;
                                    double diffpeaktt0allshwr            =   shwrhits[itr]->PeakTime() - (T0_value0/1000)*2; //calculate the peak time in ticks
                                    fshwrallhits_xpos[fshwrallhits]      =   detProp.ConvertTicksToX(diffpeaktt0allshwr,shwrhits[itr]->WireID().Plane,shwrhits[itr]->WireID().TPC,shwrhits[itr]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
                                    double Wirestartallshwr[3], Wireendallshwr[3];
                                    fGeom->WireEndPoints(shwrhits[itr]->WireID());
                                    // fGeom->WireEndPoints(shwrhits[itr]->WireID().Cryostat, shwrhits[itr]->WireID().TPC, shwrhits[itr]->WireID().Plane, shwrhits[itr]->WireID().Wire, Wirestartallshwr, Wireendallshwr);
                                    //			   fshwrallhits_ypos[fshwrallhits]      =   Wirestartallshwr[1];
                                    fshwrallhits_zpos[fshwrallhits]      =   Wirestartallshwr[2];
                                    if(hits_spacept.at(shwrhits[itr].key()).size())
                                    fshwrallhits_ypos[fshwrallhits]      =   hits_spacept.at(shwrhits[itr].key())[0]->XYZ()[1];
                                    if(!hits_spacept.at(shwrhits[itr].key()).size())
                                    fshwrallhits_ypos[fshwrallhits]      =   303.5;
                                    fshwrallhits_mult[fshwrallhits]      =   shwrhits[itr]->Multiplicity();
                                    fshwrallhits_sigptime[fshwrallhits]  =   shwrhits[itr]->SigmaPeakTime();
                                    fshwrallhits_sigchrg[fshwrallhits]   =   shwrhits[itr]->SigmaIntegral();
                                    fshwrallhits_sigpamp[fshwrallhits]	=   shwrhits[itr]->SigmaPeakAmplitude();
                                    fshwrallhits_dof[fshwrallhits]       =   shwrhits[itr]->DegreesOfFreedom();
                                    fshwrallhits_gof[fshwrallhits]       =   shwrhits[itr]->GoodnessOfFit();
                                    fshwrallhits_ptminusRMS[fshwrallhits]= shwrhits[itr]->PeakTimeMinusRMS(5.0);
                                    fshwrallhits_ptplusRMS[fshwrallhits]	 = shwrhits[itr]->PeakTimePlusRMS(5.0);

                                    fshwrallhits++;
                              
                                    // looping over the collection plane  hits of the selected shower
                                    if(shwrhits[itr]->WireID().Plane == 2)
                                    {
                                      fshwrhits_chno[fshwrcolhits]      =   shwrhits[itr]->Channel();
                                      fshwrhits_peakT[fshwrcolhits]     =   shwrhits[itr]->PeakTime();
                                      fshwrhits_charge[fshwrcolhits]    =   shwrhits[itr]->Integral();///(TMath::Exp(-(hits_peakT[longcolhits]-800)*500/tau)); //multiplied by 500nsec to convert time ticks to actual generation time
                                      fshwrhits_wire[fshwrcolhits]      =   shwrhits[itr]->WireID().Wire;
                                      fshwrhits_plane[fshwrcolhits]     =   shwrhits[itr]->WireID().Plane;
                                      fshwrhits_TPC[fshwrcolhits]       =   shwrhits[itr]->WireID().TPC;
                                      double diffpeaktt0shwr           =    shwrhits[itr]->PeakTime() - (T0_value0/1000)*2; //calculate the peak time in ticks
                                      fshwrhits_xpos[fshwrcolhits]      =   detProp.ConvertTicksToX(diffpeaktt0shwr,shwrhits[itr]->WireID().Plane,shwrhits[itr]->WireID().TPC,shwrhits[itr]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
                                      double Wirestartshwr[3], Wireendshwr[3];
                                      fGeom->WireEndPoints(shwrhits[itr]->WireID());
                                      // fGeom->WireEndPoints(shwrhits[itr]->WireID().Cryostat, shwrhits[itr]->WireID().TPC, shwrhits[itr]->WireID().Plane, shwrhits[itr]->WireID().Wire, Wirestartshwr, Wireendshwr);
                                      fshwrhits_zpos[fshwrcolhits]      =   Wirestartshwr[2];
                                      if(hits_spacept.at(shwrhits[itr].key()).size())
                                      fshwrhits_ypos[fshwrcolhits]      =   hits_spacept.at(shwrhits[itr].key())[0]->XYZ()[1];
                                      if(!hits_spacept.at(shwrhits[itr].key()).size())
                                      fshwrhits_ypos[fshwrcolhits]      =   303.5;
                                      fshwrhits_mult[fshwrcolhits]      =   shwrhits[itr]->Multiplicity();
                                      fshwrhits_sigptime[fshwrcolhits]  =   shwrhits[itr]->SigmaPeakTime();
                                      fshwrhits_sigchrg[fshwrcolhits]   =   shwrhits[itr]->SigmaIntegral();
                                      fshwrhits_sigpamp[fshwrcolhits]    =   shwrhits[itr]->SigmaPeakAmplitude();
                                      fshwrhits_dof[fshwrcolhits]       =   shwrhits[itr]->DegreesOfFreedom();
                                      fshwrhits_gof[fshwrcolhits]       =   shwrhits[itr]->GoodnessOfFit();
                                      fshwrhits_ptminusRMS[fshwrcolhits]	  = shwrhits[itr]->PeakTimeMinusRMS(5.0);
                                      fshwrhits_ptplusRMS[fshwrcolhits]	 = shwrhits[itr]->PeakTimePlusRMS(5.0);

                                      fshwrcolhits++;
                                    }// if(shwrhits[itr]->WireID().Plane == 2)
                                  }// for(itr = shwrhits.begin(); itr < shwrhits.end(); itr++)
                              
                                  //-------------Filling number of nearby hits--------------------//
                                  fnearhitct = nearhitct;
            
                                  for(int lp=0; lp<fnearhitct; lp++)
                                  {
                                    fnearhits_key[lp]     = nearhits_key[lp];
                                    fnearhits_peakT[lp]   = nearhits_peakT[lp];
                                    fnearhits_charge[lp]  = nearhits_charge[lp];
                                    fnearhits_wire[lp]    = nearhits_wire[lp];
                                    fnearhits_chno[lp]    = nearhits_chno[lp];
                                    fnearhits_TPC[lp]     = nearhits_TPC[lp];
                                    fnearhits_plane[lp]   = nearhits_plane[lp];
                                    fnearhits_xpos[lp]    = nearhits_xpos[lp];
                                    fnearhits_ypos[lp]    = nearhits_ypos[lp];			   
                                    fnearhits_zpos[lp]    = nearhits_zpos[lp];
                                    fnearhits_mult[lp]    = nearhits_mult[lp];
                                    fnearhits_sigptime[lp]= nearhits_sigptime[lp];
                                    fnearhits_sigchrg[lp] = nearhits_sigchrg[lp];
                                    fnearhits_sigpamp[lp]  = nearhits_sigpamp[lp];
                                    fnearhits_dof[lp]     = nearhits_dof[lp];
                                    fnearhits_gof[lp]     = nearhits_gof[lp];
                                    fnearhits_ptminusRMS[lp]= nearhits_ptminusRMS[lp];
                                    fnearhits_ptplusRMS[lp] = nearhits_ptplusRMS[lp];
                                    fnearhits_cnnMichel[lp] = nearhits_cnnMichel[lp];
                                    fnearhits_cnnEM[lp]     = nearhits_cnnEM[lp];
                                    fnearhits_cnnTrack[lp]  = nearhits_cnnTrack[lp];
                                  }// for(int lp=0; lp<fnearhitct; lp++)
                                }// if (abs(shwr_dist0)< _shwr_dist && MinHitPeakTimeshwrcol > _minhitpeakcut  && MaxHitPeakTimeshwrcol < _maxhitpeakcut)  
                              
                                //=======================================================================//
                                //------------------------ CALORIMETRY TIME -----------------------------//
                                //=======================================================================//
                                std::vector<art::Ptr<anab::Calorimetry>> calos=tracks_caloSCE.at(NTracks);                         
                                std::cout<<"calos.size() "<<calos.size()<<std::endl;
                                fhitsU = 0;
                                fhitsV = 0;
                                fhitsY = 0;

                                for(size_t ical = 0; ical<calos.size(); ++ical)
                                {
                                  if(!calos[ical]) continue;
                                  if(!calos[ical]->PlaneID().isValid) continue;

                                  int planenum = calos[ical]->PlaneID().Plane;
                                  if(planenum<0 || planenum>2) continue;
                                  const size_t NHits = calos[ical] -> dEdx().size();
                                  fntrkhits=int(NHits);
                                  for(size_t iHit = 0; iHit < NHits; ++iHit)
                                  {
                                    const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
                                    if(planenum == 0)
                                    {
                                      ftrkdqdxU[fhitsU]=(calos[ical] -> dQdx())[iHit];
                                      ftrkdedxU[fhitsU]=(calos[ical] -> dEdx())[iHit];
                                      ftrkresrangeU[fhitsU]=(calos[ical]->ResidualRange())[iHit];
                                      ftrkhitxU[fhitsU]=TrkPos.X();
                                      ftrkhityU[fhitsU]=TrkPos.Y();
                                      ftrkhitzU[fhitsU]=TrkPos.Z();
                                      ftrkpitchU[fhitsU]=(calos[ical]->TrkPitchVec())[iHit];

                                      fhitsU++;
                                    }// if(planenum == 0)
                                    if(planenum == 1)
                                    {
                                      ftrkdqdxV[fhitsV]=(calos[ical] -> dQdx())[iHit];
                                      ftrkdedxV[fhitsV]=(calos[ical] -> dEdx())[iHit];
                                      ftrkresrangeV[fhitsV]=(calos[ical]->ResidualRange())[iHit];
                                      ftrkhitxV[fhitsV]=TrkPos.X();
                                      ftrkhityV[fhitsV]=TrkPos.Y();
                                      ftrkhitzV[fhitsV]=TrkPos.Z();
                                      ftrkpitchV[fhitsV]=(calos[ical]->TrkPitchVec())[iHit];
                                      fhitsV++;
                                    }// if(planenum == 1)
                                    if(planenum == 2)
                                    {
                                      ftrkdqdxY[fhitsY]=(calos[ical] -> dQdx())[iHit];
                                      ftrkdedxY[fhitsY]=(calos[ical] -> dEdx())[iHit];
                                      ftrkresrangeY[fhitsY]=(calos[ical]->ResidualRange())[iHit];
                                      ftrkhitxY[fhitsY]=TrkPos.X();
                                      ftrkhityY[fhitsY]=TrkPos.Y();
                                      ftrkhitzY[fhitsY]=TrkPos.Z();
                                      ftrkpitchY[fhitsY]=(calos[ical]->TrkPitchVec())[iHit];
                                      fhitsY++;
                                    }// if(planenum == 2)
                                  } // loop over iHit.. //for(size_t iHit = 0; iHit < NHits; ++iHit)
                                }//for(size_t ical = 0; ical<calos.size(); ++ical)
                                
                                //---------------------------Taking the midpoint of the linear nearby hits----------//
                                int counter =0;
                                int zpos=0, zpossq=0;
                                double xposzpos=0, xpos=0;
                                double slope=-999, yintercept=-999;
                                for(int i=0; i <fnearhitct; i++)
                                {
                                  counter++;
                                  xposzpos += fnearhits_xpos[NTracks]*fnearhits_zpos[NTracks];
                                  xpos += fnearhits_xpos[NTracks];
                                  zpos += fnearhits_zpos[NTracks];
                                  zpossq += fnearhits_zpos[NTracks]*fnearhits_zpos[NTracks];  
                                } // for(int i=0; i <fnearhitct; i++) 
                                slope = (xposzpos*counter - xpos*zpos)/(zpossq*counter - zpos*zpos);
                                yintercept = (zpossq*xpos - xposzpos*zpos)/(zpossq*counter - zpos*zpos);
                                
                                float xposend = slope*fnearhits_zpos[int(fnearhitct/2)]+yintercept;

                                /********************* For candidate muon truth information *******************************/
                                std::vector <double> trueHitsKey;
                                if(!evt.isRealData())
                                {
                                  // Get true MCParticle associated with recob::Track
                                  const simb::MCParticle *particleP = truthUtil.GetMCParticleFromRecoTrack(clockData, track,evt,label_tracks);
                                  if(!particleP) continue;
                                  const art::Ptr<simb::MCTruth> mcP=pi_serv->TrackIdToMCTruth_P(particleP->TrackId());
                                  if(!mcP) continue;

                                  TLorentzVector mcstart, mcend, mcstartdrifted, mcenddrifted;
                                  unsigned int pstarti, pendi, pstartdriftedi, penddriftedi; //mcparticle indices for starts and ends in tpc or drifted volumes
                                  double plen = utils->length(*particleP, mcstart, mcend, pstarti, pendi);
                                  double plendrifted = utils->driftedLength(detProp, *particleP, mcstartdrifted, mcenddrifted, pstartdriftedi, penddriftedi);

                                  // bool isActive = plen != 0;
                                  bool isDrifted = plendrifted!= 0;

                                  unsigned int iniPointInAV = utils->TrueParticleIniPointInAV(fiducialBounds,*particleP);
                                  unsigned int endPointInAV = utils->TrueParticleEndPointInAV(fiducialBounds,*particleP);
                                  std::cout<<"track id "<<NTracks<<" bt sim trackid "<<particleP->TrackId()<<std::endl;

                                  fmccut_trkid      = particleP->TrackId();
                                  fmccut_vx         = particleP->Vx(iniPointInAV);
                                  fmccut_vy         = particleP->Vy(iniPointInAV);
                                  fmccut_vz         = particleP->Vz(iniPointInAV);
                                  fmccut_t          = particleP->T(iniPointInAV);
                                  fmccut_endx       = particleP->Vx(endPointInAV);
                                  fmccut_endy       = particleP->Vy(endPointInAV);
                                  fmccut_endz       = particleP->Vz(endPointInAV);
                                  fmccut_endt       = particleP->T(endPointInAV);
                                  fmccut_px         = particleP->Px();
                                  fmccut_py         = particleP->Py();
                                  fmccut_pz         = particleP->Pz();
                                  fmccut_momentum   = particleP->P();
                                  fmccut_energy     = particleP->E();
                                  fmccut_endpx      = particleP->EndPx();
                                  fmccut_endpy      = particleP->EndPy();
                                  fmccut_endpz      = particleP->EndPz();
                                  fmccut_endenergy  = particleP->EndE();
                                  fmccut_pathlen    = plen;
                                  fmccut_length     = truthUtil.GetMCParticleLengthInTPCActiveVolume(*particleP,fiducialBounds[0],fiducialBounds[1],fiducialBounds[2],fiducialBounds[3],fiducialBounds[4],fiducialBounds[5]);

                                  if (isDrifted)
                                  {
                                    fmccut_vxdrifted         = particleP->Vx(iniPointInAV);
                                    fmccut_vydrifted         = particleP->Vy(iniPointInAV);
                                    fmccut_vzdrifted         = particleP->Vz(iniPointInAV);
                                    fmccut_tdrifted          = particleP->T(iniPointInAV);
                                    fmccut_endxdrifted       = particleP->Vx(endPointInAV);
                                    fmccut_endydrifted       = particleP->Vy(endPointInAV);
                                    fmccut_endzdrifted       = particleP->Vz(endPointInAV);
                                    fmccut_endtdrifted       = particleP->T(endPointInAV);
                                    fmccut_pxdrifted         = particleP->Px();
                                    fmccut_pydrifted         = particleP->Py();
                                    fmccut_pzdrifted         = particleP->Pz();
                                    fmccut_momentumdrifted   = particleP->P();
                                    fmccut_energydrifted     = particleP->E();
                                    fmccut_endpxdrifted      = particleP->EndPx();
                                    fmccut_endpydrifted      = particleP->EndPy();
                                    fmccut_endpzdrifted      = particleP->EndPz();
                                    fmccut_endenergydrifted  = particleP->EndE();
                                    fmccut_pathlendrifted    = plendrifted;
                                  }// if (isDrifted)
                                  fmccut_endprocess = int(particleP->Trajectory().ProcessToKey(particleP->EndProcess()));
                                  fmccut_theta      = particleP->Momentum().Theta();
                                  fmccut_phi        = particleP->Momentum().Phi();
                                  fmccut_pdg        = particleP->PdgCode();
                                  fmccut_status_code= particleP->StatusCode();
                                  fmccut_mass       = particleP->Mass();
                                  fmccut_ND         = particleP->NumberDaughters();
                                  fmccut_mother     = particleP->Mother();
                                  fmccut_origin     = mcP->Origin();
                                  fmccut_process    = int(particleP->Trajectory().ProcessToKey(particleP->Process()));
                                  fmccut_rescatter  = particleP->Rescatter();
                            
                                  std::cout<<"true cut mu: Pdg code "<<fmccut_pdg<<" track ID "<<fmccut_trkid<<" no of daughters "<<fmccut_ND<<" and origin_type "<<fmccut_origin<<" fmccut_endx "<<fmccut_endx<<" fmccut_endy "<<fmccut_endy<<" fmccut_endz "<<fmccut_endz<<std::endl;
                          
                                  fmcd_trkid      = mcd_trkid;
                                  fmcd_vx         = mcd_vx;
                                  fmcd_vy         = mcd_vy;
                                  fmcd_vz         = mcd_vz;
                                  fmcd_t          = mcd_t;
                                  fmcd_endx       = mcd_endx;
                                  fmcd_endy       = mcd_endy;
                                  fmcd_endz       = mcd_endz;
                                  fmcd_endt	      = mcd_endt;
                                  fmcd_px         = mcd_px;
                                  fmcd_py         = mcd_py;
                                  fmcd_pz         = mcd_pz;
                                  fmcd_momentum   = mcd_momentum;
                                  fmcd_energy     = mcd_energy;
                                  fmcd_endpx      = mcd_endpx;
                                  fmcd_endpy      = mcd_endpy;
                                  fmcd_endpz      = mcd_endpz;
                                  fmcd_endenergy  = mcd_endenergy;
                                  fmcd_pathlen    = mcd_pathlen;
                                  fmcd_vxdrifted         = mcd_vxdrifted;
                                  fmcd_vydrifted         = mcd_vydrifted;
                                  fmcd_vzdrifted         = mcd_vzdrifted;
                                  fmcd_tdrifted          = mcd_tdrifted;
                                  fmcd_endxdrifted       = mcd_endxdrifted;
                                  fmcd_endydrifted       = mcd_endydrifted;
                                  fmcd_endzdrifted       = mcd_endzdrifted;
                                  fmcd_endtdrifted       = mcd_endtdrifted;
                                  fmcd_pxdrifted         = mcd_pxdrifted;
                                  fmcd_pydrifted	       = mcd_pydrifted;
                                  fmcd_pzdrifted         = mcd_pzdrifted;
                                  fmcd_momentumdrifted   = mcd_momentumdrifted;
                                  fmcd_energydrifted     = mcd_energydrifted;
                                  fmcd_endpxdrifted      = mcd_endpxdrifted;
                                  fmcd_endpydrifted      = mcd_endpydrifted;
                                  fmcd_endpzdrifted      = mcd_endpzdrifted;
                                  fmcd_endenergydrifted  = mcd_endenergydrifted;
                                  fmcd_pathlendrifted    = mcd_pathlendrifted;
                                  fmcd_endprocess = mcd_endprocess;
                                  fmcd_theta      = mcd_theta;
                                  fmcd_phi        = mcd_phi;
                                  fmcd_pdg        = mcd_pdg;
                                  fmcd_status_code= mcd_status_code;
                                  fmcd_mass       = mcd_mass;
                                  fmcd_ND         = mcd_ND;
                                  fmcd_mother     = mcd_mother;
                                  fmcd_origin     = mcd_origin;
                                  fmcd_process    = mcd_process;
                                  fmcd_rescatter  = mcd_rescatter;
                                      
                                  //-----------Getting the average y position of the Michel track-----------------//
                                  double yposMi = 0; int countMi = 0;
                                  int trueparhitscol_key[500] = {0};

                                  std::vector<art::Ptr<recob::Hit>> hitsfromMi = bt_serv->TrackIdToHits_Ps(clockData, fmcd_trkid, hitVect);
                                  for(size_t e=0;e<hitsfromMi.size(); ++e)
                                  {			    
                                    if(hitsfromMi[e]->WireID().Plane==2)
                                    {
                                      if(!hits_spacept.at(hitsfromMi[e].key()).size()) continue;
                                      
                                      trueparhitscol_key[countMi]   = hitsfromMi[e].key();
                                      yposMi += hits_spacept.at(hitsfromMi[e].key())[0]->XYZ()[1];
                                      countMi ++;
                                    }// if(hitsfromMi[e]->WireID().Plane==2)
                                  }// for(size_t e=0;e<hitsfromMi.size(); ++e) 

                                  std::vector<art::Ptr<recob::Hit>> hitsfromMi1 = bt_serv->TrackIdToHits_Ps(clockData, -(fmcd_trkid), hitVect);
                                  for(size_t e=0;e<hitsfromMi1.size(); ++e)
                                  {			    
                                    int found_dupMi = 0;
                                    if(hitsfromMi1[e]->WireID().Plane==2)
                                    {
                                      for(int kk=0; kk<countMi; kk++)
                                      {
                                        trueparhitscol_key[countMi] = hitsfromMi1[e].key();
                                        if(trueparhitscol_key[kk]==trueparhitscol_key[countMi]) found_dupMi = 1;
                                      }// for(int kk=0; kk<countMi; kk++)
                                      
                                      if(found_dupMi == 0)
                                      {
                                        if(!hits_spacept.at(hitsfromMi1[e].key()).size()) continue;
                                        yposMi += hits_spacept.at(hitsfromMi1[e].key())[0]->XYZ()[1];
                                        countMi ++;
                                      }// if(found_dupMi == 0)
                                    }// if(hitsfromMi1[e]->WireID().Plane==2)
                                  }// for(size_t e=0;e<hitsfromMi1.size(); ++e)
                                  
                                  float avgYpos = float(yposMi)/int(countMi) ;   
                                  std::cout<<"avgMichelYpos "<<avgYpos<<std::endl;

                                  //-----------Tyring to solve overlapping hits issues -----------------//
                                  int hitmultptimeminus[20000] = {0};
                                  int hitmultptimeplus[20000] = {0};
                                  int prevptimeminus = 6000;
                                  int prevptimeplus = 0;
                                  int prevchanno = 0;
                                  int countii = 0;
                                  std::vector<int> mult_key; 
                                  
                                  int MitrackID_key[500] = {0}; int ctMi = 0;
                                    
                                  std::vector<art::Ptr<recob::Hit>> hitsfromPar = bt_serv->TrackIdToHits_Ps(clockData, fmcd_trkid, hitVect);
                                  for(size_t e=0;e<hitsfromPar.size(); ++e)
                                  {			    
                                    if(hitsfromPar[e]->WireID().Plane==2 && hitsfromPar[e]->Multiplicity()>1)
                                    {
                                      MitrackID_key[ctMi]   = hitsfromPar[e].key();
                                      ctMi++;
                                  
                                      int currchanno = hitsfromPar[e]->Channel();
                                      if(countii==0 || currchanno !=prevchanno)
                                      {
                                        //  cout<<"inside loop"<<std::endl;
                                        prevptimeminus = hitsfromPar[e]->PeakTimeMinusRMS(3);
                                        prevptimeplus = hitsfromPar[e]->PeakTimePlusRMS(3);
                                        prevchanno = hitsfromPar[e]->Channel();
                                        hitmultptimeminus[hitsfromPar[e]->Channel()] = int(hitsfromPar[e]->PeakTimeMinusRMS(3));
                                        hitmultptimeplus[hitsfromPar[e]->Channel()] = int(hitsfromPar[e]->PeakTimePlusRMS(3));
                                    
                                      } // if(countii==0 || currchanno !=prevchanno)
                                    
                                      if(hitsfromPar[e]->PeakTimeMinusRMS(3) < prevptimeminus)				           				         hitmultptimeminus[hitsfromPar[e]->Channel()] = int(hitsfromPar[e]->PeakTimeMinusRMS(3));
                                      
                                      if(hitsfromPar[e]->PeakTimePlusRMS(3) > prevptimeplus)
                                      hitmultptimeplus[hitsfromPar[e]->Channel()] = int(hitsfromPar[e]->PeakTimePlusRMS(3));
                                      
                                      prevptimeminus = hitsfromPar[e]->PeakTimeMinusRMS(3);
                                      prevptimeplus = hitsfromPar[e]->PeakTimePlusRMS(3);
                                      prevchanno = hitsfromPar[e]->Channel();
                                      
                                      mult_key.push_back(hitsfromPar[e].key());
                                  
                                      // cout<<"hitkey "<<hitsfromPar[e].key()<<" prevptimeminus "<<prevptimeminus<<" prevptimeplus "<<prevptimeplus<<" prevchanno "<<prevchanno<<" hitmultptimeminus["<<hitsfromPar[e]->Channel()<<"] "<<hitmultptimeminus[hitsfromPar[e]->Channel()]<<" hitmultptimeplus["<<hitsfromPar[e]->Channel()<<"] "<<hitmultptimeplus[hitsfromPar[e]->Channel()]<<std::endl;
                                  
                                      countii++;
                                    }// if (hitsfromPar[e]->WireID().Plane==2 && hitsfromPar[e]->Multiplicity()>1)
                                  }// for (size_t e=0;e<hitsfromPar.size(); ++e)
                                    
                                  //Now from negative particle ID 
                                    
                                  std::vector<art::Ptr<recob::Hit>> hitsfromPar1 = bt_serv->TrackIdToHits_Ps(clockData, -(fmcd_trkid), hitVect);
                                    
                                  for(size_t e=0;e<hitsfromPar1.size(); ++e)
                                  {
                                    int found_dupMi = 0;		    
                                    if(hitsfromPar1[e]->WireID().Plane==2 && hitsfromPar1[e]->Multiplicity()>1)
                                    {
                                      for(int kk=0; kk<ctMi; kk++)
                                      {
                                        MitrackID_key[ctMi] = hitsfromPar1[e].key();
                                        if(MitrackID_key[kk]==MitrackID_key[ctMi]) found_dupMi = 1;
                                      }
                                      
                                      if(found_dupMi == 0)
                                      {
                                        int currchanno = hitsfromPar1[e]->Channel();
                                        if(countii==0 || currchanno !=prevchanno)
                                        {
                                          // cout<<"inside neg track ID loop"<<std::endl;
                                          prevptimeminus = hitsfromPar1[e]->PeakTimeMinusRMS(3);
                                          prevptimeplus = hitsfromPar1[e]->PeakTimePlusRMS(3);
                                          prevchanno = hitsfromPar1[e]->Channel();
                                          hitmultptimeminus[hitsfromPar1[e]->Channel()] = int(hitsfromPar1[e]->PeakTimeMinusRMS(3));
                                          hitmultptimeplus[hitsfromPar1[e]->Channel()] = int(hitsfromPar1[e]->PeakTimePlusRMS(3));
                                          // ctMi++;
                                        }  //if(countii==0 || currchanno !=prevchanno)
                                    
                                        if(hitsfromPar1[e]->PeakTimeMinusRMS(3) < prevptimeminus)				           				 hitmultptimeminus[hitsfromPar1[e]->Channel()] = int(hitsfromPar1[e]->PeakTimeMinusRMS(3));
                                        if(hitsfromPar1[e]->PeakTimePlusRMS(3) > prevptimeplus)
                                        hitmultptimeplus[hitsfromPar1[e]->Channel()] = int(hitsfromPar1[e]->PeakTimePlusRMS(3));
                                        prevptimeminus = hitsfromPar1[e]->PeakTimeMinusRMS(3);
                                        prevptimeplus = hitsfromPar1[e]->PeakTimePlusRMS(3);
                                        prevchanno = hitsfromPar1[e]->Channel();
                                        mult_key.push_back(hitsfromPar1[e].key());
                                  
                                        // cout<<"hitkey "<<hitsfromPar1[e].key()<<" prevptimeminus "<<prevptimeminus<<" prevptimeplus "<<prevptimeplus<<" prevchanno "<<prevchanno<<" hitmultptimeminus["<<hitsfromPar1[e]->Channel()<<"] "<<hitmultptimeminus[hitsfromPar1[e]->Channel()]<<" hitmultptimeplus["<<hitsfromPar1[e]->Channel()<<"] "<<hitmultptimeplus[hitsfromPar1[e]->Channel()]<<std::endl;
                                  
                                        countii++;
                                      }//if(found_dupMi == 0)
                                    }//if(hitsfromPar1[e]->WireID().Plane==2...)
                                  }//for(size_t e=0;e<hitsfromPar1.size(); ++e)  
                                  
                                  std::vector<int> mult_chan;  
                                  std::vector<int> mult_ptimeminus;
                                  std::vector<int> mult_ptimeplus;
                                  std::vector<float> mult_trueE;
                                  
                                  //----Now counting energy of hits having higher multiplicity-------//
                                  for(int ii = 0; ii< 20000; ii++)
                                  {
                                    if(hitmultptimeminus[ii]!=0)
                                    {
                                      std::cout<<"multhit: channal "<<ii<<" ptimeminus "<<hitmultptimeminus[ii]<<" ptimeplus "<<hitmultptimeplus[ii]<<std::endl;
                                      mult_chan.push_back(ii);
                                      mult_ptimeminus.push_back(hitmultptimeminus[ii]);
                                      mult_ptimeplus.push_back(hitmultptimeplus[ii]);

                                      fReadOutWindowSize = detProp.ReadOutWindowSize();
                                      fNumberTimeSamples = detProp.NumberTimeSamples(); 	       
                                      art::Handle< std::vector<sim::SimChannel> > simchannelHandle;
                                      float totEmult = 0;	   

                                      if(evt.getByLabel("largeant", simchannelHandle))
                                      {
                                        //Loop over simchannels 
                                        for(auto const& simchannel : (*simchannelHandle))
                                        {
                                          if(fGeom->View(simchannel.Channel()) != 2) continue;	
                                          auto const& alltimeslices = simchannel.TDCIDEMap();
                              
                                          // Loop over ticks 6000
                                          for(auto const& tslice : alltimeslices)
                                          {
                                            int spill = 1;
                                            if (fReadOutWindowSize != fNumberTimeSamples) 
                                            {
                                              if (tslice.first < spill * fReadOutWindowSize || tslice.first > (spill + 1) * fReadOutWindowSize) continue;
                                            }
                                            else if (tslice.first < 0 || tslice.first > fReadOutWindowSize) continue;		 
                                            auto const& simide = tslice.second;
                                            int sim_channel = int(simchannel.Channel());
                                            // Loop over energy deposits
                                            for(auto const& eDep : simide)
                                            {
                                              if(sim_channel == ii && (eDep.trackID == fmcd_trkid || eDep.trackID == -fmcd_trkid))
                                              if(tslice.first >= hitmultptimeminus[ii] && tslice.first <= hitmultptimeplus[ii])
                                              {	
                                                totEmult += eDep.energy;
                                                std::cout<<"totEmult "<<totEmult<<" tslice.first "<<tslice.first<<" MeV"<<std::endl;	
                                              }
                                            }  //for(auto const& eDep : simide)
                                          } //for(auto const& tslice : alltimeslices)  
                                        }//for(auto const& simchannel : (*simchannelHandle))
                                      }//if(evt.getByLabel("largeant", simchannelHandle))
                                      
                                      mult_trueE.push_back(totEmult);

                                    }//if(hitmultptimeminus[ii]!=0)
                                  }//for(int ii = 0; ii< 20000; ii++)    


                                  //-----------getting a subset of the reco hits that are matched to MC particles listed in trkIDs-----------------//
                                  std::vector<float> hit_key_fromHits;
                                  std::vector<float> chan_charge_fromHits;
                                  std::vector<int> chan_no_fromHits;
                                  std::vector<float> energy_fromHits;
                                  std::vector<int> ptime_minus_rms_fromHits;
                                  std::vector<int> ptime_plus_rms_fromHits;
                                  std::vector<float> ptime_startTick_fromHits;
                                  std::vector<float> ptime_endTick_fromHits;

                                  // std::vector<art::Ptr<recob::Hit>> hitsfromPar = bt_serv->TrackIdToHits_Ps(clockData, fmcd_trkid, hitVect);
                                  std::map<int,double> mtrkide1;
                                  int mtrackid1=-1;

                                  ftrueparhitallcount = 0, ftrueparhitcolcount = 0;
                                  for(size_t e=0;e<hitsfromPar.size(); ++e)
                                  { 
                                    double truepardiffpeaktt0 = hitsfromPar[e]->PeakTime() - (T0_value0/1000)*2;
                                    double trueparhitXpos = detProp.ConvertTicksToX(truepardiffpeaktt0,hitsfromPar[e]->WireID().Plane,hitsfromPar[e]->WireID().TPC,hitsfromPar[e]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
                                    double trueparWirestart[3], trueparWireend[3];
                                    fGeom->WireEndPoints(hitsfromPar[e]->WireID());
                                    // fGeom->WireEndPoints(hitsfromPar[e]->WireID().Cryostat, hitsfromPar[e]->WireID().TPC, hitsfromPar[e]->WireID().Plane, hitsfromPar[e]->WireID().Wire, trueparWirestart, trueparWireend);
                                    double trueparhitZpos = trueparWirestart[2];
                                    //                             if(!hits_spacept.at(hitsfromPar[e].key()).size()) continue;
                                    //			     double trueparhitYpos      =   hits_spacept.at(hitsfromPar[e].key())[0]->XYZ()[1];
                                    double trueparhitYpos = 0;
                                    if(hits_spacept.at(hitsfromPar[e].key()).size())
                                    trueparhitYpos      =   hits_spacept.at(hitsfromPar[e].key())[0]->XYZ()[1];
                                    if(!hits_spacept.at(hitsfromPar[e].key()).size()) trueparhitYpos = avgYpos;


                                    double magnearhitvec = sqrt(pow(xposend-fcut_endhitx,2)+pow(fnearhits_zpos[int(fnearhitct/2)]-fcut_endhitz,2));
                                    double maghitvec = sqrt(pow(trueparhitXpos-fcut_endhitx,2)+pow(trueparhitZpos-fcut_endhitz,2));
                                    double angle_theta = acos(((xposend-fcut_endhitx)*(trueparhitXpos-fcut_endhitx) + (fnearhits_zpos[int(fnearhitct/2)]-fcut_endhitz)*(trueparhitZpos-fcut_endhitz))/(magnearhitvec*maghitvec));
                                    double angle_theta_deg = (angle_theta*180)/3.14;
                                    double maghitveccostheta = maghitvec* cos(angle_theta);
                                    double distance = sqrt(pow(trueparhitXpos-fcut_trkendx,2)+pow(trueparhitZpos-fcut_trkendz,2));
                                
                                    ftrueparhitsall_key[ftrueparhitallcount]   = hitsfromPar[e].key();
                                    ftrueparhitsall_peakT[ftrueparhitallcount] = hitsfromPar[e]->PeakTime();
                                    ftrueparhitsall_charge[ftrueparhitallcount]= hitsfromPar[e]->Integral();
                                    ftrueparhitsall_wire[ftrueparhitallcount]  = hitsfromPar[e]->WireID().Wire;
                                    ftrueparhitsall_chno[ftrueparhitallcount]  = hitsfromPar[e]->Channel();
                                    ftrueparhitsall_TPC[ftrueparhitallcount]   = hitsfromPar[e]->WireID().TPC;
                                    ftrueparhitsall_plane[ftrueparhitallcount] = hitsfromPar[e]->WireID().Plane;
                                    ftrueparhitsall_xpos[ftrueparhitallcount]  = trueparhitXpos;
                                    ftrueparhitsall_ypos[ftrueparhitallcount]  = trueparhitYpos;
                                    ftrueparhitsall_zpos[ftrueparhitallcount]  = trueparhitZpos;
                                    ftrueparhitsall_mult[ftrueparhitallcount]      = hitsfromPar[e]->Multiplicity();
                                    ftrueparhitsall_sigptime[ftrueparhitallcount]  = hitsfromPar[e]->SigmaPeakTime();
                                    ftrueparhitsall_sigchrg[ftrueparhitallcount]   = hitsfromPar[e]->SigmaIntegral();
                                    ftrueparhitsall_sigpamp[ftrueparhitallcount]    = hitsfromPar[e]->SigmaPeakAmplitude();
                                    ftrueparhitsall_dof[ftrueparhitallcount]       = hitsfromPar[e]->DegreesOfFreedom();
                                    ftrueparhitsall_gof[ftrueparhitallcount]       = hitsfromPar[e]->GoodnessOfFit();
                                    ftrueparhitsall_ptminusRMS[ftrueparhitallcount]	  = hitsfromPar[e]->PeakTimeMinusRMS(5.0);
                                    ftrueparhitsall_ptplusRMS[ftrueparhitallcount]	 = hitsfromPar[e]->PeakTimePlusRMS(5.0);

                                    ftrueparhitallcount++;
                                    if(hitsfromPar[e]->WireID().Plane==2)// && hitsfromPar[e]->Multiplicity()<2)
                                    {
                                    
                                      //--------look at true info of this hit--------//
                                      
                                      art::Ptr<recob::Hit> mhit1=hitVect[hitsfromPar[e].key()];
                                                  std::vector<sim::TrackIDE> meveIDs1 = bt_serv->HitToTrackIDEs(clockData,mhit1);
                                                        
                                      std::cout<<std::endl;
                                      // cout<<"key "<<hitsfromPar[e].key()<<" angle_theta_deg "<<angle_theta_deg<<" maghitveccostheta "<<maghitveccostheta<<" mult "<<hitsfromPar[e]->Multiplicity()<<" charge "<<hitsfromPar[e]->Integral()<<std::endl;
                                      float maxen1 = -1; int maxenid1  =0;  float maxenfrac1=-1;
                                      std::map<int,int> mtrkid11;
                                      std::map<int,double> mtrkide11;
                                      std::map<int,double> mtrkidefrac11;			    

                                      for(size_t e=0;e<meveIDs1.size(); ++e)
                                      { 
                                        // mtrkide1[abs(meveIDs1[e].trackID)] += meveIDs1[e].energy;
                                        mtrkid11[abs(meveIDs1[e].trackID)] = abs(meveIDs1[e].trackID);
                                        mtrkide11[abs(meveIDs1[e].trackID)] += meveIDs1[e].energy;
                                        mtrkidefrac11[abs(meveIDs1[e].trackID)] += meveIDs1[e].energyFrac;			
                                        std::cout<<"trackID "<<meveIDs1[e].trackID<<" energy "<<meveIDs1[e].energy<<" energyFrac "<<meveIDs1[e].energyFrac<<std::endl;
                                        if(mtrkide11[abs(meveIDs1[e].trackID)]>maxen1)
                                        {
                                          maxen1 = mtrkide11[abs(meveIDs1[e].trackID)];
                                          maxenid1 = abs(meveIDs1[e].trackID);
                                          maxenfrac1 = mtrkidefrac11[abs(meveIDs1[e].trackID)];
                                          std::cout<<"maxen1 "<<maxen1<<" maxenid1 "<<maxenid1<<" RealID "<<meveIDs1[e].trackID<<" maxenfrac1 "<<maxenfrac1<<std::endl;
                                        }// if(mtrkide11[abs(meveIDs1[e].trackID)]>maxen1)
                                      } //for(size_t e=0;e<meveIDs1.size(); ++e)
                                
                                      //see what other particles take away the energy fraction;
                                      for(std::map<int,double>::iterator ii = mtrkidefrac11.begin(); ii!=mtrkidefrac11.end(); ++ii)
                                      {              
                                        const simb::MCParticle *mparticleP1 = pi_serv->TrackIdToParticle_P(abs(ii->first));
                                              // if(!mparticleP) continue;
                                              // const art::Ptr<simb::MCTruth> mmcP=pi_serv->TrackIdToMCTruth_P(mtrackid);
                                              // if(!mmcP) continue;
                                        std::cout<<"TrackID "<<ii->first<<" energyFrac "<<ii->second<<" pdg "<<mparticleP1->PdgCode()<<" E "<<mparticleP1->E()<<std::endl;
                                      }// for(std::map<int,double>::iterator ii = mtrkidefrac11.begin(); ii!=mtrkidefrac11.end(); ++ii)
                              
                                      ftrueMiEFrac[ftrueparhitcolcount] = mtrkidefrac11[abs(fmcd_trkid)];

                                      if(hitsfromPar[e]->Multiplicity()<2)
                                      {
                                        //for counting true MC energy		    
                                        for(size_t e=0;e<meveIDs1.size(); ++e)
                                        { 
                                          mtrkide1[abs(meveIDs1[e].trackID)] += meveIDs1[e].energy;
                                          // std::cout<<"abs(meveIDs1[e].trackID) "<<abs(meveIDs1[e].trackID)<<" mtrkide1[abs(meveIDs1[e].trackID)] "<<mtrkide1[abs(meveIDs1[e].trackID)]<<std::endl;
                                        }// for(size_t e=0;e<meveIDs1.size(); ++e)
                                        ////////////////////////////////
                                        hit_key_fromHits.push_back(hitsfromPar[e].key());
                                        chan_charge_fromHits.push_back(hitsfromPar[e]->Integral());
                                        chan_no_fromHits.push_back(hitsfromPar[e]->Channel());
                                        energy_fromHits.push_back(maxen1);
                                        ptime_minus_rms_fromHits.push_back(hitsfromPar[e]->PeakTimeMinusRMS(3));
                                        ptime_plus_rms_fromHits.push_back(hitsfromPar[e]->PeakTimePlusRMS(3));
                                        ptime_startTick_fromHits.push_back(hitsfromPar[e]->StartTick());
                                        ptime_endTick_fromHits.push_back(hitsfromPar[e]->EndTick());
                                      }//if(hitsfromPar[e]->Multiplicity()<2)   

                                      if(abs(maxenid1) != abs(fmcd_trkid)) continue; // to make sure that the maximum energy deposition is coming from the true Michel 
                                
                                      ftrueparhitscol_key[ftrueparhitcolcount]   = hitsfromPar[e].key();
                                      ftrueparhitscol_peakT[ftrueparhitcolcount] = hitsfromPar[e]->PeakTime();
                                      ftrueparhitscol_charge[ftrueparhitcolcount]= hitsfromPar[e]->Integral();
                                      ftrueparhitscol_wire[ftrueparhitcolcount]  = hitsfromPar[e]->WireID().Wire;
                                      ftrueparhitscol_chno[ftrueparhitcolcount]  = hitsfromPar[e]->Channel();
                                      ftrueparhitscol_TPC[ftrueparhitcolcount]   = hitsfromPar[e]->WireID().TPC;
                                      ftrueparhitscol_plane[ftrueparhitcolcount] = hitsfromPar[e]->WireID().Plane;
                                      ftrueparhitscol_xpos[ftrueparhitcolcount]  = trueparhitXpos;
                                      ftrueparhitscol_ypos[ftrueparhitcolcount]  = trueparhitYpos;
                                      ftrueparhitscol_zpos[ftrueparhitcolcount]  = trueparhitZpos;
                                      ftrueparhitscol_angledeg[ftrueparhitcolcount] = angle_theta_deg;
                                      ftrueparhitscol_maghitveccostheta[ftrueparhitcolcount] = maghitveccostheta;
                                      ftrueparhitscol_distance[ftrueparhitcolcount] = distance;
                                      ftrueparhitscol_mult[ftrueparhitcolcount]      = hitsfromPar[e]->Multiplicity();
                                      ftrueparhitscol_sigptime[ftrueparhitcolcount]  = hitsfromPar[e]->SigmaPeakTime();
                                      ftrueparhitscol_sigchrg[ftrueparhitcolcount]   = hitsfromPar[e]->SigmaIntegral();
                                      ftrueparhitscol_sigpamp[ftrueparhitcolcount]    = hitsfromPar[e]->SigmaPeakAmplitude();
                                      ftrueparhitscol_dof[ftrueparhitcolcount]       = hitsfromPar[e]->DegreesOfFreedom();
                                      ftrueparhitscol_gof[ftrueparhitcolcount]       = hitsfromPar[e]->GoodnessOfFit();
                                      ftrueparhitscol_ptminusRMS[ftrueparhitcolcount]     = hitsfromPar[e]->PeakTimeMinusRMS(5.0);
                                      ftrueparhitscol_ptplusRMS[ftrueparhitcolcount]     = hitsfromPar[e]->PeakTimePlusRMS(5.0);

                                      trueHitsKey.push_back(hitsfromPar[e].key());

                                      ftrueparhitcolcount++;
                                    }// if(hitsfromPar[e]->WireID().Plane==2)
                                    // std::cout<<"hit key from true par "<<hitsfromPar[e].key()<<std::endl;
                                    //			     } //if(same_trk == 0) 
                                  } //for(size_t e=0;e<hitsfromPar.size(); ++e)

                                  std::cout<<"\n Now from -(fmcd_trkid)"<<std::endl;
                                  // 	                    std::vector<art::Ptr<recob::Hit>> hitsfromPar1 = bt_serv->TrackIdToHits_Ps(clockData, -(fmcd_trkid), hitVect);
                                  for(size_t e=0;e<hitsfromPar1.size(); ++e)
                                  { 
                                    double truepardiffpeaktt0 = hitsfromPar1[e]->PeakTime() - (T0_value0/1000)*2;
                                    double trueparhitXpos = detProp.ConvertTicksToX(truepardiffpeaktt0,hitsfromPar1[e]->WireID().Plane,hitsfromPar1[e]->WireID().TPC,hitsfromPar1[e]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
                                    double trueparWirestart[3], trueparWireend[3];
                                    fGeom->WireEndPoints(hitsfromPar1[e]->WireID());
                                    // fGeom->WireEndPoints(hitsfromPar1[e]->WireID().Cryostat, hitsfromPar1[e]->WireID().TPC, hitsfromPar1[e]->WireID().Plane, hitsfromPar1[e]->WireID().Wire, trueparWirestart, trueparWireend);
                                    double trueparhitZpos = trueparWirestart[2];
                                    //                             if(!hits_spacept.at(hitsfromPar1[e].key()).size()) continue;
                                    //			     double trueparhitYpos      =   hits_spacept.at(hitsfromPar1[e].key())[0]->XYZ()[1];
                                    double trueparhitYpos = 0;
                                    if(hits_spacept.at(hitsfromPar1[e].key()).size())
                                    trueparhitYpos      =   hits_spacept.at(hitsfromPar1[e].key())[0]->XYZ()[1];
                                    if(!hits_spacept.at(hitsfromPar1[e].key()).size()) trueparhitYpos = avgYpos;


                                    double magnearhitvec = sqrt(pow(xposend-fcut_endhitx,2)+pow(fnearhits_zpos[int(fnearhitct/2)]-fcut_endhitz,2));
                                    double maghitvec = sqrt(pow(trueparhitXpos-fcut_endhitx,2)+pow(trueparhitZpos-fcut_endhitz,2));
                                    double angle_theta = acos(((xposend-fcut_endhitx)*(trueparhitXpos-fcut_endhitx) + (fnearhits_zpos[int(fnearhitct/2)]-fcut_endhitz)*(trueparhitZpos-fcut_endhitz))/(magnearhitvec*maghitvec));
                                    double angle_theta_deg = (angle_theta*180)/3.14;
                                    double maghitveccostheta = maghitvec* cos(angle_theta);
                                    double distance = sqrt(pow(trueparhitXpos-fcut_trkendx,2)+pow(trueparhitZpos-fcut_trkendz,2));
                                
                                    ftrueparhitsall_key[ftrueparhitallcount]   = hitsfromPar1[e].key();

                                    int found_dup = 0;
                                    for(int kk=0; kk<ftrueparhitallcount; kk++)
                                    {
                                      if(ftrueparhitsall_key[ftrueparhitallcount]==ftrueparhitsall_key[kk]) found_dup = 1;
                                    }
                                      
                                    if(found_dup == 0)
                                    {
                                      ftrueparhitsall_peakT[ftrueparhitallcount] = hitsfromPar1[e]->PeakTime();
                                      ftrueparhitsall_charge[ftrueparhitallcount]= hitsfromPar1[e]->Integral();
                                      ftrueparhitsall_wire[ftrueparhitallcount]  = hitsfromPar1[e]->WireID().Wire;
                                      ftrueparhitsall_chno[ftrueparhitallcount]  = hitsfromPar1[e]->Channel();
                                      ftrueparhitsall_TPC[ftrueparhitallcount]   = hitsfromPar1[e]->WireID().TPC;
                                      ftrueparhitsall_plane[ftrueparhitallcount] = hitsfromPar1[e]->WireID().Plane;
                                      ftrueparhitsall_xpos[ftrueparhitallcount]  = trueparhitXpos;
                                      ftrueparhitsall_ypos[ftrueparhitallcount]  = trueparhitYpos;
                                      ftrueparhitsall_zpos[ftrueparhitallcount]  = trueparhitZpos;
                                      ftrueparhitsall_mult[ftrueparhitallcount]      = hitsfromPar1[e]->Multiplicity();
                                      ftrueparhitsall_sigptime[ftrueparhitallcount]  = hitsfromPar1[e]->SigmaPeakTime();
                                      ftrueparhitsall_sigchrg[ftrueparhitallcount]   = hitsfromPar1[e]->SigmaIntegral();
                                      ftrueparhitsall_sigpamp[ftrueparhitallcount]    = hitsfromPar1[e]->SigmaPeakAmplitude();
                                      ftrueparhitsall_dof[ftrueparhitallcount]       = hitsfromPar1[e]->DegreesOfFreedom();
                                      ftrueparhitsall_gof[ftrueparhitallcount]       = hitsfromPar1[e]->GoodnessOfFit();
                                      ftrueparhitsall_ptminusRMS[ftrueparhitallcount]	  = hitsfromPar1[e]->PeakTimeMinusRMS(5.0);
                                      ftrueparhitsall_ptplusRMS[ftrueparhitallcount]	 = hitsfromPar1[e]->PeakTimePlusRMS(5.0);
                                      
                                      ftrueparhitallcount++;
                                    }//if (found_dup == 0)


                                    if(hitsfromPar1[e]->WireID().Plane==2)// && hitsfromPar1[e]->Multiplicity()<2)
                                    {
                                      ftrueparhitscol_key[ftrueparhitcolcount]   = hitsfromPar1[e].key();
                                      int found_dup1 = 0;
                                      for(int kk=0; kk<ftrueparhitcolcount; kk++)
                                      {
                                        if(ftrueparhitscol_key[ftrueparhitcolcount]==ftrueparhitscol_key[kk]) found_dup1 = 1;
                                      }
                                      
                                      if(found_dup1 == 0)
                                      {
                                        //--------look at true info of this hit--------//
                                        art::Ptr<recob::Hit> mhit2=hitVect[hitsfromPar1[e].key()];
                                        std::vector<sim::TrackIDE> meveIDs2 = bt_serv->HitToTrackIDEs(clockData,mhit2);
                                        std::cout<<std::endl;
                                        // cout<<"key "<<hitsfromPar1[e].key()<<" angle_theta_deg "<<angle_theta_deg<<" maghitveccostheta "<<maghitveccostheta<<" mult "<<hitsfromPar1[e]->Multiplicity()<<" charge "<<hitsfromPar1[e]->Integral()<<std::endl;
                                                  
                                        float maxen2 = -1; int maxenid2  =0; float maxenfrac2=-1;
                                        std::map<int,double> mtrkide12;
                                        std::map<int,double> mtrkidefrac12;

                                        for(size_t e=0;e<meveIDs2.size(); ++e)
                                        { 
                                          // mtrkide1[abs(meveIDs2[e].trackID)] += meveIDs2[e].energy;
                                          mtrkide12[abs(meveIDs2[e].trackID)] += meveIDs2[e].energy;
                                          mtrkidefrac12[abs(meveIDs2[e].trackID)] += meveIDs2[e].energyFrac;			
                                          std::cout<<"trackID "<<meveIDs2[e].trackID<<" energy "<<meveIDs2[e].energy<<" energyFrac "<<meveIDs2[e].energyFrac<<std::endl;
                                          //		                cout<<"abs(meveIDs2[e].trackID) "<<abs(meveIDs2[e].trackID)<<" mtrkide12[abs(meveIDs2[e].trackID)] "<<mtrkide12[abs(meveIDs2[e].trackID)]<<std::endl;
                                          if(mtrkide12[abs(meveIDs2[e].trackID)]>maxen2)
                                          {
                                            maxen2 = mtrkide12[abs(meveIDs2[e].trackID)];
                                            maxenid2 = abs(meveIDs2[e].trackID);
                                            maxenfrac2 = mtrkidefrac12[abs(meveIDs2[e].trackID)];
                                            std::cout<<"maxen2 "<<maxen2<<" maxenid2 "<<maxenid2<<" RealID "<<meveIDs2[e].trackID<<" maxenfrac2 "<<maxenfrac2<<std::endl;
                                          }// if(mtrkide12[abs(meveIDs2[e].trackID)]>maxen2)
                                        }//for(size_t e=0;e<meveIDs2.size(); ++e)

                                        //see what other particles take away the energy fraction;
                                        for(std::map<int,double>::iterator ii = mtrkidefrac12.begin(); ii!=mtrkidefrac12.end(); ++ii)
                                        {              
                                          const simb::MCParticle *mparticleP2 = pi_serv->TrackIdToParticle_P(abs(ii->first));
                                          // if(!mparticleP) continue;
                                          // const art::Ptr<simb::MCTruth> mmcP=pi_serv->TrackIdToMCTruth_P(mtrackid);
                                          // if(!mmcP) continue;
                                          std::cout<<"TrackID "<<ii->first<<" energyFrac "<<ii->second<<" pdg "<<mparticleP2->PdgCode()<<" E "<<mparticleP2->E()<<std::endl;
                                        } //for(std::map<int,double>::iterator ii = mtrkidefrac12.begin(); ii!=mtrkidefrac12.end(); ++ii)

                                        ftrueMiEFrac[ftrueparhitcolcount] = mtrkidefrac12[abs(fmcd_trkid)];

                                        if(hitsfromPar1[e]->Multiplicity()<2)
                                        {
                                          //for counting true MC energy
                                          for(size_t e=0;e<meveIDs2.size(); ++e)
                                          { 
                                            mtrkide1[abs(meveIDs2[e].trackID)] += meveIDs2[e].energy;
                                          //		                cout<<"abs(meveIDs2[e].trackID) "<<abs(meveIDs2[e].trackID)<<" mtrkide1[abs(meveIDs2[e].trackID)] "<<mtrkide12[abs(meveIDs2[e].trackID)]<<std::endl;
                                          }// for(size_t e=0;e<meveIDs2.size(); ++e)
                                          ////////////////////////////////
                                          hit_key_fromHits.push_back(hitsfromPar1[e].key());
                                          chan_charge_fromHits.push_back(hitsfromPar1[e]->Integral());
                                          chan_no_fromHits.push_back(hitsfromPar1[e]->Channel());
                                          energy_fromHits.push_back(maxen2);
                                          ptime_minus_rms_fromHits.push_back(hitsfromPar1[e]->PeakTimeMinusRMS(3));
                                          ptime_plus_rms_fromHits.push_back(hitsfromPar1[e]->PeakTimePlusRMS(3));
                                          ptime_startTick_fromHits.push_back(hitsfromPar1[e]->StartTick());
                                          ptime_endTick_fromHits.push_back(hitsfromPar1[e]->EndTick());
                                        } // if(hitsfromPar1[e]->Multiplicity()<2)  

                                        if(abs(maxenid2) != abs(fmcd_trkid)) continue; // to make sure that the maximum energy deposition is coming from the true Michel 
                                        ftrueparhitscol_peakT[ftrueparhitcolcount] = hitsfromPar1[e]->PeakTime();
                                        ftrueparhitscol_charge[ftrueparhitcolcount]= hitsfromPar1[e]->Integral();
                                        ftrueparhitscol_wire[ftrueparhitcolcount]  = hitsfromPar1[e]->WireID().Wire;
                                        ftrueparhitscol_chno[ftrueparhitcolcount]  = hitsfromPar1[e]->Channel();
                                        ftrueparhitscol_TPC[ftrueparhitcolcount]   = hitsfromPar1[e]->WireID().TPC;
                                        ftrueparhitscol_plane[ftrueparhitcolcount] = hitsfromPar1[e]->WireID().Plane;
                                        ftrueparhitscol_xpos[ftrueparhitcolcount]  = trueparhitXpos;
                                        ftrueparhitscol_ypos[ftrueparhitcolcount]  = trueparhitYpos;
                                        ftrueparhitscol_zpos[ftrueparhitcolcount]  = trueparhitZpos;
                                        ftrueparhitscol_angledeg[ftrueparhitcolcount] = angle_theta_deg;
                                        ftrueparhitscol_maghitveccostheta[ftrueparhitcolcount] = maghitveccostheta;
                                        ftrueparhitscol_distance[ftrueparhitcolcount] = distance;
                                        ftrueparhitscol_mult[ftrueparhitcolcount]      = hitsfromPar1[e]->Multiplicity();
                                        ftrueparhitscol_sigptime[ftrueparhitcolcount]  = hitsfromPar1[e]->SigmaPeakTime();
                                        ftrueparhitscol_sigchrg[ftrueparhitcolcount]   = hitsfromPar1[e]->SigmaIntegral();
                                        ftrueparhitscol_sigpamp[ftrueparhitcolcount]    = hitsfromPar1[e]->SigmaPeakAmplitude();
                                        ftrueparhitscol_dof[ftrueparhitcolcount]       = hitsfromPar1[e]->DegreesOfFreedom();
                                        ftrueparhitscol_gof[ftrueparhitcolcount]       = hitsfromPar1[e]->GoodnessOfFit();
                                        ftrueparhitscol_ptminusRMS[ftrueparhitcolcount]     = hitsfromPar1[e]->PeakTimeMinusRMS(5.0);
                                        ftrueparhitscol_ptplusRMS[ftrueparhitcolcount]     = hitsfromPar1[e]->PeakTimePlusRMS(5.0);

                                        trueHitsKey.push_back(hitsfromPar1[e].key());

                                        ftrueparhitcolcount++;
                                      }//if(found_dup1 == 0)
                                    }//if(hitsfromPar1[e]->WireID().Plane==2)
                                    // std::cout<<"hit key from true par "<<hitsfromPar[e].key()<<std::endl;
                                    // } //if(same_trk == 0) 
                                  }//for(size_t e=0;e<hitsfromPar1.size(); ++e)
                                  std::cout<<"ftrueparhitcolcount "<<ftrueparhitcolcount<<std::endl;
                                  double  mmaxe1 = -1;
                                  //double mtote1 = 0;
                                  for(std::map<int,double>::iterator ii = mtrkide1.begin(); ii!=mtrkide1.end(); ++ii)
                                  {
                                    if(abs(ii->first) == fmcd_trkid)
                                    {
                                      std::cout<<"\ntrkid "<<ii->first<<" energy deposited from hitMult=1 "<<ii->second<<std::endl;
                                      //	                    if((ii->second)>mmaxe1)
                                      //	                    {
                                      mmaxe1 = ii->second;
                                      mtrackid1 = abs(ii->first);
                                      //	                    }
                                    }//if(abs(ii->first) == fmcd_trkid)
                                  }//for(std::map<int,double>::iterator ii = mtrkide1.begin(); ii!=mtrkide1.end(); ++ii)
                                  
                                  //--------Total Michel energy deposited on hits including hit multiplicity >= 2--------------//
                                  double tothitmultE = 0;
                                  for(size_t nn=0; nn<mult_chan.size(); nn++)
                                  {
                                    tothitmultE += mult_trueE.at(nn);
                                    chan_no_fromHits.push_back(mult_chan.at(nn));
                                    ptime_minus_rms_fromHits.push_back(mult_ptimeminus.at(nn)); 
                                    ptime_plus_rms_fromHits.push_back(mult_ptimeplus.at(nn)); 
                                    energy_fromHits.push_back(mult_trueE.at(nn)); 
                                    std::cout<<"tothitmultE "<<tothitmultE<<std::endl;
                                  }//for(size_t nn=0; nn<mult_chan.size(); nn++)
                                
                                  std::cout<<" Total true Michel energy deposited on hits "<<(tothitmultE+mmaxe1)<<std::endl;

                                  const simb::MCParticle *mparticleP1 = pi_serv->TrackIdToParticle_P(mtrackid1);
                                  // if(!mparticleP) continue;

                                  // const art::Ptr<simb::MCTruth> mmcP=pi_serv->TrackIdToMCTruth_P(mtrackid);
                                  // if(!mmcP) continue;

                                  std::cout<<"mtrackid1 "<<mtrackid1<<" energy "<<(tothitmultE+mmaxe1)<<" pdg "<<mparticleP1->PdgCode()<<" E "<<mparticleP1->E()<<std::endl;
                          
                                  fmcd_truecuthitsE     = (tothitmultE+mmaxe1);
                                  ftrueEdepo = truthUtil.GetDepEnergyMC(evt, fGeom, mtrackid1, 2 );
                                  std::cout<<"pdg "<<mparticleP1->PdgCode()<<" trueEdepo "<<ftrueEdepo<<std::endl;
                                  
                                  //-----------getting the MC true information of the selected nearby michel hits-----------------//

                                  std::map<int,double> mtrkide;
                                  int mtrackid=-1;
                                  for(int qq =0; qq < fnearhitct; qq++)
                                  {
                                    //std::cout<<"hitKey[qq] "<<hitKey[qq]<<std::endl;
                                    art::Ptr<recob::Hit> mhit=hitVect[fnearhits_key[qq]];    
                                    std::cout<<"michel hit key "<<mhit.key()<<" peakT "<<mhit->PeakTime()<<" WireID "<<mhit->WireID().Wire<<std::endl;
                                    std::vector<sim::TrackIDE> meveIDs = bt_serv->HitToTrackIDEs(clockData,mhit);

                                    for(size_t e=0;e<meveIDs.size(); ++e)
                                    { 
                                      if(abs(meveIDs[e].trackID)!=fmccut_trkid) //the hits should not belong to the same mother track (muon)
                                      {
                                        mtrkide[meveIDs[e].trackID] += meveIDs[e].energy;
                                        std::cout<<"meveIDs[e].trackID "<<meveIDs[e].trackID<<" mtrkide[meveIDs[e].trackID] "<<mtrkide[meveIDs[e].trackID]<<std::endl;
                                      }//if(abs(meveIDs[e].trackID)!=fmccut_trkid)
                                    }//for(size_t e=0;e<meveIDs.size(); ++e)
                                  }//for(int qq =0; qq < fnearhitct; qq++)
                                  double  mmaxe = -1;
                                  double mtote = 0;
                                  for(std::map<int,double>::iterator ii = mtrkide.begin(); ii!=mtrkide.end(); ++ii)
                                  {
                                    //	std::cout<<" trkid "<<ii->first<<"  energy deposited = "<<ii->second<<std::endl;
                                    mtote += ii->second;
                                    if((ii->second)>mmaxe)
                                    {
                                      mmaxe = ii->second;
                                      mtrackid = abs(ii->first);
                                    }//if((ii->second)>mmaxe)
                                  }//for(std::map<int,double>::iterator ii = mtrkide.begin(); ii!=mtrkide.end(); ++ii)

                                  const simb::MCParticle *mparticleP = pi_serv->TrackIdToParticle_P(mtrackid);
                                  if(!mparticleP) continue;

                                  const art::Ptr<simb::MCTruth> mmcP=pi_serv->TrackIdToMCTruth_P(mtrackid);
                                  if(!mmcP) continue;

                                  TLorentzVector mmcstart, mmcend, mmcstartdrifted, mmcenddrifted;
                                  unsigned int mpstarti, mpendi, mpstartdriftedi, mpenddriftedi; //mcparticle indices for starts and ends in tpc or drifted volumes
                                  double mplen = utils->length(*mparticleP, mmcstart, mmcend, mpstarti, mpendi);
                                  double mplendrifted = utils->driftedLength(detProp, *mparticleP, mmcstartdrifted, mmcenddrifted, mpstartdriftedi, mpenddriftedi);
                                  // bool isActive = plen != 0;
                                  bool misDrifted = mplendrifted!= 0;

                                  std::cout<<"Michel hits: bt:trackid "<<mtrackid<<std::endl;
                                  fmchits_trkid      = mparticleP->TrackId();
                                  fmchits_vx         = mparticleP->Vx();
                                  fmchits_vy         = mparticleP->Vy();
                                  fmchits_vz         = mparticleP->Vz();
                                  fmchits_t          = mparticleP->T();
                                  fmchits_endx       = mparticleP->EndX();
                                  fmchits_endy       = mparticleP->EndY();
                                  fmchits_endz       = mparticleP->EndZ();
                                  fmchits_endt       = mparticleP->EndT();
                                  fmchits_px         = mparticleP->Px();
                                  fmchits_py         = mparticleP->Py();
                                  fmchits_pz         = mparticleP->Pz();
                                  fmchits_momentum   = mparticleP->P();
                                  fmchits_energy     = mparticleP->E();
                                  fmchits_endpx      = mparticleP->EndPx();
                                  fmchits_endpy      = mparticleP->EndPy();
                                  fmchits_endpz      = mparticleP->EndPz();
                                  fmchits_endenergy  = mparticleP->EndE();
                                  fmchits_pathlen    = mplen;
                                  if (misDrifted)
                                  {
                                    fmchits_vxdrifted         = mparticleP->Vx();
                                    fmchits_vydrifted         = mparticleP->Vy();
                                    fmchits_vzdrifted         = mparticleP->Vz();
                                    fmchits_tdrifted          = mparticleP->T();
                                    fmchits_endxdrifted       = mparticleP->EndX();
                                    fmchits_endydrifted       = mparticleP->EndY();
                                    fmchits_endzdrifted       = mparticleP->EndZ();
                                    fmchits_endtdrifted       = mparticleP->EndT();
                                    fmchits_pxdrifted         = mparticleP->Px();
                                    fmchits_pydrifted         = mparticleP->Py();
                                    fmchits_pzdrifted         = mparticleP->Pz();
                                    fmchits_momentumdrifted   = mparticleP->P();
                                    fmchits_energydrifted     = mparticleP->E();
                                    fmchits_endpxdrifted      = mparticleP->EndPx();
                                    fmchits_endpydrifted      = mparticleP->EndPy();
                                    fmchits_endpzdrifted      = mparticleP->EndPz();
                                    fmchits_endenergydrifted  = mparticleP->EndE();
                                    fmchits_pathlendrifted    = mplendrifted;
                                  }//if (misDrifted)
                                  fmchits_endprocess = int(mparticleP->Trajectory().ProcessToKey(mparticleP->EndProcess()));
                                  fmchits_theta      = mparticleP->Momentum().Theta();
                                  fmchits_phi        = mparticleP->Momentum().Phi();
                                  fmchits_pdg        = mparticleP->PdgCode();
                                  fmchits_status_code= mparticleP->StatusCode();
                                  fmchits_mass       = mparticleP->Mass();
                                  fmchits_ND         = mparticleP->NumberDaughters();
                                  fmchits_mother     = mparticleP->Mother();
                                  fmchits_origin     = mmcP->Origin();
                                  fmchits_process    = int(mparticleP->Trajectory().ProcessToKey(mparticleP->Process()));
                                  fmchits_rescatter  = mparticleP->Rescatter();

                                  std::cout<<"MC true Michel hits: Pdg code "<<fmchits_pdg<<" track ID "<<fmchits_trkid<<" no of daughters "<<fmchits_ND<<" and origin_type "<<fmchits_origin<<std::endl;

                                  //--------Calculate the average missing energy-----------//
                                  fReadOutWindowSize = detProp.ReadOutWindowSize();
                                  fNumberTimeSamples = detProp.NumberTimeSamples();
                                  std::cout<<"fReadOutWindowSize "<<fReadOutWindowSize<<" fNumberTimeSamples "<<fNumberTimeSamples<<std::endl;

                                  double totEMiDepo =0;
                                  double totNelectrons = 0;
                                  float avg_missing_energy = 0, avg_missing_numelec=0;
                                  favg_missing_energy = 0, favg_missing_numelec = 0;
                                  
                                  art::Handle< std::vector<sim::SimChannel> > simchannelHandle;
                                  if(evt.getByLabel("largeant", simchannelHandle))
                                  {
                                    //Loop over simchannels 
                                    for(auto const& simchannel : (*simchannelHandle))
                                    {
                                      if(fGeom->View(simchannel.Channel()) != 2) continue;
                                      auto const& alltimeslices = simchannel.TDCIDEMap();
                                      // bool savethis = false;
                                      // std::cout<<"Simchan_charge["<<simchannel.Channel()<<"] "<<chan_charge[simchannel.Channel()]<<std::endl;
                                      // if( chan_charge[simchannel.Channel()] ==  0) savethis = true;
                                
                                      double MichelEperCh = 0;
                                      double MissingEperCh = 0;
                                
                                      // Loop over ticks 60000
                                      for(auto const& tslice : alltimeslices)
                                      {
                                        int spill = 1;
                                        if (fReadOutWindowSize != fNumberTimeSamples) 
                                        {
                                          if (tslice.first < spill * fReadOutWindowSize || tslice.first > (spill + 1) * fReadOutWindowSize) continue;
                                        }
                                        else if (tslice.first < 0 || tslice.first > fReadOutWindowSize) continue;		 
                                        auto const& simide = tslice.second;
                                        float simETick = 0;
                                  
                                        //---For calculation of missing energy, only look for TDCs where we dont see a signal------//
                                        bool missed_signal=true;
                                        int sim_channel = int(simchannel.Channel());
                                        // std::cout<<"chan_charge_fromHits.size() "<<chan_charge_fromHits.size()<<std::endl;

                                        //---------- recob:hit signal--------------//
                                        for(size_t nn=0; nn<chan_no_fromHits.size(); nn++)
                                        {
                                          if(chan_no_fromHits.at(nn) == sim_channel && ptime_minus_rms_fromHits.at(nn)<=tslice.first && ptime_plus_rms_fromHits.at(nn)>=tslice.first) missed_signal=false;
                                          
                                        }

                                        // Loop over energy deposits
                                        for(auto const& eDep : simide)
                                        {
                                          if(eDep.trackID == fmcd_trkid || eDep.trackID == -fmcd_trkid)
                                          {
                                            MichelEperCh += eDep.energy;  
                                            totEMiDepo += eDep.energy;
                                            totNelectrons += eDep.numElectrons;
                                            // std::cout<<"eDep.energy "<<eDep.energy<<std::endl;
                                            // if( savethis ) 
                                            if( missed_signal ) 
                                            {
                                              simETick += eDep.energy;
                                              avg_missing_energy += eDep.energy;
                                              avg_missing_numelec += eDep.numElectrons;
                                              MissingEperCh += eDep.energy;
                                              // std::cout<<"eDep.energy "<<eDep.energy<<" avg_missing_energy "<<avg_missing_energy<<" MeV\n";
                                            }// if(missed_signal)
                                          }//if(eDep.trackID == fmcd_trkid || eDep.trackID == -fmcd_trkid)
                                        }//for(auto const& eDep : simide
                                      }//for(auto const& tslice : alltimeslices)
                                    } //for(auto const& simchannel : (*simchannelHandle))  
                                  }//if(evt.getByLabel("largeant", simchannelHandle))
                                  
                                  std::cout<<"chan_no_fromHits.size() "<<chan_no_fromHits.size()<<" chan_charge_fromHits.size() "<<chan_charge_fromHits.size()<<std::endl;

                                  favg_missing_energy = avg_missing_energy;
                                  favg_missing_numelec = avg_missing_numelec;
                                  std::cout<<"\ntotEMiDepo "<<totEMiDepo<<" avg_missing_E "<<favg_missing_energy<<" MeV\n";
                                  std::cout<<"avg_missing_numelec "<<favg_missing_numelec<<std::endl;

                                  std::cout<<"\nGEANTMiE "<<fmcd_energy*1000<<" totEMiDepo "<<totEMiDepo<<" fmcd_truecuthitsE "<<fmcd_truecuthitsE<<" avg_missing_E "<<favg_missing_energy<<" cut+missingE "<<(fmcd_truecuthitsE+favg_missing_energy)<<" MeV\n"<<std::endl;
                        
                                } //if(!evt.isRealData())  
                                
                                //------------Looping over all hits and checking if they belong to the michel cone -------------//
                                fmhitcount = 0;
                                for(size_t ll=0; ll<hitVect.size();++ll)
                                {
                                  if(hitVect[ll]->WireID().Plane==2)
                                  {
                                    //----Make sure that this hit does not belong to the  candidate muon track----//
                                    int same_trk = 0;
                                    for(size_t mm=0;mm<longtrk_hitkey.size();mm++)
                                    {
                                      if(longtrk_hitkey.at(mm)==hitVect[ll].key())
                                      { 
                                        same_trk = 1;
                                      }//if(longtrk_hitkey.at(mm)==hitVect[ll].key())
                                    }//for(size_t mm=0;mm<longtrk_hitkey.size();mm++)
                              
                                    //---To make sure that these hits dont belong to another long track-----//
                                    int long_trk = 0;
                                    for(size_t o1=0;o1<HitsTracksKey.size();o1++)
                                    {
                                      if(HitsTracksKey.at(o1)==hitVect[ll].key())
                                      {
                                        long_trk = 1;
                                      }//if(HitsTracksKey.at(o1)==hitVect[ll].key())
                                    }//for(size_t o1=0;o1<HitsTracksKey.size();o1++)
                                    int corr_hit = 0;
                                    for(size_t o1=0;o1<trueHitsKey.size();o1++)
                                    {
                                      if(trueHitsKey.at(o1)==hitVect[ll].key())
                                      {
                                        corr_hit = 1;
                                        std::cout<<"yes corr hit "<<std::endl;

                                      }//if(HitsTracksKey.at(o1)==hitVect[ll].key())
                                    }//for(size_t o1=0;o1<HitsTracksKey.size();o1++)

                                    double mdiffpeaktt0 = hitVect[ll]->PeakTime() - (T0_value0/1000)*2;
                                    double mhitXposcone = detProp.ConvertTicksToX(mdiffpeaktt0,hitVect[ll]->WireID().Plane,hitVect[ll]->WireID().TPC,hitVect[ll]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
                                    double mWirestart[3], mWireend[3];
                                    fGeom->WireEndPoints(hitVect[ll]->WireID());
                                    // fGeom->WireEndPoints(hitVect[ll]->WireID().Cryostat, hitVect[ll]->WireID().TPC, hitVect[ll]->WireID().Plane, hitVect[ll]->WireID().Wire, mWirestart, mWireend);
                                    double mhitZposcone = mWirestart[2];
                                    double mhitYposcone = 0;
                                    if(hits_spacept.at(ll).size())
                                    mhitYposcone      =   hits_spacept.at(ll)[0]->XYZ()[1];
                                    if(!hits_spacept.at(ll).size()) mhitYposcone = 303.5; //average y position in the detecor FV
                                
                                    double magnearhitvec = sqrt(pow(xposend-fcut_endhitx,2)+pow(fnearhits_zpos[int(fnearhitct/2)]-fcut_endhitz,2));
                                    double maghitvec = sqrt(pow(mhitXposcone-fcut_endhitx,2)+pow(mhitZposcone-fcut_endhitz,2));
                                    double angle_theta = acos(((xposend-fcut_endhitx)*(mhitXposcone-fcut_endhitx) + (fnearhits_zpos[int(fnearhitct/2)]-fcut_endhitz)*(mhitZposcone-fcut_endhitz))/(magnearhitvec*maghitvec));
                                    double angle_theta_deg = (angle_theta*180)/3.14;
                                    double maghitveccostheta = maghitvec* cos(angle_theta);
                                    double distance = sqrt(pow(mhitXposcone-fcut_trkendx,2)+pow(mhitZposcone-fcut_trkendz,2));
                                        
                                    if(corr_hit==1)
                                    {
                                      std::cout<<"angle_theta_deg "<<angle_theta_deg<<" maghitveccostheta "<<maghitveccostheta<<" distance "<<distance<<std::endl;
                                    } 

                                    if(long_trk ==0 && same_trk==0 /*&& abs(maghitveccostheta) < 80 && angle_theta_deg<30 && maghitvec > 0 && maghitvec < 80*/)
                                    {
                                      fmhits_key[fmhitcount]    = hitVect[ll].key();
                                      fmhits_peakT[fmhitcount] = hitVect[ll]->PeakTime();
                                      fmhits_charge[fmhitcount]= hitVect[ll]->Integral();
                                      fmhits_wire[fmhitcount]  = hitVect[ll]->WireID().Wire;
                                      fmhits_chno[fmhitcount]  = hitVect[ll]->Channel();
                                      fmhits_TPC[fmhitcount]   = hitVect[ll]->WireID().TPC;
                                      fmhits_plane[fmhitcount] = hitVect[ll]->WireID().Plane;
                                      fmhits_xpos[fmhitcount]  = mhitXposcone;
                                      fmhits_ypos[fmhitcount]  = mhitYposcone;
                                      fmhits_zpos[fmhitcount]  = mhitZposcone;
                                      fmhits_mult[fmhitcount]    = hitVect[ll]->Multiplicity();
                                      fmhits_sigptime[fmhitcount]= hitVect[ll]->SigmaPeakTime();
                                      fmhits_sigchrg[fmhitcount] = hitVect[ll]->SigmaIntegral();
                                      fmhits_sigpamp[fmhitcount]  = hitVect[ll]->SigmaPeakAmplitude();
                                      fmhits_dof[fmhitcount]     = hitVect[ll]->DegreesOfFreedom();
                                      fmhits_gof[fmhitcount]     = hitVect[ll]->GoodnessOfFit();
                                      fmhits_ptminusRMS[fmhitcount]     = hitVect[ll]->PeakTimeMinusRMS(5.0);
                                      fmhits_ptplusRMS[fmhitcount]     = hitVect[ll]->PeakTimePlusRMS(5.0);
                                      fmhits_angledeg[fmhitcount] = angle_theta_deg;
                                      fmhits_maghitveccostheta[fmhitcount] = maghitveccostheta;
                                      fmhits_distance[fmhitcount] = distance;
                                      fmhits_longtrk[fmhitcount] = long_trk;
                                      fmhits_sametrk[fmhitcount] = same_trk;			       
                                      fmhits_corrhit[fmhitcount] = corr_hit;
                                          
                                      std::array<float,4> cnn_out = hitResults.getOutput( hitVect[ll].key() );
                                      fmhits_cnnMichel[fmhitcount]     = cnn_out[ hitResults.getIndex("michel") ];
                                      fmhits_cnnEM[fmhitcount]         = cnn_out[ hitResults.getIndex("em") ];
                                      fmhits_cnnTrack[fmhitcount]      = cnn_out[ hitResults.getIndex("track") ];

                                      if(fmhits_corrhit[fmhitcount]==1)
                                      {
                                        std::cout<<"hitID "<<hitVect[ll].key()<<" long_trk "<<fmhits_longtrk[fmhitcount]<<" same_trk "<<fmhits_sametrk[fmhitcount]<<" corr_hit "<<fmhits_corrhit[fmhitcount]<<" charge "<<fmhits_charge[fmhitcount]<<" michel "<<cnn_out[ hitResults.getIndex("michel") ]<<" trk "<<cnn_out[ hitResults.getIndex("track") ]<<" em "<<cnn_out[ hitResults.getIndex("em") ]<<std::endl;
                                      }
                                      fmhitcount++;
                                    }//if(long_trk ==0 && same_trk==0)
                                  }//if(hitVect[ll]->WireID().Plane==2)
                                }//for(size_t ll=0; ll<hitVect.size();++ll)

                                trueHitsKey.clear();
                                longtrk_hitkey.clear();
                                
                                //-----------recover some michel hits from the last 10 muon hits----------------------//
                                float trunc_charge[10000] = {-999};
                                double mean_charge = 0;
                                double max_truncchrg = 0;
                                int   nnewcolhits = 0; 
                                for(int kk=1; kk<ftrkcolhits; kk++)
                                {
                                  if(kk<ftrkcolhits-1)
                                  mean_charge = (fhits_charge[kk-1]+fhits_charge[kk]+fhits_charge[kk+1])/ 3;
                                  if(kk==ftrkcolhits-1)
                                  mean_charge = (fhits_charge[kk-1]+fhits_charge[kk])/ 2;
                                  
                                  if(fhits_charge[kk]>0.2*mean_charge && fhits_charge[kk]<2*mean_charge)
                                  trunc_charge[kk] = mean_charge;
                                  else
                                  trunc_charge[kk] = fhits_charge[kk];
                                  if(kk>=ftrkcolhits-10 && kk<ftrkcolhits)
                                  {
                                    if(trunc_charge[kk]>max_truncchrg)
                                    {
                                      max_truncchrg = trunc_charge[kk];
                                      nnewcolhits = kk;
                                    }//if(trunc_charge[kk]>max_truncchrg)
                                  }//if(kk>=ftrkcolhits-10 && kk<ftrkcolhits)
                                }//for(int kk=1; kk<ftrkcolhits; kk++)
                                // cout<<"extraHits "<<ftrkcolhits-1-nnewcolhits<<"\n";

                                std::vector <int>    cone_wire;
                                std::vector <int>    cone_tpc;
                                std::vector <int>    cone_mult;
                                std::vector <double> cone_key;
                                std::vector <double> cone_charge;
                                std::vector <double> cone_ptime;
                                std::vector <double> addmu_key;
                                std::vector <double> cone_xpos;
                                std::vector <double> cone_ypos;
                                std::vector <double> cone_zpos;
                                std::vector <double> cone_angledeg;
                                std::vector <double> cone_maghitveccostheta;
                                std::vector <double> cone_distance;
                                std::vector <double> cone_sigptime;
                                std::vector <double> cone_sigchrg;
                                std::vector <double> cone_sigpamp;
                                std::vector <double> cone_dof;
                                std::vector <double> cone_gof;
                                std::vector <double> cone_ptminusRMS;
                                std::vector <double> cone_ptplusRMS;
                                std::vector <double> cone_status;
                                std::vector <double> cone_cnnMichel;
                                std::vector <double> cone_cnnEM;
                                std::vector <double> cone_cnnTrack;

                                int add_hits = 0;
                                for(int kk=0; kk<ftrkcolhits; kk++)
                                {
                                  double magnearhitvec = sqrt(pow(xposend-fcut_endhitx,2)+pow(fnearhits_zpos[int(fnearhitct/2)]-fcut_endhitz,2));
                                  double maghitvec = sqrt(pow(fhits_xpos[kk]-fcut_endhitx,2)+pow(fhits_zpos[kk]-fcut_endhitz,2));
                                  double angle_theta = 0;
                                  if(fhits_xpos[kk]!=fcut_endhitx)
                                  angle_theta = acos(((xposend-fcut_endhitx)*(fhits_xpos[kk]-fcut_endhitx) + (fnearhits_zpos[int(fnearhitct/2)]-fcut_endhitz)*(fhits_zpos[kk]-fcut_endhitz))/(magnearhitvec*maghitvec));
                                  double angle_theta_deg = (angle_theta*180)/3.14;
                                  double maghitveccostheta = maghitvec* cos(angle_theta);
                                  double distance = sqrt(pow(fhits_xpos[kk]-fcut_trkendx,2)+pow(fhits_ypos[kk]-fcut_trkendy,2)+pow(fhits_zpos[kk]-fcut_trkendz,2));

                                  if(kk>nnewcolhits && kk<ftrkcolhits)
                                  if(distance<=1.5 || (angle_theta_deg>5 && angle_theta_deg<40))
                                  {
                                    add_hits++;
                                    cone_key.push_back(fhits_key[kk]);
                                    cone_charge.push_back(fhits_charge[kk]);
                                    cone_wire.push_back(fhits_wire[kk]);
                                    cone_tpc.push_back(fhits_TPC[kk]);
                                    cone_ptime.push_back(fhits_peakT[kk]);
                                    cone_xpos.push_back(fhits_xpos[kk]);
                                    cone_ypos.push_back(fhits_ypos[kk]);
                                    cone_zpos.push_back(fhits_zpos[kk]);
                                    cone_angledeg.push_back(angle_theta_deg);
                                    cone_maghitveccostheta.push_back(maghitveccostheta);
                                    cone_distance.push_back(distance);
                                    cone_mult.push_back(fhits_mult[kk]);
                                    cone_sigptime.push_back(fhits_sigptime[kk]);
                                    cone_sigchrg.push_back(fhits_sigchrg[kk]);
                                    cone_sigpamp.push_back(fhits_sigpamp[kk]);
                                    cone_dof.push_back(fhits_dof[kk]);
                                    cone_gof.push_back(fhits_gof[kk]);
                                    cone_ptminusRMS.push_back(fhits_ptminusRMS[kk]);
                                    cone_ptplusRMS.push_back(fhits_ptplusRMS[kk]);
                                    cone_status.push_back(2);
                                    cone_cnnMichel.push_back(fhits_cnnMichel[kk]);
                                    cone_cnnEM.push_back(fhits_cnnEM[kk]);
                                    cone_cnnTrack.push_back(fhits_cnnTrack[kk]);
                                  }//if(distance<=1.5 || (angle_theta_deg>5 && angle_theta_deg<40))
                                }//for(int kk=0; kk<ftrkcolhits; kk++)
                    
                                //--------analyzing michel cone hits MC-----------------//
                                int inicone_hits = 0;
                                for(int kk=0; kk<fmhitcount; kk++)
                                {
                                  double magnearhitvec = sqrt(pow(xposend-fcut_endhitx,2)+pow(fnearhits_zpos[int(fnearhitct/2)]-fcut_endhitz,2));
                                  double maghitvec = sqrt(pow(fmhits_xpos[kk]-fcut_endhitx,2)+pow(fmhits_zpos[kk]-fcut_endhitz,2));
                                  double angle_theta = acos(((xposend-fcut_endhitx)*(fmhits_xpos[kk]-fcut_endhitx) + (fnearhits_zpos[int(fnearhitct/2)]-fcut_endhitz)*(fmhits_zpos[kk]-fcut_endhitz))/(magnearhitvec*maghitvec));
                                  double angle_theta_deg = (angle_theta*180)/3.14;
                                  double maghitveccostheta = maghitvec* cos(angle_theta);
                                  double distance = sqrt(pow(fmhits_xpos[kk]-fcut_trkendx,2)+pow(fmhits_ypos[kk]-fcut_trkendy,2)+pow(fmhits_zpos[kk]-fcut_trkendz,2));

                                  float y0 = (-(_cone_length/_cone_angle1)*angle_theta_deg)+_cone_length;
                            
                                  if(maghitveccostheta <= y0 && angle_theta_deg<_cone_angle1)
                                  {
                                    inicone_hits++;
                                    cone_key.push_back(fmhits_key[kk]); 
                                    cone_charge.push_back(fmhits_charge[kk]);
                                    cone_wire.push_back(fmhits_wire[kk]);
                                    cone_tpc.push_back(fmhits_TPC[kk]);
                                    cone_ptime.push_back(fmhits_peakT[kk]);
                                    cone_xpos.push_back(fmhits_xpos[kk]);
                                    cone_ypos.push_back(fmhits_ypos[kk]);
                                    cone_zpos.push_back(fmhits_zpos[kk]);
                                    cone_angledeg.push_back(angle_theta_deg);
                                    cone_maghitveccostheta.push_back(maghitveccostheta);
                                    cone_distance.push_back(distance);
                                    cone_mult.push_back(fmhits_mult[kk]);
                                    cone_sigptime.push_back(fmhits_sigptime[kk]);
                                    cone_sigchrg.push_back(fmhits_sigchrg[kk]);
                                    cone_sigpamp.push_back(fmhits_sigpamp[kk]);
                                    cone_dof.push_back(fmhits_dof[kk]);
                                    cone_gof.push_back(fmhits_gof[kk]);
                                    cone_ptminusRMS.push_back(fmhits_ptminusRMS[kk]);
                                    cone_ptplusRMS.push_back(fmhits_ptplusRMS[kk]);
                                    cone_status.push_back(1);
                                    cone_cnnMichel.push_back(fmhits_cnnMichel[kk]);
                                    cone_cnnEM.push_back(fmhits_cnnEM[kk]);
                                    cone_cnnTrack.push_back(fmhits_cnnTrack[kk]);
                                  } //if(maghitveccostheta <= y0 && angle_theta_deg<_cone_angle1)
                                }//for(int kk=0; kk<fmhitcount; kk++)

                                //TODO: optimize this loop                                
                                // std::vector<Cone> cones(cone_key.size());
                                // // Populate the cones vector...
                                // for(size_t i = 0; i < cone_key.size(); i++) {
                                //     cones[NTracks].key = cone_key[NTracks];
                                //     cones[NTracks].charge = cone_charge[NTracks];
                                //     cones[NTracks].wire = cone_wire[NTracks];
                                //     cones[NTracks].tpc = cone_tpc[NTracks];
                                //     cones[NTracks].ptime = cone_ptime[NTracks];
                                //     // ... set other properties ...
                                //     cones[NTracks].cnnTrack = cone_cnnTrack[NTracks];
                                // }

                                // std::sort(cones.begin(), cones.end(), [](const Cone& a, const Cone& b) {
                                //     return a.ptime < b.ptime;
                                // });

                                float o3=0,o4=0,o5=0,o6=0,o7=0,o8=0,o9=0,o10=0,o11=0,o12=0,o13=0,o14=0,o15=0,o16=0,o17=0,o18=0,o19=0,o20=0,o21=0,o22=0,o23=0,o24=0,o25=0;
                                for(size_t o1=0;o1<cone_key.size();o1++)
                                for(size_t o2=o1+1; o2 < cone_key.size(); o2++)
                                {
                                  if(cone_ptime.at(o1) > cone_ptime.at(o2))
                                  {
                                    o3 = cone_key.at(o1);
                                    cone_key.at(o1) = cone_key.at(o2);
                                    cone_key.at(o2) = o3;

                                    o4 = cone_charge.at(o1);
                                    cone_charge.at(o1) = cone_charge.at(o2);
                                    cone_charge.at(o2) = o4;

                                    o5 = cone_wire.at(o1);
                                    cone_wire.at(o1) = cone_wire.at(o2);
                                    cone_wire.at(o2) = o5;

                                    o6 = cone_tpc.at(o1);
                                    cone_tpc.at(o1) = cone_tpc.at(o2);
                                    cone_tpc.at(o2) = o6;
                                    
                                    o7 = cone_ptime.at(o1);
                                    cone_ptime.at(o1) = cone_ptime.at(o2);
                                    cone_ptime.at(o2) = o7;

                                    o8 = cone_xpos.at(o1);
                                    cone_xpos.at(o1) = cone_xpos.at(o2);
                                    cone_xpos.at(o2) = o8;

                                    o9 = cone_ypos.at(o1);
                                    cone_ypos.at(o1) = cone_ypos.at(o2);
                                    cone_ypos.at(o2) = o9;

                                    o10 = cone_zpos.at(o1);
                                    cone_zpos.at(o1) = cone_zpos.at(o2);
                                    cone_zpos.at(o2) = o10;

                                    o11 = cone_angledeg.at(o1);
                                    cone_angledeg.at(o1) = cone_angledeg.at(o2);
                                    cone_angledeg.at(o2) = o11;

                                    o12 = cone_maghitveccostheta.at(o1);
                                    cone_maghitveccostheta.at(o1) = cone_maghitveccostheta.at(o2);
                                    cone_maghitveccostheta.at(o2) = o12;

                                    o13 = cone_distance.at(o1);
                                    cone_distance.at(o1) = cone_distance.at(o2);
                                    cone_distance.at(o2) = o13;

                                    o14 = cone_mult.at(o1);
                                    cone_mult.at(o1) = cone_mult.at(o2);
                                    cone_mult.at(o2) = o14;

                                    o15 = cone_sigptime.at(o1);
                                    cone_sigptime.at(o1) = cone_sigptime.at(o2);
                                    cone_sigptime.at(o2) = o15;

                                    o16 = cone_sigchrg.at(o1);
                                    cone_sigchrg.at(o1) = cone_sigchrg.at(o2);
                                    cone_sigchrg.at(o2) = o16;

                                    o17 = cone_sigpamp.at(o1);
                                    cone_sigpamp.at(o1) = cone_sigpamp.at(o2);
                                    cone_sigpamp.at(o2) = o17;

                                    o18 = cone_dof.at(o1);
                                    cone_dof.at(o1) = cone_dof.at(o2);
                                    cone_dof.at(o2) = o18;

                                    o19 = cone_gof.at(o1);
                                    cone_gof.at(o1) = cone_gof.at(o2);
                                    cone_gof.at(o2) = o19;

                                    o20 = cone_ptminusRMS.at(o1);
                                    cone_ptminusRMS.at(o1) = cone_ptminusRMS.at(o2);
                                    cone_ptminusRMS.at(o2) = o20;

                                    o21 = cone_ptplusRMS.at(o1);
                                    cone_ptplusRMS.at(o1) = cone_ptplusRMS.at(o2);
                                    cone_ptplusRMS.at(o2) = o21;

                                    o22 = cone_status.at(o1);
                                    cone_status.at(o1) = cone_status.at(o2);
                                    cone_status.at(o2) = o22;

                                    o23 = cone_cnnMichel.at(o1);
                                    cone_cnnMichel.at(o1) = cone_cnnMichel.at(o2);
                                    cone_cnnMichel.at(o2) = o23;

                                    o24 = cone_cnnEM.at(o1);
                                    cone_cnnEM.at(o1) = cone_cnnEM.at(o2);
                                    cone_cnnEM.at(o2) = o24;

                                    o25 = cone_cnnTrack.at(o1);
                                    cone_cnnTrack.at(o1) = cone_cnnTrack.at(o2);
                                    cone_cnnTrack.at(o2) = o25;
                                  }//if(cone_ptime.at(o1) > cone_ptime.at(o2))
                                }//for(size_t o2=o1+1; o2 < cone_key.size(); o2++)
                              
                                //-----Stack up all cone + additional mu hits
                                fmichel_conesize = cone_key.size();
                                std::map<int,std::vector<double>> hit_t1;
                                std::map<int,std::vector<double>> hit_t2;
                                std::map<int,std::vector<double>> hit_y;
                                std::map<int,std::vector<double>> hit_x;
                                std::map<int,std::vector<double>> hit_key;
                                std::map<int,std::vector<double>> hit_wire;
                                std::map<int,std::vector<double>> hit_tpc;
                                std::map<int,std::vector<double>> hit_chargehit;
                                std::map<int,std::vector<double>> hit_ptime;
                                std::map<int,std::vector<double>> hit_angledeg;
                                std::map<int,std::vector<double>> hit_maghitveccostheta;
                                std::map<int,std::vector<double>> hit_distance;
                                std::map<int,std::vector<double>> hit_mult;
                                std::map<int,std::vector<double>> hit_sigptime;
                                std::map<int,std::vector<double>> hit_sigchrg;
                                std::map<int,std::vector<double>> hit_sigpamp;
                                std::map<int,std::vector<double>> hit_dof;
                                std::map<int,std::vector<double>> hit_gof;
                                std::map<int,std::vector<double>> hit_ptminusRMS;
                                std::map<int,std::vector<double>> hit_ptplusRMS;
                                std::map<int,std::vector<double>> hit_status;
                                std::map<int,std::vector<double>> hit_cnnMichel;
                                std::map<int,std::vector<double>> hit_cnnEM;
                                std::map<int,std::vector<double>> hit_cnnTrack;
                            
                                fmichel_wcount = 0;
                          
                                for(size_t o1=0;o1<cone_key.size();o1++)
                                {
                                  std::vector<art::Ptr<recob::Wire>> wirescone = wire_hits.at(cone_key.at(o1));
                                  hit_t1[wirescone[0]->Channel()].push_back(cone_ptminusRMS.at(o1));
                                  hit_t2[wirescone[0]->Channel()].push_back(cone_ptplusRMS.at(o1));
                                  hit_key[wirescone[0]->Channel()].push_back(cone_key.at(o1));
                                  hit_wire[wirescone[0]->Channel()].push_back(cone_wire.at(o1));
                                  hit_tpc[wirescone[0]->Channel()].push_back(cone_tpc.at(o1));
                                  hit_chargehit[wirescone[0]->Channel()].push_back(cone_charge.at(o1));
                                  hit_ptime[wirescone[0]->Channel()].push_back(cone_ptime.at(o1));
                                  hit_angledeg[wirescone[0]->Channel()].push_back(cone_angledeg.at(o1));
                                  hit_maghitveccostheta[wirescone[0]->Channel()].push_back(cone_maghitveccostheta.at(o1));
                                  hit_distance[wirescone[0]->Channel()].push_back(cone_distance.at(o1));
                                  hit_mult[wirescone[0]->Channel()].push_back(cone_mult.at(o1));
                                  hit_sigptime[wirescone[0]->Channel()].push_back(cone_sigptime.at(o1));
                                  hit_sigchrg[wirescone[0]->Channel()].push_back(cone_sigchrg.at(o1));
                                  hit_sigpamp[wirescone[0]->Channel()].push_back(cone_sigpamp.at(o1));
                                  hit_dof[wirescone[0]->Channel()].push_back(cone_dof.at(o1));
                                  hit_gof[wirescone[0]->Channel()].push_back(cone_gof.at(o1));
                                  hit_ptminusRMS[wirescone[0]->Channel()].push_back(cone_ptminusRMS.at(o1));
                                  hit_ptplusRMS[wirescone[0]->Channel()].push_back(cone_ptplusRMS.at(o1));
                                  hit_status[wirescone[0]->Channel()].push_back(cone_status.at(o1));
                                  hit_cnnMichel[wirescone[0]->Channel()].push_back(cone_cnnMichel.at(o1));
                                  hit_cnnEM[wirescone[0]->Channel()].push_back(cone_cnnEM.at(o1));
                                  hit_cnnTrack[wirescone[0]->Channel()].push_back(cone_cnnTrack.at(o1));
                                  hit_x[wirescone[0]->Channel()].push_back(cone_xpos.at(o1));	 	  
                                  hit_y[wirescone[0]->Channel()].push_back(cone_ypos.at(o1)); 
                                }//for(size_t o1=0;o1<cone_key.size();o1++)  


                                //one wire can have various hits so lets remove duplicate wires
                                auto const& wires = evt.getValidHandle<std::vector<recob::Wire> >(label_wires);
                                auto w1 = hit_t1.begin();
                                auto w2 = hit_t2.begin();
                                auto x  = hit_x.begin();
                                auto y  = hit_y.begin();
                                auto key = hit_key.begin(); 
                                auto wi = hit_wire.begin();
                                auto chargehit = hit_chargehit.begin();
                                auto tpc = hit_tpc.begin();
                                auto ptime = hit_ptime.begin();
                                auto angledeg = hit_angledeg.begin();
                                auto maghitveccostheta = hit_maghitveccostheta.begin();
                                auto distance = hit_distance.begin();
                                auto mult = hit_mult.begin();
                                auto sigptime = hit_sigptime.begin();
                                auto sigchrg = hit_sigchrg.begin();
                                auto sigpamp = hit_sigpamp.begin();
                                auto dof = hit_dof.begin();
                                auto gof = hit_gof.begin();
                                auto ptminusRMS = hit_ptminusRMS.begin();
                                auto ptplusRMS = hit_ptplusRMS.begin();
                                auto status = hit_status.begin();
                                auto cnnMichel = hit_cnnMichel.begin();
                                auto cnnEM = hit_cnnEM.begin();
                                auto cnnTrack = hit_cnnTrack.begin();
                            
                                while( w1 != hit_t1.end())
                                {
                                  int it_w = w1->first;
                                  int n_hits = w1->second.size();
                                  std::cout<<"it_w "<<it_w<<" n_hits "<<n_hits<<" w1->second[0] "<<w1->second[0]<<"  w2->second[0] "<<w2->second[0]<<" w1->second[n_hits-1] "<<w1->second[n_hits-1]<<" w2->second[n_hits-1] "<<w2->second[n_hits-1]<<std::endl;
                                  double t1 =  w1->second[0]; //first hit
                                  double t2 =  w2->second[n_hits-1]; //last hit
                                  double x_w, y_w, key_w, wi_w, chargehit_w, tpc_w, ptime_w, angledeg_w, maghitveccostheta_w, distance_w, mult_w, sigptime_w, 
                                  sigchrg_w, sigpamp_w, dof_w, gof_w, ptminusRMS_w, ptplusRMS_w, status_w, cnnMichel_w, cnnEM_w, cnnTrack_w;
                                  x_w = x->second[0];
                                  y_w = y->second[0];
                                  key_w = key->second[0];
                                  wi_w = wi->second[0];
                                  tpc_w = tpc->second[0];
                                  chargehit_w = chargehit->second[0];
                                  ptime_w = ptime->second[0];
                                  angledeg_w = angledeg->second[0];
                                  maghitveccostheta_w = maghitveccostheta->second[0];
                                  distance_w = distance->second[0];
                                  mult_w = mult->second[0];
                                  sigptime_w = sigptime->second[0];
                                  sigchrg_w = sigchrg->second[0];
                                  sigpamp_w = sigpamp->second[0];
                                  dof_w = dof->second[0];
                                  gof_w = gof->second[0];
                                  ptminusRMS_w = ptminusRMS->second[0];
                                  ptplusRMS_w = ptplusRMS->second[0];
                                  status_w = status->second[0];
                                  cnnMichel_w = cnnMichel->second[0];
                                  cnnEM_w = cnnEM->second[0];
                                  cnnTrack_w = cnnTrack->second[0];
                              
                                  for(auto & wire : * wires)
                                  {
                                    int channel_no = wire.Channel();
                                    int plane = fGeom->View(wire.Channel()); 
                                    if( plane != 2 ) continue;
                                    std::vector< geo::WireID > wireID= fGeom->ChannelToWire(channel_no);
                                    const geo::WireGeo* pwire = fGeom->WirePtr(wireID.at(0)); //for collection plane there is one wire per channel
                                    geo::Point_t xyzWire = pwire->GetCenter();
                                    if(it_w == channel_no)
                                    {  
                                      double charge =0.0;
                                      for(size_t i = 0; i < wire.Signal().size(); ++i)
                                      {
                                        if( i > t1 && i < t2 ) charge += wire.Signal()[NTracks];
                                      }   
                                      fmichel_zpos[fmichel_wcount] = (xyzWire.Z()); 
                                      fmichel_ypos[fmichel_wcount] = (y_w); 
                                      fmichel_xpos[fmichel_wcount] = (x_w); 
                                      fmichel_chrg[fmichel_wcount] = (charge);
                                      fmichel_chno[fmichel_wcount] = (channel_no);
                                      fmichel_key[fmichel_wcount] = (key_w);
                                      fmichel_wire[fmichel_wcount] = (wi_w);
                                      fmichel_chargehit[fmichel_wcount] = (chargehit_w);
                                      fmichel_tpc[fmichel_wcount]  = (tpc_w);
                                      fmichel_ptime[fmichel_wcount]= (ptime_w);
                                      fmichel_angledeg[fmichel_wcount]= (angledeg_w);
                                      fmichel_maghitveccostheta[fmichel_wcount]= (maghitveccostheta_w);
                                      fmichel_distance[fmichel_wcount]= (distance_w);
                                      fmichel_mult[fmichel_wcount]= (mult_w);
                                      fmichel_sigptime[fmichel_wcount]= (sigptime_w);
                                      fmichel_sigchrg[fmichel_wcount]= (sigchrg_w);
                                      fmichel_sigpamp[fmichel_wcount]= (sigpamp_w);
                                      fmichel_dof[fmichel_wcount]= (dof_w);
                                      fmichel_gof[fmichel_wcount]= (gof_w);
                                      fmichel_ptminusRMS[fmichel_wcount]= (ptminusRMS_w);
                                      fmichel_ptplusRMS[fmichel_wcount]= (ptplusRMS_w);
                                      fmichel_status[fmichel_wcount]= (status_w);
                                      fmichel_cnnMichel[fmichel_wcount]= (cnnMichel_w);
                                      fmichel_cnnEM[fmichel_wcount]= (cnnEM_w);
                                      fmichel_cnnTrack[fmichel_wcount]= (cnnTrack_w);
                                      fmichel_wcount++;
                                      std::cout<<"wi_w "<<wi_w<<" tpc_w "<<tpc_w<<" ptime_w "<<ptime_w<<" mult_w "<<mult_w<<" michel_z "<<xyzWire.Z()<<" y "<<y_w<<" x "<<x_w<<" charge "<<charge<<" hit_chrge "<<chargehit_w<<" t1 "<<t1<<" t2 "<<t2<<" status "<<status_w<<"\n"<<std::endl;
                                      break;
                                    }//if( it_w == channel_no )   
                                  }// all recob::wire
                                  w1 ++; w2 ++;
                                  x ++; y ++;
                                  key++; wi++; chargehit++; tpc++; ptime++; angledeg++; maghitveccostheta++; distance++; mult++; sigptime++; sigchrg++; sigpamp++;
                                  dof++; gof++; ptminusRMS++; ptplusRMS++, status++;
                                }//wire from cone hits 

                                cone_key.clear(); 
                                cone_charge.clear();
                                cone_wire.clear();
                                cone_tpc.clear();
                                cone_ptime.clear();
                                cone_xpos.clear();
                                cone_ypos.clear();
                                cone_zpos.clear();
                                cone_angledeg.clear();
                                cone_maghitveccostheta.clear();
                                cone_distance.clear();
                                cone_mult.clear();
                                cone_sigptime.clear();
                                cone_sigchrg.clear();
                                cone_sigpamp.clear();
                                cone_dof.clear();
                                cone_gof.clear();
                                cone_ptminusRMS.clear();
                                cone_ptplusRMS.clear();
                                cone_status.clear();
                                cone_cnnMichel.clear();
                                cone_cnnEM.clear();
                                cone_cnnTrack.clear();
                                
                                //-----------getting the MC true information of the selected cone michel hits-----------------//
                                if(!evt.isRealData()) 
                                { 
                                  std::map<int,double> mconetrkide;
                                  int mconetrackid=-1;
                                  for(int qq =0; qq < fmhitcount; qq++)
                                  {
                                    art::Ptr<recob::Hit> mconehit=hitVect[fmhits_key[qq]];    
                                    std::vector<sim::TrackIDE> meveIDs = bt_serv->HitToTrackIDEs(clockData,mconehit);

                                    for(size_t e=0;e<meveIDs.size(); ++e)
                                    { 
                                      if(abs(meveIDs[e].trackID)!=fmccut_trkid) //the hits should not belong to the same mother track (muon)
                                      {
                                        mconetrkide[meveIDs[e].trackID] += meveIDs[e].energy;
                                        // cout<<"meveIDs[e].trackID "<<meveIDs[e].trackID<<" mconetrkide[meveIDs[e].trackID] "<<mconetrkide[meveIDs[e].trackID]<<"\n";
                                      }//if(abs(meveIDs[e].trackID)!=fmccut_trkid)
                                    }//for(size_t e=0;e<meveIDs.size(); ++e)
                                  }//for(int qq =0; qq < fmhitcount; qq++)
                                  double  mconemaxe = -1;
                                  double mconetote = 0;
                                  for(std::map<int,double>::iterator ii = mconetrkide.begin(); ii!=mconetrkide.end(); ++ii)
                                  {
                                    mconetote += ii->second;
                                    if((ii->second)>mconemaxe)
                                    {
                                      mconemaxe = ii->second;
                                      mconetrackid = abs(ii->first);
                                    }
                                  }//for(std::map<int,double>::iterator ii = mconetrkide.begin(); ii!=mconetrkide.end(); ++ii)
                                  
                                  const simb::MCParticle *mconeparticleP = pi_serv->TrackIdToParticle_P(mconetrackid);
                                  if(!mconeparticleP) continue;

                                  const art::Ptr<simb::MCTruth> mconemcP=pi_serv->TrackIdToMCTruth_P(mconetrackid);
                                  if(!mconemcP) continue;

                                  TLorentzVector mconemcstart, mconemcend, mconemcstartdrifted, mconemcenddrifted;
                                  unsigned int mconepstarti, mconependi, mconepstartdriftedi, mconependdriftedi; //mcparticle indices for starts and ends in tpc or drifted volumes
                                  double mconeplen = utils->length(*mconeparticleP, mconemcstart, mconemcend, mconepstarti, mconependi);
                                  double mconeplendrifted = utils->driftedLength(detProp, *mconeparticleP, mconemcstartdrifted, mconemcenddrifted, mconepstartdriftedi, mconependdriftedi);
                                  bool mconeisDrifted = mconeplendrifted!= 0;

                                  std::cout<<"Michel cone hits: bt:trackid "<<mconetrackid<<std::endl;
                                  fmcconehits_trkid      = mconeparticleP->TrackId();
                                  fmcconehits_vx         = mconeparticleP->Vx();
                                  fmcconehits_vy         = mconeparticleP->Vy();
                                  fmcconehits_vz         = mconeparticleP->Vz();
                                  fmcconehits_t          = mconeparticleP->T();
                                  fmcconehits_endx       = mconeparticleP->EndX();
                                  fmcconehits_endy       = mconeparticleP->EndY();
                                  fmcconehits_endz       = mconeparticleP->EndZ();
                                  fmcconehits_endt       = mconeparticleP->EndT();
                                  fmcconehits_px         = mconeparticleP->Px();
                                  fmcconehits_py         = mconeparticleP->Py();
                                  fmcconehits_pz         = mconeparticleP->Pz();
                                  fmcconehits_momentum   = mconeparticleP->P();
                                  fmcconehits_energy     = mconeparticleP->E();
                                  fmcconehits_endpx      = mconeparticleP->EndPx();
                                  fmcconehits_endpy      = mconeparticleP->EndPy();
                                  fmcconehits_endpz      = mconeparticleP->EndPz();
                                  fmcconehits_endenergy  = mconeparticleP->EndE();
                                  fmcconehits_pathlen    = mconeplen;
                                  
                                  if (mconeisDrifted)
                                  {
                                    fmcconehits_vxdrifted         = mconeparticleP->Vx();
                                    fmcconehits_vydrifted         = mconeparticleP->Vy();
                                    fmcconehits_vzdrifted         = mconeparticleP->Vz();
                                    fmcconehits_tdrifted          = mconeparticleP->T();
                                    fmcconehits_endxdrifted       = mconeparticleP->EndX();
                                    fmcconehits_endydrifted       = mconeparticleP->EndY();
                                    fmcconehits_endzdrifted       = mconeparticleP->EndZ();
                                    fmcconehits_endtdrifted       = mconeparticleP->EndT();
                                    fmcconehits_pxdrifted         = mconeparticleP->Px();
                                    fmcconehits_pydrifted         = mconeparticleP->Py();
                                    fmcconehits_pzdrifted         = mconeparticleP->Pz();
                                    fmcconehits_momentumdrifted   = mconeparticleP->P();
                                    fmcconehits_energydrifted     = mconeparticleP->E();
                                    fmcconehits_endpxdrifted      = mconeparticleP->EndPx();
                                    fmcconehits_endpydrifted      = mconeparticleP->EndPy();
                                    fmcconehits_endpzdrifted      = mconeparticleP->EndPz();
                                    fmcconehits_endenergydrifted  = mconeparticleP->EndE();
                                    fmcconehits_pathlendrifted    = mconeplendrifted;
                                  }

                                  fmcconehits_endprocess = int(mconeparticleP->Trajectory().ProcessToKey(mconeparticleP->EndProcess()));
                                  fmcconehits_theta      = mconeparticleP->Momentum().Theta();
                                  fmcconehits_phi        = mconeparticleP->Momentum().Phi();
                                  fmcconehits_pdg        = mconeparticleP->PdgCode();
                                  fmcconehits_status_code= mconeparticleP->StatusCode();
                                  fmcconehits_mass       = mconeparticleP->Mass();
                                  fmcconehits_ND         = mconeparticleP->NumberDaughters();
                                  fmcconehits_mother     = mconeparticleP->Mother();
                                  fmcconehits_origin     = mconemcP->Origin();
                                  fmcconehits_process    = int(mconeparticleP->Trajectory().ProcessToKey(mconeparticleP->Process()));
                                  fmcconehits_rescatter  = mconeparticleP->Rescatter();

                                  std::cout<<"MC true Michel cone hits: Pdg code "<<fmcconehits_pdg<<" track ID "<<fmcconehits_trkid<<" no of daughters "<<fmcconehits_ND<<" and origin_type "<<fmcconehits_origin<<std::endl;
                                
                                  std::vector<art::Ptr<recob::Hit>> shwrHits= shower_hits_assn.at(fcutshwr_key); //storing hits for ith track
                                  std::map<int,double> shwride;
                                  int shwrid=-1;

                                  for(size_t h=0; h<shwrHits.size();h++)
                                  {
                                    art::Ptr<recob::Hit> shhit=shwrHits[h];
                                    std::vector<sim::TrackIDE> sheveIDs = bt_serv->HitToTrackIDEs(clockData,shhit);
                                    for(size_t e=0;e<sheveIDs.size(); ++e)
                                    {
                                      if(abs(sheveIDs[e].trackID)!=fmccut_trkid) //the hits should not belong to the same mother particle (muon)
                                      {
                                        shwride[sheveIDs[e].trackID] += sheveIDs[e].energy;
                                      }//if(abs(sheveIDs[e].trackID)!=fmccut_trkid)
                                    }//for(size_t e=0;e<sheveIDs.size(); ++e)
                                  }//for(size_t h=0; h<shwrHits.size();h++)
                                  double  shmaxe = -1;
                                  double shtote = 0;
                                  for(std::map<int,double>::iterator ii = shwride.begin(); ii!=shwride.end(); ++ii)
                                  {
                                    //	cout<<" trkid "<<ii->first<<"  energy deposited = "<<ii->second<<std::endl;
                                    shtote += ii->second;
                                    if((ii->second)>shmaxe)
                                    {
                                      shmaxe = ii->second;
                                      shwrid = ii->first;
                                    }
                                  }//for(std::map<int,double>::iterator ii = shwride.begin(); ii!=shwride.end(); ++ii)
                            
                                  const simb::MCParticle *shparticleP = pi_serv->TrackIdToParticle_P(shwrid);
                                  if(!shparticleP) continue;
                                  const art::Ptr<simb::MCTruth> shmcP=pi_serv->TrackIdToMCTruth_P(shwrid);
                                  if(!shmcP) continue;

                                  TLorentzVector shmcstart, shmcend, shmcstartdrifted, shmcenddrifted;
                                  unsigned int shpstarti, shpendi, shpstartdriftedi, shpenddriftedi; //mcparticle indices for starts and ends in tpc or drifted volumes
                                  double shplen = utils->length(*shparticleP, shmcstart, shmcend, shpstarti, shpendi);
                                  double shplendrifted = utils->driftedLength(detProp, *shparticleP, shmcstartdrifted, shmcenddrifted, shpstartdriftedi, shpenddriftedi);
                                  bool isshDrifted = shplendrifted!= 0;

                                  fmcshwr_trkid     = shparticleP->TrackId();
                                  fmcshwr_vx        = shparticleP->Vx();
                                  fmcshwr_vy        = shparticleP->Vy();
                                  fmcshwr_vz        = shparticleP->Vz();
                                  fmcshwr_t         = shparticleP->T();
                                  fmcshwr_endx      = shparticleP->EndX();
                                  fmcshwr_endy      = shparticleP->EndY();
                                  fmcshwr_endz      = shparticleP->EndZ();
                                  fmcshwr_endt      = shparticleP->EndT();
                                  fmcshwr_px        = shparticleP->Px();
                                  fmcshwr_py        = shparticleP->Py();
                                  fmcshwr_pz        = shparticleP->Pz();
                                  fmcshwr_momentum  = shparticleP->P();
                                  fmcshwr_energy    = shparticleP->E();
                                  fmcshwr_endpx     = shparticleP->EndPx();
                                  fmcshwr_endpy     = shparticleP->EndPy();
                                  fmcshwr_endpz     = shparticleP->EndPz();
                                  fmcshwr_endenergy = shparticleP->EndE();
                                  fmcshwr_pathlen   = shplen;
                                  if (isshDrifted)
                                  {
                                    fmcshwr_vxdrifted        = shparticleP->Vx();
                                    fmcshwr_vydrifted        = shparticleP->Vy();
                                    fmcshwr_vzdrifted        = shparticleP->Vz();
                                    fmcshwr_tdrifted         = shparticleP->T();
                                    fmcshwr_endxdrifted      = shparticleP->EndX();
                                    fmcshwr_endydrifted      = shparticleP->EndY();
                                    fmcshwr_endzdrifted      = shparticleP->EndZ();
                                    fmcshwr_endtdrifted      = shparticleP->EndT();
                                    fmcshwr_pxdrifted        = shparticleP->Px();
                                    fmcshwr_pydrifted        = shparticleP->Py();
                                    fmcshwr_pzdrifted        = shparticleP->Pz();
                                    fmcshwr_momentumdrifted  = shparticleP->P();
                                    fmcshwr_energydrifted    = shparticleP->E();
                                    fmcshwr_endpxdrifted     = shparticleP->EndPx();
                                    fmcshwr_endpydrifted     = shparticleP->EndPy();
                                    fmcshwr_endpzdrifted     = shparticleP->EndPz();
                                    fmcshwr_endenergydrifted = shparticleP->EndE();
                                    fmcshwr_pathlendrifted   = shplendrifted;
                                  }//if (isshDrifted)
                                  fmcshwr_endprocess = int(shparticleP->Trajectory().ProcessToKey(shparticleP->EndProcess()));
                                  fmcshwr_theta      = shparticleP->Momentum().Theta();
                                  fmcshwr_phi        = shparticleP->Momentum().Phi();
                                  fmcshwr_pdg        = shparticleP->PdgCode();
                                  fmcshwr_status_code= shparticleP->StatusCode();
                                  fmcshwr_mass       = shparticleP->Mass();
                                  fmcshwr_ND         = shparticleP->NumberDaughters();
                                  fmcshwr_mother     = shparticleP->Mother();
                                  fmcshwr_origin     = shmcP->Origin();
                                  fmcshwr_process    = int(shparticleP->Trajectory().ProcessToKey(shparticleP->Process()));
                                  fmcshwr_rescatter  = shparticleP->Rescatter();

                                  std::cout<<"MC shower: Pdg code "<<fmcshwr_pdg<<" track ID "<<fmcshwr_trkid<<" no of daughters "<<fmcshwr_ND<<" and origin_type "<<fmcshwr_origin<<" energy "<<fmcshwr_energy<<std::endl;
                                }//if(!evt.isRealData())  
                      

                                double TPC_trigger_offset = 0.0;
                                // Get trigger to TPC Offset
                                // auto const* detclock = lar::providerFrom<detinfo::DetectorClocksService>();
                                TPC_trigger_offset = clockData.TriggerOffsetTPC();
                                std::cout << "TPC time offset from trigger: " << TPC_trigger_offset << "\n";

                                //No more external/internal optical flashes TODO: rewrite this
                                // float flsh_trk_time = 99999;
                                // totflash = 0;
                          
                                // for (size_t nn = 0; nn < flashlistInt.size() && nn < kMaxFlashes ; ++nn)
                                // {
                                //   art::Ptr<recob::OpFlash> pflashInt(flashListHandleInt, nn);
                                //   const recob::OpFlash& flashInt = *pflashInt;

                                //   // Find the closest in time external flash
                                //   double minflashtime = 9999999999999, minflashtime1 = 999; 
                                //   fext_trigger_time = 999;

                                //   for (size_t nm = 0; nm < flashlistExt.size() && nm < kMaxFlashes ; ++nm)
                                //   {
                                //     //loop over extrnal flashes
                                //     art::Ptr<recob::OpFlash> pmflashExt(flashListHandleExt,nm);
                                //     const recob::OpFlash& mflashExt = *pmflashExt;

                                //     minflashtime1 = abs(flashInt.Time() - mflashExt.Time());
                                //     if(minflashtime1 < minflashtime )
                                //     {
                                //       minflashtime = minflashtime1;
                                //       fext_trigger_time = mflashExt.Time();
                                //     }//if(minflashtime1 < minflashtime )
                                //   }//for (size_t nm = 0; nm < flashlistExt.size() && nm < kMaxFlashes ; ++nm)
                                                                                          
                                //   float fall_flash_time_diffs = (flashInt.Time()-fext_trigger_time) - (fcut_trkrecotime); // Time diff relative to trigger			  
                                //   if(abs(fall_flash_time_diffs)< flsh_trk_time && flashInt.TotalPE()>20)
                                //   {
                                //     flsh_trk_time = abs(fall_flash_time_diffs);
                                //     fflashcut_track_dist_int  = flsh_trk_time;
                                //     fflashcut_time_int        = flashInt.Time(); // Time relative to trigger
                                //     fflashcut_pe_int          = flashInt.TotalPE();
                                //     fflashcut_ycenter_int     = flashInt.YCenter(); // Geometric center in y
                                //     fflashcut_zcenter_int     = flashInt.ZCenter(); // Geometric center in z
                                //     fflashcut_ywidth_int      = flashInt.YWidth(); // Geometric width in y
                                //     fflashcut_zwidth_int      = flashInt.ZWidth(); // Geometric width in z
                                //     fflashcut_timewidth_int   = flashInt.TimeWidth(); // Width of the flash in time
                                //     fflashcut_abstime_int     = flashInt.AbsTime(); // Time by PMT readout clock (?)
                                //     fflashcut_frame_int       = flashInt.Frame(); // Frame number 
                                //   }//if(abs(fall_flash_time_diffs)< flsh_trk_time && flashInt.TotalPE()>20)
                                //   totflash++;
                                // }// loop over flashes
                                // std::cout<<"\n"<<"cut flash track time diff "<<fflashcut_track_dist_int<<" cut flash PE "<<fflashcut_pe_int<<"\n"<<std::endl;

                              }//if(nearhitct>= _minhitcountmichel && nearhitct< _maxhitcountmichel)
                            }//if(MinHitPeakTime > fMaxHitCountMichel)
                          }//if(MinHitPeakTime > fMinHitCountMichel)
                        }//if(track.Length()> fMuonTrackLengthCut )
                      }//if(bcount==0) //unbroken tracks
                    }//if((trkendz<APAnegbound1 || trkendz>APAposbound1) && (trkendz<APAnegbound2 || trkendz>APAposbound2))
                  }//if(ccrosser == 1)

                  CutReco->Fill();
                  fcut_mu++;

                  //loop over uncontained tracks
                  if((pos_ini.X()<Xnegbound || pos_ini.X()>Xposbound || pos_ini.Y()<Ynegbound || pos_ini.Y()>Yposbound || pos_ini.Z()<Znegbound || pos_ini.Z()>Zposbound) && 
                  (end.X()<Xnegbound || pos_end.X()<Xposbound || pos_end.Y()<Ynegbound || pos_end.Y()>Yposbound || pos_end.Z()<Znegbound || pos_end.Z()>Zposbound) && track.Length()>100)
                  {
                    funcont_trks++;
                  }
                        
                  //if track cross two z boundaries
                  if((pos_ini.Z()<Znegbound || pos_ini.Z()>Zposbound) && (pos_end.Z()<Znegbound || pos_end.Z()>Zposbound) && track.Length()>100)
                  {
                    fcrossZ_trks++;
                  }	

                }//if (endInFidVol)

                // // reset hit counts
                // for (int plane = 0; plane < 3; plane++) {fNCloseHitsIni[plane] = 0; fNCloseHitsEnd[plane] = 0;}
                // // loop over hits --> store the cutected
                // std::vector< recob::Hit > taggedHitsIni;
                // std::vector< recob::Hit > taggedHitsEnd;
                // for (size_t hdx = 0; hdx < hitResults.size(); ++hdx) 
                // {
                //   // cnn output vector: [em,trk,none,michel]
                //   std::array<float,MVA_LENGTH> cnn_out = hitResults.getOutput(hdx);
                //   // Looking at hits with large CNN values classified as michel
                //   // std::cout << "Your cutected fCNNThreshold is: " << fCNNThreshold << std::endl; //TO_DO: Find a way to save all the values given in the configuration file and print them just once
                //   if (cnn_out[hitResults.getIndex("michel")] > fCNNThreshold) 
                //   {
                //     // need to check at both ends of the track
                //     if (iniInFidVol) 
                //     {
                //       finiinFV_trks++;  // start in bounds tracks per event count
                //       // _fendinFV_trks++; // total start in bounds track count

                //       for (int plane = 0; plane < 3; plane++) // for each plane check if michel tagged hits are close to the track end
                //       {
                //         // project the track onto the plane
                //         geo::GeometryCore const & geom = *art::ServiceHandle<geo::Geometry>();
                //         geo::TPCID tpcid = geom.FindTPCAtPosition(ini_pt);
                //         if (!tpcid.isValid) continue;
                //         auto const & trackIni2D = pma::GetProjectionToPlane(ini,
                //                                                             plane,
                //                                                             geom.HasTPC(tpcid),
                //                                                             geom.PositionToCryostatID(ini_pt).Cryostat);
                //         double ini2D[2] = {trackIni2D.X(), trackIni2D.Y()};
                //         const recob::Hit & hit = hitResults.item(hdx); // get the hit itself
                //         if (hit.View() == plane) 
                //         {
                //             // std::cout << "hit.View() " << std::endl;
                //             // std::cout << ini_pt.X() << " " << ini_pt.Y() << " " << ini_pt.Z() << std::endl;
                //             // std::cout << ini2D[0] << " " << ini2D[1] << std::endl;

                //           // check that the hit is close to the track endpoint
                //           // and add to tagged hits if so
                //           if ( utils->hitCloseToTrackEnd( detProp, fRadiusThreshold, ini_pt, ini2D, hit, geom ) ) 
                //           {
                //             fNCloseHitsIni[plane] += 1;
                //             if (plane == 2) { taggedHitsIni.push_back( hit ); };
                //           }//if ( utils->hitCloseToTrackEnd( detProp, fRadiusThreshold, ini_pt, ini2D, hit, geom ) )
                //         }//if (hit.View() == plane)
                //       }//for ( int plane = 0; plane < 3; plane++ )
                //     }//if ( iniInFidVol )
                //     if (endInFidVol) 
                //     {
                //       // In each plane check if the michel tagged hits are close to the end of the track
                //       for (int plane = 0; plane < 3; plane++) 
                //       {
                //         // project the track onto the plane 
                //         geo::GeometryCore const & geom = *art::ServiceHandle<geo::Geometry>();
                //         geo::TPCID tpcid = geom.FindTPCAtPosition(end_pt);
                //         if (!tpcid.isValid) continue;
                //         auto const & trackEnd2D = pma::GetProjectionToPlane(end,
                //                                                             plane,
                //                                                             geom.HasTPC(tpcid),
                //                                                             geom.PositionToCryostatID(end_pt).Cryostat);
                //         double end2D[2] = {trackEnd2D.X(), trackEnd2D.Y()};
                //         const recob::Hit & hit = hitResults.item(hdx); // get the hit itself
                //         if (hit.View() == plane) {

                //           // check that the hit is close to the track endpoint
                //           // and add to tagged hits if so
                //           if ( utils->hitCloseToTrackEnd( detProp, fRadiusThreshold, end_pt, end2D, hit, geom ) ) 
                //           {
                //             fNCloseHitsEnd[plane] += 1;
                //             if (plane == 2) { taggedHitsEnd.push_back(hit); }
                //           }//if ( utils->hitCloseToTrackEnd( detProp, fRadiusThreshold, end_pt, end2D, hit, geom ) )
                //         }//if (hit.View() == plane)
                //       }//for ( int plane = 0; plane < 3; plane++ )
                //     }//if (endInFidVol)
                //   }//if (cnn_out[hitResults.getIndex("michel")] > fCNNThreshold)
                // } // end of loop over hits

                // std::cout << "fNCloseHitsIni[0] " << fNCloseHitsIni[0] << std::endl;
                // std::cout << "fNCloseHitsIni[1] " << fNCloseHitsIni[1] << std::endl;
                // std::cout << "fNCloseHitsIni[2] " << fNCloseHitsIni[2] << std::endl;
                // std::cout << "fNCloseHitsEnd[0] " << fNCloseHitsEnd[0] << std::endl;
                // std::cout << "fNCloseHitsEnd[1] " << fNCloseHitsEnd[1] << std::endl;
                // std::cout << "fNCloseHitsEnd[2] " << fNCloseHitsEnd[2] << std::endl;
                // // event selection decision: clusters of hits in all planes at one end of the track
                // bool iniSelected = ( fNCloseHitsIni[0] > fCloseHitsThreshold && fNCloseHitsIni[1] > fCloseHitsThreshold && fNCloseHitsIni[2] > fCloseHitsThreshold );
                // bool endSelected = ( fNCloseHitsEnd[0] > fCloseHitsThreshold && fNCloseHitsEnd[1] > fCloseHitsThreshold && fNCloseHitsEnd[2] > fCloseHitsThreshold );
                // bool evtSelected = ( iniSelected || endSelected );
                // if ( iniSelected && endSelected ) { evtSelected = false; }
                
                // // check if the event was correclty tagged
                // //  i.e. do the tagged hits correspond to a michel 
                // bool particleIsMichel(false);
                // // std::cout << "evtSelected " << evtSelected << std::endl;
                // // std::cout << "taggedHitsIni.size() " << taggedHitsIni.size() << std::endl;
                // // std::cout << "taggedHitsEnd.size() " << taggedHitsEnd.size() << std::endl;
                // if (evtSelected) 
                // {
                //   if ( iniSelected ) { particleIsMichel = utils->areHitsMichel(clockData, taggedHitsIni); }
                //   if ( endSelected ) { particleIsMichel = utils->areHitsMichel(clockData, taggedHitsEnd); }
                // }//if (evtSelected)
                // // std::cout << "particleIsMichel " << particleIsMichel << std::endl;
                // // update purity and efficiecny numbers and draw example events
                // if (  particleIsMichel && evtSelected ) { fYesSelected += 1; }
                // if ( !particleIsMichel && evtSelected ) { fNoSelected  += 1; }
              
              }//if (iniInFidVol || endInFidVol)
              
              NTracks++;
            
            }//if(T0_pfpsFromTrack_vect.size())
          }//if (pfpsFromTrack.size())
        }//for(auto const & track : tracks)

        HitsTracksKey.clear();
        AllReco->Fill();
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
