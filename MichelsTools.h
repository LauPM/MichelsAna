// ========================================================================================
// MichelsUtils.h
// This library is used as an auxiliar file for MichelsAna_module.cc.
// Here we define all the useful functions that may be applied in the module.
// 
// @authors     : Laura Pérez-Molina
// @created     : Apr, 2024 
//=========================================================================================

#ifndef MichelsAna_TOOLS_H
#define MichelsAna_TOOLS_H

#include <iostream>
#include <vector>
#include <fcntl.h>
#include <memory>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <sstream>
#include <fstream>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCParticle.h" // simb::MCParticle
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h" 
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"


#include "TH1I.h"
#include "TH1F.h"

extern const double ActiveBounds[6]; // Cryostat boundaries ( neg x, pos x, neg y, pos y, neg z, pos z )
extern const int Xnegbound;
extern const int Xposbound;
extern const int Ynegbound;
extern const int Yposbound;
extern const int Znegbound;
extern const int Zposbound;
extern const int APAnegbound1;
extern const int APAposbound1;
extern const int APAnegbound2;
extern const int APAposbound2;
// Limits in the analysis output
//TODO: remove unused and repeated variables
//TODO: Define from fcl
extern const float _longtrklen     ;
extern const float hits_track_length;
extern const float _minhitpeakcut  ;
extern const float _maxhitpeakcut  ;
extern const float _muendycut      ;
extern const float _muendzmincut   ;
extern const float _muendzmaxcut   ;
extern const float _absxdiff          ;
extern const int   _abszdiff          ;
extern const int   _hitdist           ;
extern const int   _minhitcountmichel ;
extern const int   _maxhitcountmichel ;
extern const int   _shwr_dist         ;
//Hit cone constants
extern const float _disthitend ;
extern const float _anglehitend;
extern const float _cone_length;
extern const float _cone_angle1;

struct ParticleData 
{
  double iniPosX = -999, iniPosY = -999, iniPosZ = -999, iniPosT = -999;
  double endPosX = -999, endPosY = -999, endPosZ = -999, endPosT = -999;
  double Energy  = -999, Moment  = -999;
  double iniMomX = -999, iniMomY = -999, iniMomZ = -999;
  double endMomX = -999, endMomY = -999, endMomZ = -999;
  int PdgCode = -999, TrackId = -999, MumTID = -999;
};
struct FilteredOutput 
{
  // int TrID = -999;
  // float Length = -999;
  // float Energy = -999;
  // float Moment = -999;
  // ParticleData iniPart;
  // float iniDist = -999;
  // float iniAngle = -999;
  // float iniDistEnd = -999;
  // ParticleData endPart;
  // std::string Type = "None";
}; 

namespace michels
{
    class MichelsBasicUtils
    {
        public:
            MichelsBasicUtils(fhicl::ParameterSet const &p);
            // ------------------- PERSONALIZED FUNCTIONS ------------------- //
            ParticleData extractParticleData(const simb::MCParticle* particle);
            std::tuple<int, float, float, float, float, ParticleData, ParticleData, std::string> MuonCandidateAnalysis(const art::Event & evt, int fNMichel, int i_part, art::ValidHandle<std::vector<simb::MCTruth>> AllParticles, const simb::MCTruth particle);
            
            //.......................... SELECTION ..........................//
            bool hitCloseToTrackEnd(detinfo::DetectorPropertiesData const& detProp,
                                    double radius, geo::Point_t const& end,
                                    double end2D[2], recob::Hit hit, 
                                    geo::GeometryCore const& geom); 
            bool IsPointInBounds(double *v, TVector3 const & p);
            int trackMatching( 	detinfo::DetectorClocksData const &  	clockData, int trackIndex, art::FindManyP<recob::Hit> hitsFromTracks );

            bool isMuonDecaying(const simb::MCParticle* particle, std::vector<art::Ptr<simb::MCParticle>, std::allocator<art::Ptr<simb::MCParticle> > > daughters);
            bool isMuonDecaying(const simb::MCParticle* particle);
            bool areHitsMichel(detinfo::DetectorClocksData const& clockData,
                     const std::vector<recob::Hit> & hits );
            //.......................... COMPUTATION ..........................//
            double length(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, 
                          unsigned int &starti, unsigned int &endi);
            // double driftedLength(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi);
            double driftedLength(detinfo::DetectorPropertiesData const& detProp, 
                                 const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, 
                                 unsigned int &starti, unsigned int &endi);
            void reset(bool deepClean);
            bool InMyMap(int TrID, std::map<int, simb::MCParticle> ParMap);
            unsigned int TrueParticleIniPointInAV(int *v, simb::MCParticle const & p);
            unsigned int TrueParticleEndPointInAV(int *v, simb::MCParticle const & p);
            void orderRecoStartEnd(TVector3 &start, TVector3 &end, 
                                   TVector3 &start_dir, TVector3 &end_dir);
            void SwitchEndPoints(TVector3 &start, TVector3 &end, 
                                 TVector3 &start_dir, TVector3 &end_dir);
            bool insideFidVol(geo::Point_t const &pos);
            //.......................... STYLE ..........................//
            void resume_stdout(int fd); // resumes the standard output
            void PrintInColor(std::string MyString, int Color); // prints a string in a given color
            int GetColor(std::string ColorName); // returns an integer that corresponds to a given color name
            // Conversions to string
            std::string str(int i);
            std::string str(unsigned int i);
            std::string str(double i);
            std::string str(float i);
            std::string str(std::vector<int> i);
            std::string str(std::vector<double> i);
            std::string str(std::vector<float> i);

            // --------------------- VARIABLES ------------------------ //
            double fFidVolCut; // Size of cut to select a fiducial volume
          
        // private:
        //     // From fhicl configuration
        //     const double ftry;
    };

    class MichelsCutManager
    {
      MichelsCutManager(fhicl::ParameterSet const &p);

      // void addCut(CutFunction cut);
      // void removeCut(int index);
      // std::vector<std::optional<CutData>> applyCuts(const EventData& event);
    };
}

#endif // MichelsAna_BASICTOOL_H