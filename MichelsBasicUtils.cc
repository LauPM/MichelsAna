#include "MichelsTools.h"

const double ActiveBounds[6] = {DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX};

// Geometrical variables
const int Xnegbound = -330;
const int Xposbound =  330;
const int Ynegbound =  50;
const int Yposbound =  550;
const int Znegbound =  50;
const int Zposbound =  645;
const int APAnegbound1 = 226-10;
const int APAposbound1 = 236+10;
const int APAnegbound2 = 456-10;    
const int APAposbound2 = 472+10;
const double fThicknessStartVolume = 30;

// Limits in the analysis output
//TODO: remove unused and repeated variables
//TODO: Define from fcl
const float _longtrklen        = 75;  //75
const float hits_track_length  = 10;  //75
const float _minhitpeakcut     = 200;
const float _maxhitpeakcut     = 5800;
const float _muendycut         = 80; //80 cm
const float _muendzmincut      = 80; //80 cm
const float _muendzmaxcut      = 610; //610 cm
const float _absxdiff          = 12; //in cm
const int   _abszdiff          = 10;  //in cm
const int   _hitdist           = 10; //in cm
const int   _minhitcountmichel = 5;  //5
const int   _maxhitcountmichel = 40;  //40
const int   _shwr_dist         = 10; //10
//Hit cone constants
const float _disthitend        = 80;
const float _anglehitend       = 0.523599; // 30 deg in radian
const float _cone_length       = 15;
const float _cone_angle1       = 60;


namespace michels

{

  MichelsBasicUtils::MichelsBasicUtils(fhicl::ParameterSet const &p) //: ftry(p.get<double>("try"))
  {};

  // TODO: Clean up the code and remove unused variables
  // TODO: Make sure all the functions are giving the expected output
  
  // ------------ YOUR BASIC FUNCTIONS DEFINITIONS HERE ------------ //  
  
  // Compute relevant information of your particle
  ParticleData MichelsBasicUtils::extractParticleData(const simb::MCParticle* particle)
  {
    ParticleData data;
    data.PdgCode = particle->PdgCode();
    data.TrackId = particle->TrackId();
    data.MumTID  = particle->Mother();
    data.iniPosX = particle->Position().X();
    data.iniPosY = particle->Position().Y();
    data.iniPosZ = particle->Position().Z();
    data.iniPosT = particle->T();
    data.endPosX = particle->EndX();
    data.endPosY = particle->EndY();
    data.endPosZ = particle->EndZ();
    data.endPosT = particle->EndT();
    data.Energy  = particle->E();
    data.Moment  = particle->P(); //std::sqrt(std::pow(Momentum(i).E(),2.) - std::pow(fmass,2.)
    data.iniMomX = particle->Px();
    data.iniMomY = particle->Py();
    data.iniMomZ = particle->Pz();
    data.endMomX = particle->EndPx();
    data.endMomY = particle->EndPy();
    data.endMomZ = particle->EndPz();

                
    return data;
  }

  std::tuple<int, float, float, float, float, ParticleData, ParticleData, std::string> MichelsBasicUtils::MuonCandidateAnalysis(const art::Event & evt, int fNMichel, int i_part, art::ValidHandle<std::vector<simb::MCTruth>> AllParticles, const simb::MCTruth particle)
  {
    ParticleData mumData, dauData;

    int Origin    = particle.Origin();     // Origin of the particles (0=unknown, 1=beam, 2=cosmic, 3=decay, 4=nucl. capture, 5=external process)
    int Daughters = particle.NParticles(); // Number of descendant particles in this MCTruth (fPartList=n particles in the event)
    const simb::MCParticle* mum_particle = &particle.GetParticle(i_part);
    float PDG = mum_particle->PdgCode();
    float TID = mum_particle->TrackId(); // Saving the TrackID of the generator particle
    float E   = mum_particle->E();
    float X   = mum_particle->Vx();
    float Y   = mum_particle->Vy();
    float Z   = mum_particle->Vz();

		std::string sInfo = "";
    sInfo = sInfo + "\nPrimary beam/cosmic/Ncapt [1/2/4]: " + str(Origin) +", (PDG,TID): (" + str(PDG) + "," + str(TID) + "), E: " + str(E) + " GeV";
    sInfo = sInfo + "\nVertex position (" + str(X) + ", " + str(Y) + ", " + str(Z) + ") cm";
    
    art::FindManyP<simb::MCParticle> mums(AllParticles, evt, "largeant"); //HARD CODED
    auto daughters = mums.at(i_part); // Get the list of MCParticle daughters ("largeant" label) of i_gen primary particle
    int N_all = daughters.size();     // Number of descendant particles in this MCParticle
    sInfo = sInfo + "\nPrimary Particle " + str(i_part) + " has " + str(Daughters) + " daughters. (Total particles: " + str(N_all) +")";
    
    geo::Point_t const mcEnd{ mum_particle->EndX(), mum_particle->EndY(), mum_particle->EndZ() };
    if (abs(PDG) == 13 && (mum_particle->NumberDaughters())>0) // Selecting muons
    {
      sInfo = sInfo + "\n";
      sInfo = sInfo + "\nMuon!";
      if ( isMuonDecaying( mum_particle, daughters ) && insideFidVol( mcEnd ) ) // Only for decaying muons we save more details
      {
        mumData = extractParticleData(mum_particle);
        // --- Loop over the MCParticle ALL descendants stored in trueParticle --- //
        for (int iDau = 1; iDau<N_all; iDau++) // all largeant objects (possible daughters)
        {
          auto daughtersPtr = daughters.at(iDau);
          const simb::MCParticle* daughter  = daughtersPtr.get();
          if (daughter->Mother() == TID) //genTrackID[i_part]
          {
            sInfo = sInfo + "\n";
            sInfo = sInfo + "\nPossible MC Michel candidate found!";
            fNMichel += 1;

            dauData = extractParticleData(daughter);
            // gen_g4IniMomX.push_back(iniMomX);
            // gen_g4IniMomY.push_back(iniMomY);
            // gen_g4IniMomZ.push_back(iniMomZ);
            // gen_g4EndMomX.push_back(endMomX);
            // gen_g4EndMomY.push_back(endMomY);
            // gen_g4EndMomZ.push_back(endMomZ);
          } //if (trueParticle.Mother() == 0)
        } //for loop filling MCInfo  
      } // if we have a muon decaying
    } // if the particle is a muon
    else
    {
      sInfo = sInfo + "\n";
      sInfo = sInfo + "\nNot selected event: PDG code is not muon. Not saving more details of this event on the tree.";
    }

    return std::make_tuple(fNMichel, PDG, TID, Origin, Daughters, mumData, dauData, sInfo);
  }
  

  //........................... SELECTION ...........................//
  // checks if a position is inside the fiducial volume
  bool MichelsBasicUtils::insideFidVol( geo::Point_t const &pos ) 
  {

    bool inside = true;
    // bool inside = false;
    // To_DO:solve the problem of the CryostatGeo
    // geo::GeometryCore const& fGeom = *art::ServiceHandle<geo::Geometry>();

    // geo::TPCID tpcID(fGeom.FindTPCAtPosition(pos));
    // if (fGeom.HasTPC(tpcID)) 
    // {
    //   const geo::TPCGeo& tpcgeo = fGeom.GetElement(tpcID);
    //   double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
    //   double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
    //   double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();
      
    //   for (size_t iCryo=0; iCryo<fGeom.Ncryostats(); iCryo++) 
    //   {
    //     const geo::CryostatGeo& cryostat = fGeom.Cryostat(iCryo);
    //     for (size_t iTPC=0; iTPC<cryostat.NTPC(); iTPC++) 
    //     {
    //       const geo::TPCGeo& tpc = cryostat.TPC(iTPC);
    //       if (tpc.MinX() < minx) minx = tpc.MinX();
    //       if (tpc.MaxX() > maxx) maxx = tpc.MaxX();
    //       if (tpc.MinY() < miny) miny = tpc.MinY();
    //       if (tpc.MaxY() > maxy) maxy = tpc.MaxY();
    //       if (tpc.MinZ() < minz) minz = tpc.MinZ();
    //       if (tpc.MaxZ() > maxz) maxz = tpc.MaxZ();
    //     }
    //   }

    //   double dista = fabs(minx - pos.X());
    //   double distb = fabs(pos.X() - maxx);
    //   if ((pos.X() > minx) && (pos.X() < maxx) &&
    //       (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
        
    //   dista = fabs(miny - pos.Y());
    //   distb = fabs(pos.Y()-maxy);
    //   if (inside && (pos.Y() > miny) && (pos.Y() < maxy) &&
    //       (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
    //   else inside = false;
      
    //   dista = fabs(minz - pos.Z());
    //   distb = fabs(pos.Z() - maxz);
    //   if (inside && (pos.Z() > minz) && (pos.Z() < maxz) &&
    //       (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
    //   else inside = false;
    // }
    return inside;
  }

  bool MichelsBasicUtils::IsPointInBounds(double *v, TVector3 const & p)  
  {
    return ((p.Y()>=(v[3]-fThicknessStartVolume) && p.Y()<=v[3]) || (p.X()>=v[0] && p.X()<=(v[0]+fThicknessStartVolume)) || (p.X()<=v[1] && p.X()>=(v[1]-fThicknessStartVolume)) || (p.Z()>=v[4] && p.Z()<=(v[4]+fThicknessStartVolume)) || (p.Z()<=v[5] && p.Z()>=(v[5]-fThicknessStartVolume)));
  }
  // Checks if the hits from a given track match the track with a given index
  int MichelsBasicUtils::trackMatching(detinfo::DetectorClocksData const &  clockData, 
                                       int trackIndex, art::FindManyP<recob::Hit> hitsFromTracks)
  {
    std::map<int,double> trackID_E;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    for (size_t h = 0; h < hitsFromTracks.at(trackIndex).size(); ++h)
    {
      for (auto const & id : bt_serv->HitToTrackIDEs(clockData,hitsFromTracks.at(trackIndex)[h]))
      {
        trackID_E[id.trackID] += id.energy;
      }
    }

    double max_e = 0.0; double tot_e = 0.0;
    int best_id = 0;
    for (std::map<int,double>::iterator it = trackID_E.begin(); it != trackID_E.end(); ++it)
    {
      tot_e += it->second;
      if (it->second > max_e)
      {
        max_e = it->second;
        best_id = it->first;
      }
    }

    if ((max_e > 0.0) && (true) && (tot_e > 0.0)) {return best_id;}
    else {return -999;}

  }

  // checks if the particle responsible for a track has an associated muon decay
  bool MichelsBasicUtils::isMuonDecaying(const simb::MCParticle* particle, std::vector<art::Ptr<simb::MCParticle>, std::allocator<art::Ptr<simb::MCParticle> > > daughters) 
  { 
    bool hasElectron = false, hasNuMu = false, hasNuE = false;
    for (auto const& daughtersPtr : daughters)
    {
      const simb::MCParticle* daughter  = daughtersPtr.get();
      if (daughter->Mother() == particle->TrackId()) 
      {
        int d_pdg = abs(daughter->PdgCode());
        if      (d_pdg == 11) hasElectron = true;
        else if (d_pdg == 12) hasNuE      = true;
        else if (d_pdg == 14) hasNuMu     = true;
      }
    }
    return (hasElectron && hasNuMu && hasNuE);
  }

  // checks if the particle responsible for a track has an associated muon decay
  bool MichelsBasicUtils::isMuonDecaying(const simb::MCParticle* particle) 
  { 
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    bool hasElectron = false; bool hasNuE = false; bool hasNuMu = false;
    for (int i = 0; i < particle->NumberDaughters(); ++i) 
    {
      const simb::MCParticle* daughter = pi_serv->TrackIdToParticle_P((particle->Daughter(i)));
      if (daughter && daughter->Mother() == particle->TrackId()) 
      {
        int d_pdg = abs(daughter->PdgCode());
        if      (d_pdg == 11) hasElectron = true;
        else if (d_pdg == 12) hasNuE      = true;
        else if (d_pdg == 14) hasNuMu     = true;
      }
    }
    return (hasElectron && hasNuMu && hasNuE);
  }

  // checks if a 2d hit is close to a given position
  bool MichelsBasicUtils::hitCloseToTrackEnd(detinfo::DetectorPropertiesData const& detProp,
                                      double radius, geo::Point_t const& end, double end2D[2], recob::Hit hit, geo::GeometryCore const & geom ) 
  {
    bool close = false;
    auto const & hitLocation = pma::WireDriftToCm(detProp,
                                                  hit.WireID().Wire, hit.PeakTime(), hit.View(), geom.FindTPCAtPosition(end).TPC, geom.PositionToCryostatID(end).Cryostat);

    double deltaXToHit = hitLocation.X() - end2D[0];
    double deltaYToHit = hitLocation.Y() - end2D[1];
    double displacementToHit = pow(pow(deltaXToHit,2)+pow(deltaYToHit,2),0.5);
    // std::cout << "Displacement to hit: " << displacementToHit << std::endl;
    if (displacementToHit < radius) close = true; 

    return close;
  }

  // checks if vector of hits come mostly from michel charge
  bool MichelsBasicUtils::areHitsMichel(detinfo::DetectorClocksData const& clockData,
                                const std::vector< recob::Hit > & hits ) 
  {
    const simb::MCParticle *  mcParticle=0;
    art::ServiceHandle< cheat::BackTrackerService > bt_serv;
    art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
    std::unordered_map< int, double > trkIDE;
    for ( auto const & hit : hits ) 
    {
      for ( auto const & ide : bt_serv->HitToTrackIDEs(clockData, hit) ) { trkIDE[ide.trackID] += ide.energy; }
    }
    
    int best_id(0);
    double tot_e(0), max_e(0);
    for ( auto const & contrib : trkIDE ) 
    {

      tot_e += contrib.second;
      if ( contrib.second > max_e ) 
      {
        max_e = contrib.second;
        best_id = contrib.first;
      }
    }

    if ( (max_e > 0) && (tot_e > 0) ) 
    {
      if ( best_id < 0 ) 
      {
        best_id = -best_id;
      }
      mcParticle = pi_serv->TrackIdToParticle_P( best_id );
    }
  
    if (mcParticle != 0) 
    { 
      return ( abs(mcParticle->PdgCode()) && ( mcParticle->Process() == "Decay" || mcParticle->Process() == "muMinusCaptureAtRest") ) ;
    }
    else 
    {
      std::cout << "No match found" << hits.size() << std::endl;
      return false; 
    }

  }

  //......................................................
  // This function get the length of the MC particle
  // Computed trajectory by trajectory (without the manual shifting for x correction)
  double MichelsBasicUtils::length(const simb::MCParticle& p, 
                                 TLorentzVector& start, 
                                 TLorentzVector& end, 
                                 unsigned int &starti, 
                                 unsigned int &endi)
  {
    double result = 0.;
    TVector3 disp;
    bool first = true;

    for(unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i) 
    {
      // check if the particle is inside a TPC
      if (p.Vx(i) >= Xnegbound && p.Vx(i) <= Xposbound && p.Vy(i) >= Ynegbound && p.Vy(i) <= Yposbound && p.Vz(i) >= Znegbound && p.Vz(i) <= Zposbound)
      {
        if(first)
        {
          start = p.Position(i);
          first = false;
          starti = i;
        }
        else
        {
          disp -= p.Position(i).Vect();
          result += disp.Mag();
        }
        disp = p.Position(i).Vect();
        end = p.Position(i);
        endi = i;
      }
    }
    return result;
  } // end of length function
  
  //......................................................
  // This function drifts the MCParticle trajectory for fair
  // comparisson between Reco and Truth i.e. not account for 
  // an interaction occuring with the beam dump

  // double MichelsBasicUtils::driftedLength(detinfo::DetectorPropertiesData const& detProp,
  //                                         const sim::MCTrack& mctrack, 
  //                                         TLorentzVector& tpcstart, 
  //                                         TLorentzVector& tpcend, 
  //                                         TLorentzVector& tpcmom){
  //     auto const* geom = lar::providerFrom<geo::Geometry>();
  
  //     //compute the drift x range
  //     double vDrift = detProp.DriftVelocity()*1e-3; //cm/ns
  //     double xrange[2] = {DBL_MAX, -DBL_MAX };
  //     for (unsigned int c=0; c<geom->Ncryostats(); ++c) {
  //       for (unsigned int t=0; t<geom->NTPC(c); ++t) {
  //         double Xat0 = detProp.ConvertTicksToX(0,0,t,c);
  //         double XatT = detProp.ConvertTicksToX(detProp.NumberTimeSamples(),0,t,c);
  //         xrange[0] = std::min(std::min(Xat0, XatT), xrange[0]);
  //         xrange[1] = std::max(std::max(Xat0, XatT), xrange[1]);
  //       }
  //     }
  
  //     double result = 0.;
  //     TVector3 disp;
  //     bool first = true;
  
  //     for(auto step: mctrack) {
  //       // check if the particle is inside a TPC
  //       if (step.X() >= ActiveBounds[0] && step.X() <= ActiveBounds[1] &&
  //           step.Y() >= ActiveBounds[2] && step.Y() <= ActiveBounds[3] &&
  //           step.Z() >= ActiveBounds[4] && step.Z() <= ActiveBounds[5] ){
  //         // Doing some manual shifting to account for
  //         // an interaction not occuring with the beam dump
  //         // we will reconstruct an x distance different from
  //         // where the particle actually passed to to the time
  //         // being different from in-spill interactions
  //         double newX = step.X()+(step.T()*vDrift);
  //         if (newX < xrange[0] || newX > xrange[1]) continue;
  
  //        TLorentzVector pos(newX,step.Y(),step.Z(),step.T());
  //         if(first){
  //           tpcstart = pos;
  //           tpcmom = step.Momentum();
  //           first = false;
  //         }
  //         else {
  //           disp -= pos.Vect();
  //           result += disp.Mag();
  //         }
  //         disp = pos.Vect();
  //         tpcend = pos;
  //       }
  //     }
  //     return result;
  //   }

  // Length of MC particle, trajectory by trajectory (with the manual shifting for x correction)
  double MichelsBasicUtils::driftedLength(detinfo::DetectorPropertiesData const& detProp,
                                          const simb::MCParticle& p, 
                                          TLorentzVector& start, 
                                          TLorentzVector& end, 
                                          unsigned int &starti, 
                                          unsigned int &endi)
  {
    const geo::Geometry* geom = &*art::ServiceHandle<geo::Geometry>(); //Geometry
    // auto const* geom = lar::providerFrom<geo::Geometry>();

    //compute the drift x range
    double vDrift = detProp.DriftVelocity()*1e-3; //cm/ns
    double xrange[2] = {DBL_MAX, -DBL_MAX };
    for (unsigned int c=0; c<geom->Ncryostats(); ++c) {
      for (unsigned int t=0; t<geom->NTPC(); ++t) {
        double Xat0 = detProp.ConvertTicksToX(0,0,t,c);
        double XatT = detProp.ConvertTicksToX(detProp.NumberTimeSamples(),0,t,c);
        xrange[0] = std::min(std::min(Xat0, XatT), xrange[0]);
        xrange[1] = std::max(std::max(Xat0, XatT), xrange[1]);
      }
    }

    double result = 0.;
    TVector3 disp;
    bool first = true;

    for(unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i) 
    {
      if (p.Vx(i) >= Xnegbound && p.Vx(i) <= Xposbound &&
      p.Vy(i) >= Ynegbound && p.Vy(i) <= Yposbound &&
      p.Vz(i) >= Znegbound && p.Vz(i) <= Zposbound)
      {
      // check if the particle is inside a TPC
      // if (p.Vx(i) >= ActiveBounds[0] && p.Vx(i) <= ActiveBounds[1] &&
      //     p.Vy(i) >= ActiveBounds[2] && p.Vy(i) <= ActiveBounds[3] &&
      //     p.Vz(i) >= ActiveBounds[4] && p.Vz(i) <= ActiveBounds[5]){
        // Doing some manual shifting to account for
        // an interaction not occuring with the beam dump
        // we will reconstruct an x distance different from
        // where the particle actually passed to to the time
        // being different from in-spill interactions
        double newX = p.Vx(i)+(p.T(i)*vDrift);
        if (newX < xrange[0] || newX > xrange[1]) continue;
        TLorentzVector pos(newX,p.Vy(i),p.Vz(i),p.T());
        if(first){
          start = pos;
          starti=i;
          first = false;
        }
        else {
          disp -= pos.Vect();
          result += disp.Mag();
        }
        disp = pos.Vect();
        end = pos;
        endi = i;
      }
    }
    return result;
  }

  void MichelsBasicUtils::reset(bool deepClean)
  {

    // fMCParticleParentTrackID.clear();
    return;
  }

    //......................................................
    // This function checks if a given TrackID is in a given map
    bool MichelsBasicUtils::InMyMap(int TrID, std::map<int, simb::MCParticle> ParMap)
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
    // This function
    unsigned int MichelsBasicUtils::TrueParticleIniPointInAV(int *v, simb::MCParticle const & p)  
    {
      for (unsigned int t = 0; t < p.NumberTrajectoryPoints(); t++)
      {
        if (p.Vx(t) >= v[0] && p.Vx(t) <= v[1] && p.Vy(t) >= v[2] && p.Vy(t) <= v[3] && p.Vz(t) >= v[4] && p.Vz(t) <= v[5])//
        return t;
      }
      return 999;
    }

    //......................................................
    // This function
    unsigned int MichelsBasicUtils::TrueParticleEndPointInAV(int *v, simb::MCParticle const & p) 
    {
      for (unsigned int t = p.NumberTrajectoryPoints()-1; t >= 0; t--)
      {
        if (p.Vx(t) >= v[0] && p.Vx(t) <= v[1] && p.Vy(t) >= v[2] && p.Vy(t) <= v[3] && p.Vz(t) >= v[4] && p.Vz(t) <= v[5])//
        return t;
      }
      return 999;
    }  

    //......................................................
    // This function 
    void MichelsBasicUtils::orderRecoStartEnd(TVector3 &start, 
                                              TVector3 &end, 
                                              TVector3 &start_dir, 
                                              TVector3 &end_dir) 
    {
      double prov_x, prov_y, prov_z;
      if (end.Y() > start.Y()) 
      {
        prov_x = start.X();
        prov_y = start.Y();
        prov_z = start.Z();
        start.SetXYZ(end.X(), end.Y(), end.Z());
        end.SetXYZ(prov_x, prov_y, prov_z);
        for (int k = 0; k<3; k++) 
        {
          TVector3 prov;
          prov = start_dir;
          start_dir = end_dir;
          end_dir = prov;
        }
        return;
      }
      else return;
    }
      
    //......................................................
    // Switch reco start and end point
    void MichelsBasicUtils::SwitchEndPoints(TVector3 &start, 
                                    TVector3 &end, 
                                    TVector3 &start_dir, 
                                    TVector3 &end_dir)
    {
      double prov_x, prov_y, prov_z;
      prov_x = start.X();
      prov_y = start.Y();
      prov_z = start.Z();
      start.SetXYZ(end.X(), end.Y(), end.Z());
      end.SetXYZ(prov_x, prov_y, prov_z);
      for (int k = 0; k<3; k++) 
      {
        TVector3 prov;
        prov = start_dir;
        start_dir = end_dir;
        end_dir = prov;
      }
      return;
    }

    //......................................................
    // UTILS
    // Initialize a vector which contains the limits of the Active Volume (taken from MicroBooNE AnalysisTree module)
    // //
    // //

    /*  //========================================================================
      //Modified backtracker 
      std::vector< art::Ptr< recob::Hit > > MichelsBasicUtils::TrackIdToHits_ps(detinfo::DetectorClocksData const & clockData,int tkId, std::vector< art::Ptr< recob::Hit >> const & hitsIn ) const
    {
        // returns a subset of the hits in the hitsIn collection that are matched
        // to the given track
    
        // temporary vector of TrackIds and Ptrs to hits so only one
        // loop through the (possibly large) hitsIn collection is needed
        std::vector<art::Ptr<recob::Hit>> hitList;
        std::vector<sim::TrackIDE> trackIDE;
        for (auto itr = hitsIn.begin(); itr != hitsIn.end(); ++itr) 
        {
          trackIDE.clear();
          art::Ptr<recob::Hit> const& hit = *itr;
          trackIDE = this.ChannelToTrackIDEs_ps(
            clockData, hit.Channel(), hit.PeakTimeMinusRMS(), hit.PeakTimePlusRMS());
          for (auto itr_trakIDE = trackIDE.begin(); itr_trakIDE != trackIDE.end(); ++itr_trakIDE)
          {
            if (itr_trakIDE.trackID == tkId && itr_trakIDE.energyFrac > fMinHitEnergyFraction)
              hitList.push_back(hit);
          } // itr_trakIDE
        }   // itr
        return hitList;
      }
      //---------------------------  
      std::vector<sim::TrackIDE>
      MichelsBasicUtils::ChannelToTrackIDEs_ps(detinfo::DetectorClocksData const& clockData,
                                      raw::ChannelID_t channel,
                                      const double hit_start_time,
                                      const double hit_end_time) const
      {
        std::vector<sim::TrackIDE> trackIDEs;
        double totalE = 0.;
        try {
          art::Ptr<sim::SimChannel> schannel = this.FindSimChannel_ps(channel);
    
          // loop over the electrons in the channel and grab those that are in time
          // with the identified hit start and stop times
          int start_tdc = clockData.TPCTick2TDC(hit_start_time);
          int end_tdc = clockData.TPCTick2TDC(hit_end_time);
          if (start_tdc < 0) start_tdc = 0;
          if (end_tdc < 0) end_tdc = 0;
          std::vector<sim::IDE> simides = schannel.TrackIDsAndEnergies(start_tdc, end_tdc);
    
          // first get the total energy represented by all track ids for
          // this channel and range of tdc values
          for (size_t e = 0; e < simides.size(); ++e)
            totalE += simides[e].energy;
    
          // protect against a divide by zero below
          if (totalE < 1.e-5) totalE = 1.;
    
          // loop over the entries in the map and fill the input vectors
    
          for (size_t e = 0; e < simides.size(); ++e) {
    
            if (simides[e].trackID == sim::NoParticleId) continue;
    
            sim::TrackIDE info;
            info.trackID = simides[e].trackID;
            info.energyFrac = simides[e].energy / totalE;
            info.energy = simides[e].energy;
            info.numElectrons = simides[e].numElectrons;
    
            trackIDEs.push_back(info);
          }
        } // end try
        catch (cet::exception const& e) {
          mf::LogWarning("BackTracker") << "caught exception \n" << e;
        }
    
        return trackIDEs;
      }   //--------------------------

      //-----------------------------------------------------------------------
      art::Ptr<sim::SimChannel>
      MichelsBasicUtilsy::FindSimChannel_ps(raw::ChannelID_t channel) const
      {
        art::Ptr<sim::SimChannel> chan;
        auto ilb = std::lower_bound(fSimChannels.begin(),
                                    fSimChannels.end(),
                                    channel,
                                    [](art::Ptr<sim::SimChannel> a, raw::ChannelID_t channel) {
                                      return (a.Channel() < channel);
                                    });
        if (ilb != fSimChannels.end())
          if ((*ilb).Channel() == channel) { chan = *ilb; }
        if (!chan)
          throw cet::exception("BackTracker") << "No sim::SimChannel corresponding "
                                              << "to channel: " << channel << "\n";
        return chan;
      }
      //-----------------   
    */  

    /*
    void MichelsBasicUtils::initActiveVol(double *fActiveBounds)  {
      fActiveBounds[0] = fActiveBounds[2] = fActiveBounds[4] = DBL_MAX;
      fActiveBounds[1] = fActiveBounds[3] = fActiveBounds[5] = -DBL_MAX;
      double abs_X_collection = 0;
      auto const* geom = lar::providerFrom<geo::Geometry>();
      for (geo::TPCGeo const& TPC: geom.IterateTPCs())  {
        double origin[3] = {0.};
        double center[3] = {0.};
        //geo::BoxBoundedGeo const& box = TPC.ActiveBoundingBox(); // It says the function does no exist...
        //double center[3] = {box.CenterX(), box.CenterY(), box.CenterZ()};
        TPC.LocalToWorld(origin, center); // had to modify CMakeLists.txt to make this work
        double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length()};

        //std::cout << "Width: " << 2*TPC.HalfWidth()
        //          << "\nHeight: " << 2*TPC.HalfHeight()
        //          << "\nLength: " << TPC.Length()
        //          << std::endl;

        if( center[0] - tpcDim[0] < fActiveBounds[0]) fActiveBounds[0] = center[0] - tpcDim[0];
        if( center[0] - tpcDim[0] > fActiveBounds[1]) fActiveBounds[1] = center[0] + tpcDim[0];
        if( center[1] - tpcDim[1] < fActiveBounds[2]) fActiveBounds[2] = center[1] - tpcDim[1];
        if( center[1] - tpcDim[1] > fActiveBounds[3]) fActiveBounds[3] = center[1] + tpcDim[1];
        if( center[2] - tpcDim[2] < fActiveBounds[4]) fActiveBounds[4] = center[2] - tpcDim[2];
        if( center[2] - tpcDim[2] > fActiveBounds[5]) fActiveBounds[5] = center[2] + tpcDim[2];

        //check coordinates of collection plane
        geo::PlaneGeo collectionPlane = TPC.LastPlane();
        double planeOrigin[3] = {0.};
        double planeCenter[3] = {0.};
        collectionPlane.LocalToWorld(planeOrigin, planeCenter);
        //std::cout << "++++++++++++++++++++" << std::setprecision(10) << std::endl
        //          << "Drift distance: " << TPC.DriftDistance() << std::endl
        //          << "x Collection Plane: " << planeCenter[0] << std::endl
        //          << "++++++++++++++++++++" << std::endl;
        if (TPC.DriftDistance() > 25)
          abs_X_collection = planeCenter[0];

        std::cout << "TPCs' drift lengths" << std::endl;
        std::cout << "TPC ID: " << TPC.ID() << " Drift distance: " << TPC.DriftDistance() << " Drift direction: " << TPC.DriftDirection() << std::endl;

      } // for all TPC

      std::cout << "Active Boundaries:  " << std::setprecision(10)
          << "\n\tx:  " << fActiveBounds[0] << " to " << fActiveBounds[1]
          << "\n\ty:  " << fActiveBounds[2] << " to " << fActiveBounds[3]
          << "\n\tz:  " << fActiveBounds[4] << " to " << fActiveBounds[5]
          << std::endl;
      std::cout << abs_X_collection << std::endl;
      fActiveBounds_eff[0] = -abs(abs_X_collection);
      fActiveBounds_eff[1] = abs(abs_X_collection);
      fActiveBounds_eff[2] = fActiveBounds[2];
      fActiveBounds_eff[3] = fActiveBounds[3];
      fActiveBounds_eff[4] = fActiveBounds[4];
      fActiveBounds_eff[5] = fActiveBounds[5];

    } 

    */



  // ------------ YOUR STYLE FUNCTIONS HERE ------------ //  

  //......................................................
  // This function creates a terminal color printout
  void MichelsBasicUtils::PrintInColor(std::string MyString, int Color)
  {
      mf::LogInfo("MichelsBasicUtils") << "\033[" << Color << "m" << MyString << "\033[0m";
    return;
  }

  // ......................................................
  // This function returns an integer that corresponds to a given color name
  int MichelsBasicUtils::GetColor(std::string ColorName)
  {
    if (ColorName == "black")               return 30;
    else if (ColorName == "red")            return 31;
    else if (ColorName == "green")          return 32;
    else if (ColorName == "yellow")         return 33;
    else if (ColorName == "blue")           return 34;
    else if (ColorName == "magenta")        return 35;
    else if (ColorName == "cyan")           return 36;
    else if (ColorName == "white")          return 37;
    else if (ColorName == "bright_black")   return 90;
    else if (ColorName == "bright_red")     return 91;
    else if (ColorName == "bright_green")   return 92;
    else if (ColorName == "bright_yellow")  return 93;
    else if (ColorName == "bright_blue")    return 94;
    else if (ColorName == "bright_magenta") return 95;
    else if (ColorName == "bright_cyan")    return 96;
    else if (ColorName == "bright_white")   return 97;
    else
    {
      mf::LogError("MichelsBasicUtils") << "Color " << ColorName << " not recognized. Returning white.";
      return 37;
    }
    return 0;
  }

  std::string MichelsBasicUtils::str(int i)
  {
    std::stringstream ss;
    ss << i;
    return ss.str();
  }
  std::string MichelsBasicUtils::str(unsigned int i)
  {
    std::stringstream ss;
    ss << i;
    return ss.str();
  }
  std::string MichelsBasicUtils::str(double i)
  {
    std::stringstream ss;
    ss << i;
    return ss.str();
  }
  std::string MichelsBasicUtils::str(float i)
  {
    std::stringstream ss;
    ss << i;
    return ss.str();
  }
  std::string MichelsBasicUtils::str(std::vector<int> i)
  {
    std::stringstream ss;
    for (int j = 0; j < int(i.size()); j++)
    {
      ss << i[j] << " ";
    }
    return ss.str();
  }
  std::string MichelsBasicUtils::str(std::vector<double> i)
  {
    std::stringstream ss;
    for (int j = 0; j < int(i.size()); j++)
    {
      ss << i[j] << " ";
    }
    return ss.str();
  }
  std::string MichelsBasicUtils::str(std::vector<float> i)
  {
    std::stringstream ss;
    for (int j = 0; j < int(i.size()); j++)
    {
      ss << i[j] << " ";
    }
    return ss.str();
  }

  void MichelsBasicUtils::resume_stdout(int fd)
  {
    std::fflush(stdout);
    dup2(fd, 1);
    close(fd);
  }


} // end of namespace michels