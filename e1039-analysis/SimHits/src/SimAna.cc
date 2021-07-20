#include <TFile.h>
#include <TTree.h>
#include <iomanip>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <interface_main/SQDimuonVector_v1.h>
#include <interface_main/SQHit.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQTrackVector_v1.h>
#include <interface_main/SQDimuonVector_v1.h>

#include <phool/getClass.h>

#include <ktracker/SRecEvent.h>

#include <interface_main/SQEvent_v1.h>
#include <interface_main/SQMCEvent.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include "SimAna.h"

SimAna::SimAna(const std::string& name): SubsysReco(name), legacyContainer(true)
{}

SimAna::~SimAna() {}

int SimAna::Init(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int SimAna::InitRun(PHCompositeNode *topNode) {
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
    return ret;

  eventID = 0;
  ResetEvalVars();
  MakeTree();

  return Fun4AllReturnCodes::EVENT_OK;
}

// find one g4 hit at station
PHG4Hit* SimAna::FindG4HitAtStation(const int trk_id, const PHG4HitContainer* g4hc) {
  PHG4Hit* hit = nullptr;
  PHG4HitContainer::ConstRange range = g4hc->getHits();
  for (PHG4HitContainer::ConstIterator it = range.first; it != range.second;
       it++) {
    PHG4Hit *tmphit = it->second;
    if (tmphit->get_trkid() == trk_id) {
      hit = tmphit;
      break;
    }
  }
  return hit;
}

// find several g4 hits at one station (e.g for emcal)
std::vector<PHG4Hit*> SimAna::FindG4HitsAtStation(const int trk_id, const PHG4HitContainer* g4hc) {
  std::vector<PHG4Hit*> vhit;
  PHG4HitContainer::ConstRange range = g4hc->getHits();
  for (PHG4HitContainer::ConstIterator it = range.first; it != range.second;
       it++) {
    PHG4Hit *tmphit = it->second;
    if (tmphit->get_trkid() == trk_id) {
      vhit.push_back(tmphit);
    }
  }
  return vhit;
}

// find g4shower
PHG4Shower *SimAna::get_primary_shower(PHG4Particle *primary) {
  PHG4Shower *shower = nullptr;
  for (auto iter = _truth->GetPrimaryShowerRange().first;
       iter != _truth->GetPrimaryShowerRange().second; ++iter) {
    PHG4Shower *tmpshower = iter->second;
    if (tmpshower->get_parent_particle_id() == primary->get_track_id()) {
      shower = tmpshower;
      break;
    }
  }
  return shower;
}

// find best reco track
SRecTrack* SimAna::FindBestMomRecTrack(SRecEvent *recEvent,  const float true_TargetP) {
  double dP = 100.; // delta(momentum)
  double hold_dP = 99999.;
  
  SRecTrack* Best_recTrack =  NULL;
  for(int itrack=0; itrack<recEvent->getNTracks(); ++itrack){   
    if (hold_dP>dP) hold_dP = dP;
    SRecTrack *recTrack = &recEvent->getTrack(itrack);
    dP = fabs(true_TargetP - recTrack->getTargetMom().Mag());
   
    //Finding out best match track in terms of energy
    if(dP-hold_dP<0.) Best_recTrack = recTrack;  
  }
  return Best_recTrack;
}

// find common ids for reco and truth tracks
int SimAna::FindCommonHitIDs(std::vector<int>& hitidvec1, std::vector<int>& hitidvec2) {
  //This function assumes the input vectors have been sorted
  auto iter = hitidvec1.begin();
  auto jter = hitidvec2.begin();

  int nCommon = 0;
  while(iter != hitidvec1.end() && jter != hitidvec2.end()) {
    if(*iter < *jter) {
      ++iter;
    } else {
      if(!(*jter < *iter)) {
        ++nCommon;
        ++iter;
      }
      ++jter;
    }
  }

  return nCommon;
}

void SimAna::set_legacy_rec_container(bool b) { 
  legacyContainer = b; 
}

int SimAna::ResetEvalVars() {
  // hits
  for(int i=0; i<1000; ++i) {
    hit_detid[i]        = std::numeric_limits<short>::max();
    hit_elmid[i]        = std::numeric_limits<short>::max();
    hit_trkid[i]        = std::numeric_limits<short>::max();
    hit_driftdis[i]     = std::numeric_limits<float>::max();
    hit_pos[i]          = std::numeric_limits<float>::max();
    hit_edep[i]         = std::numeric_limits<float>::max();
    hit_truthx[i]       = std::numeric_limits<float>::max();
    hit_truthy[i]       = std::numeric_limits<float>::max();
    hit_truthz[i]       = std::numeric_limits<float>::max();
    hit_truthpx[i]      = std::numeric_limits<float>::max();
    hit_truthpy[i]      = std::numeric_limits<float>::max();
    hit_truthpz[i]      = std::numeric_limits<float>::max();
  }

  // truth tracks and rec tracks
  n_truthtracks = 0;
  n_tracks = 0;

  for(int i=0; i<100; ++i) {
    truthtrack_charge[i] = std::numeric_limits<int>::max();
    truthtrack_x_st1[i] = std::numeric_limits<int>::max();
    truthtrack_y_st1[i] = std::numeric_limits<int>::max();
    truthtrack_z_st1[i] = std::numeric_limits<int>::max();
    truthtrack_px_st1[i] = std::numeric_limits<int>::max();
    truthtrack_py_st1[i] = std::numeric_limits<int>::max();
    truthtrack_pz_st1[i] = std::numeric_limits<int>::max();
    truthtrack_x_st3[i] = std::numeric_limits<int>::max();
    truthtrack_y_st3[i] = std::numeric_limits<int>::max();
    truthtrack_z_st3[i] = std::numeric_limits<int>::max();
    truthtrack_px_st3[i] = std::numeric_limits<int>::max();
    truthtrack_py_st3[i] = std::numeric_limits<int>::max();
    truthtrack_pz_st3[i] = std::numeric_limits<int>::max();
    truthtrack_x_vtx[i] = std::numeric_limits<int>::max();
    truthtrack_y_vtx[i] = std::numeric_limits<int>::max();
    truthtrack_z_vtx[i] = std::numeric_limits<int>::max();
    truthtrack_px_vtx[i] = std::numeric_limits<int>::max();
    truthtrack_py_vtx[i] = std::numeric_limits<int>::max();
    truthtrack_pz_vtx[i] = std::numeric_limits<int>::max();
  }

  for(int i=0; i<100; ++i) {
    track_charge[i]      = std::numeric_limits<int>::max();
    track_nhits[i]       = std::numeric_limits<int>::max();
    track_x_target[i]    = std::numeric_limits<float>::max();
    track_y_target[i]    = std::numeric_limits<float>::max();
    track_z_target[i]    = std::numeric_limits<float>::max();
    track_px_target[i]   = std::numeric_limits<float>::max();
    track_py_target[i]   = std::numeric_limits<float>::max();
    track_pz_target[i]   = std::numeric_limits<float>::max();
    track_x_st1[i]       = std::numeric_limits<float>::max();
    track_y_st1[i]       = std::numeric_limits<float>::max();
    track_z_st1[i]       = std::numeric_limits<float>::max();
    track_px_st1[i]      = std::numeric_limits<float>::max();
    track_py_st1[i]      = std::numeric_limits<float>::max();
    track_pz_st1[i]      = std::numeric_limits<float>::max();
    track_x_vtx[i]       = std::numeric_limits<float>::max();
    track_y_vtx[i]       = std::numeric_limits<float>::max();
    track_z_vtx[i]       = std::numeric_limits<float>::max();
    track_px_vtx[i]      = std::numeric_limits<float>::max();
    track_py_vtx[i]      = std::numeric_limits<float>::max();
    track_pz_vtx[i]      = std::numeric_limits<float>::max();
    track_m[i]           = std::numeric_limits<float>::max();
    track_chisq[i]       = std::numeric_limits<float>::max();
    track_prob[i]        = std::numeric_limits<float>::max();
    track_quality[i]     = std::numeric_limits<float>::max();
    track_nhits_st1[i]   = std::numeric_limits<float>::max();
    track_nhits_st2[i]   = std::numeric_limits<float>::max();
    track_nhits_st3[i]   = std::numeric_limits<float>::max();
  }

  // dimuons
  n_truthdimuons = 0;
  for(int i=0; i<100; ++i) {
    truthdimuon_mass[i]        = std::numeric_limits<float>::max();
    truthdimuon_x_vtx[i]       = std::numeric_limits<float>::max();
    truthdimuon_y_vtx[i]       = std::numeric_limits<float>::max();
    truthdimuon_z_vtx[i]       = std::numeric_limits<float>::max();
    truthdimuon_px[i]          = std::numeric_limits<float>::max();
    truthdimuon_py[i]          = std::numeric_limits<float>::max();
    truthdimuon_pz[i]          = std::numeric_limits<float>::max();
    truthdimuon_pmom_x[i]      = std::numeric_limits<float>::max();
    truthdimuon_pmom_y[i]      = std::numeric_limits<float>::max();
    truthdimuon_pmom_z[i]      = std::numeric_limits<float>::max();
    truthdimuon_nmom_x[i]      = std::numeric_limits<float>::max();
    truthdimuon_nmom_y[i]      = std::numeric_limits<float>::max();
    truthdimuon_nmom_z[i]      = std::numeric_limits<float>::max();
  }

  n_dimuons = 0;
  for(int i=0; i<100; ++i) {
    dimuon_mass[i]       = std::numeric_limits<float>::max();
    dimuon_chisq[i]      = std::numeric_limits<float>::max();
    dimuon_x_vtx[i]      = std::numeric_limits<float>::max();
    dimuon_y_vtx[i]      = std::numeric_limits<float>::max();
    dimuon_z_vtx[i]      = std::numeric_limits<float>::max();
    dimuon_px[i]         = std::numeric_limits<float>::max();
    dimuon_py[i]         = std::numeric_limits<float>::max();
    dimuon_pz[i]         = std::numeric_limits<float>::max();
    dimuon_pmom_x[i]     = std::numeric_limits<float>::max();
    dimuon_pmom_y[i]     = std::numeric_limits<float>::max();
    dimuon_pmom_z[i]     = std::numeric_limits<float>::max();
    dimuon_nmom_x[i]     = std::numeric_limits<float>::max();
    dimuon_nmom_y[i]     = std::numeric_limits<float>::max();
    dimuon_nmom_z[i]     = std::numeric_limits<float>::max();
    dimuon_ppos_x[i]     = std::numeric_limits<float>::max();
    dimuon_ppos_y[i]     = std::numeric_limits<float>::max();
    dimuon_ppos_z[i]     = std::numeric_limits<float>::max();
    dimuon_npos_x[i]     = std::numeric_limits<float>::max();
    dimuon_npos_y[i]     = std::numeric_limits<float>::max();
    dimuon_npos_z[i]     = std::numeric_limits<float>::max();
  }

  // truth information
  // g4 showers
  n_showers =0;
  for(int i=0; i<1000; ++i) {
    sx_ecal[i]          = std::numeric_limits<float>::max();
    sy_ecal[i]          = std::numeric_limits<float>::max();
    sz_ecal[i]          = std::numeric_limits<float>::max();
    sedep_ecal[i]       = std::numeric_limits<float>::max();
  }

  // primary hits
  n_primaries = 0;
  for (int i = 0; i < 1000; ++i) {
    gtrkid[i] = std::numeric_limits<int>::max();
    gpid[i] = std::numeric_limits<int>::max();
    gvx[i] = std::numeric_limits<float>::max();
    gvy[i] = std::numeric_limits<float>::max();
    gvz[i] = std::numeric_limits<float>::max();
    gpx[i] = std::numeric_limits<float>::max();
    gpy[i] = std::numeric_limits<float>::max();
    gpz[i] = std::numeric_limits<float>::max();
    gpt[i] = std::numeric_limits<float>::max();
    geta[i] = std::numeric_limits<float>::max();
    gphi[i] = std::numeric_limits<float>::max();
    ge[i] = std::numeric_limits<float>::max();

    nhits_ecal[i] = 0;
    for (int j = 0; j < 100; ++j) {
      gx_ecal[i][j] = std::numeric_limits<int>::max();
      gy_ecal[i][j] = std::numeric_limits<int>::max();
      gz_ecal[i][j] = std::numeric_limits<int>::max();
      gpx_ecal[i][j] = std::numeric_limits<int>::max();
      gpy_ecal[i][j] = std::numeric_limits<int>::max();
      gpz_ecal[i][j] = std::numeric_limits<int>::max();
      gedep_ecal[i][j] = std::numeric_limits<int>::max();
    }

    gx_st1[i] = std::numeric_limits<float>::max();
    gy_st1[i] = std::numeric_limits<float>::max();
    gz_st1[i] = std::numeric_limits<float>::max();
    gpx_st1[i] = std::numeric_limits<float>::max();
    gpy_st1[i] = std::numeric_limits<float>::max();
    gpz_st1[i] = std::numeric_limits<float>::max();
    gx_st2[i] = std::numeric_limits<float>::max();
    gy_st2[i] = std::numeric_limits<float>::max();
    gz_st2[i] = std::numeric_limits<float>::max();
    gpx_st2[i] = std::numeric_limits<float>::max();
    gpy_st2[i] = std::numeric_limits<float>::max();
    gpz_st2[i] = std::numeric_limits<float>::max();

    gx_h1[i] = std::numeric_limits<float>::max();
    gy_h1[i] = std::numeric_limits<float>::max();
    gz_h1[i] = std::numeric_limits<float>::max();
    gpx_h1[i] = std::numeric_limits<float>::max();
    gpy_h1[i] = std::numeric_limits<float>::max();
    gpz_h1[i] = std::numeric_limits<float>::max();
    gx_h2[i] = std::numeric_limits<float>::max();
    gy_h2[i] = std::numeric_limits<float>::max();
    gz_h2[i] = std::numeric_limits<float>::max();
    gpx_h2[i] = std::numeric_limits<float>::max();
    gpy_h2[i] = std::numeric_limits<float>::max();
    gpz_h2[i] = std::numeric_limits<float>::max();
    gx_h3[i] = std::numeric_limits<float>::max();
    gy_h3[i] = std::numeric_limits<float>::max();
    gz_h3[i] = std::numeric_limits<float>::max();
    gpx_h3[i] = std::numeric_limits<float>::max();
    gpy_h3[i] = std::numeric_limits<float>::max();
    gpz_h3[i] = std::numeric_limits<float>::max();
    gx_h4[i] = std::numeric_limits<float>::max();
    gy_h4[i] = std::numeric_limits<float>::max();
    gz_h4[i] = std::numeric_limits<float>::max();
    gpx_h4[i] = std::numeric_limits<float>::max();
    gpy_h4[i] = std::numeric_limits<float>::max();
    gpz_h4[i] = std::numeric_limits<float>::max();
  }

  for (int i = 0; i < 5; ++i) {
    fpga_trigger[i] = 0;
  }

  return 0;
}

int SimAna::process_event(PHCompositeNode* topNode) {

  ResetEvalVars();

  // filling arrays
  n_hits = 0;
  for(int ihit=0; ihit<_hitVector->size(); ++ihit) {
    SQHit *hit = _hitVector->at(ihit);
    int hitID = hit->get_hit_id();
    hit_detid[n_hits]      = hit->get_detector_id();
    hit_elmid[n_hits]      = hit->get_element_id();
    hit_trkid[n_hits]      = hit->get_track_id();
    hit_driftdis[n_hits]   = hit->get_drift_distance();
    hit_pos[n_hits]        = hit->get_pos();
    hit_edep[n_hits]       = hit->get_edep();
    hit_truthx[n_hits]     = hit->get_truth_x();
    hit_truthy[n_hits]     = hit->get_truth_y();
    hit_truthz[n_hits]     = hit->get_truth_z();
    hit_truthpx[n_hits]    = hit->get_truth_px();
    hit_truthpy[n_hits]    = hit->get_truth_py();
    hit_truthpz[n_hits]    = hit->get_truth_pz();
    ++n_hits;
    if(n_hits>=1000) break;
  }

  // tracks
  int n_recTracks = legacyContainer ? _recEvent->getNTracks() : _recTrackVector->size();
  n_truthtracks = 0;
  n_tracks = 0;
  for(int itrk = 0; itrk < _trackVector->size(); ++itrk) {
    SQTrack* track = _trackVector->at(itrk);

    truthtrack_charge[n_truthtracks] = track->get_charge();
    truthtrack_x_st1[n_truthtracks] = (track->get_pos_st1()).X();
    truthtrack_y_st1[n_truthtracks] = (track->get_pos_st1()).Y();
    truthtrack_z_st1[n_truthtracks] = (track->get_pos_st1()).Z();
    truthtrack_px_st1[n_truthtracks] = (track->get_mom_st1()).Px();
    truthtrack_py_st1[n_truthtracks] = (track->get_mom_st1()).Py();
    truthtrack_pz_st1[n_truthtracks] = (track->get_mom_st1()).Pz();
    truthtrack_x_st3[n_truthtracks] = (track->get_pos_st3()).X();
    truthtrack_y_st3[n_truthtracks] = (track->get_pos_st3()).Y();
    truthtrack_z_st3[n_truthtracks] = (track->get_pos_st3()).Z();
    truthtrack_px_st3[n_truthtracks] = (track->get_mom_st3()).Px();
    truthtrack_py_st3[n_truthtracks] = (track->get_mom_st3()).Py();
    truthtrack_pz_st3[n_truthtracks] = (track->get_mom_st3()).Pz();
    truthtrack_x_vtx[n_truthtracks] = (track->get_pos_vtx()).X();
    truthtrack_y_vtx[n_truthtracks] = (track->get_pos_vtx()).Y();
    truthtrack_z_vtx[n_truthtracks] = (track->get_pos_vtx()).Z();
    truthtrack_px_vtx[n_truthtracks] = (track->get_mom_vtx()).Px();
    truthtrack_py_vtx[n_truthtracks] = (track->get_mom_vtx()).Py();
    truthtrack_pz_vtx[n_truthtracks] = (track->get_mom_vtx()).Pz();

    int recid = track->get_rec_track_id();
    if(recid >= 0 && recid < n_recTracks) {
      SRecTrack* recTrack = legacyContainer ? &(_recEvent->getTrack(recid)) : dynamic_cast<SRecTrack*>(_recTrackVector->at(recid));
      //std::cout << "******************** (recTrack->getTargetMom()).Px() " << (recTrack->getTargetMom()).Px() << std::endl;
      track_charge[n_tracks] = recTrack->getCharge();
      track_nhits[n_tracks] = recTrack->getNHits();
      track_x_target[n_tracks] = (recTrack->getTargetPos()).X();
      track_y_target[n_tracks] = (recTrack->getTargetPos()).Y();
      track_z_target[n_tracks] = (recTrack->getTargetPos()).Z();
      track_px_target[n_tracks] = (recTrack->getTargetMom()).Px();
      track_py_target[n_tracks] = (recTrack->getTargetMom()).Py();
      track_pz_target[n_tracks] = (recTrack->getTargetMom()).Pz();
      track_x_st1[n_tracks] = (recTrack->getPositionVecSt1()).X();
      track_y_st1[n_tracks] = (recTrack->getPositionVecSt1()).Y();
      track_z_st1[n_tracks] = (recTrack->getPositionVecSt1()).Z();
      track_px_st1[n_tracks] = (recTrack->getMomentumVecSt1()).Px();
      track_py_st1[n_tracks] = (recTrack->getMomentumVecSt1()).Py();
      track_pz_st1[n_tracks] = (recTrack->getMomentumVecSt1()).Pz();
      track_x_vtx[n_tracks] = (recTrack->getVertexPos()).X();
      track_y_vtx[n_tracks] = (recTrack->getVertexPos()).Y();
      track_z_vtx[n_tracks] = (recTrack->getVertexPos()).Z();
      track_px_vtx[n_tracks] = (recTrack->getVertexMom()).X();
      track_py_vtx[n_tracks] = (recTrack->getVertexMom()).Y();
      track_pz_vtx[n_tracks] = (recTrack->getVertexMom()).Z();
      track_chisq[n_tracks] = recTrack->getChisq();
      track_prob[n_tracks] = recTrack->getProb();
      track_quality[n_tracks] = recTrack->getQuality();
      track_nhits_st1[n_tracks] = recTrack->getNHitsInStation(1);
      track_nhits_st2[n_tracks] = recTrack->getNHitsInStation(2);
      track_nhits_st3[n_tracks] = recTrack->getNHitsInStation(3);
      ++n_tracks;
      if(n_tracks >= 100) break;
    }
    ++n_truthtracks;
    if (n_truthtracks >= 100)
      break;
  }

  // vertices
  n_truthdimuons = 0;
  int nDimuons = _dimuonVector->size();
  for(int i = 0; i < nDimuons; ++i) {
    SQDimuon* dimuon = _dimuonVector->at(i);
    // truth dimuon
    // from https://github.com/E1039-Collaboration/e1039-core/blob/master/simulation/g4dst/TruthNodeMaker.cc#L133-L155
    truthdimuon_mass[n_truthdimuons] = dimuon->get_mom().M();
    truthdimuon_x_vtx[n_truthdimuons] = (dimuon->get_pos()).X();
    truthdimuon_y_vtx[n_truthdimuons] = (dimuon->get_pos()).Y();
    truthdimuon_z_vtx[n_truthdimuons] = (dimuon->get_pos()).Z();
    truthdimuon_px[n_truthdimuons] = (dimuon->get_mom()).X();
    truthdimuon_py[n_truthdimuons] = (dimuon->get_mom()).Y();
    truthdimuon_pz[n_truthdimuons] = (dimuon->get_mom()).Z();
    truthdimuon_pmom_x[n_truthdimuons] = (dimuon->get_mom_pos()).Px();
    truthdimuon_pmom_y[n_truthdimuons] = (dimuon->get_mom_pos()).Py();
    truthdimuon_pmom_z[n_truthdimuons] = (dimuon->get_mom_pos()).Pz();
    truthdimuon_nmom_x[n_truthdimuons] = (dimuon->get_mom_neg()).Px();
    truthdimuon_nmom_y[n_truthdimuons] = (dimuon->get_mom_neg()).Py();
    truthdimuon_nmom_z[n_truthdimuons] = (dimuon->get_mom_neg()).Pz();
    ++n_truthdimuons;
    if(n_truthdimuons >= 100) break;
  }

  n_dimuons = 0;
  int nRecDimuons = legacyContainer ? _recEvent->getNDimuons() : (_recDimuonVector ? _recDimuonVector->size() : -1);
  for(int i = 0; i < nRecDimuons; ++i) {
    SRecDimuon* recDimuon = legacyContainer ? &(_recEvent->getDimuon(i)) : dynamic_cast<SRecDimuon*>(_recDimuonVector->at(i));
    dimuon_mass[n_dimuons] = recDimuon->mass;
    dimuon_chisq[n_dimuons] = recDimuon->get_chisq();
    dimuon_x_vtx[n_dimuons] = (recDimuon->vtx).X();
    dimuon_y_vtx[n_dimuons] = (recDimuon->vtx).Y();
    dimuon_z_vtx[n_dimuons] = (recDimuon->vtx).Z();
    dimuon_px[n_dimuons] = (recDimuon->get_mom()).X();
    dimuon_py[n_dimuons] = (recDimuon->get_mom()).Y();
    dimuon_pz[n_dimuons] = (recDimuon->get_mom()).Z();
    dimuon_pmom_x[n_dimuons] = (recDimuon->p_pos).Px(); //4-momentum of the muon tracks after vertex fit
    dimuon_pmom_y[n_dimuons] = (recDimuon->p_pos).Py();
    dimuon_pmom_z[n_dimuons] = (recDimuon->p_pos).Pz();
    dimuon_nmom_x[n_dimuons] = (recDimuon->p_neg).Px();
    dimuon_nmom_y[n_dimuons] = (recDimuon->p_neg).Py();
    dimuon_nmom_z[n_dimuons] = (recDimuon->p_neg).Pz();
    dimuon_ppos_x[n_dimuons] = (recDimuon->vtx_pos).X(); // vertex position
    dimuon_ppos_y[n_dimuons] = (recDimuon->vtx_pos).Y();
    dimuon_ppos_z[n_dimuons] = (recDimuon->vtx_pos).Z();
    dimuon_npos_x[n_dimuons] = (recDimuon->vtx_neg).X(); 
    dimuon_npos_y[n_dimuons] = (recDimuon->vtx_neg).Y();
    dimuon_npos_z[n_dimuons] = (recDimuon->vtx_neg).Z();

    ++n_dimuons;
    if(n_dimuons >= 100) break;
  }
  
  n_showers = 0;
  for (auto iter = _truth->GetPrimaryParticleRange().first;
       iter != _truth->GetPrimaryParticleRange().second; ++iter) {
    PHG4Particle *primary = iter->second;
    int ECAL_volume = PHG4HitDefs::get_volume_id("G4HIT_EMCal");
    PHG4Shower *shower = get_primary_shower(primary);
    if (shower != 0) {
      sx_ecal[n_showers] = shower->get_x();
      sy_ecal[n_showers] = shower->get_y();
      sz_ecal[n_showers] = shower->get_z();
      sedep_ecal[n_showers] = shower->get_edep(ECAL_volume);
      n_showers++;
    }
    if (n_showers >= 1000)
      break;
  }

  /*
  _truth->identify();
  for (auto iterp = _truth->GetSecondaryParticleRange().first;
       iterp != _truth->GetSecondaryParticleRange().second; ++iterp) {
    PHG4Particle *secondary = iterp->second;
    std::cout << " secondary particle " << secondary->get_pid() << " e" << secondary->get_e() << std::endl;
  }
  */

  n_primaries = 0;
  for (auto iterp = _truth->GetPrimaryParticleRange().first;
       iterp != _truth->GetPrimaryParticleRange().second; ++iterp) {
    PHG4Particle *primary = iterp->second;
    gpid[n_primaries] = primary->get_pid();

    int vtx_id = primary->get_vtx_id();
    PHG4VtxPoint *vtx = _truth->GetVtx(vtx_id);
    gvx[n_primaries] = vtx->get_x();
    gvy[n_primaries] = vtx->get_y();
    gvz[n_primaries] = vtx->get_z();

    TVector3 mom(primary->get_px(), primary->get_py(), primary->get_pz());
    gpx[n_primaries] = primary->get_px();
    gpy[n_primaries] = primary->get_py();
    gpz[n_primaries] = primary->get_pz();
    gpt[n_primaries] = mom.Pt();
    geta[n_primaries] = mom.Eta();
    gphi[n_primaries] = mom.Phi();
    ge[n_primaries] = primary->get_e();

    int trkID = primary->get_track_id();
    gtrkid[n_primaries] = trkID;


    // G4Hits at different stations                                                                                                                                                                       
    if(g4hc_ecal){
      std::vector<PHG4Hit*> g4hits = FindG4HitsAtStation(trkID, g4hc_ecal);
      nhits_ecal[n_primaries] =0;
      for(int iecal=0; iecal<g4hits.size(); ++iecal) {
	PHG4Hit* g4hit = g4hits[iecal];
	gx_ecal[n_primaries][iecal] = g4hit->get_x(0);
	gy_ecal[n_primaries][iecal] = g4hit->get_y(0);
	gz_ecal[n_primaries][iecal] = g4hit->get_z(0);
	gpx_ecal[n_primaries][iecal] = g4hit->get_px(0);
	gpy_ecal[n_primaries][iecal] = g4hit->get_py(0);
	gpz_ecal[n_primaries][iecal] = g4hit->get_pz(0);
	gedep_ecal[n_primaries][iecal] = g4hit->get_edep();
	++nhits_ecal[n_primaries];
	if(iecal>=100){
	  std::cout << "More than 100 hits in EMCAL " << std::endl;
	  break;
	}
      }
    }

    PHG4Hit *st1hit = FindG4HitAtStation(trkID, g4hc_d1x);
    if (st1hit) {
      gx_st1[n_primaries] = st1hit->get_x(0);
      gy_st1[n_primaries] = st1hit->get_y(0);
      gz_st1[n_primaries] = st1hit->get_z(0);
      gpx_st1[n_primaries] = st1hit->get_px(0);
      gpy_st1[n_primaries] = st1hit->get_py(0);
      gpz_st1[n_primaries] = st1hit->get_pz(0);
    }
    PHG4Hit *st2hit = FindG4HitAtStation(trkID, g4hc_d2xp);
    if (st2hit) {
      gx_st2[n_primaries] = st2hit->get_x(0);
      gy_st2[n_primaries] = st2hit->get_y(0);
      gz_st2[n_primaries] = st2hit->get_z(0);
      gpx_st2[n_primaries] = st2hit->get_px(0);
      gpy_st2[n_primaries] = st2hit->get_py(0);
      gpz_st2[n_primaries] = st2hit->get_pz(0);
    }
    PHG4Hit *st3hit = FindG4HitAtStation(trkID, g4hc_d3px);
    if (!st3hit)
      PHG4Hit *st3hit = FindG4HitAtStation(trkID, g4hc_d3mx);
    if (st3hit) {
      gx_st3[n_primaries] = st3hit->get_x(0);
      gy_st3[n_primaries] = st3hit->get_y(0);
      gz_st3[n_primaries] = st3hit->get_z(0);
      gpx_st3[n_primaries] = st3hit->get_px(0);
      gpy_st3[n_primaries] = st3hit->get_py(0);
      gpz_st3[n_primaries] = st3hit->get_pz(0);
    }

    PHG4Hit *h1hit = FindG4HitAtStation(trkID, g4hc_h1t);
    if (!h1hit)
      PHG4Hit *h1hit = FindG4HitAtStation(trkID, g4hc_h1b);
    if (h1hit) {
      gx_h1[n_primaries] = h1hit->get_x(0);
      gy_h1[n_primaries] = h1hit->get_y(0);
      gz_h1[n_primaries] = h1hit->get_z(0);
      gpx_h1[n_primaries] = h1hit->get_px(0);
      gpy_h1[n_primaries] = h1hit->get_py(0);
      gpz_h1[n_primaries] = h1hit->get_pz(0);
    }
    PHG4Hit *h2hit = FindG4HitAtStation(trkID, g4hc_h2t);
    if (!h2hit)
      PHG4Hit *h2hit = FindG4HitAtStation(trkID, g4hc_h2b);
    if (h2hit) {
      gx_h2[n_primaries] = h2hit->get_x(0);
      gy_h2[n_primaries] = h2hit->get_y(0);
      gz_h2[n_primaries] = h2hit->get_z(0);
      gpx_h2[n_primaries] = h2hit->get_px(0);
      gpy_h2[n_primaries] = h2hit->get_py(0);
      gpz_h2[n_primaries] = h2hit->get_pz(0);
    }
    PHG4Hit *h3hit = FindG4HitAtStation(trkID, g4hc_h3t);
    if (!h3hit)
      PHG4Hit *h3hit = FindG4HitAtStation(trkID, g4hc_h3b);
    if (h3hit) {
      gx_h3[n_primaries] = h3hit->get_x(0);
      gy_h3[n_primaries] = h3hit->get_y(0);
      gz_h3[n_primaries] = h3hit->get_z(0);
      gpx_h3[n_primaries] = h3hit->get_px(0);
      gpy_h3[n_primaries] = h3hit->get_py(0);
      gpz_h3[n_primaries] = h3hit->get_pz(0);
    }

    PHG4Hit *h4hit = FindG4HitAtStation(trkID, g4hc_h4t);
    if (!h4hit)
      PHG4Hit *h4hit = FindG4HitAtStation(trkID, g4hc_h4b);
    if (h4hit) {
      gx_h4[n_primaries] = h4hit->get_x(0);
      gy_h4[n_primaries] = h4hit->get_y(0);
      gz_h4[n_primaries] = h4hit->get_z(0);
      gpx_h4[n_primaries] = h4hit->get_px(0);
      gpy_h4[n_primaries] = h4hit->get_py(0);
      gpz_h4[n_primaries] = h4hit->get_pz(0);
    }
    ++n_primaries;
  }

  //dptrigger
  //std::cout<<"DS in SimAna: "<<_sqEvent->get_trigger(SQEvent::MATRIX1)<<" "<<_sqEvent->get_trigger(SQEvent::MATRIX2)<<" "<<_sqEvent->get_trigger(SQEvent::MATRIX3)<<" "<<_sqEvent->get_trigger(SQEvent::MATRIX4)<<" "<<_sqEvent->get_trigger(SQEvent::MATRIX5)<<" "<<std::endl;

  fpga_trigger[0] = _sqEvent->get_trigger(SQEvent::MATRIX1);
  fpga_trigger[1] = _sqEvent->get_trigger(SQEvent::MATRIX2);
  fpga_trigger[2] = _sqEvent->get_trigger(SQEvent::MATRIX3);
  fpga_trigger[3] = _sqEvent->get_trigger(SQEvent::MATRIX4);
  fpga_trigger[4] = _sqEvent->get_trigger(SQEvent::MATRIX5);

  weight = _sqMCEvent->get_weight();

  ++eventID;
  saveTree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

int SimAna::End(PHCompositeNode *topNode) {
  saveFile->cd();
  saveTree->Write();
  saveFile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

int SimAna::GetNodes(PHCompositeNode* topNode)
{
  _sqEvent   = findNode::getClass<SQEvent_v1>(topNode, "SQEvent");
  if(!_sqEvent) return Fun4AllReturnCodes::ABORTEVENT;

  _sqMCEvent = findNode::getClass<SQMCEvent     >(topNode, "SQMCEvent");
  if(!_sqMCEvent) return Fun4AllReturnCodes::ABORTEVENT;

  _hitVector    = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if(!_hitVector) return Fun4AllReturnCodes::ABORTEVENT;

  _trackVector  = findNode::getClass<SQTrackVector>(topNode, "SQTruthTrackVector");
  if(!_trackVector) {
    std::cout << "ERROR:: did not find SQTruthTrackVector" << std::endl;
  }

  _dimuonVector = findNode::getClass<SQDimuonVector>(topNode, "SQTruthDimuonVector");
  if(!_dimuonVector) {
    std::cout << "ERROR:: did not find SQTruthDimuonVector" << std::endl;
  }

  if(legacyContainer) {
    _recEvent = findNode::getClass<SRecEvent>(topNode, "SRecEvent");
    if(!_recEvent) {
      _recEvent = nullptr;
      std::cout << "ERROR:: no RecEvent " << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  else { 
    _recTrackVector  = findNode::getClass<SQTrackVector>(topNode, "SQRecTrackVector");
    _recDimuonVector = findNode::getClass<SQDimuonVector>(topNode, "SQRecDimuonVector");
    if(!_recTrackVector) {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  _truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth) return Fun4AllReturnCodes::ABORTEVENT;

  g4hc_d1x  = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D1X");
  g4hc_d2xp = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D2Xp");
  g4hc_d3px = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3pXp");
  g4hc_d3mx = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3mXp");
  if (! g4hc_d1x) 
    g4hc_d1x = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D0X"); // D0X is considered as station 1
  if ( !g4hc_d1x || !g4hc_d3px || !g4hc_d3mx) {
    std::cout << "SimAna::GetNode No drift chamber node " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  g4hc_h1t = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1T");
  g4hc_h1b = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1B");
  g4hc_h1l = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1L");
  g4hc_h1r = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1R");
  g4hc_h2t = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2T");
  g4hc_h2b = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2B");
  g4hc_h2l = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2L");
  g4hc_h2r = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2R");
  g4hc_h3t = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H3T");
  g4hc_h3b = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H3B");
  g4hc_h4t = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4T");
  g4hc_h4b = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4B");
  if (!g4hc_h1t || !g4hc_h1b || !g4hc_h2t || !g4hc_h2b || !g4hc_h3t ||
      !g4hc_h3b || !g4hc_h4t || !g4hc_h4b) {
    std::cout << "SimAna::GetNode No nodoscope node, abort " <<std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  g4hc_p1y1 = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1Y1");
  g4hc_p1y2 = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1Y2");
  g4hc_p1x1 = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1X1");
  g4hc_p1x2 = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1X2");
  g4hc_p2x1 = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2X1");
  g4hc_p2x2 = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2X2");
  g4hc_p2y1 = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2Y1");
  g4hc_p2y2 = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2Y2");
  if (!g4hc_p1y1 || !g4hc_p1y2 || !g4hc_p1x1 || !g4hc_p1x2 || !g4hc_p2x1 ||
      !g4hc_p2x2 || !g4hc_p2y1 || !g4hc_p2y2) {
    std::cout << "SimAna::GetNode No prototube node, abort " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  g4hc_ecal = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_EMCal");
  if (!g4hc_ecal) {
    std::cout << "SimAna::GetNode No EMcal node" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void SimAna::MakeTree()
{
  saveFile= new TFile(saveNameOut, "RECREATE");
  saveTree = new TTree("Events", "Tree Created by SimAna");
  saveTree->Branch("eventID", &eventID, "eventID/I");

  saveTree->Branch("n_hits",        &n_hits,          "n_hits/I");
  saveTree->Branch("hit_detID",     hit_detid,        "hit_detID[n_hits]/I");
  saveTree->Branch("hit_elmID",     hit_elmid,        "hit_elmID[n_hits]/I");
  saveTree->Branch("hit_trkID",     hit_trkid,        "hit_trkID[n_hits]/I");
  saveTree->Branch("hit_driftdis",  hit_driftdis,     "hit_driftdis[n_hits]/F");
  saveTree->Branch("hit_pos",       hit_pos,          "hit_pos[n_hits]/F");
  saveTree->Branch("hit_edep",      hit_edep,         "hit_edep[n_hits]/F");
  saveTree->Branch("hit_truthx",    hit_truthx,       "hit_truthx[n_hits]/F");
  saveTree->Branch("hit_truthy",    hit_truthy,       "hit_truthy[n_hits]/F");
  saveTree->Branch("hit_truthz",    hit_truthz,       "hit_truthz[n_hits]/F");
  saveTree->Branch("hit_truthpx",   hit_truthpx,      "hit_truthpx[n_hits]/F");
  saveTree->Branch("hit_truthpy",   hit_truthpy,      "hit_truthpy[n_hits]/F");
  saveTree->Branch("hit_truthpz",   hit_truthpz,      "hit_truthpz[n_hits]/F");

  saveTree->Branch("n_truthtracks",             &n_truthtracks,           "n_truthtracks/I");
  saveTree->Branch("truthtrack_charge",         truthtrack_charge,        "truthtrack_charge[n_truthtracks]/I");
  saveTree->Branch("truthtrack_x_st1",          truthtrack_x_st1,         "truthtrack_x_st1[n_truthtracks]/F");
  saveTree->Branch("truthtrack_y_st1",          truthtrack_y_st1,         "truthtrack_y_st1[n_truthtracks]/F");
  saveTree->Branch("truthtrack_z_st1",          truthtrack_z_st1,         "truthtrack_z_st1[n_truthtracks]/F");
  saveTree->Branch("truthtrack_px_st1",         truthtrack_px_st1,        "truthtrack_px_st1[n_truthtracks]/F");
  saveTree->Branch("truthtrack_py_st1",         truthtrack_py_st1,        "truthtrack_py_st1[n_truthtracks]/F");
  saveTree->Branch("truthtrack_pz_st1",         truthtrack_pz_st1,        "truthtrack_pz_st1[n_truthtracks]/F");
  saveTree->Branch("truthtrack_x_st3",          truthtrack_x_st3,         "truthtrack_x_st3[n_truthtracks]/F");
  saveTree->Branch("truthtrack_y_st3",          truthtrack_y_st3,         "truthtrack_y_st3[n_truthtracks]/F");
  saveTree->Branch("truthtrack_z_st3",          truthtrack_z_st3,         "truthtrack_z_st3[n_truthtracks]/F");
  saveTree->Branch("truthtrack_px_st3",         truthtrack_px_st3,        "truthtrack_px_st3[n_truthtracks]/F");
  saveTree->Branch("truthtrack_py_st3",         truthtrack_py_st3,        "truthtrack_py_st3[n_truthtracks]/F");
  saveTree->Branch("truthtrack_pz_st3",         truthtrack_pz_st3,        "truthtrack_pz_st3[n_truthtracks]/F");
  saveTree->Branch("truthtrack_x_vtx",          truthtrack_x_vtx,         "truthtrack_x_vtx[n_truthtracks]/F");
  saveTree->Branch("truthtrack_y_vtx",          truthtrack_y_vtx,         "truthtrack_y_vtx[n_truthtracks]/F");
  saveTree->Branch("truthtrack_z_vtx",          truthtrack_z_vtx,         "truthtrack_z_vtx[n_truthtracks]/F");
  saveTree->Branch("truthtrack_px_vtx",         truthtrack_px_vtx,        "truthtrack_px_vtx[n_truthtracks]/F");
  saveTree->Branch("truthtrack_py_vtx",         truthtrack_py_vtx,        "truthtrack_py_vtx[n_truthtracks]/F");
  saveTree->Branch("truthtrack_pz_vtx",         truthtrack_pz_vtx,        "truthtrack_pz_vtx[n_truthtracks]/F");

  saveTree->Branch("n_tracks",              &n_tracks,            "n_tracks/I");
  saveTree->Branch("track_charge",         track_charge,        "track_charge[n_tracks]/I");
  saveTree->Branch("track_nhits",          track_nhits,         "track_nhits[n_tracks]/I");
  saveTree->Branch("track_x_target",       track_x_target,      "track_x_target[n_tracks]/F");
  saveTree->Branch("track_y_target",       track_y_target,      "track_y_target[n_tracks]/F");
  saveTree->Branch("track_z_target",       track_z_target,      "track_z_target[n_tracks]/F");
  saveTree->Branch("track_px_target",      track_px_target,     "track_px_target[n_tracks]/F");
  saveTree->Branch("track_py_target",      track_py_target,     "track_py_target[n_tracks]/F");
  saveTree->Branch("track_pz_target",      track_pz_target,     "track_pz_target[n_tracks]/F");
  saveTree->Branch("track_x_st1",          track_x_st1,         "track_x_st1[n_tracks]/F");
  saveTree->Branch("track_y_st1",          track_y_st1,         "track_y_st1[n_tracks]/F");
  saveTree->Branch("track_z_st1",          track_z_st1,         "track_z_st1[n_tracks]/F");
  saveTree->Branch("track_px_st1",         track_px_st1,        "track_px_st1[n_tracks]/F");
  saveTree->Branch("track_py_st1",         track_py_st1,        "track_py_st1[n_tracks]/F");
  saveTree->Branch("track_pz_st1",         track_pz_st1,        "track_pz_st1[n_tracks]/F");
  saveTree->Branch("track_x_vtx",          track_x_vtx,         "track_x_vtx[n_tracks]/F");
  saveTree->Branch("track_y_vtx",          track_y_vtx,         "track_y_vtx[n_tracks]/F");
  saveTree->Branch("track_z_vtx",          track_z_vtx,         "track_z_vtx[n_tracks]/F");
  saveTree->Branch("track_px_vtx",         track_px_vtx,        "track_px_vtx[n_tracks]/F");
  saveTree->Branch("track_py_vtx",         track_py_vtx,        "track_py_vtx[n_tracks]/F");
  saveTree->Branch("track_pz_vtx",         track_pz_vtx,        "track_pz_vtx[n_tracks]/F");
  saveTree->Branch("track_m",              track_m,             "track_m[n_tracks]/F");
  saveTree->Branch("track_chisq",          track_chisq,         "track_chisq[n_tracks]/F");
  saveTree->Branch("track_prob",           track_prob,          "track_prob[n_tracks]/F");
  saveTree->Branch("track_quality",        track_quality,       "track_quality[n_tracks]/F");
  saveTree->Branch("track_nhits_st1",      track_nhits_st1,     "track_nhits_st1[n_tracks]/I");
  saveTree->Branch("track_nhits_st2",      track_nhits_st2,     "track_nhits_st2[n_tracks]/I");
  saveTree->Branch("track_nhits_st3",      track_nhits_st3,     "track_nhits_st3[n_tracks]/I");

  saveTree->Branch("n_truthdimuons",     &n_truthdimuons,   "n_truthdimuons/I");
  saveTree->Branch("truthdimuon_mass",   truthdimuon_mass,  "truthdimuon_mass[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_x_vtx",  truthdimuon_x_vtx, "truthdimuon_x_vtx[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_y_vtx",  truthdimuon_y_vtx, "truthdimuon_y_vtx[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_z_vtx",  truthdimuon_z_vtx, "truthdimuon_z_vtx[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_px",     truthdimuon_px,    "truthdimuon_px[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_py",     truthdimuon_py,    "truthdimuon_py[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_pz",     truthdimuon_pz,    "truthdimuon_pz[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_pmom_x", truthdimuon_pmom_x, "truthdimuon_pmom_x[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_pmom_y", truthdimuon_pmom_y, "truthdimuon_pmom_y[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_pmom_z", truthdimuon_pmom_z, "truthdimuon_pmom_z[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_nmom_x", truthdimuon_nmom_x, "truthdimuon_nmom_x[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_nmom_y", truthdimuon_nmom_y, "truthdimuon_nmom_y[n_truthdimuons]/F");
  saveTree->Branch("truthdimuon_nmom_z", truthdimuon_nmom_z, "truthdimuon_nmom_z[n_truthdimuons]/F");

  saveTree->Branch("n_dimuons",     &n_dimuons,    "n_dimuons/I");
  saveTree->Branch("dimuon_mass",   dimuon_mass,   "dimuon_mass[n_dimuons]/F");
  saveTree->Branch("dimuon_chisq",  dimuon_chisq,  "dimuon_chisq[n_dimuons]/F");
  saveTree->Branch("dimuon_x_vtx",  dimuon_x_vtx,  "dimuon_x_vtx[n_dimuons]/F");
  saveTree->Branch("dimuon_y_vtx",  dimuon_y_vtx,  "dimuon_y_vtx[n_dimuons]/F");
  saveTree->Branch("dimuon_z_vtx",  dimuon_z_vtx,  "dimuon_z_vtx[n_dimuons]/F");
  saveTree->Branch("dimuon_px",     dimuon_px,     "dimuon_px[n_dimuons]/F");
  saveTree->Branch("dimuon_py",     dimuon_py,     "dimuon_py[n_dimuons]/F");
  saveTree->Branch("dimuon_pz",     dimuon_pz,     "dimuon_pz[n_dimuons]/F");
  saveTree->Branch("dimuon_pmom_x", dimuon_pmom_x, "dimuon_pmom_x[n_dimuons]/F");
  saveTree->Branch("dimuon_pmom_y", dimuon_pmom_y, "dimuon_pmom_y[n_dimuons]/F");
  saveTree->Branch("dimuon_pmom_z", dimuon_pmom_z, "dimuon_pmom_z[n_dimuons]/F");
  saveTree->Branch("dimuon_nmom_x", dimuon_nmom_x, "dimuon_nmom_x[n_dimuons]/F");
  saveTree->Branch("dimuon_nmom_y", dimuon_nmom_y, "dimuon_nmom_y[n_dimuons]/F");
  saveTree->Branch("dimuon_nmom_z", dimuon_nmom_z, "dimuon_nmom_z[n_dimuons]/F");
  saveTree->Branch("dimuon_ppos_x", dimuon_ppos_x, "dimuon_ppos_x[n_dimuons]/F");
  saveTree->Branch("dimuon_ppos_y", dimuon_ppos_y, "dimuon_ppos_y[n_dimuons]/F");
  saveTree->Branch("dimuon_ppos_z", dimuon_ppos_z, "dimuon_ppos_z[n_dimuons]/F");
  saveTree->Branch("dimuon_npos_x", dimuon_npos_x, "dimuon_npos_x[n_dimuons]/F");
  saveTree->Branch("dimuon_npos_y", dimuon_npos_y, "dimuon_npos_y[n_dimuons]/F");
  saveTree->Branch("dimuon_npos_z", dimuon_npos_z, "dimuon_npos_z[n_dimuons]/F");

  saveTree->Branch("n_showers",     &n_showers,       "n_showers/I");
  saveTree->Branch("sx_ecal",       &sx_ecal,         "sx_ecal[n_showers]/F");
  saveTree->Branch("sy_ecal",       &sy_ecal,         "sy_ecal[n_showers]/F");
  saveTree->Branch("sz_ecal",       &sz_ecal,         "sz_ecal[n_showers]/F");
  saveTree->Branch("sedep_ecal",    &sedep_ecal,      "sedep_ecal[n_showers]/F");

  saveTree->Branch("n_primaries",   &n_primaries,        "n_primaries/I");
  saveTree->Branch("gtrkid",        gtrkid,              "gtrkid[n_primaries]/I");
  saveTree->Branch("gpid",          gpid,                "gpid[n_primaries]/I");
  saveTree->Branch("gvx",           gvx,                 "gvx[n_primaries]/F");
  saveTree->Branch("gvy",           gvy,                 "gvy[n_primaries]/F");
  saveTree->Branch("gvz",           gvz,                 "gvz[n_primaries]/F");
  saveTree->Branch("gpx",           gpx,                 "gpx[n_primaries]/F");
  saveTree->Branch("gpy",           gpy,                 "gpy[n_primaries]/F");
  saveTree->Branch("gpz",           gpz,                 "gpz[n_primaries]/F");
  saveTree->Branch("gpt",           gpt,                 "gpt[n_primaries]/F");
  saveTree->Branch("geta",          geta,                "geta[n_primaries]/F");
  saveTree->Branch("gphi",          gphi,                "gphi[n_primaries]/F");
  saveTree->Branch("ge",            ge,                  "ge[n_primaries]/F");

  saveTree->Branch("nhits_ecal",    nhits_ecal,          "nhits_ecal[n_primaries]/I");
  saveTree->Branch("gx_ecal",       gx_ecal,             "gx_ecal[n_primaries][100]/F"); // not sure how to make this with the right size                                                                
  saveTree->Branch("gy_ecal",       gy_ecal,             "gy_ecal[n_primaries][100]/F");
  saveTree->Branch("gz_ecal",       gz_ecal,             "gz_ecal[n_primaries][100]/F");
  saveTree->Branch("gpx_ecal",      gpx_ecal,            "gpx_ecal[n_primaries][100]/F");
  saveTree->Branch("gpy_ecal",      gpy_ecal,            "gpy_ecal[n_primaries][100]/F");
  saveTree->Branch("gpz_ecal",      gpz_ecal,            "gpz_ecal[n_primaries][100]/F");
  saveTree->Branch("gedep_ecal",    gedep_ecal,          "gedep_ecal[n_primaries][100]/F");

  saveTree->Branch("gx_st1",        gx_st1,              "gx_st1[n_primaries]/F");
  saveTree->Branch("gy_st1",        gy_st1,              "gy_st1[n_primaries]/F");
  saveTree->Branch("gz_st1",        gz_st1,              "gz_st1[n_primaries]/F");
  saveTree->Branch("gpx_st1",       gpx_st1,             "gpx_st1[n_primaries]/F");
  saveTree->Branch("gpy_st1",       gpy_st1,             "gpy_st1[n_primaries]/F");
  saveTree->Branch("gpz_st1",       gpz_st1,             "gpz_st1[n_primaries]/F");

  saveTree->Branch("gx_st2",        gx_st2,              "gx_st2[n_primaries]/F");
  saveTree->Branch("gy_st2",        gy_st2,              "gy_st2[n_primaries]/F");
  saveTree->Branch("gz_st2",        gz_st2,              "gz_st2[n_primaries]/F");
  saveTree->Branch("gpx_st2",       gpx_st2,             "gpx_st2[n_primaries]/F");
  saveTree->Branch("gpy_st2",       gpy_st2,             "gpy_st2[n_primaries]/F");
  saveTree->Branch("gpz_st2",       gpz_st2,             "gpz_st2[n_primaries]/F");

  saveTree->Branch("gx_st3",        gx_st3,              "gx_st3[n_primaries]/F");
  saveTree->Branch("gy_st3",        gy_st3,              "gy_st3[n_primaries]/F");
  saveTree->Branch("gz_st3",        gz_st3,              "gz_st3[n_primaries]/F");
  saveTree->Branch("gpx_st3",       gpx_st3,             "gpx_st3[n_primaries]/F");
  saveTree->Branch("gpy_st3",       gpy_st3,             "gpy_st3[n_primaries]/F");
  saveTree->Branch("gpz_st3",       gpz_st3,             "gpz_st3[n_primaries]/F");

  saveTree->Branch("gx_h1",         gx_h1,               "gx_h1[n_primaries]/F");
  saveTree->Branch("gy_h1",         gy_h1,               "gy_h1[n_primaries]/F");
  saveTree->Branch("gz_h1",         gz_h1,               "gz_h1[n_primaries]/F");
  saveTree->Branch("gpx_h1",        gpx_h1,              "gpx_h1[n_primaries]/F");
  saveTree->Branch("gpy_h1",        gpy_h1,              "gpy_h1[n_primaries]/F");
  saveTree->Branch("gpz_h1",        gpz_h1,              "gpz_h1[n_primaries]/F");

  saveTree->Branch("gx_h2",         gx_h2,               "gx_h2[n_primaries]/F");
  saveTree->Branch("gy_h2",         gy_h2,               "gy_h2[n_primaries]/F");
  saveTree->Branch("gz_h2",         gz_h2,               "gz_h2[n_primaries]/F");
  saveTree->Branch("gpx_h2",        gpx_h2,              "gpx_h2[n_primaries]/F");
  saveTree->Branch("gpy_h2",        gpy_h2,              "gpy_h2[n_primaries]/F");
  saveTree->Branch("gpz_h2",        gpz_h2,              "gpz_h2[n_primaries]/F");

  saveTree->Branch("gx_h3",         gx_h3,               "gx_h3[n_primaries]/F");
  saveTree->Branch("gy_h3",         gy_h3,               "gy_h3[n_primaries]/F");
  saveTree->Branch("gz_h3",         gz_h3,               "gz_h3[n_primaries]/F");
  saveTree->Branch("gpx_h3",        gpx_h3,              "gpx_h3[n_primaries]/F");
  saveTree->Branch("gpy_h3",        gpy_h3,              "gpy_h3[n_primaries]/F");
  saveTree->Branch("gpz_h3",        gpz_h3,              "gpz_h3[n_primaries]/F");

  saveTree->Branch("gx_h4",         gx_h4,               "gx_h4[n_primaries]/F");
  saveTree->Branch("gy_h4",         gy_h4,               "gy_h4[n_primaries]/F");
  saveTree->Branch("gz_h4",         gz_h4,               "gz_h4[n_primaries]/F");
  saveTree->Branch("gpx_h4",        gpx_h4,              "gpx_h4[n_primaries]/F");
  saveTree->Branch("gpy_h4",        gpy_h4,              "gpy_h4[n_primaries]/F");
  saveTree->Branch("gpz_h4",        gpz_h4,              "gpz_h4[n_primaries]/F");

  saveTree->Branch("gx_p1",         gx_p1,               "gx_p1[n_primaries]/F");
  saveTree->Branch("gy_p1",         gy_p1,               "gy_p1[n_primaries]/F");
  saveTree->Branch("gz_p1",         gz_p1,               "gz_p1[n_primaries]/F");
  saveTree->Branch("gpx_p1",        gpx_p1,              "gpx_p1[n_primaries]/F");
  saveTree->Branch("gpy_p1",        gpy_p1,              "gpy_p1[n_primaries]/F");
  saveTree->Branch("gpz_p1",        gpz_p1,              "gpz_p1[n_primaries]/F");

  saveTree->Branch("fpga_trigger",  fpga_trigger,        "fpga_trigger[5]/O");

  saveTree->Branch("weight",        &weight,              "weight/F");

}
