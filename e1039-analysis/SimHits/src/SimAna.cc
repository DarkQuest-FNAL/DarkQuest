#include <iomanip>
#include <TFile.h>
#include <TTree.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/getClass.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQHit.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4VtxPoint.h>

#include "SimAna.h"

SimAna::SimAna(const std::string& name): SubsysReco(name)
{}

SimAna::~SimAna()
{}

int SimAna::Init(PHCompositeNode* topNode)
{
  n_hits = 0;
  for(int i=0; i<10000; ++i) {
    hit_detid[i]        = std::numeric_limits<short>::max();
    hit_elmid[i]        = std::numeric_limits<short>::max();
    hit_driftdis[i]     = std::numeric_limits<float>::max();
    hit_pos[i]          = std::numeric_limits<float>::max();
    hit_edep[i]         = std::numeric_limits<float>::max();
    hit_truthx[i]       = std::numeric_limits<float>::max();
    hit_truthy[i]       = std::numeric_limits<float>::max();
    hit_truthz[i]       = std::numeric_limits<float>::max();
  }

  n_showers =0;
  for(int i=0; i<1000; ++i) {
    sx_ecal[i]          = std::numeric_limits<float>::max();
    sy_ecal[i]          = std::numeric_limits<float>::max();
    sz_ecal[i]          = std::numeric_limits<float>::max();
    sedep_ecal[i]       = std::numeric_limits<float>::max();
  }

  n_primaries = 0;
  for(int i=0; i<1000; ++i) {
    gtrkid[i]     = std::numeric_limits<int>::max();
    gpid[i]       = std::numeric_limits<int>::max();
    gvx[i]        = std::numeric_limits<float>::max();
    gvy[i]        = std::numeric_limits<float>::max();
    gvz[i]        = std::numeric_limits<float>::max();
    gpx[i]        = std::numeric_limits<float>::max();
    gpy[i]        = std::numeric_limits<float>::max();
    gpz[i]        = std::numeric_limits<float>::max();
    gpt[i]        = std::numeric_limits<float>::max();
    geta[i]       = std::numeric_limits<float>::max();
    gphi[i]       = std::numeric_limits<float>::max();
    ge[i]         = std::numeric_limits<float>::max();

    nhits_ecal[i] = 0;
    for(int j=0; j<100; ++j) {
      gx_ecal[i][j]    = std::numeric_limits<int>::max();
      gy_ecal[i][j]    = std::numeric_limits<int>::max();
      gz_ecal[i][j]    = std::numeric_limits<int>::max();
      gpx_ecal[i][j]   = std::numeric_limits<int>::max();
      gpy_ecal[i][j]   = std::numeric_limits<int>::max();
      gpz_ecal[i][j]   = std::numeric_limits<int>::max();
      gedep_ecal[i][j] = std::numeric_limits<int>::max();
    }

    gx_st1[i]     = std::numeric_limits<float>::max();
    gy_st1[i]     = std::numeric_limits<float>::max();
    gz_st1[i]     = std::numeric_limits<float>::max();
    gpx_st1[i]    = std::numeric_limits<float>::max();
    gpy_st1[i]    = std::numeric_limits<float>::max();
    gpz_st1[i]    = std::numeric_limits<float>::max();
    gx_st2[i]     = std::numeric_limits<float>::max();
    gy_st2[i]     = std::numeric_limits<float>::max();
    gz_st2[i]     = std::numeric_limits<float>::max();
    gpx_st2[i]    = std::numeric_limits<float>::max();
    gpy_st2[i]    = std::numeric_limits<float>::max();
    gpz_st2[i]    = std::numeric_limits<float>::max();
    gx_st3[i]     = std::numeric_limits<float>::max();
    gy_st3[i]     = std::numeric_limits<float>::max();
    gz_st3[i]     = std::numeric_limits<float>::max();
    gpx_st3[i]    = std::numeric_limits<float>::max();
    gpy_st3[i]    = std::numeric_limits<float>::max();
    gpz_st3[i]    = std::numeric_limits<float>::max();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SimAna::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  eventID = 0;
  MakeTree();
  return Fun4AllReturnCodes::EVENT_OK;
}

PHG4Hit* SimAna::FindG4HitAtStation(const int trk_id, const PHG4HitContainer* g4hc) {
  PHG4Hit* hit = nullptr;
  PHG4HitContainer::ConstRange range = g4hc->getHits();
  for (PHG4HitContainer::ConstIterator it = range.first; it != range.second; it++) {
    PHG4Hit* tmphit = it->second;
    if (tmphit->get_trkid() == trk_id) {
      hit = tmphit;
      break;
    }
  }
  return hit;
}


std::vector<PHG4Hit*> SimAna::FindG4HitsAtStation(const int trk_id, const PHG4HitContainer* g4hc) {
  std::vector<PHG4Hit*> vhit;
  PHG4HitContainer::ConstRange range = g4hc->getHits();
  for (PHG4HitContainer::ConstIterator it = range.first; it != range.second; it++) {
    PHG4Hit* tmphit = it->second;
    if (tmphit->get_trkid() == trk_id) {
      vhit.push_back(tmphit);
    }
  }
  return vhit;
}

// find g4shower                                                                                                                                                                                         
PHG4Shower* SimAna::get_primary_shower(PHG4Particle* primary) {
  PHG4Shower* shower = nullptr;

  PHG4TruthInfoContainer::ShowerRange range = _truth->GetSecondaryShowerRange();
  for (PHG4TruthInfoContainer::ShowerIterator iter = range.first;
       iter != range.second;
       ++iter)
    {
      PHG4Shower* tmpshower = iter->second;
      std::cout << "secshower z " << tmpshower->get_z() << " parent id " << tmpshower->get_parent_particle_id() << " shower parent " << tmpshower->get_parent_shower_id()<<  std::endl;
    }

  for (auto iter=_truth->GetPrimaryShowerRange().first; iter!=_truth->GetPrimaryShowerRange().second; ++iter) {
    PHG4Shower* tmpshower = iter->second;
    if (tmpshower->get_parent_particle_id() == primary->get_track_id()) {
      shower = tmpshower;
      int ECAL_volume = PHG4HitDefs::get_volume_id("G4HIT_EMCal");

      int h1t_volume = PHG4HitDefs::get_volume_id("G4HIT_H1T");
      int h1b_volume = PHG4HitDefs::get_volume_id("G4HIT_H1B");
      int h1l_volume = PHG4HitDefs::get_volume_id("G4HIT_H1L");
      int h1r_volume = PHG4HitDefs::get_volume_id("G4HIT_H1R");
      int h2t_volume = PHG4HitDefs::get_volume_id("G4HIT_H2T");
      int h2b_volume = PHG4HitDefs::get_volume_id("G4HIT_H2B");
      int h2l_volume = PHG4HitDefs::get_volume_id("G4HIT_H2L");
      int h2r_volume = PHG4HitDefs::get_volume_id("G4HIT_H2R");
      int h3t_volume = PHG4HitDefs::get_volume_id("G4HIT_H3T");
      int h3b_volume = PHG4HitDefs::get_volume_id("G4HIT_H3B");
      int h3l_volume = PHG4HitDefs::get_volume_id("G4HIT_H3L");
      int h3r_volume = PHG4HitDefs::get_volume_id("G4HIT_H3R");
      int h4t_volume = PHG4HitDefs::get_volume_id("G4HIT_H4T");
      int h4b_volume = PHG4HitDefs::get_volume_id("G4HIT_H4B");

      int d1x_volume = PHG4HitDefs::get_volume_id("G4HIT_D1X");
      int d0x_volume = PHG4HitDefs::get_volume_id("G4HIT_D0X");
      int d2xp_volume = PHG4HitDefs::get_volume_id("G4HIT_D2Xp");
      int d3pxp_volume = PHG4HitDefs::get_volume_id("G4HIT_D3pXp");
      int d3mxp_volume = PHG4HitDefs::get_volume_id("G4HIT_D3mXp");

      int abs_volume = PHG4HitDefs::get_volume_id("MUID_absorber");

      std::cout << "primary shower z " << shower->get_z() << " shower id " << tmpshower->get_id() << " pz " << primary->get_pz() << " e " << primary->get_e() << std::endl;
      std::cout << "edep emcal " <<  shower->get_edep(ECAL_volume) << std::endl;
      std::cout << " h1t " << shower->get_edep(h1t_volume) << " h1b " << shower->get_edep(h1b_volume) << " h1r " <<  shower->get_edep(h1r_volume) << " h1l " << shower->get_edep(h1l_volume);
      std::cout << " h2t " << shower->get_edep(h2t_volume) << " h2b " << shower->get_edep(h2b_volume) << " h2r " <<  shower->get_edep(h2r_volume) << " h2l " << shower->get_edep(h2l_volume);
      std::cout << " h3t " << shower->get_edep(h3t_volume) << " h3b " << shower->get_edep(h3b_volume) << std::endl;
      std::cout << " d0x " << shower->get_edep(d0x_volume) << " d1x " << shower->get_edep(d1x_volume);
      std::cout << " d2xp " << shower->get_edep(d2xp_volume) << " d3pxp " << shower->get_edep(d3pxp_volume) << " d3mxp " << shower->get_edep(d3mxp_volume) << std::endl;
      std::cout << " absorber " << shower->get_edep(abs_volume) << std::endl;

      double fmag_maxkick = 2.9; //GeV/c
      double fmag_minz = 0.0;
      double fmag_maxz = 500.0;
      double kmag_maxkick = 0.414; //GeV/c
      //center of KMag is at z=1041.8
      double kmag_minz = 891.8;
      double kmag_maxz = 1191.8;

      double fmag_kick, kmag_kick;
      double fmag_center, kmag_center;
      if (vx[2]>fmag_maxz) {
	fmag_kick = 0.0;
	fmag_center = 0.0;
      } else if (vx[2]>fmag_minz) {
	fmag_kick = fmag_maxkick*(fmag_maxz-vx[2])/(fmag_maxz-fmag_minz);
	fmag_center = (fmag_maxz+vx[2])/2.0;
      } else {
	fmag_kick = fmag_maxkick;
	fmag_center = (fmag_maxz+fmag_minz)/2.0;
      }

      int vtx_id =  primary->get_vtx_id();
      PHG4VtxPoint* vtx = _truth->GetVtx(vtx_id);
      double gvx = vtx->get_x();
      double gvy = vtx->get_y();
      double gvz = vtx->get_z();

      PHG4Hit* st1hit = FindG4HitAtStation(trkID, g4hc_d1x);
      if(st1hit){
	std::cout << "st1 x " << st1hit->get_x(0) << " y " << st1hit->get_y(0) <<" z " << st1hit->get_z(0) << std::endl;                                                                                 
	std::cout << "    px " << st1hit->get_px(0) << " py " << st1hit->get_py(0) << " pz " << st1hit->get_pz(0) << std::endl;

	double tx1_st1 = (primary->get_px() + fmag_kick)/primary->get_pz();
	double x1_st1 = gvx + (fmag_center-gvz)*(primary->get_px()/primary->get_pz()) - fmag_center*tx1_st1;

	std::cout << "cal px " <<   
      }
      
      PHG4Hit* h2hit = FindG4HitAtStation(primary->get_track_id(),g4hc_h2t);
      if(!h2hit)
	PHG4Hit* h2hit = FindG4HitAtStation(primary->get_track_id(),g4hc_h2b);
      if(h2hit){
	std::cout << "h2y t/b x " << h2hit->get_x(0) << " y " << h2hit->get_y(0) <<" z " << h2hit->get_z(0) << " pz " << h2hit->get_pz(0) << std::endl;
      }

      PHG4Hit* h3hit = FindG4HitAtStation(primary->get_track_id(),g4hc_h3t);
      if(!h3hit)
	PHG4Hit* h3hit = FindG4HitAtStation(primary->get_track_id(),g4hc_h3b);
      if(h3hit){
	std::cout << "h3y t/b x " << h3hit->get_x(0) << " y " << h3hit->get_y(0) <<" z " << h3hit->get_z(0) << std::endl;
      }

      PHG4Hit* h4hit = FindG4HitAtStation(primary->get_track_id(),g4hc_h4t);
      if(!h4hit)
	PHG4Hit* h4hit = FindG4HitAtStation(primary->get_track_id(),g4hc_h4b);
      if(h4hit){
	std::cout << "h4y t/b x " << h4hit->get_x(0) << " y " << h4hit->get_y(0) <<" z " << h4hit->get_z(0) << std::endl; 
      }
      break;
    }
  }
  return shower;
}


int SimAna::process_event(PHCompositeNode* topNode)
{
  //std::cout << "new event " <<std::endl;
  n_hits = 0;
  for(int ihit=0; ihit<hitVector->size(); ++ihit) {
    SQHit *hit = hitVector->at(ihit);
    int hitID = hit->get_hit_id();
    hit_detid[n_hits]      = hit->get_detector_id();
    hit_elmid[n_hits]      = hit->get_element_id();
    hit_driftdis[n_hits]   = hit->get_drift_distance();
    hit_pos[n_hits]        = hit->get_pos();
    hit_edep[n_hits]       = hit->get_edep();
    hit_truthx[n_hits] = hit->get_truth_x();
    hit_truthy[n_hits] = hit->get_truth_y();
    hit_truthz[n_hits] = hit->get_truth_z();
    ++n_hits;
    if(n_hits>=10000) break;
    if(hit->get_detector_id() == 100){
      std::cout << "emcal hit edep " << hit->get_edep() << std::endl;
    }
  }

  n_showers = 0;
  for(auto iter=_truth->GetPrimaryParticleRange().first; iter!=_truth->GetPrimaryParticleRange().second; ++iter) {
    PHG4Particle * primary = iter->second;
    int ECAL_volume = PHG4HitDefs::get_volume_id("G4HIT_EMCal");
    PHG4Shower* shower = get_primary_shower(primary);
    if(shower !=0){
      sx_ecal[n_showers] = shower->get_x();
      sy_ecal[n_showers] = shower->get_y();
      sz_ecal[n_showers] = shower->get_z();
      sedep_ecal[n_showers] = shower->get_edep(ECAL_volume);
      //std::cout << "shower " << shower->get_x() << " y " << shower->get_y() <<" z " << shower->get_z() << " edep " << shower->get_edep(ECAL_volume) << std::endl;
      n_showers++;
    }
    if(n_showers>=1000) break;
  }
  //std::cout << "nshowers " << n_showers << std::endl;

  n_primaries = 0;
  for(auto iterp=_truth->GetPrimaryParticleRange().first; iterp!=_truth->GetPrimaryParticleRange().second; ++iterp) {
    PHG4Particle * primary = iterp->second;
    gpid[n_primaries] = primary->get_pid();
    //std::cout<<"primary " << n_primaries << " id " << primary->get_pid() << std::endl;

    int vtx_id =  primary->get_vtx_id();
    PHG4VtxPoint* vtx = _truth->GetVtx(vtx_id);
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
    ge[n_primaries] =  primary->get_e();

    int trkID = primary->get_track_id();
    gtrkid[n_primaries] = trkID;

    // G4Hits at different stations                                                                                                                                                                       
    if(g4hc_ecal){
      std::vector<PHG4Hit*> g4hits = FindG4HitsAtStation(trkID, g4hc_ecal);
      nhits_ecal[n_primaries] =0;
      //std::cout << "nhits in EMCAL " << g4hits.size() << std::endl;
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

    PHG4Hit* st1hit = FindG4HitAtStation(trkID, g4hc_d1x);
    if(st1hit){
      gx_st1[n_primaries] = st1hit->get_x(0);
      gy_st1[n_primaries] = st1hit->get_y(0);
      gz_st1[n_primaries] = st1hit->get_z(0);
      //std::cout << "st1 x " << st1hit->get_x(0) << " y " << st1hit->get_y(0) <<" z " << st1hit->get_z(0) << std::endl;
      gpx_st1[n_primaries] = st1hit->get_px(0);
      gpy_st1[n_primaries] = st1hit->get_py(0);
      gpz_st1[n_primaries] = st1hit->get_pz(0);
    }
    PHG4Hit* st2hit = FindG4HitAtStation(trkID, g4hc_d2xp);
    if(st2hit){
      gx_st2[n_primaries] = st2hit->get_x(0);
      gy_st2[n_primaries] = st2hit->get_y(0);
      gz_st2[n_primaries] = st2hit->get_z(0);
      //std::cout << "st2 x " << st2hit->get_x(0) << " y " << st2hit->get_y(0) <<" z " << st2hit->get_z(0) << std::endl;
      gpx_st2[n_primaries] = st2hit->get_px(0);
      gpy_st2[n_primaries] = st2hit->get_py(0);
      gpz_st2[n_primaries] = st2hit->get_pz(0);
    }
    PHG4Hit* st3hit = FindG4HitAtStation(trkID, g4hc_d3px);
    if(!st3hit)
      PHG4Hit* st3hit = FindG4HitAtStation(trkID, g4hc_d3mx);
    if(st3hit){
      gx_st3[n_primaries] = st3hit->get_x(0);
      gy_st3[n_primaries] = st3hit->get_y(0);
      gz_st3[n_primaries] = st3hit->get_z(0);
      //std::cout << "st3 x " << st3hit->get_x(0) << " y " << st3hit->get_y(0) <<" z " << st3hit->get_z(0) << std::endl;
      gpx_st3[n_primaries] = st3hit->get_px(0);
      gpy_st3[n_primaries] = st3hit->get_py(0);
      gpz_st3[n_primaries] = st3hit->get_pz(0);
    }

    PHG4Hit* h1hit = FindG4HitAtStation(trkID,g4hc_h1t);
    if(!h1hit)
      PHG4Hit* h1hit = FindG4HitAtStation(trkID,g4hc_h1b);
    if(h1hit){
      gx_h1[n_primaries] = h1hit->get_x(0);
      gy_h1[n_primaries] = h1hit->get_y(0);
      gz_h1[n_primaries] = h1hit->get_z(0);
      //std::cout << "h1 x " << h1hit->get_x(0) << " y " << h1hit->get_y(0) <<" z " << h1hit->get_z(0) << std::endl;
      gpx_h1[n_primaries] = h1hit->get_px(0);
      gpy_h1[n_primaries] = h1hit->get_py(0);
      gpz_h1[n_primaries] = h1hit->get_pz(0);
    }
    PHG4Hit* h2hit = FindG4HitAtStation(trkID,g4hc_h2t);
    if(!h2hit)
      PHG4Hit* h2hit = FindG4HitAtStation(trkID,g4hc_h2b);
    if(h2hit){
      gx_h2[n_primaries] = h2hit->get_x(0);
      gy_h2[n_primaries] = h2hit->get_y(0);
      gz_h2[n_primaries] = h2hit->get_z(0);
      //std::cout << "h2 x " << h2hit->get_x(0) << " y " << h2hit->get_y(0) <<" z " << h2hit->get_z(0) << std::endl;
      gpx_h2[n_primaries] = h2hit->get_px(0);
      gpy_h2[n_primaries] = h2hit->get_py(0);
      gpz_h2[n_primaries] = h2hit->get_pz(0);
    }
    PHG4Hit* h3hit = FindG4HitAtStation(trkID,g4hc_h3t);
    if(!h3hit)
      PHG4Hit* h3hit = FindG4HitAtStation(trkID,g4hc_h3b);
    if(h3hit){
      gx_h3[n_primaries] = h3hit->get_x(0);
      gy_h3[n_primaries] = h3hit->get_y(0);
      gz_h3[n_primaries] = h3hit->get_z(0);
      //std::cout << "h3 x " << h3hit->get_x(0) << " y " << h3hit->get_y(0) <<" z " << h3hit->get_z(0) << std::endl;
      gpx_h3[n_primaries] = h3hit->get_px(0);
      gpy_h3[n_primaries] = h3hit->get_py(0);
      gpz_h3[n_primaries] = h3hit->get_pz(0);
    }
    PHG4Hit* h4hit = FindG4HitAtStation(trkID,g4hc_h4t);
    if(!h4hit)
      PHG4Hit* h4hit = FindG4HitAtStation(trkID,g4hc_h4b);
    if(h4hit){
      gx_h4[n_primaries] = h4hit->get_x(0);
      gy_h4[n_primaries] = h4hit->get_y(0);
      gz_h4[n_primaries] = h4hit->get_z(0);
      //std::cout << "h4y t/b x " << h4hit->get_x(0) << " y " << h4hit->get_y(0) <<" z " << h4hit->get_z(0) << std::endl;
      gpx_h4[n_primaries] = h4hit->get_px(0);
      gpy_h4[n_primaries] = h4hit->get_py(0);
      gpz_h4[n_primaries] = h4hit->get_pz(0);
    }
    ++n_primaries;
  }

  ++eventID;
  saveTree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

int SimAna::End(PHCompositeNode* topNode)
{
  saveFile->cd();
  saveTree->Write();
  saveFile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

int SimAna::GetNodes(PHCompositeNode* topNode)
{
  hitVector    = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if(!hitVector) return Fun4AllReturnCodes::ABORTEVENT;

  _truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth) return Fun4AllReturnCodes::ABORTEVENT;

  g4hc_d1x  = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D1X");
  g4hc_d2xp = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D2Xp");
  g4hc_d3px = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3pXp");
  g4hc_d3mx = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3mXp");
  if (! g4hc_d1x) g4hc_d1x = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D0X");
  if ( !g4hc_d1x || !g4hc_d3px || !g4hc_d3mx) {
    std::cout << "didnt get d node " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  g4hc_h1t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1T");
  g4hc_h1b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1B");
  g4hc_h2t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2T");
  g4hc_h2b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2B");
  g4hc_h3t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H3T");
  g4hc_h3b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H3B");
  g4hc_h4t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4T");
  g4hc_h4b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4B");
  if (!g4hc_h1t || !g4hc_h1b || !g4hc_h2t || !g4hc_h2b ||
      !g4hc_h3t || !g4hc_h3b || !g4hc_h4t || !g4hc_h4b ) {
    std::cout << "didnt get h node  " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  g4hc_p1y1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1Y1");
  g4hc_p1y2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1Y2");
  g4hc_p1x1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1X1");
  g4hc_p1x2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1X2");
  g4hc_p2x1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2X1");
  g4hc_p2x2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2X2");
  g4hc_p2y1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2Y1");
  g4hc_p2y2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2Y2");
  if (!g4hc_p1y1 || !g4hc_p1y2 || !g4hc_p1x1 || !g4hc_p1x2 ||
      !g4hc_p2x1 || !g4hc_p2x2 || !g4hc_p2y1 || !g4hc_p2y2   ) {
    std::cout << "didnt get p node " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  g4hc_ecal  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_EMCal");
  if (!g4hc_ecal) {
    std::cout << "didnt get ecal node " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void SimAna::MakeTree()
{
  saveFile= new TFile("output.root", "RECREATE");
  saveTree = new TTree("Truth", "Truth Tree Created by SimAna");
  saveTree->Branch("eventID", &eventID, "eventID/I");
  saveTree->Branch("n_hits",        &n_hits,          "n_hits/I");
  saveTree->Branch("hit_detID",     hit_detid,        "hit_detID[n_hits]/I");
  saveTree->Branch("hit_elmID",     hit_elmid,        "hit_elmID[n_hits]/I");
  saveTree->Branch("hit_driftdis",  hit_driftdis,     "hit_driftdis[n_hits]/F");
  saveTree->Branch("hit_pos",       hit_pos,          "hit_pos[n_hits]/F");
  saveTree->Branch("hit_edep",      hit_edep,         "hit_edep[n_hits]/F");
  saveTree->Branch("hit_truthx",    hit_truthx,       "hit_truthx[n_hits]/F");
  saveTree->Branch("hit_truthy",    hit_truthy,       "hit_truthy[n_hits]/F");
  saveTree->Branch("hit_truthz",    hit_truthz,       "hit_truthz[n_hits]/F");

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

  saveTree->Branch("gx_p1",         gx_p1,               "gx_p1[n_primaries]/F");
  saveTree->Branch("gy_p1",         gy_p1,               "gy_p1[n_primaries]/F");
  saveTree->Branch("gz_p1",         gz_p1,               "gz_p1[n_primaries]/F");
  saveTree->Branch("gpx_p1",        gpx_p1,              "gpx_p1[n_primaries]/F");
  saveTree->Branch("gpy_p1",        gpy_p1,              "gpy_p1[n_primaries]/F");
  saveTree->Branch("gpz_p1",        gpz_p1,              "gpz_p1[n_primaries]/F");

  saveTree->Branch("gx_h4",         gx_h4,              "gx_h4[n_primaries]/F");
  saveTree->Branch("gy_h4",         gy_h4,              "gy_h4[n_primaries]/F");
  saveTree->Branch("gz_h4",         gz_h4,              "gz_h4[n_primaries]/F");
  saveTree->Branch("gpx_h4",        gpx_h4,             "gpx_h4[n_primaries]/F");
  saveTree->Branch("gpy_h4",        gpy_h4,             "gpy_h4[n_primaries]/F");
  saveTree->Branch("gpz_h4",        gpz_h4,             "gpz_h4[n_primaries]/F");

}
