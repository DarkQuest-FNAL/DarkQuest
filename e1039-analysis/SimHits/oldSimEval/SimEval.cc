#include "SimEval.h"

#include <interface_main/SQHit.h>
#include <interface_main/SQHit_v1.h>
#include <interface_main/SQMCHit_v1.h>
#include <interface_main/SQHitMap_v1.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQEvent_v1.h>
#include <interface_main/SQRun_v1.h>

#include <geom_svc/GeomSvc.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4VtxPoint.h>

#include <TFile.h>
#include <TTree.h>

#include <cstring>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <limits>
#include <tuple>
#include <typeinfo>

#include <boost/lexical_cast.hpp>


using namespace std;

SimEval::SimEval(const std::string& name) :
  SubsysReco(name),
  _event(0),
  _event_header(nullptr),
  _hit_map(nullptr),
  _hit_vector(nullptr),
  _out_name("SimEval.root")
{
}

int SimEval::Init(PHCompositeNode* topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int SimEval::InitRun(PHCompositeNode* topNode) {
  ResetEvalVars();
  InitEvalTree();
  p_geomSvc = GeomSvc::instance();
  
  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int SimEval::process_event(PHCompositeNode* topNode) {
  int ret = Fun4AllReturnCodes::ABORTRUN;
  
  ret = TruthEval(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  ++_event;
  
  return ret;
}

// find g4shower
PHG4Shower* SimEval::get_primary_shower(PHG4Particle* primary) {
  PHG4Shower* shower = nullptr;
  for (auto iter=_truth->GetPrimaryShowerRange().first; iter!=_truth->GetPrimaryShowerRange().second; ++iter) {
    PHG4Shower* tmpshower = iter->second;
    if (tmpshower->get_parent_particle_id() == primary->get_track_id()) {
      shower = tmpshower;
      break;
    }
    std::cout << " not finding shower " << std::endl;
  }
  return shower;
}

// find g4hit
PHG4Hit* SimEval::FindG4HitAtStation(const int trk_id, const PHG4HitContainer* g4hc) {
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

// for multiple hits
std::vector<PHG4Hit*> SimEval::FindG4HitsAtStation(const int trk_id, const PHG4HitContainer* g4hc) {
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

// per event
int SimEval::TruthEval(PHCompositeNode* topNode)
{
  ResetEvalVars();
  if(_event_header) {  
    event_id    = _event_header->get_event_id();    
  }

  if(_hit_vector) {
    n_hits = 0;

    // SQHit truth information =/= the information that you give the simulation.
    for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
      SQHit *hit = _hit_vector->at(ihit);
      
      if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) {
	hit->identify();
      }
      
      // get info about hit
      int hitID = hit->get_hit_id();
      hit_detid[n_hits]      = hit->get_detector_id();
      hit_elmid[n_hits]      = hit->get_element_id();
      hit_driftdis[n_hits]   = hit->get_drift_distance();
      hit_pos[n_hits]        = hit->get_pos();
      hit_detz[n_hits]       = p_geomSvc->getPlanePosition(hit->get_detector_id());
      hit_edep[n_hits]       = hit->get_edep();

      if(_truth) {
	int track_id = hit->get_track_id(); 
	int det_id = hit->get_detector_id();

	hit_truthx[n_hits] = hit->get_truth_x();
	hit_truthy[n_hits] = hit->get_truth_y();
	hit_truthz[n_hits] = hit->get_truth_z();	
	if(hit->get_detector_id()<100){ // something weird is happening here for emcal
	  double uVec[3] = {p_geomSvc->getPlane(hit->get_detector_id()).uVec[0],
			    p_geomSvc->getPlane(hit->get_detector_id()).uVec[1],
			    p_geomSvc->getPlane(hit->get_detector_id()).uVec[2]
	  };
	  hit_truthpos[n_hits] = (hit_truthx[n_hits])*uVec[0] + (hit_truthy[n_hits])*uVec[1] + (hit_truthz[n_hits]-p_geomSvc->getPlane(hit->get_detector_id()).zc)*uVec[2];
	}
      }
      ++n_hits;
      if(n_hits>=10000) break;
    }
  } // end SQ Hit_vector

  std::cout << "new loop over truth " << std::endl;

  // Info from PHG4Particle
  if(_truth) {

    // Showers
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
	std::cout << "shower " << shower->get_x() << " y " << shower->get_y() <<" z " << shower->get_z() << std::endl;
	n_showers++;
      }
      if(n_showers>=1000) break;
    }

    // Primaries
    n_primaries = 0;
    for(auto iterp=_truth->GetPrimaryParticleRange().first; iterp!=_truth->GetPrimaryParticleRange().second; ++iterp) {
      PHG4Particle * primary = iterp->second;
      gpid[n_primaries] = primary->get_pid();

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
	for(int iecal=0; iecal<g4hits.size(); ++iecal) {     
	  PHG4Hit* g4hit = g4hits[iecal];
	  gx_ecal[n_primaries][iecal] = g4hit->get_x(0);
          gy_ecal[n_primaries][iecal] = g4hit->get_y(0);
          gz_ecal[n_primaries][iecal] = g4hit->get_z(0);
          gpx_ecal[n_primaries][iecal] = g4hit->get_px(0)/1000.; 
          gpy_ecal[n_primaries][iecal] = g4hit->get_py(0)/1000.; 
          gpz_ecal[n_primaries][iecal] = g4hit->get_pz(0)/1000.; 
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
	gpx_st1[n_primaries] = st1hit->get_px(0)/1000.;
	gpy_st1[n_primaries] = st1hit->get_py(0)/1000.;
	gpz_st1[n_primaries] = st1hit->get_pz(0)/1000.;
      }
      PHG4Hit* st2hit = FindG4HitAtStation(trkID, g4hc_d2xp);
      if(st2hit){
	gx_st2[n_primaries] = st2hit->get_x(0);
	gy_st2[n_primaries] = st2hit->get_y(0);
        gz_st2[n_primaries] = st2hit->get_z(0);
	gpx_st2[n_primaries] = st2hit->get_px(0)/1000.;
	gpy_st2[n_primaries] = st2hit->get_py(0)/1000.;
	gpz_st2[n_primaries] = st2hit->get_pz(0)/1000.;
      }
      PHG4Hit* st3hit = FindG4HitAtStation(trkID, g4hc_d3px);
      if(!st3hit)
	PHG4Hit* st3hit = FindG4HitAtStation(trkID, g4hc_d3mx);
      if(st3hit){
	gx_st3[n_primaries] = st3hit->get_x(0);
	gy_st3[n_primaries] = st3hit->get_y(0);
	gz_st3[n_primaries] = st3hit->get_z(0);
	gpx_st3[n_primaries] = st3hit->get_px(0)/1000.;
	gpy_st3[n_primaries] = st3hit->get_py(0)/1000.;
	gpz_st3[n_primaries] = st3hit->get_pz(0)/1000.;
      }

      PHG4Hit* h1hit = FindG4HitAtStation(trkID,g4hc_h1t);
      if(!h1hit)
	PHG4Hit* h1hit = FindG4HitAtStation(trkID,g4hc_h1b);
      if(h1hit){
        gx_h1[n_primaries] = h1hit->get_x(0);
        gy_h1[n_primaries] = h1hit->get_y(0);
        gz_h1[n_primaries] = h1hit->get_z(0);
        gpx_h1[n_primaries] = h1hit->get_px(0)/1000.;
        gpy_h1[n_primaries] = h1hit->get_py(0)/1000.;
        gpz_h1[n_primaries] = h1hit->get_pz(0)/1000.;
      }
      PHG4Hit* h2hit = FindG4HitAtStation(trkID,g4hc_h2t);
      if(!h2hit)
        PHG4Hit* h2hit = FindG4HitAtStation(trkID,g4hc_h2b);
      if(h2hit){
        gx_h2[n_primaries] = h2hit->get_x(0);
        gy_h2[n_primaries] = h2hit->get_y(0);
        gz_h2[n_primaries] = h2hit->get_z(0);
	std::cout << "h2 x " << h2hit->get_x(0) << " y " << h2hit->get_y(0) <<" z " << h2hit->get_z(0) << std::endl;
        gpx_h2[n_primaries] = h2hit->get_px(0)/1000.;
        gpy_h2[n_primaries] = h2hit->get_py(0)/1000.;
        gpz_h2[n_primaries] = h2hit->get_pz(0)/1000.;
      }
      PHG4Hit* h3hit = FindG4HitAtStation(trkID,g4hc_h3t);
      if(!h3hit)
        PHG4Hit* h3hit = FindG4HitAtStation(trkID,g4hc_h3b);
      if(h3hit){
        gx_h3[n_primaries] = h3hit->get_x(0);
        gy_h3[n_primaries] = h3hit->get_y(0);
        gz_h3[n_primaries] = h3hit->get_z(0);
        gpx_h3[n_primaries] = h3hit->get_px(0)/1000.;
        gpy_h3[n_primaries] = h3hit->get_py(0)/1000.;
        gpz_h3[n_primaries] = h3hit->get_pz(0)/1000.;
      }
      PHG4Hit* h4hit = FindG4HitAtStation(trkID,g4hc_h4t);
      if(!h4hit)
        PHG4Hit* h4hit = FindG4HitAtStation(trkID,g4hc_h4b);
      if(h4hit){
        gx_h4[n_primaries] = h4hit->get_x(0);
        gy_h4[n_primaries] = h4hit->get_y(0);
        gz_h4[n_primaries] = h4hit->get_z(0);
	std::cout << "h4y t/b x " << h4hit->get_x(0) << " y " << h4hit->get_y(0) <<" z " << h4hit->get_z(0) << std::endl;
        gpx_h4[n_primaries] = h4hit->get_px(0)/1000.;
        gpy_h4[n_primaries] = h4hit->get_py(0)/1000.;
        gpz_h4[n_primaries] = h4hit->get_pz(0)/1000.;
      }
      PHG4Hit* h4y2lhit = FindG4HitAtStation(trkID,g4hc_h4y2l);
      if(h4y2lhit){
	gx_h4y2l[n_primaries] = h4y2lhit->get_x(0);
        gy_h4y2l[n_primaries] = h4y2lhit->get_y(0);
        gz_h4y2l[n_primaries] = h4y2lhit->get_z(0);
	std::cout << "h4y2l x " << h4y2lhit->get_x(0) << " y " << h4y2lhit->get_y(0) << " z " << h4y2lhit->get_z(0) << std::endl;
        gpx_h4y2l[n_primaries] = h4y2lhit->get_px(0)/1000.;
        gpy_h4y2l[n_primaries] = h4y2lhit->get_py(0)/1000.;
        gpz_h4y2l[n_primaries] = h4y2lhit->get_pz(0)/1000.;
      }
      PHG4Hit* h4y2rhit = FindG4HitAtStation(trkID,g4hc_h4y2r);
      if(h4y2rhit){
        gx_h4y2r[n_primaries] = h4y2rhit->get_x(0);
        gy_h4y2r[n_primaries] = h4y2rhit->get_y(0);
        gz_h4y2r[n_primaries] = h4y2rhit->get_z(0);
        gpx_h4y2r[n_primaries] = h4y2rhit->get_px(0)/1000.;
        gpy_h4y2r[n_primaries] = h4y2rhit->get_py(0)/1000.;
        gpz_h4y2r[n_primaries] = h4y2rhit->get_pz(0)/1000.;
      }
      ++n_primaries;
    } 
  }
  
  _tout_truth->Fill();

  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
    std::cout << "Leaving SimEval::TruthEval: " << _event << std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int SimEval::End(PHCompositeNode* topNode) {
  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
    std::cout << "SimEval::End" << std::endl;
  
  PHTFileServer::get().cd(_out_name.c_str());
  _tout_truth->Write();
  return Fun4AllReturnCodes::EVENT_OK;
}

int SimEval::InitEvalTree() {
  PHTFileServer::get().open(_out_name.c_str(), "RECREATE");

  _tout_truth = new TTree("Truth", "Truth Eval");
  _tout_truth->Branch("eventID",       &event_id,        "eventID/I");

  _tout_truth->Branch("n_hits",        &n_hits,          "n_hits/I");
  _tout_truth->Branch("hit_detID",     hit_detid,        "hit_detID[n_hits]/I");
  _tout_truth->Branch("hit_elmID",     hit_elmid,        "hit_elmID[n_hits]/I");
  _tout_truth->Branch("hit_driftdis",  hit_driftdis,     "hit_driftdis[n_hits]/F");
  _tout_truth->Branch("hit_pos",       hit_pos,          "hit_pos[n_hits]/F");
  _tout_truth->Branch("hit_detZ",      hit_detz,         "hit_detZ[n_hits]/F");
  _tout_truth->Branch("hit_edep",      hit_edep,         "hit_edep[n_hits]/F");
  _tout_truth->Branch("hit_truthx",    hit_truthx,       "hit_truthx[n_hits]/F");
  _tout_truth->Branch("hit_truthy",    hit_truthy,       "hit_truthy[n_hits]/F");
  _tout_truth->Branch("hit_truthz",    hit_truthz,       "hit_truthz[n_hits]/F");
  _tout_truth->Branch("hit_truthpos",  hit_truthpos,     "hit_truthpos[n_hits]/F");

  _tout_truth->Branch("n_showers",     &n_showers,       "n_showers/I");
  _tout_truth->Branch("sx_ecal",       &sx_ecal,         "sx_ecal[n_showers]/F");
  _tout_truth->Branch("sy_ecal",       &sy_ecal,         "sy_ecal[n_showers]/F");
  _tout_truth->Branch("sz_ecal",       &sz_ecal,         "sz_ecal[n_showers]/F");
  _tout_truth->Branch("sedep_ecal",    &sedep_ecal,      "sedep_ecal[n_showers]/F");

  _tout_truth->Branch("n_primaries",   &n_primaries,        "n_primaries/I");
  _tout_truth->Branch("gtrkid",        gtrkid,              "gtrkid[n_primaries]/I");
  _tout_truth->Branch("gpid",          gpid,                "gpid[n_primaries]/I");
  _tout_truth->Branch("gvx",           gvx,                 "gvx[n_primaries]/F");
  _tout_truth->Branch("gvy",           gvy,                 "gvy[n_primaries]/F");
  _tout_truth->Branch("gvz",           gvz,                 "gvz[n_primaries]/F");
  _tout_truth->Branch("gpx",           gpx,                 "gpx[n_primaries]/F");
  _tout_truth->Branch("gpy",           gpy,                 "gpy[n_primaries]/F");
  _tout_truth->Branch("gpz",           gpz,                 "gpz[n_primaries]/F");
  _tout_truth->Branch("gpt",           gpt,                 "gpt[n_primaries]/F");
  _tout_truth->Branch("geta",          geta,                "geta[n_primaries]/F");
  _tout_truth->Branch("gphi",          gphi,                "gphi[n_primaries]/F");
  _tout_truth->Branch("ge",            ge,                  "ge[n_primaries]/F");

  _tout_truth->Branch("nhits_ecal",    nhits_ecal,          "nhits_ecal[n_primaries]/I");
  _tout_truth->Branch("gx_ecal",       gx_ecal,             "gx_ecal[n_primaries][100]/F"); // not sure how to make this with the right size
  _tout_truth->Branch("gy_ecal",       gy_ecal,             "gy_ecal[n_primaries][100]/F");
  _tout_truth->Branch("gz_ecal",       gz_ecal,             "gz_ecal[n_primaries][100]/F");
  _tout_truth->Branch("gpx_ecal",      gpx_ecal,            "gpx_ecal[n_primaries][100]/F");
  _tout_truth->Branch("gpy_ecal",      gpy_ecal,            "gpy_ecal[n_primaries][100]/F");
  _tout_truth->Branch("gpz_ecal",      gpz_ecal,            "gpz_ecal[n_primaries][100]/F");
  _tout_truth->Branch("gedep_ecal",    gedep_ecal,          "gedep_ecal[n_primaries][100]/F");

  _tout_truth->Branch("gx_st1",        gx_st1,              "gx_st1[n_primaries]/F");
  _tout_truth->Branch("gy_st1",        gy_st1,              "gy_st1[n_primaries]/F");
  _tout_truth->Branch("gz_st1",        gz_st1,              "gz_st1[n_primaries]/F");
  _tout_truth->Branch("gpx_st1",       gpx_st1,             "gpx_st1[n_primaries]/F");
  _tout_truth->Branch("gpy_st1",       gpy_st1,             "gpy_st1[n_primaries]/F");
  _tout_truth->Branch("gpz_st1",       gpz_st1,             "gpz_st1[n_primaries]/F");
  _tout_truth->Branch("gx_st2",        gx_st2,              "gx_st2[n_primaries]/F");
  _tout_truth->Branch("gy_st2",        gy_st2,              "gy_st2[n_primaries]/F");
  _tout_truth->Branch("gz_st2",        gz_st2,              "gz_st2[n_primaries]/F");
  _tout_truth->Branch("gpx_st2",       gpx_st2,             "gpx_st2[n_primaries]/F");
  _tout_truth->Branch("gpy_st2",       gpy_st2,             "gpy_st2[n_primaries]/F");
  _tout_truth->Branch("gpz_st2",       gpz_st2,             "gpz_st2[n_primaries]/F");
  _tout_truth->Branch("gx_st3",        gx_st3,              "gx_st3[n_primaries]/F");
  _tout_truth->Branch("gy_st3",        gy_st3,              "gy_st3[n_primaries]/F");
  _tout_truth->Branch("gz_st3",        gz_st3,              "gz_st3[n_primaries]/F");
  _tout_truth->Branch("gpx_st3",       gpx_st3,             "gpx_st3[n_primaries]/F");
  _tout_truth->Branch("gpy_st3",       gpy_st3,             "gpy_st3[n_primaries]/F");
  _tout_truth->Branch("gpz_st3",       gpz_st3,             "gpz_st3[n_primaries]/F");

  _tout_truth->Branch("gx_h1",         gx_h1,               "gx_h1[n_primaries]/F");
  _tout_truth->Branch("gy_h1",         gy_h1,               "gy_h1[n_primaries]/F");
  _tout_truth->Branch("gz_h1",         gz_h1,               "gz_h1[n_primaries]/F");
  _tout_truth->Branch("gpx_h1",        gpx_h1,              "gpx_h1[n_primaries]/F");
  _tout_truth->Branch("gpy_h1",        gpy_h1,              "gpy_h1[n_primaries]/F");
  _tout_truth->Branch("gpz_h1",        gpz_h1,              "gpz_h1[n_primaries]/F");
  _tout_truth->Branch("gx_h2",         gx_h2,               "gx_h2[n_primaries]/F");
  _tout_truth->Branch("gy_h2",         gy_h2,               "gy_h2[n_primaries]/F");
  _tout_truth->Branch("gz_h2",         gz_h2,               "gz_h2[n_primaries]/F");
  _tout_truth->Branch("gpx_h2",        gpx_h2,              "gpx_h2[n_primaries]/F");
  _tout_truth->Branch("gpy_h2",        gpy_h2,              "gpy_h2[n_primaries]/F");
  _tout_truth->Branch("gpz_h2",        gpz_h2,              "gpz_h2[n_primaries]/F");

  _tout_truth->Branch("gx_p1",         gx_p1,               "gx_p1[n_primaries]/F");
  _tout_truth->Branch("gy_p1",         gy_p1,               "gy_p1[n_primaries]/F");
  _tout_truth->Branch("gz_p1",         gz_p1,               "gz_p1[n_primaries]/F");
  _tout_truth->Branch("gpx_p1",        gpx_p1,              "gpx_p1[n_primaries]/F");
  _tout_truth->Branch("gpy_p1",        gpy_p1,              "gpy_p1[n_primaries]/F");
  _tout_truth->Branch("gpz_p1",        gpz_p1,              "gpz_p1[n_primaries]/F");

  _tout_truth->Branch("gx_h4",         gx_h4,              "gx_h4[n_primaries]/F");
  _tout_truth->Branch("gy_h4",         gy_h4,              "gy_h4[n_primaries]/F");
  _tout_truth->Branch("gz_h4",         gz_h4,              "gz_h4[n_primaries]/F");
  _tout_truth->Branch("gpx_h4",        gpx_h4,             "gpx_h4[n_primaries]/F");
  _tout_truth->Branch("gpy_h4",        gpy_h4,             "gpy_h4[n_primaries]/F");
  _tout_truth->Branch("gpz_h4",        gpz_h4,             "gpz_h4[n_primaries]/F");

  _tout_truth->Branch("gx_h4y2l",      gx_h4y2l,            "gx_h4y2l[n_primaries]/F");
  _tout_truth->Branch("gy_h4y2l",      gy_h4y2l,            "gy_h4y2l[n_primaries]/F");
  _tout_truth->Branch("gz_h4y2l",      gz_h4y2l,            "gz_h4y2l[n_primaries]/F");
  _tout_truth->Branch("gpx_h4y2l",     gpx_h4y2l,           "gpx_h4y2l[n_primaries]/F");
  _tout_truth->Branch("gpy_h4y2l",     gpy_h4y2l,           "gpy_h4y2l[n_primaries]/F");
  _tout_truth->Branch("gpz_h4y2l",     gpz_h4y2l,           "gpz_h4y2l[n_primaries]/F");
  _tout_truth->Branch("gx_h4y2r",      gx_h4y2r,            "gx_h4y2r[n_primaries]/F");
  _tout_truth->Branch("gy_h4y2r",      gy_h4y2r,            "gy_h4y2r[n_primaries]/F");
  _tout_truth->Branch("gz_h4y2r",      gz_h4y2r,            "gz_h4y2r[n_primaries]/F");
  _tout_truth->Branch("gpx_h4y2r",     gpx_h4y2r,           "gpx_h4y2r[n_primaries]/F");
  _tout_truth->Branch("gpy_h4y2r",     gpy_h4y2r,           "gpy_h4y2r[n_primaries]/F");
  _tout_truth->Branch("gpz_h4y2r",     gpz_h4y2r,           "gpz_h4y2r[n_primaries]/F");

  _tout_truth->Branch("gx_dp1",        gx_dp1,              "gx_dp1[n_primaries]/F");
  _tout_truth->Branch("gy_dp1",        gy_dp1,              "gy_dp1[n_primaries]/F");
  _tout_truth->Branch("gz_dp1",        gz_dp1,              "gz_dp1[n_primaries]/F");
  _tout_truth->Branch("gpx_dp1",       gpx_dp1,             "gpx_dp1[n_primaries]/F");
  _tout_truth->Branch("gpy_dp1",       gpy_dp1,             "gpy_dp1[n_primaries]/F");
  _tout_truth->Branch("gpz_dp1",       gpz_dp1,             "gpz_dp1[n_primaries]/F");

  _tout_truth->Branch("gx_dp2",        gx_dp2,              "gx_dp2[n_primaries]/F");
  _tout_truth->Branch("gy_dp2",        gy_dp2,              "gy_dp2[n_primaries]/F");
  _tout_truth->Branch("gz_dp2",        gz_dp2,              "gz_dp2[n_primaries]/F");
  _tout_truth->Branch("gpx_dp2",       gpx_dp2,             "gpx_dp2[n_primaries]/F");
  _tout_truth->Branch("gpy_dp2",       gpy_dp2,             "gpy_dp2[n_primaries]/F");
  _tout_truth->Branch("gpz_dp2",       gpz_dp2,             "gpz_dp2[n_primaries]/F");

  return 0;
}

int SimEval::ResetEvalVars() {
  event_id = std::numeric_limits<int>::max();
 
  n_hits = 0;
  for(int i=0; i<10000; ++i) {
    hit_detid[i]        = std::numeric_limits<short>::max();
    hit_elmid[i]        = std::numeric_limits<short>::max();
    hit_driftdis[i]     = std::numeric_limits<float>::max();
    hit_pos[i]          = std::numeric_limits<float>::max();
    hit_detz[i]         = std::numeric_limits<float>::max();
    hit_edep[i]         = std::numeric_limits<float>::max();

    hit_truthx[i]       = std::numeric_limits<float>::max();
    hit_truthy[i]       = std::numeric_limits<float>::max();
    hit_truthz[i]       = std::numeric_limits<float>::max();
    hit_truthpos[i]     = std::numeric_limits<float>::max();
  }

  n_showers =0;
  for(int i=0; i<1000; ++i) {                                                                                                                                                           sx_ecal[i]          = std::numeric_limits<float>::max();
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

    gx_h1[i]      = std::numeric_limits<float>::max();
    gy_h1[i]      = std::numeric_limits<float>::max();
    gz_h1[i]      = std::numeric_limits<float>::max();
    gpx_h1[i]     = std::numeric_limits<float>::max();
    gpy_h1[i]     = std::numeric_limits<float>::max();
    gpz_h1[i]     = std::numeric_limits<float>::max();
    gx_h2[i]      = std::numeric_limits<float>::max();
    gy_h2[i]      = std::numeric_limits<float>::max();
    gz_h2[i]      = std::numeric_limits<float>::max();
    gpx_h2[i]     = std::numeric_limits<float>::max();
    gpy_h2[i]     = std::numeric_limits<float>::max();
    gpz_h2[i]     = std::numeric_limits<float>::max();

    gx_p1[i]      = std::numeric_limits<float>::max();
    gy_p1[i]      = std::numeric_limits<float>::max();
    gz_p1[i]      = std::numeric_limits<float>::max();
    gpx_p1[i]     = std::numeric_limits<float>::max();
    gpy_p1[i]     = std::numeric_limits<float>::max();
    gpz_p1[i]     = std::numeric_limits<float>::max();

    gx_dp1[i]     = std::numeric_limits<float>::max();
    gy_dp1[i]     = std::numeric_limits<float>::max();
    gz_dp1[i]     = std::numeric_limits<float>::max();
    gpx_dp1[i]     = std::numeric_limits<float>::max();
    gpy_dp1[i]     = std::numeric_limits<float>::max();
    gpz_dp1[i]     = std::numeric_limits<float>::max();
    gx_dp2[i]     = std::numeric_limits<float>::max();
    gy_dp2[i]     = std::numeric_limits<float>::max();
    gz_dp2[i]     = std::numeric_limits<float>::max();
    gpx_dp2[i]     = std::numeric_limits<float>::max();
    gpy_dp2[i]     = std::numeric_limits<float>::max();
    gpz_dp2[i]     = std::numeric_limits<float>::max();

    gx_h4y2l[i]     = std::numeric_limits<float>::max();
    gy_h4y2l[i]     = std::numeric_limits<float>::max();
    gz_h4y2l[i]     = std::numeric_limits<float>::max();
    gpx_h4y2l[i]     = std::numeric_limits<float>::max();
    gpy_h4y2l[i]     = std::numeric_limits<float>::max();
    gpz_h4y2l[i]     = std::numeric_limits<float>::max();
    gx_h4y2r[i]     = std::numeric_limits<float>::max();
    gy_h4y2r[i]     = std::numeric_limits<float>::max();
    gz_h4y2r[i]     = std::numeric_limits<float>::max();
    gpx_h4y2r[i]     = std::numeric_limits<float>::max();
    gpy_h4y2r[i]     = std::numeric_limits<float>::max();
    gpz_h4y2r[i]     = std::numeric_limits<float>::max();

  }

  return 0;
}

int SimEval::GetNodes(PHCompositeNode* topNode) {

  _event_header = findNode::getClass<SQEvent>(topNode, "SQEvent");
  if (!_event_header) {
    cout << "!_event_header" << endl;
  }
  /*
  if(_hit_container_type.find("Map") != std::string::npos) {
    _hit_map = findNode::getClass<SQHitMap>(topNode, "SQHitMap");
    if (!_hit_map) {
      cout << "!_hit_map" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    }
  */
  //if(_hit_container_type.find("Vector") != std::string::npos) {
  _hit_vector = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!_hit_vector) {
    cout << "!_hit_vector" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  //}

  _truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth) {
    cout << "!_truth" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Drift chambers
  g4hc_d1x  = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D1X");
  g4hc_d2xp = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D2Xp");
  g4hc_d3px = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3pXp");
  g4hc_d3mx = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3mXp");
  if (! g4hc_d1x) g4hc_d1x = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D0X");

  if ( !g4hc_d1x || !g4hc_d3px || !g4hc_d3mx) {
    cout << "Failed at getting nodes dc "<< endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Hodoscopes
  g4hc_h1t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1T");
  g4hc_h1b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1B");
  g4hc_h2t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2T");
  g4hc_h2b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2B");
  g4hc_h3t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H3T");
  g4hc_h3b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H3B");
  g4hc_h4t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4T");
  g4hc_h4b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4B");
  g4hc_h4y2l= findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4Y2L");
  g4hc_h4y2r= findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4Y2R");
  if (!g4hc_h1t || !g4hc_h1b || !g4hc_h2t || !g4hc_h2b ||
      !g4hc_h3t || !g4hc_h3b || !g4hc_h4t || !g4hc_h4b || !g4hc_h4y2l || !g4hc_h4y2r  ) {
    cout << "Failed at getting nodes hodos" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Prop tubes 
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
    cout << "Failed at getting nodes prop" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Ecal
  g4hc_ecal  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_EMCal");
  if (!g4hc_ecal) {
    cout << "Failed at getting nodes ecal " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
