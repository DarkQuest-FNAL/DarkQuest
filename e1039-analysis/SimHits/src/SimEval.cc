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

#define LogDebug(exp) std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogError(exp) std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogWarning(exp) std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl

using namespace std;

SimEval::SimEval(const std::string& name) :
  SubsysReco(name),
  _hit_container_type("Vector"),
  _event(0),
  _run_header(nullptr),
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

// These Eval functions run per event.
// They limit on how many hits/tracks are in one single event
int SimEval::TruthEval(PHCompositeNode* topNode)
{
  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
    std::cout << "Entering SimEval::TruthEval: " << _event << std::endl;
  
  ResetEvalVars();

  if(_event_header) {  
    event_id    = _event_header->get_event_id();    
  }

  // These maps are just counting the number of hits
  std::map<int, int> parID_nhits_dc;
  std::map<int, int> parID_nhits_hodo;
  std::map<int, int> parID_nhits_prop;
  std::map<int, int> parID_nhits_dp;
  std::map<int, int> parID_nhits_H4Y;

  std::map<int, std::map<int, int> > parID_detid_elmid;

  typedef std::tuple<int, int> ParDetPair;
  std::map<ParDetPair, int> parID_detID_ihit;
  std::map<int, int> hitID_ihit;
  
  if(_hit_vector) {
    n_hits = 0;

    // this is not the PHG4 truth information ("true truth") information.
    // SQHit truth information =/= the information that you give the simulation.
    // For example, if the information you give says hit 5 should be at 1,1,1
    // the "truth" information from SQHit might not give 1,1,1 below.

    for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
      SQHit *hit = _hit_vector->at(ihit);
      
      if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) {
	LogInfo(hit->get_detector_id());
	hit->identify();
      }
      
      // get info about hit
      int hitID = hit->get_hit_id();
      hit_detid[n_hits]      = hit->get_detector_id();
      hit_elmid[n_hits]      = hit->get_element_id();
      hit_driftdis[n_hits]   = hit->get_drift_distance();
      hit_pos[n_hits]        = hit->get_pos();
      hit_detz[n_hits]       = p_geomSvc->getPlanePosition(hit->get_detector_id());
      
      if(_truth) {
	int track_id = hit->get_track_id();
	int det_id = hit->get_detector_id();

	// NOTE HERE:
	// parID is not the particle ID, but the trackID.
	// It indexes by n_hits, starting from 1
	// Anything indexed by n_tracks starts at 0
	// Make sure you pay attention what is indexed where. 
	
	parID_detID_ihit[std::make_tuple(track_id, det_id)] = ihit;
	
	auto detid_elmid_iter = parID_detid_elmid.find(track_id);
	if(detid_elmid_iter != parID_detid_elmid.end()) {
	  detid_elmid_iter->second.insert(std::pair<int, int>(det_id, hit->get_element_id()));
	} 
	else {
	  std::map<int, int> detid_elmid;
	  detid_elmid.insert(std::pair<int, int>(det_id, hit->get_element_id()));
	  parID_detid_elmid[track_id] = detid_elmid;
	}
	
	// drift chamber hits
	if(hit->get_detector_id() >= 1 and hit->get_detector_id() <=30) {
	  if(parID_nhits_dc.find(track_id)!=parID_nhits_dc.end())
	    parID_nhits_dc[track_id] = parID_nhits_dc[track_id]+1;
	  else
	    parID_nhits_dc[track_id] = 1;
	}

	// hodoscope hits
	if(hit->get_detector_id() >= 31 and hit->get_detector_id() <=46) {
	  if(parID_nhits_hodo.find(track_id)!=parID_nhits_hodo.end())
	    parID_nhits_hodo[track_id] = parID_nhits_hodo[track_id]+1;
	  else
	    parID_nhits_hodo[track_id] = 1;
	}
	
	// prop tube hits
	if(hit->get_detector_id() >= 47 and hit->get_detector_id() <=54) {
	  if(parID_nhits_prop.find(track_id)!=parID_nhits_prop.end())
	    parID_nhits_prop[track_id] = parID_nhits_prop[track_id]+1;
	  else
	    parID_nhits_prop[track_id] = 1;
	}
	
	// dp hits
	if((hit->get_detector_id() >= 55 and hit->get_detector_id() <=62)) {
	  if(parID_nhits_dp.find(track_id)!=parID_nhits_dp.end())
	    parID_nhits_dp[track_id] = parID_nhits_dp[track_id]+1;
	  else
	    parID_nhits_dp[track_id] = 1;
	}
	
	// h4y hits
	if((hit->get_detector_id() == 43 or hit->get_detector_id() ==44)) {
	  //it is possible that a hit will hit both det 43 and 44 since they over lap,
	  //but we are only interested in if it hits 1 of the h4y, so as long as it hits, take that 
	  // to be 1 hit.
	  parID_nhits_H4Y[track_id] = 1;
	}
	
	//PHG4Hit is what has the true truth information.
	hit_truthx[n_hits] = hit->get_truth_x();
	hit_truthy[n_hits] = hit->get_truth_y();
	hit_truthz[n_hits] = hit->get_truth_z();
	
	double uVec[3] = {p_geomSvc->getPlane(hit->get_detector_id()).uVec[0],
			  p_geomSvc->getPlane(hit->get_detector_id()).uVec[1],
			  p_geomSvc->getPlane(hit->get_detector_id()).uVec[2]
	};
	hit_truthpos[n_hits] = (hit_truthx[n_hits])*uVec[0] + (hit_truthy[n_hits])*uVec[1] + (hit_truthz[n_hits]-p_geomSvc->getPlane(hit->get_detector_id()).zc)*uVec[2];
	
      }
      ++n_hits;
      if(n_hits>=10000) break;
    }
  } // end SQ Hit_vector
  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) LogInfo("ghit eval finished");
  
  // PHG4Particle is the particle information from geant. This has all the information 
  // that you put into the simulation and all the information the geant created when 
  // propogating the particle
  if(_truth) {
    for(auto iter=_truth->GetPrimaryParticleRange().first; iter!=_truth->GetPrimaryParticleRange().second; ++iter) {
      if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
	std::cout << "looping over phg4particle " << std::endl;
      PHG4Particle * par = iter->second;
      gpid[n_tracks] = par->get_pid();

      int vtx_id =  par->get_vtx_id();
      PHG4VtxPoint* vtx = _truth->GetVtx(vtx_id);
      gvx[n_tracks] = vtx->get_x();
      gvy[n_tracks] = vtx->get_y();
      gvz[n_tracks] = vtx->get_z();
      
      TVector3 mom(par->get_px(), par->get_py(), par->get_pz());
      gpx[n_tracks] = par->get_px();
      gpy[n_tracks] = par->get_py();
      gpz[n_tracks] = par->get_pz();
      gpt[n_tracks] = mom.Pt();
      geta[n_tracks] = mom.Eta();
      gphi[n_tracks] = mom.Phi();
      
      // Not particle ID, trackID
      int parID = par->get_track_id();
      gparid[n_tracks] = parID;
      
      // The detector ID and names are listed in e1039-core/packages/geom_svc/GeomSvc.cxx.

      // st1, 2, 3 are the drift chambers here   
      PHG4HitContainer *D1X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D0X");
      if (!D1X_hits)
	D1X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D1X");
      
      if (!D1X_hits){
	if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
	  cout << Name() << " Could not locate g4 hit node " << "G4HIT_D0X or G4HIT_D1X" << endl;
      }
      
      // trackID + detID -> SQHit -> PHG4Hit -> momentum
      // detID 1-6 deal with D0, 7-12 D1
      for(int det_id=1; det_id<=12; ++det_id) {
	auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	if(iter != parID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and D1X_hits) {
	    PHG4Hit* g4hit =  D1X_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_st1[n_tracks]  = g4hit->get_x(0);
	      gy_st1[n_tracks]  = g4hit->get_y(0);
	      gz_st1[n_tracks]  = g4hit->get_z(0);
	      gpx_st1[n_tracks] = g4hit->get_px(0)/1000.;
	      gpy_st1[n_tracks] = g4hit->get_py(0)/1000.;
	      gpz_st1[n_tracks] = g4hit->get_pz(0)/1000.;
	      if(gpz_st1[n_tracks] <0){
		std::cout << "WARNING:: Negative z-momentum at Station 1!" << std::endl;
	      }
	      break;
	    }
	  }
	}
      }// end st1 det id loop  
      if(verbosity>=2) std::cout << "station 1 truth info done." << std::endl;
      
      PHG4HitContainer *D2X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D2X");
      if(!D2X_hits){
	D2X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D2Xp");
      }
      if (!D2X_hits)
	{
	  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
	    cout << Name() << " Could not locate g4 hit node " <<  "G4HIT_D2X" << endl;
	}
      
      // detID 13-18 are D2
      // trackID + detID -> SQHit -> PHG4Hit -> momentum
      for(int det_id=13; det_id<=18; ++det_id) {
	auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	if(iter != parID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and D2X_hits) {
	    PHG4Hit* g4hit =  D2X_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_st2[n_tracks]  = g4hit->get_x(0);
	      gy_st2[n_tracks]  = g4hit->get_y(0);
	      gz_st2[n_tracks]  = g4hit->get_z(0);
	      gpx_st2[n_tracks] = g4hit->get_px(0)/1000.;
	      gpy_st2[n_tracks] = g4hit->get_py(0)/1000.;
	      gpz_st2[n_tracks] = g4hit->get_pz(0)/1000.;
	      if(gpz_st2[n_tracks] <0){
		std::cout << "WARNING:: Negative z-momentum at Station 2!" << std::endl;
	      }
	      break;
	    }
	  }
	}
      }// end st2 det id loop  
      if(verbosity>=2) std::cout << "station 2 truth info done." << std::endl;
      
      PHG4HitContainer *D3X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D3pX");
      if (!D3X_hits){
	D3X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D3pXp");
      }
      
      // detID 19-24 are D3p, 25-30 are D3m
      for(int det_id=19; det_id<=24; ++det_id) {
	auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	if(iter != parID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and D3X_hits) {
	    PHG4Hit* g4hit =  D3X_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_st3[n_tracks]  = g4hit->get_x(0);
	      gy_st3[n_tracks]  = g4hit->get_y(0);
	      gz_st3[n_tracks]  = g4hit->get_z(0);
	      gpx_st3[n_tracks] = g4hit->get_px(0)/1000.;
	      gpy_st3[n_tracks] = g4hit->get_py(0)/1000.;
	      gpz_st3[n_tracks] = g4hit->get_pz(0)/1000.;
	      if(gpz_st3[n_tracks] <0){
		std::cout << "WARNING:: Negative z-momentum at Station 3!" << std::endl;
	      }
	      break;
	    }
	  }
	}
      } // end st3 det id loop
      
      D3X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D3mX"); 
      if (!D3X_hits){
	D3X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D3mXp");
      }
      
      if(!D3X_hits){
	std::cout << "Could not locate D3X container." << std::endl;
      }
      
      // detID  25-30 are D3m
      for(int det_id=25; det_id<=30; ++det_id) {
	auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	if(iter != parID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and D3X_hits) {
	    PHG4Hit* g4hit =  D3X_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_st3[n_tracks]  = g4hit->get_x(0);
	      gy_st3[n_tracks]  = g4hit->get_y(0);
	      gz_st3[n_tracks]  = g4hit->get_z(0);
	      gpx_st3[n_tracks] = g4hit->get_px(0)/1000.;
	      gpy_st3[n_tracks] = g4hit->get_py(0)/1000.;
	      gpz_st3[n_tracks] = g4hit->get_pz(0)/1000.;
	      if(gpz_st3[n_tracks] <0){
		std::cout << "WARNING:: Negative z-momentum at Station 3!" << std::endl;
	      }
	      break;
	    }
	  }
	}
      }
      if(verbosity>=2) std::cout << "station 3 truth info done." << std::endl;
      
      // H1_H3 are hodoscopes here
      // detID 33 and 34 are for H1L/R. This is to check possible roads, if comment
      // out if you do not want this info.
      PHG4HitContainer *H1_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1LL");
      if(!H1_hits){
	H1_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1R");
      }
      
      for(int det_id=33; det_id<=34; ++det_id) {
	auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	if(iter != parID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and H1_hits) {
	    PHG4Hit* g4hit =  H1_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      int h1barID = hit->get_element_id();
	      if(h1barID<11) gbarID_h1[n_tracks] = hit->get_element_id();
	      if(h1barID>10) gbarID_h1[n_tracks] = hit->get_element_id() - 8;
	      break;
	    }
	  }
	}
      }
      
      // detID 35 and 36 are for H2L/R. This is to check possible roads, if comment
      // out if you do not want this info.
      PHG4HitContainer *H2_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2LL");
      if(!H2_hits){
	H2_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2R");
      }
      
      for(int det_id=35; det_id<=36; ++det_id) {
	auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	if(iter != parID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and H2_hits) {
	    PHG4Hit* g4hit =  H2_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      int h2barID = hit->get_element_id();
	      if(h2barID<11) gbarID_h2[n_tracks] = hit->get_element_id();
	      if(h2barID>9) gbarID_h2[n_tracks] = hit->get_element_id() - 8;
	      break;
	    }
	  }
	}
      }
      
      //detID 43 and 44 are for H4Y2L, which is the one used in DP tracking.
      PHG4HitContainer *H4Y_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4Y2L");
      for(int det_id=43; det_id<44; ++det_id) {
	auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	if(iter != parID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and H4Y_hits) {
	    if(verbosity >= Fun4AllBase::VERBOSITY_A_LOT) {
	      LogDebug("h4y2lhit: " << iter->second);
	    }
	    PHG4Hit* g4hit =  H4Y_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      int h4ybarID = hit->get_element_id();
	      gbarID_h4y[n_tracks] = -hit->get_element_id();
	      if (h4ybarID > 8) {
		gquad_h4y[n_tracks] = 8;
	      } else {
		gquad_h4y[n_tracks] = 10;
	      }
	      break;
	    }
	  }
	}
      }
      
      H4Y_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4Y2R");
      for(int det_id=44; det_id<45; ++det_id) {
	auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	if(iter != parID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and H4Y_hits) {
	    PHG4Hit* g4hit =  H4Y_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      int h4ybarID = hit->get_element_id();
	      gbarID_h4y[n_tracks] = hit->get_element_id();
	      if (h4ybarID > 8) {
		gquad_h4y[n_tracks] = 9;
	      } else {
		gquad_h4y[n_tracks] = 11;
	      }
	      break;
	    }
	  }
	}
      }
      if(verbosity>=2) std::cout << "H4Y Truth info done." << std::endl;

      // DP1 
      PHG4HitContainer* DP1Container;
      for(int det_id = 55; det_id<=58; det_id++){
	switch(det_id){
	case 55:
	  DP1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP1TL");
	  break;
	  
	case 56:
	  DP1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP1TR");
	  break;
	  
	case 57:
	  DP1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP1BL");
	  break;
	  
	case 58:
	  DP1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP1BR");
	  break;
	  
	}
	
	auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	if(iter != parID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and DP1Container) {
	    PHG4Hit* g4hit =  DP1Container->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_dp1[n_tracks] = g4hit->get_x(0);
	      gy_dp1[n_tracks] = g4hit->get_y(0);
	      gz_dp1[n_tracks] = g4hit->get_z(0);
	      if(hit->get_detector_id()==det_id){
		switch(det_id){
		case 55:
		  gbarID_dp1[n_tracks] = hit->get_element_id();
		  gquad_dp1[n_tracks] = (80+hit->get_element_id())*(-1);
		  break;
		case 56:
		  gbarID_dp1[n_tracks] = hit->get_element_id();
		  gquad_dp1[n_tracks] = (80+hit->get_element_id());
		  break;
		case 57:
		  gbarID_dp1[n_tracks] = 81-hit->get_element_id();
		  gquad_dp1[n_tracks] = (hit->get_element_id())*(-1);
		  break;
		case 58:
		  gbarID_dp1[n_tracks] = 81-hit->get_element_id();
		  gquad_dp1[n_tracks] = hit->get_element_id();
		  break;
		}
	      }
	      break;
	    }
	  }
	}
      }
      if(verbosity>=2) std::cout << "DP1 Truth info done." << std::endl;
      
      // DP2
      PHG4HitContainer* DP2Container;
      for(int det_id = 59; det_id<=62; det_id++){
	switch(det_id){
	case 59:
	  DP2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP2TL");
	  break;
	  
	case 60:
	  DP2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP2TR");
	  break;
	  
	case 61:
	  DP2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP2BL");
	  break;
	  
	case 62:
	  DP2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP2BR");
	  break;
	}
	
	auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	if(iter != parID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and DP2Container) {
	    PHG4Hit* g4hit =  DP2Container->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_dp2[n_tracks] = g4hit->get_x(0);
	      gy_dp2[n_tracks] = g4hit->get_y(0);
	      gz_dp2[n_tracks] = g4hit->get_z(0);
	      if(hit->get_detector_id()==det_id){
		switch(det_id){
		case 59:
		  gbarID_dp2[n_tracks] = hit->get_element_id();
		  gquad_dp2[n_tracks] = (50+hit->get_element_id())*(-1);
		  break;
		case 60:
		  gbarID_dp2[n_tracks] = hit->get_element_id();
		  gquad_dp2[n_tracks] = (50+hit->get_element_id());
		  break;
		case 61:
		  gbarID_dp2[n_tracks] = 51-hit->get_element_id();
		  gquad_dp2[n_tracks] = (hit->get_element_id())*(-1);
		  break;
		case 62:
		  gbarID_dp2[n_tracks] = 51-hit->get_element_id();
		  gquad_dp2[n_tracks] = hit->get_element_id();
		  break;
		}
	      }
	      break;
	    }
	  }
	}
      }
      if(verbosity>=2) std::cout << "DP2 Truth info done." << std::endl;
	
      // total gen hits (except DP)
      gnhits[n_tracks] = parID_nhits_dc[parID] + parID_nhits_hodo[parID] + parID_nhits_prop[parID];
      gndc[n_tracks] = parID_nhits_dc[parID];
      gnhodo[n_tracks] = parID_nhits_hodo[parID];
      gnprop[n_tracks] = parID_nhits_prop[parID];
      gnDP[n_tracks] = parID_nhits_dp[parID];
      gnH4Y[n_tracks] = parID_nhits_H4Y[parID];
      
      for(auto detid_elmid : parID_detid_elmid[parID]) {
	int detID = detid_elmid.first;
	int elmID = detid_elmid.second;
	gelmid[n_tracks][detID] = elmID;
      }
      if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) LogInfo("particle eval finished");
      ++n_tracks;
      if(n_tracks>=1000) break;
    }
  	
    _tout_truth->Fill();
    
  }
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

  _tout_truth->Branch("nHits",         &n_hits,          "nHits/I");
  _tout_truth->Branch("hit_detID",     hit_detid,        "hit_detID[nHits]/I");
  _tout_truth->Branch("hit_elmID",     hit_elmid,        "hit_elmID[nHits]/I");
  _tout_truth->Branch("hit_driftdis",  hit_driftdis,     "hit_driftdis[nHits]/F");
  _tout_truth->Branch("hit_pos",       hit_pos,          "hit_pos[nHits]/F");
  _tout_truth->Branch("hit_detZ",      hit_detz,         "hit_detZ[nHits]/F");
  _tout_truth->Branch("hit_truthx",    hit_truthx,       "hit_truthx[nHits]/F");
  _tout_truth->Branch("hit_truthy",    hit_truthy,       "hit_truthy[nHits]/F");
  _tout_truth->Branch("hit_truthz",    hit_truthz,       "hit_truthz[nHits]/F");
  _tout_truth->Branch("hit_truthpos",  hit_truthpos,     "hit_truthpos[nHits]/F");

  _tout_truth->Branch("n_tracks",      &n_tracks,           "n_tracks/I");
  _tout_truth->Branch("gparid",        gparid,              "gparid[n_tracks]/I");
  _tout_truth->Branch("gpid",          gpid,                "gpid[n_tracks]/I");
  _tout_truth->Branch("gvx",           gvx,                 "gvx[n_tracks]/F");
  _tout_truth->Branch("gvy",           gvy,                 "gvy[n_tracks]/F");
  _tout_truth->Branch("gvz",           gvz,                 "gvz[n_tracks]/F");
  _tout_truth->Branch("gpx",           gpx,                 "gpx[n_tracks]/F");
  _tout_truth->Branch("gpy",           gpy,                 "gpy[n_tracks]/F");
  _tout_truth->Branch("gpz",           gpz,                 "gpz[n_tracks]/F");
  _tout_truth->Branch("gpt",           gpt,                 "gpt[n_tracks]/F");
  _tout_truth->Branch("geta",          geta,                "geta[n_tracks]/F");
  _tout_truth->Branch("gphi",          gphi,                "gphi[n_tracks]/F");
  _tout_truth->Branch("gx_st1",        gx_st1,              "gx_st1[n_tracks]/F");
  _tout_truth->Branch("gy_st1",        gy_st1,              "gy_st1[n_tracks]/F");
  _tout_truth->Branch("gz_st1",        gz_st1,              "gz_st1[n_tracks]/F");
  _tout_truth->Branch("gx_st2",        gx_st2,              "gx_st2[n_tracks]/F");
  _tout_truth->Branch("gy_st2",        gy_st2,              "gy_st2[n_tracks]/F");
  _tout_truth->Branch("gz_st2",        gz_st2,              "gz_st2[n_tracks]/F");
  _tout_truth->Branch("gx_st3",        gx_st3,              "gx_st3[n_tracks]/F");
  _tout_truth->Branch("gy_st3",        gy_st3,              "gy_st3[n_tracks]/F");
  _tout_truth->Branch("gz_st3",        gz_st3,              "gz_st3[n_tracks]/F");
  _tout_truth->Branch("gx_dp1",        gx_dp1,              "gx_dp1[n_tracks]/F");
  _tout_truth->Branch("gy_dp1",        gy_dp1,              "gy_dp1[n_tracks]/F");
  _tout_truth->Branch("gz_dp1",        gz_dp1,              "gz_dp1[n_tracks]/F");
  _tout_truth->Branch("gx_dp2",        gx_dp2,              "gx_dp2[n_tracks]/F");
  _tout_truth->Branch("gy_dp2",        gy_dp2,              "gy_dp2[n_tracks]/F");
  _tout_truth->Branch("gz_dp2",        gz_dp2,              "gz_dp2[n_tracks]/F");
  _tout_truth->Branch("gpx_st1",       gpx_st1,             "gpx_st1[n_tracks]/F");
  _tout_truth->Branch("gpy_st1",       gpy_st1,             "gpy_st1[n_tracks]/F");
  _tout_truth->Branch("gpz_st1",       gpz_st1,             "gpz_st1[n_tracks]/F");
  _tout_truth->Branch("gpx_st2",       gpx_st2,             "gpx_st2[n_tracks]/F");
  _tout_truth->Branch("gpy_st2",       gpy_st2,             "gpy_st2[n_tracks]/F");
  _tout_truth->Branch("gpz_st2",       gpz_st2,             "gpz_st2[n_tracks]/F");
  _tout_truth->Branch("gpx_st3",       gpx_st3,             "gpx_st3[n_tracks]/F");
  _tout_truth->Branch("gpy_st3",       gpy_st3,             "gpy_st3[n_tracks]/F");
  _tout_truth->Branch("gpz_st3",       gpz_st3,             "gpz_st3[n_tracks]/F");

  _tout_truth->Branch("gbarID_h1",     gbarID_h1,           "gbarID_h1[n_tracks]/I");
  _tout_truth->Branch("gbarID_h2",     gbarID_h2,           "gbarID_h2[n_tracks]/I");
  _tout_truth->Branch("gbarID_h4y",    gbarID_h4y,          "gbarID_h4y[n_tracks]/I");
  _tout_truth->Branch("gbarID_dp1",    gbarID_dp1,          "gbarID_dp1[n_tracks]/I");
  _tout_truth->Branch("gbarID_dp2",    gbarID_dp2,          "gbarID_dp2[n_tracks]/I");

  _tout_truth->Branch("gquad_h4y",     gquad_h4y,           "gquad_h4y[n_tracks]/I");
  _tout_truth->Branch("gquad_dp1",     gquad_dp1,           "gquad_dp1[n_tracks]/I");
  _tout_truth->Branch("gquad_dp2",     gquad_dp2,           "gquad_dp2[n_tracks]/I");

  _tout_truth->Branch("gnhits",        gnhits,              "gnhits[n_tracks]/I");
  _tout_truth->Branch("gndc",          gndc,                "gndc[n_tracks]/I");
  _tout_truth->Branch("gnhodo",        gnhodo,              "gnhodo[n_tracks]/I");
  _tout_truth->Branch("gnprop",        gnprop,              "gnprop[n_tracks]/I");
  _tout_truth->Branch("gnDP",          gnDP,                "gnDP[n_tracks]/F");
  _tout_truth->Branch("gnH4Y",         gnH4Y,               "gnH4Y[n_tracks]/F");
  _tout_truth->Branch("gelmid",        gelmid,              "gelmid[n_tracks][54]/I");

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

    hit_truthx[i]       = std::numeric_limits<float>::max();
    hit_truthy[i]       = std::numeric_limits<float>::max();
    hit_truthz[i]       = std::numeric_limits<float>::max();
    hit_truthpos[i]     = std::numeric_limits<float>::max();
  }

  n_tracks = 0;
  for(int i=0; i<1000; ++i) {
    gparid[i]     = std::numeric_limits<int>::max();
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

    gx_dp1[i]     = std::numeric_limits<float>::max();
    gy_dp1[i]     = std::numeric_limits<float>::max();
    gz_dp1[i]     = std::numeric_limits<float>::max();
    gx_dp2[i]     = std::numeric_limits<float>::max();
    gy_dp2[i]     = std::numeric_limits<float>::max();
    gz_dp2[i]     = std::numeric_limits<float>::max();

    gbarID_h1[i]  = std::numeric_limits<float>::max();
    gbarID_h2[i]  = std::numeric_limits<float>::max();
    gbarID_h4y[i] = std::numeric_limits<float>::max();
    gbarID_dp1[i] = std::numeric_limits<float>::max();
    gbarID_dp2[i] = std::numeric_limits<float>::max();

    gquad_dp1[i]  = std::numeric_limits<float>::max();
    gquad_dp2[i]  = std::numeric_limits<float>::max();
    gquad_h4y[i]  = std::numeric_limits<float>::max();

    gnhits[i]     = std::numeric_limits<int>::max();
    gndc[i]       = std::numeric_limits<int>::max();
    gnhodo[i]     = std::numeric_limits<int>::max();
    gnprop[i]     = std::numeric_limits<int>::max();
    gnDP[i]       = std::numeric_limits<float>::max();
    gnH4Y[i]      = std::numeric_limits<float>::max();

    for(int j=0; j<55; ++j) {
    	gelmid[i][j] = std::numeric_limits<int>::max();
    }

  }

  return 0;
}

int SimEval::GetNodes(PHCompositeNode* topNode) {

  _run_header = findNode::getClass<SQRun>(topNode, "SQRun");
  if (!_run_header) {
    LogError("!_run_header");
  }

  _event_header = findNode::getClass<SQEvent>(topNode, "SQEvent");
  if (!_event_header) {
    LogError("!_event_header");
  }

  if(_hit_container_type.find("Map") != std::string::npos) {
    _hit_map = findNode::getClass<SQHitMap>(topNode, "SQHitMap");
    if (!_hit_map) {
      LogError("!_hit_map");
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(_hit_container_type.find("Vector") != std::string::npos) {
    _hit_vector = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
    if (!_hit_vector) {
      LogError("!_hit_vector");
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  _truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth) {
    LogError("!_truth");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}







