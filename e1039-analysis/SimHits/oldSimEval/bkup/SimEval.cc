#include "SimEval.h"

#include <interface_main/SQHit.h>
#include <interface_main/SQHit_v1.h>
#include <interface_main/SQMCHit_v1.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQEvent_v1.h>

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
  _event_header(nullptr),
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

  // map of trackID, detID (detector ID) and element ID (e.g. bars of hodoscopes)
  std::map<int, std::map<int, int> > trkID_detid_elmid;

  // map of detID (detector ID) and
  typedef std::tuple<int, int> ParDetPair;
  std::map<ParDetPair, int> trkID_detID_ihit;
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
      hit_edep[n_hits]       = hit->get_edep();

      if(_truth) {
	int track_id = hit->get_track_id();
	int det_id = hit->get_detector_id();

	// NOTE HERE:
	// trkID is the trackID
	// It indexes by n_hits, starting from 1
	// Anything indexed by n_tracks starts at 0
	// Make sure you pay attention what is indexed where. 
	
	trkID_detID_ihit[std::make_tuple(track_id, det_id)] = ihit;
	
	auto detid_elmid_iter = trkID_detid_elmid.find(track_id);
	if(detid_elmid_iter != trkID_detid_elmid.end()) {
	  detid_elmid_iter->second.insert(std::pair<int, int>(det_id, hit->get_element_id()));
	} 
	else {
	  std::map<int, int> detid_elmid;
	  detid_elmid.insert(std::pair<int, int>(det_id, hit->get_element_id()));
	  trkID_detid_elmid[track_id] = detid_elmid;
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

  // trackID + detID -> SQHit -> PHG4Hit -> momentum
  
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
      ge[n_tracks] =  par->get_e();
      
      int trkID = par->get_track_id();
      gtrkid[n_tracks] = trkID;
      
      // The detector ID and names are listed in e1039-core/packages/geom_svc/GeomSvc.cxx.

      // ecal
      PHG4HitContainer *ECAL_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_EMCal");
      for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
	SQHit *hit = _hit_vector->at(ihit);
	if(hit and ECAL_hits) {
	  PHG4Hit* g4hit =  ECAL_hits->findHit(hit->get_g4hit_id());
	  if (g4hit) {
	    gx_ecal[n_tracks]  = g4hit->get_x(0);
	    gy_ecal[n_tracks]  = g4hit->get_y(0);
	    gz_ecal[n_tracks]  = g4hit->get_z(0);
	    gpx_ecal[n_tracks] = g4hit->get_px(0)/1000.;
	    gpy_ecal[n_tracks] = g4hit->get_py(0)/1000.;
	    gpz_ecal[n_tracks] = g4hit->get_pz(0)/1000.;
	    gedep_ecal[n_tracks] = g4hit->get_edep();
	  }
	}
      }

      // st1, 2, 3 are the drift chambers here   
      PHG4HitContainer *D1X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D0X");
      if (!D1X_hits)
	D1X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D1X");
      
      // detID 1-6 deal with D0, 7-12 D1
      for(int det_id=1; det_id<=12; ++det_id) {
	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
	if(iter != trkID_detID_ihit.end()) {
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
              gedep_st1[n_tracks] = g4hit->get_edep();
	      if(gpz_st1[n_tracks] <0){
		std::cout << "WARNING:: Negative z-momentum at Station 1! " << gpz_st1[n_tracks] <<std::endl;
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
	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
	if(iter != trkID_detID_ihit.end()) {
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
	      gedep_st2[n_tracks] = g4hit->get_edep();
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
	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
	if(iter != trkID_detID_ihit.end()) {
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
              gedep_st3[n_tracks] = g4hit->get_edep();
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
	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
	if(iter != trkID_detID_ihit.end()) {
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
              gedep_st3[n_tracks] = g4hit->get_edep();
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
	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
	if(iter != trkID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and H1_hits) {
	    PHG4Hit* g4hit =  H1_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_h1[n_tracks]  = g4hit->get_x(0);
              gy_h1[n_tracks]  = g4hit->get_y(0);
              gz_h1[n_tracks]  = g4hit->get_z(0);
              gpx_h1[n_tracks] = g4hit->get_px(0)/1000.;
              gpy_h1[n_tracks] = g4hit->get_py(0)/1000.;
              gpz_h1[n_tracks] = g4hit->get_pz(0)/1000.;
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
	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
	if(iter != trkID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and H2_hits) {
	    PHG4Hit* g4hit =  H2_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_h2[n_tracks]  = g4hit->get_x(0);
              gy_h2[n_tracks]  = g4hit->get_y(0);
              gz_h2[n_tracks]  = g4hit->get_z(0);
              gpx_h2[n_tracks] = g4hit->get_px(0)/1000.;
              gpy_h2[n_tracks] = g4hit->get_py(0)/1000.;
              gpz_h2[n_tracks] = g4hit->get_pz(0)/1000.;
	      int h2barID = hit->get_element_id();
	      if(h2barID<11) gbarID_h2[n_tracks] = hit->get_element_id();
	      if(h2barID>9) gbarID_h2[n_tracks] = hit->get_element_id() - 8;
	      break;
	    }
	  }
	}
      }

      // detID 47 to 54 are prop tubes
      // P1
      PHG4HitContainer* P1Container;
      for(int det_id = 47; det_id<=50; det_id++){
        switch(det_id){
        case 47:
          P1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1Y1");
          break;

        case 48:
          P1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1Y2");
          break;

        case 49:
          P1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1X1");
          break;

        case 50:
          P1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1X2");
          break;

	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
	if(iter != trkID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and P1Container) {
	    PHG4Hit* g4hit =  P1Container->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_p1[n_tracks] = g4hit->get_x(0);
	      gy_p1[n_tracks] = g4hit->get_y(0);
	      gz_p1[n_tracks] = g4hit->get_z(0);
	      gpx_p1[n_tracks] = g4hit->get_px(0)/1000.;
	      gpy_p1[n_tracks] = g4hit->get_py(0)/1000.;
	      gpz_p1[n_tracks] = g4hit->get_pz(0)/1000.;
              gedep_p1[n_tracks] = g4hit->get_edep();
              break;
	    }
	  }
	}
	}
      }

      // P2
      PHG4HitContainer* P2Container;
      for(int det_id = 51; det_id<=54; det_id++){
        switch(det_id){
        case 51:
          P2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2X1");
          break;

        case 52:
          P2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2X2");
          break;

        case 53:
          P2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2Y1");
          break;

        case 54:
          P2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2Y2");
          break;
        }
	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
        if(iter != trkID_detID_ihit.end()) {
          SQHit *hit = _hit_vector->at(iter->second);
          if(hit and P2Container) {
            PHG4Hit* g4hit =  P2Container->findHit(hit->get_g4hit_id());
            if (g4hit) {
              gx_p2[n_tracks] = g4hit->get_x(0);
              gy_p2[n_tracks] = g4hit->get_y(0);
              gz_p2[n_tracks] = g4hit->get_z(0);
              gpx_p2[n_tracks] = g4hit->get_px(0)/1000.;
              gpy_p2[n_tracks] = g4hit->get_py(0)/1000.;
              gpz_p2[n_tracks] = g4hit->get_pz(0)/1000.;
	      gedep_p2[n_tracks] = g4hit->get_edep();
              break;
            }
          }
        }
      }

      //detID 43 and 44 are for H4Y2L, which is the one used in DP tracking.
      PHG4HitContainer *H4Y_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4Y2L");
      for(int det_id=43; det_id<44; ++det_id) {
	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
	if(iter != trkID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and H4Y_hits) {
	    if(verbosity >= Fun4AllBase::VERBOSITY_A_LOT) {
	      LogDebug("h4y2lhit: " << iter->second);
	    }
	    PHG4Hit* g4hit =  H4Y_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_h4y2l[n_tracks]  = g4hit->get_x(0);
              gy_h4y2l[n_tracks]  = g4hit->get_y(0);
              gz_h4y2l[n_tracks]  = g4hit->get_z(0);
              gpx_h4y2l[n_tracks] = g4hit->get_px(0)/1000.;
              gpy_h4y2l[n_tracks] = g4hit->get_py(0)/1000.;
              gpz_h4y2l[n_tracks] = g4hit->get_pz(0)/1000.;
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
	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
	if(iter != trkID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and H4Y_hits) {
	    PHG4Hit* g4hit =  H4Y_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_h4y2r[n_tracks]  = g4hit->get_x(0);
              gy_h4y2r[n_tracks]  = g4hit->get_y(0);
              gz_h4y2r[n_tracks]  = g4hit->get_z(0);
              gpx_h4y2r[n_tracks] = g4hit->get_px(0)/1000.;
              gpy_h4y2r[n_tracks] = g4hit->get_py(0)/1000.;
              gpz_h4y2r[n_tracks] = g4hit->get_pz(0)/1000.;
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
	
	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
	if(iter != trkID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and DP1Container) {
	    PHG4Hit* g4hit =  DP1Container->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_dp1[n_tracks] = g4hit->get_x(0);
	      gy_dp1[n_tracks] = g4hit->get_y(0);
	      gz_dp1[n_tracks] = g4hit->get_z(0);
              gpx_dp1[n_tracks] = g4hit->get_px(0)/1000.;
              gpy_dp1[n_tracks] = g4hit->get_py(0)/1000.;
              gpz_dp1[n_tracks] = g4hit->get_pz(0)/1000.;
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
	
	auto iter = trkID_detID_ihit.find(std::make_tuple(trkID, det_id));
	if(iter != trkID_detID_ihit.end()) {
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(hit and DP2Container) {
	    PHG4Hit* g4hit =  DP2Container->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      gx_dp2[n_tracks] = g4hit->get_x(0);
	      gy_dp2[n_tracks] = g4hit->get_y(0);
	      gz_dp2[n_tracks] = g4hit->get_z(0);
              gpx_dp2[n_tracks] = g4hit->get_px(0)/1000.;
              gpy_dp2[n_tracks] = g4hit->get_py(0)/1000.;
              gpz_dp2[n_tracks] = g4hit->get_pz(0)/1000.;
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
  _tout_truth->Branch("hit_edep",      hit_edep,         "hit_edep[nHits]/F");
  _tout_truth->Branch("hit_truthx",    hit_truthx,       "hit_truthx[nHits]/F");
  _tout_truth->Branch("hit_truthy",    hit_truthy,       "hit_truthy[nHits]/F");
  _tout_truth->Branch("hit_truthz",    hit_truthz,       "hit_truthz[nHits]/F");
  _tout_truth->Branch("hit_truthpos",  hit_truthpos,     "hit_truthpos[nHits]/F");

  _tout_truth->Branch("n_tracks",      &n_tracks,           "n_tracks/I");
  _tout_truth->Branch("gtrkid",        gtrkid,              "gtrkid[n_tracks]/I");
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

  _tout_truth->Branch("gx_ecal",       gx_ecal,             "gx_ecal[n_tracks]/F");
  _tout_truth->Branch("gy_ecal",       gy_ecal,             "gy_ecal[n_tracks]/F");
  _tout_truth->Branch("gz_ecal",       gz_ecal,             "gz_ecal[n_tracks]/F");
  _tout_truth->Branch("gpx_ecal",      gpx_ecal,            "gpx_ecal[n_tracks]/F");
  _tout_truth->Branch("gpy_ecal",      gpy_ecal,            "gpy_ecal[n_tracks]/F");
  _tout_truth->Branch("gpz_ecal",      gpz_ecal,            "gpz_ecal[n_tracks]/F");
  _tout_truth->Branch("gedep_ecal",    gedep_ecal,          "gedep_ecal[n_tracks]/F");

  _tout_truth->Branch("gx_st1",        gx_st1,              "gx_st1[n_tracks]/F");
  _tout_truth->Branch("gy_st1",        gy_st1,              "gy_st1[n_tracks]/F");
  _tout_truth->Branch("gz_st1",        gz_st1,              "gz_st1[n_tracks]/F");
  _tout_truth->Branch("gpx_st1",       gpx_st1,             "gpx_st1[n_tracks]/F");
  _tout_truth->Branch("gpy_st1",       gpy_st1,             "gpy_st1[n_tracks]/F");
  _tout_truth->Branch("gpz_st1",       gpz_st1,             "gpz_st1[n_tracks]/F");

  _tout_truth->Branch("gx_st2",        gx_st2,              "gx_st2[n_tracks]/F");
  _tout_truth->Branch("gy_st2",        gy_st2,              "gy_st2[n_tracks]/F");
  _tout_truth->Branch("gz_st2",        gz_st2,              "gz_st2[n_tracks]/F");
  _tout_truth->Branch("gpx_st2",       gpx_st2,             "gpx_st2[n_tracks]/F");
  _tout_truth->Branch("gpy_st2",       gpy_st2,             "gpy_st2[n_tracks]/F");
  _tout_truth->Branch("gpz_st2",       gpz_st2,             "gpz_st2[n_tracks]/F");

  _tout_truth->Branch("gx_st3",        gx_st3,              "gx_st3[n_tracks]/F");
  _tout_truth->Branch("gy_st3",        gy_st3,              "gy_st3[n_tracks]/F");
  _tout_truth->Branch("gz_st3",        gz_st3,              "gz_st3[n_tracks]/F");
  _tout_truth->Branch("gpx_st3",       gpx_st3,             "gpx_st3[n_tracks]/F");
  _tout_truth->Branch("gpy_st3",       gpy_st3,             "gpy_st3[n_tracks]/F");
  _tout_truth->Branch("gpz_st3",       gpz_st3,             "gpz_st3[n_tracks]/F");

  _tout_truth->Branch("gx_h1",         gx_h1,               "gx_h1[n_tracks]/F");
  _tout_truth->Branch("gy_h1",         gy_h1,               "gy_h1[n_tracks]/F");
  _tout_truth->Branch("gz_h1",         gz_h1,               "gz_h1[n_tracks]/F");
  _tout_truth->Branch("gpx_h1",        gpx_h1,              "gpx_h1[n_tracks]/F");
  _tout_truth->Branch("gpy_h1",        gpy_h1,              "gpy_h1[n_tracks]/F");
  _tout_truth->Branch("gpz_h1",        gpz_h1,              "gpz_h1[n_tracks]/F");

  _tout_truth->Branch("gx_h2",         gx_h2,               "gx_h2[n_tracks]/F");
  _tout_truth->Branch("gy_h2",         gy_h2,               "gy_h2[n_tracks]/F");
  _tout_truth->Branch("gz_h2",         gz_h2,               "gz_h2[n_tracks]/F");
  _tout_truth->Branch("gpx_h2",        gpx_h2,              "gpx_h2[n_tracks]/F");
  _tout_truth->Branch("gpy_h2",        gpy_h2,              "gpy_h2[n_tracks]/F");
  _tout_truth->Branch("gpz_h2",        gpz_h2,              "gpz_h2[n_tracks]/F");

  _tout_truth->Branch("gx_p1",         gx_p1,               "gx_p1[n_tracks]/F");
  _tout_truth->Branch("gy_p1",         gy_p1,               "gy_p1[n_tracks]/F");
  _tout_truth->Branch("gz_p1",         gz_p1,               "gz_p1[n_tracks]/F");
  _tout_truth->Branch("gpx_p1",        gpx_p1,              "gpx_p1[n_tracks]/F");
  _tout_truth->Branch("gpy_p1",        gpy_p1,              "gpy_p1[n_tracks]/F");
  _tout_truth->Branch("gpz_p1",        gpz_p1,              "gpz_p1[n_tracks]/F");
  _tout_truth->Branch("gedep_p1",      gedep_p1,            "gedep_p1[n_tracks]/F");

  _tout_truth->Branch("gx_p2",         gx_p2,               "gx_p2[n_tracks]/F");
  _tout_truth->Branch("gy_p2",         gy_p2,               "gy_p2[n_tracks]/F");
  _tout_truth->Branch("gz_p2",         gz_p2,               "gz_p2[n_tracks]/F");
  _tout_truth->Branch("gpx_p2",        gpx_p2,              "gpx_p2[n_tracks]/F");
  _tout_truth->Branch("gpy_p2",        gpy_p2,              "gpy_p2[n_tracks]/F");
  _tout_truth->Branch("gpz_p2",        gpz_p2,              "gpz_p2[n_tracks]/F");
  _tout_truth->Branch("gedep_p2",      gedep_p2,            "gedep_p2[n_tracks]/F");

  _tout_truth->Branch("gx_h4y2l",      gx_h4y2l,            "gx_h4y2l[n_tracks]/F");
  _tout_truth->Branch("gy_h4y2l",      gy_h4y2l,            "gy_h4y2l[n_tracks]/F");
  _tout_truth->Branch("gz_h4y2l",      gz_h4y2l,            "gz_h4y2l[n_tracks]/F");
  _tout_truth->Branch("gpx_h4y2l",     gpx_h4y2l,           "gpx_h4y2l[n_tracks]/F");
  _tout_truth->Branch("gpy_h4y2l",     gpy_h4y2l,           "gpy_h4y2l[n_tracks]/F");
  _tout_truth->Branch("gpz_h4y2l",     gpz_h4y2l,           "gpz_h4y2l[n_tracks]/F");

  _tout_truth->Branch("gx_h4y2r",      gx_h4y2r,            "gx_h4y2r[n_tracks]/F");
  _tout_truth->Branch("gy_h4y2r",      gy_h4y2r,            "gy_h4y2r[n_tracks]/F");
  _tout_truth->Branch("gz_h4y2r",      gz_h4y2r,            "gz_h4y2r[n_tracks]/F");
  _tout_truth->Branch("gpx_h4y2r",     gpx_h4y2r,           "gpx_h4y2r[n_tracks]/F");
  _tout_truth->Branch("gpy_h4y2r",     gpy_h4y2r,           "gpy_h4y2r[n_tracks]/F");
  _tout_truth->Branch("gpz_h4y2r",     gpz_h4y2r,           "gpz_h4y2r[n_tracks]/F");

  _tout_truth->Branch("gx_dp1",        gx_dp1,              "gx_dp1[n_tracks]/F");
  _tout_truth->Branch("gy_dp1",        gy_dp1,              "gy_dp1[n_tracks]/F");
  _tout_truth->Branch("gz_dp1",        gz_dp1,              "gz_dp1[n_tracks]/F");
  _tout_truth->Branch("gpx_dp1",       gpx_dp1,             "gpx_dp1[n_tracks]/F");
  _tout_truth->Branch("gpy_dp1",       gpy_dp1,             "gpy_dp1[n_tracks]/F");
  _tout_truth->Branch("gpz_dp1",       gpz_dp1,             "gpz_dp1[n_tracks]/F");

  _tout_truth->Branch("gx_dp2",        gx_dp2,              "gx_dp2[n_tracks]/F");
  _tout_truth->Branch("gy_dp2",        gy_dp2,              "gy_dp2[n_tracks]/F");
  _tout_truth->Branch("gz_dp2",        gz_dp2,              "gz_dp2[n_tracks]/F");
  _tout_truth->Branch("gpx_dp2",       gpx_dp2,             "gpx_dp2[n_tracks]/F");
  _tout_truth->Branch("gpy_dp2",       gpy_dp2,             "gpy_dp2[n_tracks]/F");
  _tout_truth->Branch("gpz_dp2",       gpz_dp2,             "gpz_dp2[n_tracks]/F");

  _tout_truth->Branch("gbarID_h1",     gbarID_h1,           "gbarID_h1[n_tracks]/I");
  _tout_truth->Branch("gbarID_h2",     gbarID_h2,           "gbarID_h2[n_tracks]/I");
  _tout_truth->Branch("gbarID_h4y",    gbarID_h4y,          "gbarID_h4y[n_tracks]/I");
  _tout_truth->Branch("gbarID_dp1",    gbarID_dp1,          "gbarID_dp1[n_tracks]/I");
  _tout_truth->Branch("gbarID_dp2",    gbarID_dp2,          "gbarID_dp2[n_tracks]/I");

  _tout_truth->Branch("gquad_h4y",     gquad_h4y,           "gquad_h4y[n_tracks]/I");
  _tout_truth->Branch("gquad_dp1",     gquad_dp1,           "gquad_dp1[n_tracks]/I");
  _tout_truth->Branch("gquad_dp2",     gquad_dp2,           "gquad_dp2[n_tracks]/I");

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

  n_tracks = 0;
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

    gx_ecal[i]    = std::numeric_limits<float>::max();
    gy_ecal[i]    = std::numeric_limits<float>::max();
    gz_ecal[i]    = std::numeric_limits<float>::max();
    gpx_ecal[i]   = std::numeric_limits<float>::max();
    gpy_ecal[i]   = std::numeric_limits<float>::max();
    gpz_ecal[i]   = std::numeric_limits<float>::max();
    gedep_ecal[i] = std::numeric_limits<float>::max();

    gx_st1[i]     = std::numeric_limits<float>::max();
    gy_st1[i]     = std::numeric_limits<float>::max();
    gz_st1[i]     = std::numeric_limits<float>::max();
    gpx_st1[i]    = std::numeric_limits<float>::max();
    gpy_st1[i]    = std::numeric_limits<float>::max();
    gpz_st1[i]    = std::numeric_limits<float>::max();
    gedep_st1[i]  = std::numeric_limits<float>::max();

    gx_st2[i]     = std::numeric_limits<float>::max();
    gy_st2[i]     = std::numeric_limits<float>::max();
    gz_st2[i]     = std::numeric_limits<float>::max();
    gpx_st2[i]    = std::numeric_limits<float>::max();
    gpy_st2[i]    = std::numeric_limits<float>::max();
    gpz_st2[i]    = std::numeric_limits<float>::max();
    gedep_st2[i]  = std::numeric_limits<float>::max();

    gx_st3[i]     = std::numeric_limits<float>::max();
    gy_st3[i]     = std::numeric_limits<float>::max();
    gz_st3[i]     = std::numeric_limits<float>::max();
    gpx_st3[i]    = std::numeric_limits<float>::max();
    gpy_st3[i]    = std::numeric_limits<float>::max();
    gpz_st3[i]    = std::numeric_limits<float>::max();
    gedep_st3[i]  = std::numeric_limits<float>::max();

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
    gz_p1[i]      =  std::numeric_limits<float>::max();
    gpx_p1[i]     = std::numeric_limits<float>::max();
    gpy_p1[i]     = std::numeric_limits<float>::max();
    gpz_p1[i]     = std::numeric_limits<float>::max();
    gedep_p1[i]   = std::numeric_limits<float>::max();
    gx_p2[i]      = std::numeric_limits<float>::max();
    gy_p2[i]      = std::numeric_limits<float>::max();
    gz_p2[i]      = std::numeric_limits<float>::max();
    gpx_p2[i]     = std::numeric_limits<float>::max();
    gpy_p2[i]     = std::numeric_limits<float>::max();
    gpz_p2[i]     = std::numeric_limits<float>::max();
    gedep_p2[i]   = std::numeric_limits<float>::max();

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

    gbarID_h1[i]  = std::numeric_limits<float>::max();
    gbarID_h2[i]  = std::numeric_limits<float>::max();
    gbarID_h4y[i] = std::numeric_limits<float>::max();
    gbarID_dp1[i] = std::numeric_limits<float>::max();
    gbarID_dp2[i] = std::numeric_limits<float>::max();

    gquad_dp1[i]  = std::numeric_limits<float>::max();
    gquad_dp2[i]  = std::numeric_limits<float>::max();
    gquad_h4y[i]  = std::numeric_limits<float>::max();

  }

  return 0;
}

int SimEval::GetNodes(PHCompositeNode* topNode) {

  _event_header = findNode::getClass<SQEvent>(topNode, "SQEvent");
  if (!_event_header) {
    LogError("!_event_header");
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
