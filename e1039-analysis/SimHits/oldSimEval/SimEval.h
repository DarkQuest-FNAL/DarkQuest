#ifndef _H_SimEval_H_
#define _H_SimEval_H_

// ROOT
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

// Fun4All includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <vector>
#include <string>
#include <iostream>
#include <list>
#include <map>
//#include <algorithm>

class SQRun;
class SQSpillMap;

class SQEvent;
class SQHitMap;
class SQHitVector;

class PHG4TruthInfoContainer;
class PHG4HitContainer;
class PHG4Hit;
class PHG4Shower;
class PHG4Particle;

class SRecEvent;

class GeomSvc;

class TFile;
class TTree;

class SimEval: public SubsysReco {

public:

  SimEval(const std::string &name = "SimEval");
  virtual ~SimEval() {
  }
  
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  /// what was the primary shower that is associated with this particle?                                                                                                         
  PHG4Shower* get_primary_shower(PHG4Particle* particle);
  // what is the G4Hit at this station for this particle
  PHG4Hit* FindG4HitAtStation(const int trk_id, const PHG4HitContainer* g4hc);
  // what are the  G4Hit at this station for this particle
  std::vector<PHG4Hit*> FindG4HitsAtStation(const int trk_id, const PHG4HitContainer* g4hc);

  int InitEvalTree();
  int ResetEvalVars();
  
  const std::string& get_out_name() const {
    return _out_name;
  }
  
  void set_out_name(const std::string& outName) {
    _out_name = outName;
  }
  
private:

  int GetNodes(PHCompositeNode *topNode);
  
  int TruthEval(PHCompositeNode *topNode);
  
  //std::string _hit_container_type;

  size_t _event;
  
  SQRun* _run_header;
  SQEvent * _event_header;
  SQHitMap *_hit_map;
  SQHitVector *_hit_vector;

  PHG4TruthInfoContainer* _truth;

  PHG4HitContainer *g4hc_d1x;
  PHG4HitContainer *g4hc_d2xp;
  PHG4HitContainer *g4hc_d3px;
  PHG4HitContainer *g4hc_d3mx;

  PHG4HitContainer *g4hc_h1t;
  PHG4HitContainer *g4hc_h1b;
  PHG4HitContainer *g4hc_h2t;
  PHG4HitContainer *g4hc_h2b;
  PHG4HitContainer *g4hc_h3t;
  PHG4HitContainer *g4hc_h3b;
  PHG4HitContainer *g4hc_h4t;
  PHG4HitContainer *g4hc_h4b;
  PHG4HitContainer *g4hc_h4y2l;
  PHG4HitContainer *g4hc_h4y2r;

  PHG4HitContainer *g4hc_p1y1;
  PHG4HitContainer *g4hc_p1y2;
  PHG4HitContainer *g4hc_p1x1;
  PHG4HitContainer *g4hc_p1x2;
  PHG4HitContainer *g4hc_p2x1;
  PHG4HitContainer *g4hc_p2x2;
  PHG4HitContainer *g4hc_p2y1;
  PHG4HitContainer *g4hc_p2y2;
  
  PHG4HitContainer *g4hc_ecal;

  std::string _out_name;
  TTree* _tout_truth;
  
  int run_id;
  int event_id;
  
  // These are SQHits, info per hit
  int n_hits;
  int hit_detid[10000];
  int hit_elmid[10000];
  float hit_driftdis[10000];
  float hit_pos[10000];
  float hit_detz[10000];
  float hit_edep[10000];
  float hit_truthx[10000];
  float hit_truthy[10000];
  float hit_truthz[10000];
  float hit_truthpos[10000];
  
  int n_showers;
  float sx_ecal[1000];
  float sy_ecal[1000];
  float sz_ecal[1000];
  float sedep_ecal[1000];

  int n_primaries;
  int gtrkid[1000];
  int gpid[1000];
  float gvx[1000];
  float gvy[1000];
  float gvz[1000];
  float gpx[1000];
  float gpy[1000];
  float gpz[1000];
  float gpt[1000];
  float geta[1000];
  float gphi[1000];
  float ge[1000];

  int nhits_ecal[1000];
  float gx_ecal[1000][100];
  float gy_ecal[1000][100];
  float gz_ecal[1000][100];
  float gpx_ecal[1000][100];
  float gpy_ecal[1000][100];
  float gpz_ecal[1000][100];
  float gedep_ecal[1000][100];

  float gx_st1[1000];
  float gy_st1[1000];
  float gz_st1[1000];
  float gpx_st1[1000];
  float gpy_st1[1000];
  float gpz_st1[1000];
  float gx_st2[1000];
  float gy_st2[1000];
  float gz_st2[1000];
  float gpx_st2[1000];
  float gpy_st2[1000];
  float gpz_st2[1000];
  float gx_st3[1000];
  float gy_st3[1000];
  float gz_st3[1000];
  float gpx_st3[1000];
  float gpy_st3[1000];
  float gpz_st3[1000];

  float gx_dp1[1000];
  float gy_dp1[1000];
  float gz_dp1[1000];
  float gpx_dp1[1000];
  float gpy_dp1[1000];
  float gpz_dp1[1000];
  float gx_dp2[1000];
  float gy_dp2[1000];
  float gz_dp2[1000];
  float gpx_dp2[1000];
  float gpy_dp2[1000];
  float gpz_dp2[1000];

  float gx_h1[1000];
  float gy_h1[1000];
  float gz_h1[1000];
  float gpx_h1[1000];
  float gpy_h1[1000];
  float gpz_h1[1000];
  float gx_h2[1000];
  float gy_h2[1000];
  float gz_h2[1000];
  float gpx_h2[1000];
  float gpy_h2[1000];
  float gpz_h2[1000];
  float gx_h3[1000];
  float gy_h3[1000];
  float gz_h3[1000];
  float gpx_h3[1000];
  float gpy_h3[1000];
  float gpz_h3[1000];
  float gx_h4[1000];
  float gy_h4[1000];
  float gz_h4[1000];
  float gpx_h4[1000];
  float gpy_h4[1000];
  float gpz_h4[1000];

  float gx_p1[1000];
  float gy_p1[1000];
  float gz_p1[1000];
  float gpx_p1[1000];
  float gpy_p1[1000];
  float gpz_p1[1000];

  float gx_h4y2l[1000];
  float gy_h4y2l[1000];
  float gz_h4y2l[1000];
  float gpx_h4y2l[1000];
  float gpy_h4y2l[1000];
  float gpz_h4y2l[1000];
  float gx_h4y2r[1000];
  float gy_h4y2r[1000];
  float gz_h4y2r[1000];
  float gpx_h4y2r[1000];
  float gpy_h4y2r[1000];
  float gpz_h4y2r[1000];

  float gbarID_h1[1000];
  float gbarID_h2[1000];
  float gbarID_h4y[1000];
  float gbarID_dp1[1000];
  float gbarID_dp2[1000];

  float gquad_dp1[1000];
  float gquad_dp2[1000];
  float gquad_h4y[1000];
  
  GeomSvc *p_geomSvc;
};


#endif /* _H_SimEval_H_ */
