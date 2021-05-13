#ifndef _Sim_Ana__H_
#define _Sim_Ana__H_

#include <map>
#include <fun4all/SubsysReco.h>
#include <TString.h>
#include <TVector3.h>

class TFile;
class TTree;
class SQHitVector;
class SQTrackVector;

class PHG4TruthInfoContainer;
class PHG4HitContainer;
class SRecEvent;
class SRecTrack;
class PHG4Hit;
class PHG4Shower;
class PHG4Particle;

class SimAna: public SubsysReco 
{
public:
  SimAna(const std::string& name = "SimAna");
  virtual ~SimAna();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  PHG4Hit* FindG4HitAtStation(const int trk_id, const PHG4HitContainer* g4hc);
  std::vector<PHG4Hit*> FindG4HitsAtStation(const int trk_id, const PHG4HitContainer* g4hc);
  PHG4Shower* get_primary_shower(PHG4Particle* particle);
  void checkKinematics(PHG4Particle* primary);

  void set_out_name(const char* outName) {
    saveNameOut = outName;
  }
private:
  int GetNodes(PHCompositeNode* topNode);
  void MakeTree();

  SQHitVector* hitVector;
  SRecEvent* _recEvent;
  PHG4TruthInfoContainer* _truth;

  PHG4HitContainer *g4hc_d1x;
  PHG4HitContainer *g4hc_d2xp;
  PHG4HitContainer *g4hc_d3px;
  PHG4HitContainer *g4hc_d3mx;

  PHG4HitContainer *g4hc_h1t;
  PHG4HitContainer *g4hc_h1b;
  PHG4HitContainer *g4hc_h1l;
  PHG4HitContainer *g4hc_h1r;
  PHG4HitContainer *g4hc_h2t;
  PHG4HitContainer *g4hc_h2b;
  PHG4HitContainer *g4hc_h2l;
  PHG4HitContainer *g4hc_h2r;
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

  // Output
  TString saveName;
  const char* saveNameOut;
  TFile* saveFile;
  int eventID;
  TTree* saveTree;

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

  int n_tracks;
  int track_charge[100];
  int track_nhits[100];
  float track_x_target[100];
  float track_y_target[100];
  float track_z_target[100];
  float track_px_target[100];
  float track_py_target[100];
  float track_pz_target[100];
  float track_x_st1[100];
  float track_y_st1[100];
  float track_z_st1[100];
  float track_px_st1[100];
  float track_py_st1[100];
  float track_pz_st1[100];
  float track_x_vtx[100];
  float track_y_vtx[100];
  float track_z_vtx[100];
  float track_m[100];
  float track_chisq[100];
  float track_prob[100];
  float track_quality[100];
  int track_nhits_st1[100];
  int track_nhits_st2[100];
  int track_nhits_st3[100];

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
};

#endif
