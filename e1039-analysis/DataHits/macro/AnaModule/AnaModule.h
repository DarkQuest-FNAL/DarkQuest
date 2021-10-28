#ifndef _ANA_Module__H_
#define _ANA_Module__H_

#include <map>
#include <fun4all/SubsysReco.h>
#include <TString.h>
#include <TVector3.h>
#include <ktracker/SRawEvent.h>

class TFile;
class TTree;
class SQHitVector;
class SRecTrack;
class SRecEvent;

class AnaModule: public SubsysReco 
{
public:
  AnaModule() {;}
  virtual ~AnaModule() {;}

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  void set_output_filename(const TString& n) { saveName = n; }
  void set_reco(bool b){ saveReco = b;}

private:
  int GetNodes(PHCompositeNode* topNode);
  int ResetEvalVars();
  void MakeTree();

  // Input
  SRawEvent* rawEvent;
  SRecEvent* recEvent;

  // Output
  TString saveName;
  TFile* saveFile;
  TTree* saveTree;

  Int_t runID, spillID, eventID;
  Int_t triggerBits;
  Int_t nHits;
  Int_t detectorID[15000], elementID[15000];
  Double_t tdcTime[15000], driftDistance[15000], pos[15000];

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
  float track_px_vtx[100];
  float track_py_vtx[100];
  float track_pz_vtx[100];
  float track_m[100];
  float track_chisq[100];
  float track_prob[100];
  float track_quality[100];
  int track_nhits_st1[100];
  int track_nhits_st2[100];
  int track_nhits_st3[100];

  int n_dimuons;
  float dimuon_mass[100];
  float dimuon_chisq[100];
  float dimuon_x_vtx[100];
  float dimuon_y_vtx[100];
  float dimuon_z_vtx[100];
  float dimuon_px[100];
  float dimuon_py[100];
  float dimuon_pz[100];
  float dimuon_pmom_x[100];
  float dimuon_pmom_y[100];
  float dimuon_pmom_z[100];
  float dimuon_nmom_x[100];
  float dimuon_nmom_y[100];
  float dimuon_nmom_z[100];
  float dimuon_ppos_x[100];
  float dimuon_ppos_y[100];
  float dimuon_ppos_z[100];
  float dimuon_npos_x[100];
  float dimuon_npos_y[100];
  float dimuon_npos_z[100];

  bool saveReco=false;
};

#endif
