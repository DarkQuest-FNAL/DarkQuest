#ifndef _TREE_DATA__H_
#define _TREE_DATA__H_
#include <TLorentzVector.h>

struct EventData {
  int job_id;
  int event_id;
  int trig_bits;
  int rfp01;
  int rfp00;
  int rfm01;
  double weight;
  int rec_stat;
  int n_h1x;
  int n_h2x;
  int n_h3x;
  int n_h4x;
  int n_d1;
  int n_d2;
  int n_d3;

  EventData();
  virtual ~EventData() {;}

  ClassDef(EventData, 1);
};

struct TrackData {
  int            charge;
  TVector3       pos_vtx;
  TLorentzVector mom_vtx;

  TrackData();
  virtual ~TrackData() {;}

  ClassDef(TrackData, 1);
};

struct DimuonData {
  TVector3       pos;
  TLorentzVector mom;
  TLorentzVector mom_pos;
  TLorentzVector mom_neg;
  double mass;
  double pT;
  double x1;
  double x2;
  double xF;
  double costh_cs;
  double phi_cs;

  DimuonData();
  virtual ~DimuonData() {;}

  ClassDef(DimuonData, 1);
};

typedef std::vector<TrackData > TrackList;
typedef std::vector<DimuonData> DimuonList;

#endif // _TREE_DATA__H_
