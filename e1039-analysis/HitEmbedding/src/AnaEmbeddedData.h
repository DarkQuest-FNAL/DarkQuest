#ifndef _ANA_EMBEDDED_DATA__H_
#define _ANA_EMBEDDED_DATA__H_
#include <map>
#include <fun4all/SubsysReco.h>
#include "TreeData.h"
class TFile;
class TTree;
class SQEvent;
class SRecEvent;
class SQMCEvent;
class SQHitVector;
class SQTrackVector;
class SQDimuonVector;

/// A SubsysReco module to analyze the hit-embedded data.
/**
 * It is to check the basic quality of generated data.
 * The non-embedded data (i.e. the signal data before the hit embedding) can be analyzed also,
 * in order to be compared with the result of the embedded data.
 */
class AnaEmbeddedData: public SubsysReco {
  /// Input
  SQEvent       * mi_evt;
  SQEvent       * mi_evt_emb;
  SQMCEvent     * mi_evt_mc;
  SQMCEvent     * mi_evt_mc_emb;
  SQHitVector   * mi_vec_hit;
  SQTrackVector * mi_vec_trk;
  SQDimuonVector* mi_vec_dim;
  SRecEvent     * mi_srec;

  /// Output
  TFile*     mo_file;
  TTree*     mo_tree;
  EventData  mo_evt;
  TrackList  mo_trk_true;
  TrackList  mo_trk_reco;
  DimuonList mo_dim_true;
  DimuonList mo_dim_reco;

 public:
  AnaEmbeddedData(const std::string name="AnaEmbeddedData");
  virtual ~AnaEmbeddedData() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:
  void DivideHitVector(const SQHitVector* vec_in, std::map<int, SQHitVector*>& vec_sep, const bool in_time=false);
};

#endif /* _ANA_EMBEDDED_DATA__H_ */
