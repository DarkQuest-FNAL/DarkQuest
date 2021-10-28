#ifndef _GEN_EMBEDDING_DATA__H_
#define _GEN_EMBEDDING_DATA__H_
#include <fun4all/SubsysReco.h>
class TFile;
class TTree;
class SQEvent;
class SQMCEvent;
class SQHitVector;
class SQTrackVector;
class SQDimuonVector;

/// An example class to analyze the simulated uDST file.
class GenEmbeddingData: public SubsysReco {
  /// Input
  SQEvent       * mi_evt;
  SQHitVector   * mi_vec_hit;
  SQMCEvent     * mi_sim_evt;
  SQTrackVector * mi_sim_vec_trk;
  SQDimuonVector* mi_sim_vec_dim;

  /// Output
  std::string m_name_file;
  std::string m_name_tree;
  TFile* mo_file;
  TTree* mo_tree;

 public:
  GenEmbeddingData(const std::string name="GenEmbeddingData");
  virtual ~GenEmbeddingData() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void SetFileName(const std::string name) { m_name_file = name; }
  void SetTreeName(const std::string name) { m_name_tree = name; }
  std::string GetFileName() const { return m_name_file; }
  std::string GetTreeName() const { return m_name_tree; }
};

#endif // _GEN_EMBEDDING_DATA__H_
