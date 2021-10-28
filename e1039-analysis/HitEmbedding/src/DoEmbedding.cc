#include <fstream>
#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <interface_main/SQEvent_v1.h>
#include <interface_main/SQMCEvent_v1.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQTrackVector_v1.h>
#include <interface_main/SQDimuonVector_v1.h>
#include "DoEmbedding.h"
using namespace std;

DoEmbedding::DoEmbedding(const std::string name)
  : SubsysReco(name)
  , m_overwrite_rf_info(true)
  , m_hit_id_shift(10000)
  , m_trk_id_shift(10000)
  , m_dim_id_shift(10000)
  , m_idx_emb_file(0)
  , m_idx_emb_evt (0)
  , m_emb_data_has_sim_evt(false)
  , m_emb_data_has_sim_trk(false)
  , m_emb_data_has_sim_dim(false)
  , m_file_emb(0)
  , m_tree_emb(0)
  , m_emb_sqevt    (0)
  , m_emb_sqmcevt  (0)
  , m_emb_sqvec_hit(0)
  , m_emb_sqvec_trk(0)
  , m_emb_sqvec_dim(0)
{
  ;
}

DoEmbedding::~DoEmbedding()
{
  if (m_emb_sqevt    ) delete m_emb_sqevt    ;
  if (m_emb_sqmcevt  ) delete m_emb_sqmcevt  ;
  if (m_emb_sqvec_hit) delete m_emb_sqvec_hit;
  if (m_emb_sqvec_trk) delete m_emb_sqvec_trk;
  if (m_emb_sqvec_dim) delete m_emb_sqvec_dim;
}

int DoEmbedding::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int DoEmbedding::InitRun(PHCompositeNode* topNode)
{
  mi_evt     = findNode::getClass<SQEvent       >(topNode, "SQEvent");
  mi_vec_hit = findNode::getClass<SQHitVector   >(topNode, "SQHitVector");
  if (!mi_evt || !mi_vec_hit) {
    cout << PHWHERE << ":  Cannot find SQEvent and/or SQHitVector." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  mi_sim_evt = findNode::getClass<SQMCEvent>(topNode, "SQMCEvent");
  if (! mi_sim_evt) {
    mi_sim_evt = new SQMCEvent_v1();
    dstNode->addNode(new PHIODataNode<PHObject>(mi_sim_evt, "SQMCEvent", "PHObject"));
  }
  mi_sim_vec_trk = findNode::getClass<SQTrackVector>(topNode, "SQTruthTrackVector");
  if (! mi_sim_vec_trk) {
    mi_sim_vec_trk = new SQTrackVector_v1();
    dstNode->addNode(new PHIODataNode<PHObject>(mi_sim_vec_trk, "SQTruthTrackVector", "PHObject"));
  }
  mi_sim_vec_dim = findNode::getClass<SQDimuonVector>(topNode, "SQTruthDimuonVector");
  if (! mi_sim_vec_dim) {
    mi_sim_vec_dim = new SQDimuonVector_v1();
    dstNode->addNode(new PHIODataNode<PHObject>(mi_sim_vec_dim, "SQTruthDimuonVector", "PHObject"));
  }
  mi_evt_emb = findNode::getClass<SQEvent>(topNode, "SQEventEmb");
  if (! mi_evt_emb) {
    mi_evt_emb = new SQEvent_v1();
    dstNode->addNode(new PHIODataNode<PHObject>(mi_evt_emb, "SQEventEmb", "PHObject"));
  }
  mi_sim_evt_emb = findNode::getClass<SQMCEvent>(topNode, "SQMCEventEmb");
  if (! mi_sim_evt_emb) {
    mi_sim_evt_emb = new SQMCEvent_v1();
    dstNode->addNode(new PHIODataNode<PHObject>(mi_sim_evt_emb, "SQMCEventEmb", "PHObject"));
  }

  m_emb_sqevt     = new SQEvent_v1();
  m_emb_sqmcevt   = new SQMCEvent_v1();
  m_emb_sqvec_hit = new SQHitVector_v1();
  m_emb_sqvec_trk = new SQTrackVector_v1();
  m_emb_sqvec_dim = new SQDimuonVector_v1();

  return Fun4AllReturnCodes::EVENT_OK;
}

int DoEmbedding::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() > 9) {
    cout << "DoEmbedding::process_event(): Start.\n"
         << "  run " << m_emb_sqevt->get_run_id() << ", spill " << m_emb_sqevt->get_spill_id() << ", event " << m_emb_sqevt->get_event_id() << endl;
  }

  if (! GetNextEmbEvent()) {
    // No embedding data available.  Do nothing.
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  *mi_evt_emb = *m_emb_sqevt;
  if (m_emb_data_has_sim_evt) *mi_sim_evt_emb = *m_emb_sqmcevt;

  if (m_overwrite_rf_info) {
    mi_evt->set_qie_turn_id(m_emb_sqevt->get_qie_turn_id());
    mi_evt->set_qie_rf_id  (m_emb_sqevt->get_qie_rf_id  ());
    for (int ii = -16; ii <= 16; ii++) {
      mi_evt->set_qie_rf_intensity(ii, m_emb_sqevt->get_qie_rf_intensity(ii));
    }
  }

  if (Verbosity() > 9) cout << "  N of hits to be embedded: " << m_emb_sqvec_hit->size() << endl;
  for (SQHitVector::Iter it = m_emb_sqvec_hit->begin(); it != m_emb_sqvec_hit->end(); it++) {
    SQHit* hit_emb = *it;
    hit_emb->set_hit_id  (hit_emb->get_hit_id  () + m_hit_id_shift);
    hit_emb->set_track_id(hit_emb->get_track_id() + m_trk_id_shift);
    mi_vec_hit->push_back(hit_emb);
  }

  if (m_emb_data_has_sim_trk) {
    if (Verbosity() > 9) cout << "  N of tracks to be embedded: " << m_emb_sqvec_trk->size() << endl;
    for (SQTrackVector::Iter it = m_emb_sqvec_trk->begin(); it != m_emb_sqvec_trk->end(); it++) {
      SQTrack* trk_emb = *it;
      trk_emb->set_track_id(trk_emb->get_track_id() + m_trk_id_shift);
      mi_sim_vec_trk->push_back(trk_emb);
    }
  }

  if (m_emb_data_has_sim_dim) {
    if (Verbosity() > 9) cout << "  N of dimuons to be embedded: " << m_emb_sqvec_dim->size() << endl;
    for (SQDimuonVector::Iter it = m_emb_sqvec_dim->begin(); it != m_emb_sqvec_dim->end(); it++) {
      SQDimuon* dim_emb = *it;
      dim_emb->set_dimuon_id   (dim_emb->get_dimuon_id   () + m_dim_id_shift);
      dim_emb->set_track_id_pos(dim_emb->get_track_id_pos() + m_trk_id_shift);
      dim_emb->set_track_id_neg(dim_emb->get_track_id_neg() + m_trk_id_shift);
      mi_sim_vec_dim->push_back(dim_emb);
    }
  }

  if (Verbosity() > 9) cout << "DoEmbedding::process_event(): End." << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int DoEmbedding::End(PHCompositeNode* topNode)
{
  CloseEmbDataFile();
  return Fun4AllReturnCodes::EVENT_OK;
}

void DoEmbedding::AddEmbDataFile(const char* fn_root)
{
  m_list_emb_file.push_back(fn_root);
}

void DoEmbedding::AddEmbDataFiles(const char* fn_list)
{
  ifstream ifs(fn_list);
  string fn_root;
  while (ifs >> fn_root) AddEmbDataFile(fn_root.c_str());
  ifs.close();
}

int DoEmbedding::GetNumEmbEvents()
{
  int num = 0;
  for (vector<string>::iterator it = m_list_emb_file.begin(); it != m_list_emb_file.end(); it++) {
    OpenEmbDataFile(it->c_str());
    num += m_tree_emb->GetEntries();
    CloseEmbDataFile();
  }
  return num;
}

void DoEmbedding::OpenEmbDataFile(const char* fn_root)
{
  m_file_emb = new TFile(fn_root);
  if (! m_file_emb->IsOpen()) {
    cout << "ERROR:  Cannot open the embedding-data file, " << fn_root << ".  Abort." << endl;
    exit(1);
  }
  m_tree_emb = (TTree*)m_file_emb->Get("tree");
  if (! m_tree_emb) {
    cout << "ERROR:  Cannot get the embedding-data tree.  Abort." << endl;
    exit(1);
  }
  if (Verbosity() > 0) {
    cout << "DoEmbedding::OpenEmbDataFile(): " << fn_root << ", N = " << m_tree_emb->GetEntries() << "." << endl;
  }
}

void DoEmbedding::CloseEmbDataFile()
{
  if (! m_file_emb) return;
  m_file_emb->Close();
  m_file_emb = 0;
  m_tree_emb = 0;
}

bool DoEmbedding::GetNextEmbEvent()
{
  if (! m_file_emb) { // Need to open an embedding-data file
    if (m_idx_emb_file >= m_list_emb_file.size()) return false; // No file available.
    string fn_root = m_list_emb_file[m_idx_emb_file];
    OpenEmbDataFile(fn_root.c_str());
    m_tree_emb->SetBranchAddress("SQEvent"    , &m_emb_sqevt    );
    m_tree_emb->SetBranchAddress("SQHitVector", &m_emb_sqvec_hit);
    m_emb_data_has_sim_evt = (m_tree_emb->FindBranch("SQMCEvent"     ) != 0);
    m_emb_data_has_sim_trk = (m_tree_emb->FindBranch("SQTruthTrackVector" ) != 0);
    m_emb_data_has_sim_dim = (m_tree_emb->FindBranch("SQTruthDimuonVector") != 0);
    if (m_emb_data_has_sim_evt) m_tree_emb->SetBranchAddress("SQMCEvent"          , &m_emb_sqmcevt  );
    if (m_emb_data_has_sim_trk) m_tree_emb->SetBranchAddress("SQTruthTrackVector" , &m_emb_sqvec_trk);
    if (m_emb_data_has_sim_dim) m_tree_emb->SetBranchAddress("SQTruthDimuonVector", &m_emb_sqvec_dim);

    m_idx_emb_file++;
    m_idx_emb_evt = 0;
  }

  if (m_idx_emb_evt >= m_tree_emb->GetEntries()) { // No event available.
    CloseEmbDataFile();
    return GetNextEmbEvent();
  }

  m_tree_emb->GetEntry(m_idx_emb_evt++);

  return true;
}
