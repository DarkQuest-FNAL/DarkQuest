#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <phool/getClass.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQMCEvent.h>
#include <interface_main/SQHitVector.h>
#include <interface_main/SQTrackVector.h>
#include <interface_main/SQDimuonVector.h>
#include "GenEmbeddingData.h"
using namespace std;

GenEmbeddingData::GenEmbeddingData(const std::string name)
  : SubsysReco(name)
  , m_name_file("embedding_data.root")
  , m_name_tree("tree")
{
  ;
}

int GenEmbeddingData::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int GenEmbeddingData::InitRun(PHCompositeNode* topNode)
{
  mi_evt     = findNode::getClass<SQEvent       >(topNode, "SQEvent");
  mi_vec_hit = findNode::getClass<SQHitVector   >(topNode, "SQHitVector");
  if (!mi_evt || !mi_vec_hit) {
    cout << PHWHERE << ":  Cannot find SQEvent and/or SQHitVector." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  mi_sim_evt     = findNode::getClass<SQMCEvent     >(topNode, "SQMCEvent");
  mi_sim_vec_trk = findNode::getClass<SQTrackVector >(topNode, "SQTruthTrackVector");
  mi_sim_vec_dim = findNode::getClass<SQDimuonVector>(topNode, "SQTruthDimuonVector");

  mo_file = new TFile(m_name_file.c_str(), "RECREATE");
  mo_tree = new TTree(m_name_tree.c_str(), "Created by GenEmbeddingData");
  mo_tree->Branch("SQEvent"    , &mi_evt);
  mo_tree->Branch("SQHitVector", &mi_vec_hit);
  if (mi_sim_evt    ) mo_tree->Branch("SQMCEvent"          , &mi_sim_evt);
  if (mi_sim_vec_trk) mo_tree->Branch("SQTruthTrackVector" , &mi_sim_vec_trk);
  if (mi_sim_vec_dim) mo_tree->Branch("SQTruthDimuonVector", &mi_sim_vec_dim);

  return Fun4AllReturnCodes::EVENT_OK;
}

int GenEmbeddingData::process_event(PHCompositeNode* topNode)
{
  mo_tree->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

int GenEmbeddingData::End(PHCompositeNode* topNode)
{
  mo_file->cd();
  mo_file->Write();
  mo_file->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}
