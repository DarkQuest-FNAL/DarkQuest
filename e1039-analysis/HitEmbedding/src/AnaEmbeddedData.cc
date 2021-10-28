#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <geom_svc/GeomSvc.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQMCEvent.h>
#include <interface_main/SQHitVector.h>
#include <interface_main/SQTrackVector.h>
#include <interface_main/SQDimuonVector.h>
#include <ktracker/SRecEvent.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <UtilAna/UtilSQHit.h>
#include <UtilAna/UtilDimuon.h>
#include "AnaEmbeddedData.h"
using namespace std;

AnaEmbeddedData::AnaEmbeddedData(const std::string name)
  : SubsysReco(name)
{
  ;
}

int AnaEmbeddedData::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEmbeddedData::InitRun(PHCompositeNode* topNode)
{
  mi_evt     = findNode::getClass<SQEvent       >(topNode, "SQEvent");
  mi_vec_hit = findNode::getClass<SQHitVector   >(topNode, "SQHitVector");
  if (!mi_evt || !mi_vec_hit) {
    cout << PHWHERE << ":  Cannot find SQEvent and/or SQHitVector." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  mi_evt_emb    = findNode::getClass<SQEvent       >(topNode, "SQEventEmb");
  mi_evt_mc     = findNode::getClass<SQMCEvent     >(topNode, "SQMCEvent");
  mi_evt_mc_emb = findNode::getClass<SQMCEvent     >(topNode, "SQMCEventEmb");
  mi_vec_trk    = findNode::getClass<SQTrackVector >(topNode, "SQTruthTrackVector");
  mi_vec_dim    = findNode::getClass<SQDimuonVector>(topNode, "SQTruthDimuonVector");
  mi_srec       = findNode::getClass<SRecEvent     >(topNode, "SRecEvent");
  if (!mi_evt_emb   ) cout << "No SQEventEmb node." << endl;
  if (!mi_evt_mc    ) cout << "No SQMCEvent node." << endl;
  if (!mi_evt_mc_emb) cout << "No SQMCEventEmb node." << endl;
  if (!mi_vec_trk   ) cout << "No SQTruthTrackVector node." << endl;
  if (!mi_vec_dim   ) cout << "No SQTruthDimuonVector node." << endl;
  if (!mi_srec      ) cout << "No SRecEvent node." << endl;

  mo_file = new TFile("ana_tree.root", "RECREATE");
  mo_tree = new TTree("tree", "Created by AnaEmbeddedData");
  mo_tree->Branch("evt"     , &mo_evt);
  mo_tree->Branch("trk_true", &mo_trk_true);
  mo_tree->Branch("trk_reco", &mo_trk_reco);
  mo_tree->Branch("dim_true", &mo_dim_true);
  mo_tree->Branch("dim_reco", &mo_dim_reco);

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEmbeddedData::process_event(PHCompositeNode* topNode)
{
  ///
  /// Event info
  ///
  mo_evt.job_id    = mi_evt->get_spill_id(); // Ppill ID in simulation means job ID
  mo_evt.event_id  = mi_evt->get_event_id();
  mo_evt.trig_bits = mi_evt->get_trigger();
  mo_evt.rfp01     = mi_evt->get_qie_rf_intensity(+1);
  mo_evt.rfp00     = mi_evt->get_qie_rf_intensity( 0);
  mo_evt.rfm01     = mi_evt->get_qie_rf_intensity(-1);

  if (mi_evt_mc) {
    mo_evt.weight = mi_evt_mc->get_weight();
  }
  if (mi_srec) {
    mo_evt.rec_stat = mi_srec->getRecStatus();
  }

  ///
  /// Hit info
  ///
  map<int, int> list_cnt; // [det ID] -> hit count
  for (SQHitVector::ConstIter it = mi_vec_hit->begin(); it != mi_vec_hit->end(); it++) {
    const SQHit* hit = *it;
    //if (! hit->is_in_time()) continue;
    int det_id = hit->get_detector_id();
    list_cnt[det_id]++;
  }

  GeomSvc* geom = GeomSvc::instance();
  mo_evt.n_h1x = list_cnt[geom->getDetectorID("H1T")] + list_cnt[geom->getDetectorID("H1B")];
  mo_evt.n_h2x = list_cnt[geom->getDetectorID("H2T")] + list_cnt[geom->getDetectorID("H2B")];
  mo_evt.n_h3x = list_cnt[geom->getDetectorID("H3T")] + list_cnt[geom->getDetectorID("H3B")];
  mo_evt.n_h4x = list_cnt[geom->getDetectorID("H4T")] + list_cnt[geom->getDetectorID("H4B")];
  mo_evt.n_d1  = list_cnt[geom->getDetectorID("D0X")];
  mo_evt.n_d2  = list_cnt[geom->getDetectorID("D2X")];
  mo_evt.n_d3  = list_cnt[geom->getDetectorID("D3pX")] + list_cnt[geom->getDetectorID("D3mX")];

  ///
  /// Truth track/dimuon info
  ///
  if (mi_vec_trk) {
    mo_trk_true.clear();
    for (unsigned int ii = 0; ii < mi_vec_trk->size(); ii++) {
      SQTrack* trk = mi_vec_trk->at(ii);
      TrackData td;
      td.charge  = trk->get_charge();
      td.mom_vtx = trk->get_mom_vtx();
      mo_trk_true.push_back(td);
    }
  }
  if (mi_vec_dim) {
    mo_dim_true.clear();
    for (unsigned int ii = 0; ii < mi_vec_dim->size(); ii++) {
      SQDimuon* dim = mi_vec_dim->at(ii);
      DimuonData dd;
      dd.mom     = dim->get_mom();
      dd.mom_pos = dim->get_mom_pos();
      dd.mom_neg = dim->get_mom_neg();
      UtilDimuon::CalcVar(dim, dd.mass, dd.pT, dd.x1, dd.x2, dd.xF);
      UtilDimuon::Lab2CollinsSoper(dim, dd.costh_cs, dd.phi_cs);
      mo_dim_true.push_back(dd);
    }
  }

  ///
  /// Reconstructed track/dimuon info
  ///
  if (mi_srec) {
    mo_trk_reco.clear();
    for (int ii = 0; ii < mi_srec->getNTracks(); ii++) {
      SRecTrack* trk = &mi_srec->getTrack(ii);
      TrackData td;
      td.charge  = trk->getCharge();
      td.mom_vtx = trk->getMomentumVertex();
      mo_trk_reco.push_back(td);
    }
    mo_dim_reco.clear();
    for (int ii = 0; ii < mi_srec->getNDimuons(); ii++) {
      SRecDimuon* dim = &mi_srec->getDimuon(ii);
      DimuonData dd;
      dd.mom     = dim->p_pos + dim->p_neg;
      dd.mom_pos = dim->p_pos;
      dd.mom_neg = dim->p_neg;
      UtilDimuon::CalcVar(dd.mom_pos, dd.mom_neg, dd.mass, dd.pT, dd.x1, dd.x2, dd.xF);
      UtilDimuon::Lab2CollinsSoper(dd.mom_pos.Vect(), dd.mom_neg.Vect(), dd.costh_cs, dd.phi_cs);
      mo_dim_reco.push_back(dd);
    }
  }

  mo_tree->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEmbeddedData::End(PHCompositeNode* topNode)
{
  mo_file->cd();
  mo_file->Write();
  mo_file->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

/// Function to divide a SQHitVector object into plane-by-plane objects.
/**
 * Usage: 
 * @code
 *   map<int, SQHitVector*> vec_hit_det;
 *   DivideHitVector(mi_vec_hit, vec_hit_det);
 *   for (map<int, SQHitVector*>::iterator it = vec_hit_det.begin(); it != vec_hit_det.end(); it++) {
 *     int det_id = it->first;
 *     SQHitVector* hv = it->second;
 *     cout << " det " << det_id << ":" << hv->size();
 *   }
 * @endcode
 */
void AnaEmbeddedData::DivideHitVector(const SQHitVector* vec_in, std::map<int, SQHitVector*>& vec_sep, const bool in_time)
{
  SQHitVector* vec0 = vec_in->Clone();
  vec0->clear();
  for (SQHitVector::ConstIter it = mi_vec_hit->begin(); it != mi_vec_hit->end(); it++) {
    const SQHit* hit = *it;
    if (in_time && ! hit->is_in_time()) continue;
    int det_id = hit->get_detector_id();
    SQHitVector* hv = vec_sep[det_id];
    if (hv == 0) vec_sep[det_id] = hv = vec0->Clone();
    hv->push_back(hit);
  }
  delete vec0;
}
