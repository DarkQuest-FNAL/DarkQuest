#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include "AnaCleanAndMessyData.h"
using namespace std;

AnaCleanAndMessyData::AnaCleanAndMessyData()
  : m_verb(0)
  , m_cl_file(0)
  , m_cl_tree(0)
  , m_me_file(0)
  , m_me_tree(0)
{
  ;
}

AnaCleanAndMessyData::~AnaCleanAndMessyData()
{
  ;
}

void AnaCleanAndMessyData::Init(const char* fn_clean, const char* fn_messy)
{
  m_cl_file = new TFile(fn_clean);
  m_cl_tree = (TTree*)m_cl_file->Get("tree");
  m_cl_evt      = new EventData ();
  m_cl_trk_true = new TrackList ();
  m_cl_trk_reco = new TrackList ();
  m_cl_dim_true = new DimuonList();
  m_cl_dim_reco = new DimuonList();
  m_cl_tree->SetBranchAddress("evt"     , &m_cl_evt);
  m_cl_tree->SetBranchAddress("trk_true", &m_cl_trk_true);
  m_cl_tree->SetBranchAddress("trk_reco", &m_cl_trk_reco);
  m_cl_tree->SetBranchAddress("dim_true", &m_cl_dim_true);
  m_cl_tree->SetBranchAddress("dim_reco", &m_cl_dim_reco);

  m_me_file = new TFile(fn_messy);
  m_me_tree = (TTree*)m_me_file->Get("tree");
  m_me_evt      = new EventData ();
  m_me_trk_true = new TrackList ();
  m_me_trk_reco = new TrackList ();
  m_me_dim_true = new DimuonList();
  m_me_dim_reco = new DimuonList();
  m_me_tree->SetBranchAddress("evt"     , &m_me_evt);
  m_me_tree->SetBranchAddress("trk_true", &m_me_trk_true);
  m_me_tree->SetBranchAddress("trk_reco", &m_me_trk_reco);
  m_me_tree->SetBranchAddress("dim_true", &m_me_dim_true);
  m_me_tree->SetBranchAddress("dim_reco", &m_me_dim_reco);

  gSystem->mkdir("result", true);
  m_out_file = new TFile("result/output.root", "RECREATE");
  m_h1_trk_pos_cl = new TH1D("h1_trk_pos_cl", ";RF+00;", 20, 0, 1000);
  m_h1_trk_pos_me = new TH1D("h1_trk_pos_me", ";RF+00;", 20, 0, 1000);
  m_h1_trk_neg_cl = new TH1D("h1_trk_neg_cl", ";RF+00;", 20, 0, 1000);
  m_h1_trk_neg_me = new TH1D("h1_trk_neg_me", ";RF+00;", 20, 0, 1000);
  m_h1_dim_cl     = new TH1D("h1_dim_cl"    , ";RF+00;", 20, 0, 1000);
  m_h1_dim_me     = new TH1D("h1_dim_me"    , ";RF+00;", 20, 0, 1000);
}

/// Function to analyze a pair of non-embedded and embedded (i.e. clean and messy) data.
/**
 * The main part of this function is to match the event between the two data.
 */
void AnaCleanAndMessyData::Analyze()
{
  int n_cl_evt = m_cl_tree->GetEntries();
  int n_me_evt = m_me_tree->GetEntries();
  int i_cl_evt = 0;
  int i_me_evt = 0;

  bool no_event = false;
  while (! no_event) {
    if (i_cl_evt >= n_cl_evt || i_me_evt >= n_me_evt) return;
    m_cl_tree->GetEntry(i_cl_evt);
    m_me_tree->GetEntry(i_me_evt);
    pair<int, int> job_evt_cl(m_cl_evt->job_id, m_cl_evt->event_id);
    pair<int, int> job_evt_me(m_me_evt->job_id, m_me_evt->event_id);

    while (job_evt_cl != job_evt_me) { // job+event IDs are different
      if (job_evt_cl < job_evt_me) {
        i_cl_evt++;
        if (i_cl_evt >= n_cl_evt) return;
        m_cl_tree->GetEntry(i_cl_evt);
        job_evt_cl = pair<int, int>(m_cl_evt->job_id, m_cl_evt->event_id);
      } else { // >
        i_me_evt++;
        if (i_me_evt >= n_me_evt) return;
        m_me_tree->GetEntry(i_me_evt);
        job_evt_me = pair<int, int>(m_me_evt->job_id, m_me_evt->event_id);
      }
    }

    if (Verbosity() > 9) {
      cout << "AnaCleanAndMessyData::Analyze():  Job ID " << m_cl_evt->job_id << ", Event ID " << m_cl_evt->event_id << ": Clean " << i_cl_evt << "/" << n_cl_evt << ", Messy " << i_me_evt << "/" << n_me_evt << endl;
    }
    AnalyzeEvent();
    i_cl_evt++;
    i_me_evt++;
  }
}

void AnaCleanAndMessyData::End()
{
  if (m_cl_file) m_cl_file->Close();
  if (m_me_file) m_me_file->Close();
  if (m_out_file) DrawAndWriteOutput();
}

/// Function to analyze one event.
/**
 * This function should be modified as your analysis needs.
 */
void AnaCleanAndMessyData::AnalyzeEvent()
{
  double ww = m_cl_evt->weight;
  int rfp01 = m_me_evt->rfp01;
  int rfp00 = m_me_evt->rfp00;
  int rfm01 = m_me_evt->rfm01;
  int n_h1x = m_me_evt->n_h1x;
  int n_h2x = m_me_evt->n_h2x;
  int n_h3x = m_me_evt->n_h3x;
  int n_h4x = m_me_evt->n_h4x;

  int n_trk_pos_cl = 0;
  int n_trk_neg_cl = 0;
  for (int ii = 0; ii < m_cl_trk_reco->size(); ii++) {
    TrackData* td = &m_cl_trk_reco->at(ii);
    if (td->charge > 0) n_trk_pos_cl++;
    else                n_trk_neg_cl++;
  }

  int n_trk_pos_me = 0;
  int n_trk_neg_me = 0;
  for (int ii = 0; ii < m_me_trk_reco->size(); ii++) {
    TrackData* td = &m_me_trk_reco->at(ii);
    if (td->charge > 0) n_trk_pos_me++;
    else                n_trk_neg_me++;
  }

  int n_dim_cl = m_cl_dim_reco->size();
  int n_dim_me = m_me_dim_reco->size();

  m_h1_trk_pos_cl->Fill(rfp00, n_trk_pos_cl);
  m_h1_trk_pos_me->Fill(rfp00, n_trk_pos_me);
  m_h1_trk_neg_cl->Fill(rfp00, n_trk_neg_cl);
  m_h1_trk_neg_me->Fill(rfp00, n_trk_neg_me);
  m_h1_dim_cl    ->Fill(rfp00, n_dim_cl    );
  m_h1_dim_me    ->Fill(rfp00, n_dim_me    );
}

/// Function to be called in End() to make, draw and write output objects.
void AnaCleanAndMessyData::DrawAndWriteOutput()
{
  m_out_file->cd();
  m_out_file->Write();
  TEfficiency* teff_trk_pos = new TEfficiency(*m_h1_trk_pos_me, *m_h1_trk_pos_cl);
  TEfficiency* teff_trk_neg = new TEfficiency(*m_h1_trk_neg_me, *m_h1_trk_neg_cl);
  TEfficiency* teff_dim     = new TEfficiency(*m_h1_dim_me    , *m_h1_dim_cl    );
  teff_trk_pos->SetName("teff_trk_pos");
  teff_trk_neg->SetName("teff_trk_neg");
  teff_dim    ->SetName("teff_dim");
  teff_trk_pos->SetTitle(";RF+00;Reco. efficiency of #mu^{#plus}");
  teff_trk_neg->SetTitle(";RF+00;Reco. efficiency of #mu^{#minus}");
  teff_dim    ->SetTitle(";RF+00;Reco. efficiency of dimuon");

  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  teff_trk_pos->Draw();
  c1->SaveAs("result/teff_trk_pos.png");
  teff_trk_neg->Draw();
  c1->SaveAs("result/teff_trk_neg.png");
  teff_dim    ->Draw();
  c1->SaveAs("result/teff_dim.png");
  delete c1;

  teff_trk_pos->Write();
  teff_trk_neg->Write();
  teff_dim    ->Write();
  m_out_file->Close();
}
