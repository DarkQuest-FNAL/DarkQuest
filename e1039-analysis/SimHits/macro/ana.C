#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"

using namespace std;

double ratio_error(
    const double a,
    const double b,
    const double ea,
    const double eb
    ) {
  double r = a/b;
  double er = r*sqrt(
      (ea/a)*(ea/a) + (eb/b)*(eb/b)
      );
  return er;
}

double binom_error(
    const double a,
    const double b
    ){
  double r=a/b;
  return sqrt(r*(1-r)/b);
}

TH1D * getEffHist(
    const char* hname,
    const TH1D* h1,
    const TH1D* h2
    ) {
  TH1D *h = (TH1D*)h1->Clone(hname);
  int nbin = h->GetNbinsX();
  for(int ibin=1; ibin<=nbin; ++ibin) {
    double a = h1->GetBinContent(ibin);
    double b = h2->GetBinContent(ibin);
    double r = 0;
    double e = 0;
    if(b>0) {
      r = a/b;
      e = binom_error(a, b);
    }
    h->SetBinContent(ibin, r);
    h->SetBinError(ibin, e);
  }

  return h;
}

void ana_eval(
    const char* var_name = "dimu_gphi",
    const char* title = ";dimu gphi;",
    const int nbin = 16,
    const double xmin = -3.14,
    const double xmax = 3.14
    ) {
  gStyle->SetOptFit();

  TFile *f = TFile::Open("test-trig-emu.root","read");
  TTree *T = (TTree*) f->Get("Truth");

  TH2D *h2_trig_hodo = new TH2D("h2_trig_hodo",
      ";ntracks with nhodo>7; emu-trigger",
      3, -0.5, 2.5,
      11, -0.5, 10.5
      );

  TH1D *h1_eval_all = new TH1D("h1_eval_all",title, nbin, xmin, xmax);
  TH1D *h1_eval_geo = new TH1D("h1_eval_geo",title, nbin, xmin, xmax);
  TH1D *h1_eval_tri = new TH1D("h1_eval_tri",title, nbin, xmin, xmax);

  unsigned short emu_trigger = 0;
  int n_tracks = 0;
  int gnhodo[100];
  float dimu_gphi[100];
  float eval_var[100];

  T->SetBranchAddress("emu_trigger", &emu_trigger);
  T->SetBranchAddress("n_tracks",    &n_tracks);
  T->SetBranchAddress("gnhodo",      &gnhodo);
  T->SetBranchAddress("dimu_gphi",   &dimu_gphi);
  T->SetBranchAddress(var_name,      &eval_var);


  for(int ientry=0; ientry<T->GetEntries(); ++ientry) {
    T->GetEntry(ientry);
    //cout << emu_trigger << ", " << gnhodo[0] << endl;

    int ntrack_hodo = 0;
    for (int i=0; i<n_tracks; ++i) {
      if(gnhodo[i]>7) ++ntrack_hodo;
    }

    for (int i=0; i<10; ++i) {
      if((emu_trigger&(1<<i)) > 0) {
        h2_trig_hodo->Fill(ntrack_hodo,i);
      }
    }

    h1_eval_all->Fill(eval_var[0]);
    if(ntrack_hodo>1)
      h1_eval_geo->Fill(eval_var[0]);
    if(  (emu_trigger&(1<<5))>0
      || (emu_trigger&(1<<6))>0
      )
      h1_eval_tri->Fill(eval_var[0]);
  }

  TCanvas *c1 = new TCanvas("c1","c1"); c1->SetGrid();
  h2_trig_hodo->SetStats(0);
  h2_trig_hodo->Draw("colz,text");

  TCanvas *c2 = new TCanvas("c2","c2"); c2->SetGrid();
  h1_eval_all->SetStats(0);
  h1_eval_all->SetLineColor(kBlack);
  h1_eval_all->Draw("e");
  //h1_eval_all->SetMinimum(0.1);
  //h1_eval_geo->SetLineColor(kBlue);
  //h1_eval_geo->Draw("same");
  //h1_eval_tri->SetLineColor(kRed);
  //h1_eval_tri->Draw("same");

  TF1 *fcos = new TF1("fcos","[0]*([1]*cos(2*x)+1)");

  TH1D *h1_eval_geo_eff = getEffHist("h1_eval_geo_eff",
      h1_eval_geo,
      h1_eval_all
      );
  TCanvas *c3 = new TCanvas("c3","c3"); c3->SetGrid();
  h1_eval_geo_eff->Draw("e");
  if(var_name=="dimu_gphi")
    h1_eval_geo_eff->Fit(fcos);
  else
    h1_eval_geo_eff->SetStats(0);

  TH1D *h1_eval_tri_eff = getEffHist("h1_eval_tri_eff",
      h1_eval_tri,
      h1_eval_geo
      );
  TCanvas *c4 = new TCanvas("c4","c4"); c4->SetGrid();
  h1_eval_tri_eff->Draw("e");
  if(var_name=="dimu_gphi")
    h1_eval_tri_eff->Fit(fcos);
  else
    h1_eval_tri_eff->SetStats(0);

}





void ana(){
  //ana_eval();
  //ana_eval("dimu_gmass", ";dimu gmass;", 8, 2, 10);
  ana_eval("dimu_gpt", ";dimu gpt;", 8, 0, 4);
}







