R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libana_embedding)
R__LOAD_LIBRARY(libktracker)
using namespace std;

int Fun4Sim(const char* fn_sig, const char* fn_emb, const int n_evt_in=0)
{
  ///
  /// Global parameters
  ///
  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;
  recoConsts *rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);

  Fun4AllServer *se = Fun4AllServer::instance();

  /// Hit embedding
  DoEmbedding* do_emb = new DoEmbedding();
  do_emb->Verbosity(10);
  do_emb->AddEmbDataFile(fn_emb);
  int n_evt_emb = do_emb->GetNumEmbEvents();
  se->registerSubsystem(do_emb);

  ///
  /// Reconstruction
  ///
  SQReco* reco = new SQReco();
  //reco->Verbosity(10);
  reco->use_geom_io_node(true);
  reco->set_evt_reducer_opt("none");
  se->registerSubsystem(reco);

  VertexFit* vertexing = new VertexFit();
  //vertexing->Verbosity(1);
  vertexing->enable_fit_target_center();
  se->registerSubsystem(vertexing);

  ///
  /// Input, output and execution
  ///
  Fun4AllInputManager* man_in = new Fun4AllDstInputManager("DSTIN");
  se->registerInputManager(man_in);

  Fun4AllDstOutputManager* man_out = new Fun4AllDstOutputManager("DSTOUT", "DST.root");
  se->registerOutputManager(man_out);

  /// Get the number of events in DST.  This function will be moved into Fun4AllInputManager.
  TFile* file_sig = new TFile(fn_sig);
  TTree* tree_sig = (TTree*)file_sig->Get("T");
  int n_evt_sig = tree_sig->GetEntries();
  file_sig->Close();

  int n_evt_both = n_evt_sig < n_evt_emb  ?  n_evt_sig  :  n_evt_emb;
  int n_evt = n_evt_in == 0 || n_evt_in > n_evt_both  ?  n_evt_both  :  n_evt_in;
  cout << "N of events:\n"
       << "  In signal    data file = " << n_evt_sig << "\n"
       << "  In embedding data file = " << n_evt_emb << "\n"
       << "  To be processed        = " << n_evt << endl;

  man_in->fileopen(fn_sig);
  se->run(n_evt);
  
  se->End();
  se->PrintTimer();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
