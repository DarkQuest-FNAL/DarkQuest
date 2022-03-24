R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libdptrigger)
R__LOAD_LIBRARY(libevt_filter)
R__LOAD_LIBRARY(libsim_ana)
using namespace std;

int RunEmbedding(
    const char* fn_sig = "/seaquest/users/yfeng/DarkQuest/DarkQuest/e1039-analysis/SimHits/macro/output_DST.root",
    const char* fn_emb = "/pnfs/e1039/persistent/users/kenichi/data_emb_e906/0001/embedding_data.root", 
    const int n_evt_in=100)
{
  ///
  /// Global parameters
  ///
  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;
  recoConsts *rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  rc->set_CharFlag("fMagFile",
                   "$E1039_RESOURCE/geometry/magnetic_fields/tab.Fmag");
  rc->set_CharFlag("kMagFile",
                   "$E1039_RESOURCE/geometry/magnetic_fields/tab.Kmag");

  rc->set_BoolFlag("TRACK_DISPLACED",true);
  rc->set_BoolFlag("OLD_TRACKING",false);
  rc->set_IntFlag("MaxHitsDC0", 300);
  rc->set_IntFlag("MaxHitsDC1", 300);
  rc->set_IntFlag("MaxHitsDC2", 300);
  rc->set_IntFlag("MaxHitsDC3p", 300);
  rc->set_IntFlag("MaxHitsDC3m", 300);

  Fun4AllServer *se = Fun4AllServer::instance();

  /// Hit embedding
  //gSystem->Load("libsim_ana.so");
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

  
  // truth node maker after tracking
  TruthNodeMaker* truthMaker = new TruthNodeMaker();
  truthMaker->set_legacy_rec_container(1);
  truthMaker->Verbosity(0);
  se->registerSubsystem(truthMaker);
  

  DPTriggerAnalyzer* dptrigger = new DPTriggerAnalyzer();
  dptrigger->set_road_set_file_name("$E1039_RESOURCE/trigger/trigger_67.txt");
  dptrigger->Verbosity(0);
  se->registerSubsystem(dptrigger);

  VertexFit* vertexing = new VertexFit();
  //vertexing->Verbosity(1);
  //vertexing->enable_fit_target_center();
  se->registerSubsystem(vertexing);

  SimAna *sim_ana = new SimAna();
  std::string ofile = "test.root";
  sim_ana->set_out_name(ofile);
  sim_ana->set_legacy_rec_container(true);
  sim_ana->save_secondaries(false);
  se->registerSubsystem(sim_ana);

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
