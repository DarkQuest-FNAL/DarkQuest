#include <TSystem.h>

R__LOAD_LIBRARY(libinterface_main)
R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libanamodule)

/*
This macro takes:
1. geom.root is the standard geometry dump from running the Fun4Sim macro;
2. e906_run7.opts is provided
3. /seaquest/users/cmantill/kTracker/kTrackerDark/nimSkim_data.root
3. /seaquest/users/cmantill/DarkQuest/e1039-analysis/DataHits/macro/nimSkim_data.root (2017 data with DP hits)
*/

int RecoE906Data(const int nEvents = 20) 
{
  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;

  recoConsts* rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  // rc->set_CharFlag("AlignmentMille", "$E1039_RESOURCE/alignment/run6/align_mille.txt");
  // rc->set_CharFlag("AlignmentHodo", "$E1039_RESOURCE/alignment/run6/alignment_hodo.txt"); 
  // rc->set_CharFlag("AlignmentProp", "$E1039_RESOURCE/alignment/run6/alignment_prop.txt"); 
  // rc->set_CharFlag("Calibration", "$E1039_RESOURCE/alignment/run6/calibration.txt");
  // rc->set_CharFlag("Geometry", "user_liuk_geometry_G8_run5_2");
  // rc->set_CharFlag("MySQLURL", "mysql://e906-db1.fnal.gov:3306");
  //rc->Print();

  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(1000);
  se->Verbosity(0);
  
  GeomSvc::UseDbSvc(true);  //set to true to run E1039 style data
  //GeomSvc::UseDbSvc(false);
  GeomSvc* geom_svc = GeomSvc::instance();
  //std::cout << "print geom_svc table " << std::endl;
  //geom_svc->printTable();
  //std::cout << "done printing " << std::endl;

  SQReco* reco = new SQReco();
  reco->Verbosity(0);
  reco->set_legacy_rec_container(true);
  reco->set_geom_file_name("geom.root");
  //reco->set_enable_KF(true); //Kalman filter not needed for the track finding, disabling KF saves a lot of initialization time
  reco->setInputTy(SQReco::E1039);    //options are SQReco::E906 and SQReco::E1039
  reco->setFitterTy(SQReco::KFREF);  //not relavant for the track finding
  reco->set_evt_reducer_opt("aoce"); //if not provided, event reducer will be using JobOptsSvc to intialize; to turn off, set it to "none"
  reco->set_enable_eval(true);
  reco->set_eval_file_name("eval.root");
  //se->registerSubsystem(reco);

  //analysis module
  //gSystem->Load("libanamodule.so");
  AnaModule* ana = new AnaModule();
  ana->set_output_filename("ana.root");
  ana->set_reco(false);
  se->registerSubsystem(ana);

  //Now for a given detectorID, elementID, pos can be obtained by:
  //double pos = geom_svc->getMeasurement(detectorID, elementID);

  //For a given detectorID and tdcTime, the drift distance can be obtained by:
  //double driftdistance = geom_svc->getDriftDistance(detectorID, tdcTime);

  Fun4AllSRawEventInputManager* in = new Fun4AllSRawEventInputManager("SRawEventIM");
  //in->Verbosity(100);
  in->Verbosity(0);
  in->set_tree_name("save");
  in->set_branch_name("rawEvent");
  in->enable_E1039_translation();
  in->fileopen("/seaquest/users/cmantill/DarkQuest/e1039-analysis/DataHits/macro/nimSkim_data.root");
  se->registerInputManager(in);

  Fun4AllDstOutputManager* out = new Fun4AllDstOutputManager("DSTOUT", "result.root");
  //se->registerOutputManager(out);
  //out->Verbosity(100);
  out->Verbosity(0);
  out->AddNode("SRawEvent");
  out->AddNode("SRecEvent");

  se->run(nEvents);
  se->End();

  delete se;
  gSystem->Exit(0);
  return 0;
}
