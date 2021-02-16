#include <top/G4_Beamline.C>
#include <top/G4_Target.C>
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
#include "G4_EMCal.C"

R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libdptrigger)
R__LOAD_LIBRARY(libevt_filter)
//R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libsim_eval)

using namespace std;

#include <iostream>
#include <sstream>
using namespace std;

int Fun4SimHitMuon(
	 const int nevent = 100
    )
{
  const double target_coil_pos_z = -300;
  const double collimator_pos_z = -602.36;

  const bool do_collimator = true;
  const bool do_target = true;
  const bool do_shielding = true;
  const bool do_fmag = true;
  const bool do_kmag = true;
  const bool do_absorber = true;
  const bool do_dphodo = true;
  const bool do_station1DC = false;       //station-1 drift chamber should be turned off by default

  const double target_l = 7.9; //cm
  const double target_z = (7.9-target_l)/2.; //cm
  const int use_g4steps = 1;

  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;

  recoConsts *rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  rc->set_CharFlag("AlignmentMille", "$DIR_CMANTILL/macro/align_mille.txt");
  rc->set_CharFlag("fMagFile", "$E1039_RESOURCE/geometry/magnetic_fields/tab.Fmag");
  rc->set_CharFlag("kMagFile", "$E1039_RESOURCE/geometry/magnetic_fields/tab.Kmag");
  rc->Print();

  // Alignment
  //JobOptsSvc *jobopt_svc = JobOptsSvc::instance();
  //jobopt_svc->init("run7_sim.opts");
  GeomSvc::UseDbSvc(true);
  GeomSvc *geom_svc = GeomSvc::instance();
  geom_svc->printTable();

  // Make the Server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Generation
  PHG4ParticleGun *gun = new PHG4ParticleGun("PGUN");
  gun->set_name("mu-");
  gun->set_vtx(0, 0, 500);
  gun->set_mom(0, 0, 120);
  se->registerSubsystem(gun);

  // Fun4All G4 module
  PHG4Reco *g4Reco = new PHG4Reco();
  g4Reco->set_field_map(
			rc->get_CharFlag("fMagFile")+" "+
                        rc->get_CharFlag("kMagFile")+" "+
			Form("%f",FMAGSTR) + " " +
			Form("%f",KMAGSTR) + " " +
			"5.0",
			PHFieldConfig::RegionalConst);
  g4Reco->SetWorldSizeX(1000); // size of the world - every detector has to fit in here 
  g4Reco->SetWorldSizeY(1000);
  g4Reco->SetWorldSizeZ(5000);
  g4Reco->SetWorldShape("G4BOX");  // shape of our world - it is a box
  g4Reco->SetWorldMaterial("G4_AIR"); // this is what our world is filled with G4_Galactic, G4_AIR
  g4Reco->SetPhysicsList("FTFP_BERT"); //Geant4 Physics list to use 

  SetupBeamline(g4Reco, do_collimator, collimator_pos_z);
  if (do_target) {
    SetupTarget(g4Reco, target_coil_pos_z, target_l, target_z, use_g4steps, 0);
  }
  SetupInsensitiveVolumes(g4Reco, do_shielding, do_fmag, do_kmag, do_absorber);
  SetupSensitiveDetectors(g4Reco, do_dphodo, do_station1DC);
  SetupEMCal(g4Reco, "EMCal", 0., -110., 1930.);
  //g4Reco->Verbosity(99);
  se->registerSubsystem(g4Reco);

  // save truth info to the Node Tree
  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  // digitizer
  SQDigitizer* digitizer = new SQDigitizer("Digitizer", 0);
  digitizer->set_enable_st1dc(do_station1DC);    // these two lines need to be in sync with the parameters used
  digitizer->set_enable_dphodo(do_dphodo);       // in the SetupSensitiveVolumes() function call above
  digitizer->registerEMCal("EMCal", 100);
  //digitizer->Verbosity(99); 
  se->registerSubsystem(digitizer);

  // Trigger Emulator
  gSystem->Load("libdptrigger.so");
  DPTriggerAnalyzer* dptrigger = new DPTriggerAnalyzer();
  dptrigger->set_hit_container_choice("Vector");
  dptrigger->set_road_set_file_name(gSystem->ExpandPathName("$E1039_RESOURCE/trigger/trigger_67.txt"));
  se->registerSubsystem(dptrigger);

  // Event Filter (not requiring trigger now)
  EvtFilter *evt_filter = new EvtFilter();
  se->registerSubsystem(evt_filter);

  // tracking module
  // gSystem->Load("libktracker.so");
  // KalmanFastTrackingWrapper *ktracker = new KalmanFastTrackingWrapper();
  // ktracker->set_enable_event_reducer(true);
  // ktracker->set_DS_level(0);
  // ktracker->set_pattern_db_name(gSystem->ExpandPathName("$E1039_RESOURCE/dsearch/v1/pattern.root"));
  // se->registerSubsystem(ktracker);

  // VertexFit* vertexing = new VertexFit();
  // se->registerSubsystem(vertexing);

  // evaluation module
  gSystem->Load("libsim_eval.so");
  SimEval *sim_eval = new SimEval();
  sim_eval->set_hit_container_choice("Vector");
  stringstream ssout; ssout << "$DIR_TOP/macro/sim_eval_muon.root";
  sim_eval->set_out_name(ssout.str().c_str());
  //sim_eval->Verbosity(10);
  se->registerSubsystem(sim_eval);

  // input 
  Fun4AllInputManager *in = new Fun4AllDummyInputManager("DUMMY");
  se->registerInputManager(in);

  // output
  stringstream ssout_3; ssout_3 << "$DIR_TOP/macro/muon0_dst.root";
  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", ssout_3.str().c_str());
  se->registerOutputManager(out);
  
  // run
  se->run(nevent);
  
  // finish job - close and save output files
  se->End();
  se->PrintTimer();
  std::cout << "All done" << std::endl;

  // cleanup - delete the server and exit
  delete se;
  gSystem->Exit(0);
  return 0;
}

