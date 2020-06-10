#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <TSystem.h>
#include "G4_SensitiveDetectors.C"
#include "G4_Target.C"

R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libdptrigger)
R__LOAD_LIBRARY(libevt_filter)
//R__LOAD_LIBRARY(libktracker)
//R__LOAD_LIBRARY(libsim_eval)
#endif

#include <iostream>
#include <sstream>
using namespace std;

int Fun4Sim(
	    const int nevent = 10,
	    std::string ifile = "Brem_0.06_z500_600_eps_-6.2"
    )
{
  const double target_coil_pos_z = -300;
  const int nmu = 1;
  const int idLep = 11;

  const bool do_collimator = true;
  const bool do_target = true;
  const bool do_e1039_shielding = true;

  const double target_l = 7.9; //cm
  const double target_z = (7.9-target_l)/2.; //cm
  const int use_g4steps = 1;

  const int RUNNUMBER = 1;
  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;

  gSystem->Load("libfun4all");
  gSystem->Load("libg4detectors");
  gSystem->Load("libg4testbench");
  gSystem->Load("libg4eval");
  gSystem->Load("libg4dst");

  recoConsts *rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  //  rc->set_IntFlag("RUNNUMBER", RUNNUMBER);
  rc->Print();

  // Alignment
  JobOptsSvc *jobopt_svc = JobOptsSvc::instance();
  jobopt_svc->init("run7_sim.opts");
  GeomSvc::UseDbSvc(true);
  GeomSvc *geom_svc = GeomSvc::instance();

  // Make the Server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // HepMC reader
  HepMCNodeReader *hr = new HepMCNodeReader();
  hr->set_particle_filter_on(true);
  hr->insert_particle_filter_pid(idLep);
  hr->insert_particle_filter_pid(idLep*-1);
  se->registerSubsystem(hr);

  // Fun4All G4 module
  PHG4Reco *g4Reco = new PHG4Reco();
  //PHG4Reco::G4Seed(123);
  //g4Reco->set_field(5.);
  g4Reco->set_field_map(jobopt_svc->m_fMagFile+" "+
			jobopt_svc->m_kMagFile+" "+
			Form("%f",FMAGSTR) + " " +
			Form("%f",KMAGSTR) + " " +
			"5.0",
			PHFieldConfig::RegionalConst);
  // size of the world - every detector has to fit in here
  g4Reco->SetWorldSizeX(1000);
  g4Reco->SetWorldSizeY(1000);
  g4Reco->SetWorldSizeZ(5000);
  // shape of our world - it is a tube
  g4Reco->SetWorldShape("G4BOX");
  // this is what our world is filled with
  g4Reco->SetWorldMaterial("G4_AIR"); //G4_Galactic, G4_AIR
  // Geant4 Physics list to use
  g4Reco->SetPhysicsList("FTFP_BERT");

  // insensitive elements of the spectrometer
  PHG4E1039InsensSubsystem* insens = new PHG4E1039InsensSubsystem("Insens");
  g4Reco->registerSubsystem(insens);

  // collimator, targer and shielding between target and FMag
  gROOT->LoadMacro("G4_Target.C");
  SetupTarget(g4Reco, do_collimator, do_target, do_e1039_shielding, target_coil_pos_z, target_l, target_z, use_g4steps);

  // sensitive elements of the spectrometer
  gROOT->LoadMacro("G4_SensitiveDetectors.C");
  SetupSensitiveDetectors(g4Reco, 0);

  se->registerSubsystem(g4Reco);

  // save truth info to the Node Tree
  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  // digitizer
  DPDigitizer *digitizer = new DPDigitizer("DPDigitizer", 0);
  se->registerSubsystem(digitizer);

  // Trigger Emulator
  gSystem->Load("libdptrigger.so");
  DPTriggerAnalyzer* dptrigger = new DPTriggerAnalyzer();
  dptrigger->set_hit_container_choice("Vector");
  dptrigger->set_road_set_file_name(gSystem->ExpandPathName("$E1039_RESOURCE/trigger/trigger_67.txt"));
  se->registerSubsystem(dptrigger);

  // Event Filter
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
  //gSystem->Load("libsim_eval.so");
  //SimEval *sim_eval = new SimEval();
  //sim_eval->Verbosity(0);//(3);
  //sim_eval->set_hit_container_choice("Vector");
  //stringstream ssout; ssout << "sim_eval_" << ifile << ".root";
  //sim_eval->set_out_name(ssout.str().c_str());
  //se->registerSubsystem(sim_eval);

  // input 
  Fun4AllHepMCInputManager *in = new Fun4AllHepMCInputManager("HEPMCIN");
  //in->set_vertex_distribution_mean(0,0,target_coil_pos_z,0);
  se->registerInputManager(in);
  //stringstream ssin; ssin << "$DIR_TOP/../../lhe/displaced_Aprime_Electrons/" << ifile << ".txt";
  stringstream ssin; ssin << ifile << ".txt";
  in->fileopen(gSystem->ExpandPathName(ssin.str().c_str()));

  // DST output manager
  stringstream ssout; ssout << ifile << "0_dst.root";
  //std::cout << ssout.str().c_str() << std::endl;
  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", ssout.str().c_str());
  se->registerOutputManager(out);

  if (nevent >= 0)
  {
    se->run(nevent);

    // finish job - close and save output files
    se->End();
    //se->PrintTimer();

    // cleanup - delete the server and exit
    delete se;
    gSystem->Exit(0);
  }

  return 0;
}

