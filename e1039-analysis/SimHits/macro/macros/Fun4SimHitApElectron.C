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
R__LOAD_LIBRARY(libSQPrimaryGen)
R__LOAD_LIBRARY(libsim_eval)

#include <iostream>
#include <sstream>
using namespace std;

int Fun4SimHitApElectron(const int nevent = 10,
	    std::string ifile = "Brem_1.041330_z500_600_eps_-6.4",
			 const int idLep = 11)
{
  const double target_coil_pos_z = -300;

  const bool do_collimator = true;
  const bool do_target = true;
  const bool do_e1039_shielding = true;
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

  // Alignment
  recoConsts *rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  //rc->set_CharFlag("AlignmentMille", "$DIR_CMANTILL/macro/align_mille.txt");
  rc->Print();

  GeomSvc::UseDbSvc(true);
  GeomSvc *geom_svc = GeomSvc::instance();

  // Make the Server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Generation
  HepMCNodeReader *hr = new HepMCNodeReader();
  hr->set_particle_filter_on(true);
  hr->insert_particle_filter_pid(idLep);
  hr->insert_particle_filter_pid(idLep*-1);
  se->registerSubsystem(hr);

  PHG4Reco *g4Reco = new PHG4Reco();
  g4Reco->set_field_map(
      rc->get_CharFlag("fMagFile")+" "+
      rc->get_CharFlag("kMagFile")+" "+
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
  SetupInsensitiveVolumes(g4Reco, do_e1039_shielding, do_fmag, do_kmag, do_absorber);

  // collimator, targer and shielding between target and FMag
  SetupBeamline(g4Reco, do_collimator, target_coil_pos_z - 302.36); // Is the position correct??

  if (do_target) {
    SetupTarget(g4Reco, target_coil_pos_z, target_l, target_z, use_g4steps);
  }

  // sensitive elements of the spectrometer
  SetupSensitiveDetectors(g4Reco, do_dphodo, do_station1DC);
  SetupEMCal(g4Reco, "EMCal", 0., -110., 1930.);
  se->registerSubsystem(g4Reco);

  // save truth info to the Node Tree
  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  // digitizer
  SQDigitizer *digitizer = new SQDigitizer("DPDigitizer", 0);
  //digitizer->Verbosity(99);
  digitizer->set_enable_st1dc(do_station1DC);    // these two lines need to be in sync with the parameters used
  digitizer->set_enable_dphodo(do_dphodo);       // in the SetupSensitiveVolumes() function call above
  //digitizer->registerEMCal("EMCal", 100);
  se->registerSubsystem(digitizer);

  // Make SQ nodes for truth info
  se->registerSubsystem(new TruthNodeMaker());

  // Trigger Emulator
  DPTriggerAnalyzer* dptrigger = new DPTriggerAnalyzer();
  dptrigger->set_road_set_file_name("$E1039_RESOURCE/trigger/trigger_67.txt");
  se->registerSubsystem(dptrigger);

  // Evaluation
  gSystem->Load("libsim_eval.so");
  SimEval *sim_eval = new SimEval();
  sim_eval->set_hit_container_choice("Vector");
  stringstream ssout; ssout << "$DIR_TOP/macro/simeval_electrons_emcal/sim_eval_" << ifile << ".root";
  sim_eval->set_out_name(ssout.str().c_str());
  sim_eval->Verbosity(99);

  // Input
  Fun4AllHepMCInputManager *in = new Fun4AllHepMCInputManager("HEPMCIN");
  se->registerInputManager(in);
  stringstream ssin; ssin << "$DIR_CMANTILL/../../lhe/displaced_Aprime_Electrons/" << ifile << ".txt";
  in->fileopen(gSystem->ExpandPathName(ssin.str().c_str()));

  // Output
  stringstream ssout_3; ssout_3 << "$DIR_TOP/macro/output_electrons_emcal/" << ifile << "0_dst.root";
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
