#include <top/G4_Beamline.C>
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
#include <top/G4_Target.C>

#include "G4_EMCal.C"

R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libPHPythia8)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libdptrigger)
R__LOAD_LIBRARY(libevt_filter)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libSQPrimaryGen)
R__LOAD_LIBRARY(libsim_ana)

#include <iostream>
#include <sstream>
using namespace std;

int Fun4Sim_Muon(const int nevent = 10000) {
  const double target_coil_pos_z = -300;

  const bool do_collimator = true;
  const bool do_target = true;
  const bool do_e1039_shielding = true;
  const bool do_fmag = true;
  const bool do_kmag = true;
  const bool do_absorber = true;
  const bool do_dphodo = true;
  const bool do_station1DC =
      false;  // station-1 drift chamber should be turned off by default

  const double target_l = 7.9;                    // cm
  const double target_z = (7.9 - target_l) / 2.;  // cm
  const int use_g4steps = 1;

  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;

  const bool save_in_acc =
      false;  //< Set true to save only in-acceptance events into DST.

  recoConsts *rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  rc->set_CharFlag("AlignmentMille", "$DIR_CMANTILL/macro/align_mille.txt");
  rc->set_CharFlag("fMagFile",
                   "$E1039_RESOURCE/geometry/magnetic_fields/tab.Fmag");
  rc->set_CharFlag("kMagFile",
                   "$E1039_RESOURCE/geometry/magnetic_fields/tab.Kmag");
  rc->set_BoolFlag("COARSE_MODE", true);
  // parameters tuned for displaced tracking
  // numbers taken from
  // https://cdcvs.fnal.gov/redmine/projects/seaquest-ktracker/repository/revisions/2a8f6204797b6ce142297ea2e158756fdd151552/diff
  rc->set_DoubleFlag("TX_MAX", 0.32);
  rc->set_DoubleFlag("TY_MAX", 0.2);
  rc->set_DoubleFlag("X0_MAX", 500.0);
  rc->set_DoubleFlag("Y0_MAX", 400.0);
  rc->set_DoubleFlag("INVP_MAX", 0.5);
  rc->Print();

  GeomSvc::UseDbSvc(true);
  GeomSvc *geom_svc = GeomSvc::instance();
  std::cout << "print geometry information" << std::endl;
  geom_svc->printWirePosition();
  std::cout << " align printing " << std::endl;
  geom_svc->printAlignPar();
  std::cout << " table printing" << std::endl;
  geom_svc->printTable();
  std::cout << "done geometry printing" << std::endl;

  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
  // se->Verbosity(99);

  bool single_gun = true;

  if (single_gun)  // change the hard-coded numbers to change the initial
                   // vertex/momentum distribution
  {
    PHG4SimpleEventGenerator *genp = new PHG4SimpleEventGenerator("MUP");
    genp->add_particles("mu+", 1);  // mu+,mu-,e+,proton,pi+,Upsilon
    genp->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                           PHG4SimpleEventGenerator::Uniform,
                                           PHG4SimpleEventGenerator::Uniform);
    //genp->set_vertex_distribution_mean(0.0, 0.0, target_coil_pos_z);
    genp->set_vertex_distribution_mean(0.0, 0.0, 500.);
    genp->set_vertex_distribution_width(0.0, 0.0, 0.0);
    genp->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
    genp->set_vertex_size_parameters(0.0, 0.0);
    genp->set_pxpypz_range(-1., 1., -1., 1., 0., 100.);
    se->registerSubsystem(genp);
  }

  // Fun4All G4 module
  PHG4Reco *g4Reco = new PHG4Reco();
  g4Reco->set_field_map(
      rc->get_CharFlag("fMagFile") + " " + rc->get_CharFlag("kMagFile") + " " +
          Form("%f", FMAGSTR) + " " + Form("%f", KMAGSTR) + " " + "5.0",
      PHFieldConfig::RegionalConst);
  // size of the world - every detector has to fit in here
  g4Reco->SetWorldSizeX(1000);
  g4Reco->SetWorldSizeY(1000);
  g4Reco->SetWorldSizeZ(5000);
  // shape of our world - it is a tube
  g4Reco->SetWorldShape("G4BOX");
  // this is what our world is filled with
  g4Reco->SetWorldMaterial("G4_AIR");  // G4_Galactic, G4_AIR
  // Geant4 Physics list to use
  g4Reco->SetPhysicsList("FTFP_BERT");
  // insensitive elements of the spectrometer
  SetupInsensitiveVolumes(g4Reco, do_e1039_shielding, do_fmag, do_kmag,
                          do_absorber);
  // collimator, targer and shielding between target and FMag
  SetupBeamline(g4Reco, do_collimator,
                target_coil_pos_z - 302.36);  // Is the position correct??
  if (do_target) {
    SetupTarget(g4Reco, target_coil_pos_z, target_l, target_z, use_g4steps);
  }
  // sensitive elements of the spectrometer
  // SetupSensitiveDetectors(g4Reco, do_dphodo, do_station1DC, "SQ_ArCO2",
  // "SQ_Scintillator", 1);
  SetupSensitiveDetectors(g4Reco);
  se->registerSubsystem(g4Reco);

  // if (save_in_acc) se->registerSubsystem(new RequireParticlesInAcc());

  // save truth info to the Node Tree
  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  // digitizer
  SQDigitizer *digitizer = new SQDigitizer("DPDigitizer", 0);
  // digitizer->Verbosity(99);
  digitizer->set_enable_st1dc(do_station1DC);  // these two lines need to be in
                                               // sync with the parameters used
  digitizer->set_enable_dphodo(
      do_dphodo);  // in the SetupSensitiveVolumes() function call above
  // digitizer->registerEMCal("EMCal", 100);
  se->registerSubsystem(digitizer);

  // Make SQ nodes for truth info
  // se->registerSubsystem(new TruthNodeMaker());

  // Trigger Emulator
  DPTriggerAnalyzer *dptrigger = new DPTriggerAnalyzer();
  dptrigger->set_road_set_file_name("$E1039_RESOURCE/trigger/trigger_67.txt");
  // se->registerSubsystem(dptrigger);

  // Event Filter
  // EvtFilter *evt_filter = new EvtFilter();
  // evt_filter->Verbosity(10);
  // evt_filter->set_trigger_req(1<<5);
  // se->registerSubsystem(evt_filter);

  // Tracking module
  std::cout << "*********** Start SQReco step now..." << std::endl;
  SQReco *reco = new SQReco();
  reco->Verbosity(100);
  // reco->set_geom_file_name("support/geom.root"); //not needed as it's created
  // on the fly
  reco->set_enable_KF(true);  // Kalman filter not needed for the track finding,
                              // disabling KF saves a lot of initialization time
  reco->setInputTy(SQReco::E1039);  // options are SQReco::E906 and
                                    // SQReco::E1039
  // reco->setFitterTy(SQReco::KFREF);    //not relavant for the track finding
  reco->setFitterTy(SQReco::LEGACY);
  reco->set_evt_reducer_opt(
      "none");  // if not provided, event reducer will be using JobOptsSvc to
                // intialize; to turn off, set it to "none", for normal tracking,
                // set to something like "aoc"
  reco->set_enable_eval(true);  // set to true to generate evaluation file which
                                // includes final track candidates
  reco->set_eval_file_name("eval.root");
  reco->set_enable_eval_dst(
      false);  // set to true to include final track cnadidates in the DST tree
  // reco->add_eval_list(3);             //include back partial tracks in eval
  // tree for debuging reco->add_eval_list(2);             //include station-3+/-
  // in eval tree for debuging reco->add_eval_list(1);             //include
  // station-2 in eval tree for debugging reco->set_legacy_rec_container(false);
  se->registerSubsystem(reco);

  // VertexFit* vertexing = new VertexFit();
  // se->registerSubsystem(vertexing);

  gSystem->Load("libsim_ana.so");
  SimAna *sim_ana = new SimAna();
  se->registerSubsystem(sim_ana);

  // input - we need a dummy to drive the event loop
  Fun4AllInputManager *in = new Fun4AllDummyInputManager("DUMMY");
  se->registerInputManager(in);

  ///////////////////////////////////////////
  // Output
  ///////////////////////////////////////////

  // DST output manager
  Fun4AllDstOutputManager *out =
      new Fun4AllDstOutputManager("DSTOUT", "DST.root");
  se->registerOutputManager(out);

  se->run(nevent);

  PHGeomUtility::ExportGeomtry(se->topNode(), "geom.root");

  // finish job - close and save output files
  se->End();
  se->PrintTimer();
  std::cout << "All done" << std::endl;

  // cleanup - delete the server and exit
  delete se;
  gSystem->Exit(0);
  return 0;
}
