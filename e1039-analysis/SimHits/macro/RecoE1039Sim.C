#include "G4_EMCal.C"
#include <top/G4_Beamline.C>
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
#include <top/G4_Target.C>

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
#include <string>
using namespace std;

/*
 * Macro used to analyze simulation in SpinQuest
 * isim = 1 to run on Aprime to dimuon signal
 * isim = 2 to run on Aprime to dielectron signal
 * isim = 3 to run on single gun (default is muon, see gun options below)
 * isim = 4 to run on DY to dimuon sample generated with Pythia
 * isim = 5 to run on J/psi to dimuon sample generated with Pythia
 * isim = 6 to run on cosmic sample
 *
 * is_displaced decides to run with the displaced config or not.
 * do_analysis runs the analysis ntuple
 * for Aprime signal, always run with is_displaced to True
 */

int RecoE1039Sim(const int nevent = 200,
                const int isim = 1,
                bool is_displaced = true,
                const bool do_analysis = true,
                std::string ifile="Brem_2.750000_z500_600_eps_-6.4",
                std::string out_file = "output.root"
                )
{
  // input simulation
  bool do_aprime_muon{false},do_aprime_electron{false},do_gun{false},do_dy{false},do_jpsi{false},do_cosmic{false},do_pion{false},do_trimuon{false};
  switch(isim){
  case 1: 
    do_aprime_muon = true;
    is_displaced = true; 
    break;
  case 2:
    do_aprime_electron = true;
    is_displaced = true;
    break;
  case 3:
    do_gun = true;
    break;
  case 4:
    do_dy = true;
    break;
  case 5:
    do_jpsi = true;
    break;
  case 6:
    do_cosmic = true;
    break;
  case 7:
    do_pion = true;
    break;
  case 8:
    do_trimuon = true;
    break;
  }

  // print more information with debug mode
  const bool isDEBUG = false;
  // verbosity option
  // https://github.com/E1039-Collaboration/e1039-core/blob/master/framework/fun4all/Fun4AllBase.h#L33-L55
  // the verbosity of different modules can also be modified separately for
  // debugging
  const int verbosity = 0;

  // legacy rec container
  const bool legacy_rec_container =  true; // false is for e1039 format

  // setup detectors in SpinQuest
  const bool do_collimator = true;
  const bool do_target     = true;
  const bool do_shielding  = true;
  const bool do_fmag       = true;
  const bool do_kmag       = true;
  const bool do_absorber   = true;
  const bool do_dphodo     = true;
  const bool do_station1DC = false; // station-1 drift chamber should be turned off by default
  const bool doEMCal       = true; // emcal turned on by default

  // SpinQuest constants
  const double target_coil_pos_z = -300;
  const double target_l = 7.9;                   // cm
  const double target_z = (7.9 - target_l) / 2.; // cm
  const int use_g4steps = 1;
  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;

  // SpinQuest reco constants
  recoConsts *rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  rc->set_CharFlag(
      "AlignmentMille",
      "$DIR_CMANTILL/macro/align_mille.txt"); // alignment file needed for EMCAL
  rc->set_CharFlag("fMagFile",
                   "$E1039_RESOURCE/geometry/magnetic_fields/tab.Fmag");
  rc->set_CharFlag("kMagFile",
                   "$E1039_RESOURCE/geometry/magnetic_fields/tab.Kmag");

  if (do_cosmic) {
    rc->init("cosmic");
    rc->set_BoolFlag("COARSE_MODE", true);
    rc->set_DoubleFlag("KMAGSTR", 0.);
    rc->set_DoubleFlag("FMAGSTR", 0.);
  }

  // track quality cuts (change these parameters for displaced tracking)
  // code from
  // https://cdcvs.fnal.gov/redmine/projects/seaquest-ktracker/repository/revisions/2a8f6204797b6ce142297ea2e158756fdd151552/diff
  // these parameters affect how we reconstruct tracklets (a straight line
  // segment in a single station) and St23 (combinations of pairs of tracklets
  // from St2 and 3) we usually try to minimize tracklets  (x,y,tx,ty)
  // parameters
  if (is_displaced) {
    rc->set_DoubleFlag("TX_MAX",
                       0.32);          // maximum allowed x slope for a tracklet
    rc->set_DoubleFlag("TY_MAX", 0.2); // maximum allowed y slope for a tracklet
    rc->set_DoubleFlag("X0_MAX",
                       500.0); // maximum allowed x position for a tracklet
    rc->set_DoubleFlag("Y0_MAX",
                       400.0); // maximum allowed y position for a tracklet
    rc->set_DoubleFlag("INVP_MAX",
                       0.5); // maximum inverse momentum (invp = 1/momentum)
    rc->set_BoolFlag(
        "NOT_DISPLACED",
        false); // running on displaced mode, no fitting to the target/vertex
  }
  if (isDEBUG) {
    rc->Print();
  }

  // geometry information
  GeomSvc::UseDbSvc(true);
  GeomSvc *geom_svc = GeomSvc::instance();
  if (isDEBUG) {
    // std::cout << "print geometry information" << std::endl;
    geom_svc->printWirePosition();
    // std::cout << " align printing " << std::endl;
    geom_svc->printAlignPar();
    // std::cout << " table printing" << std::endl;
    geom_svc->printTable();
    // std::cout << "done geometry printing" << std::endl;
  }

  // make the Server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(verbosity);

  // input to the simulation
  if(do_aprime_muon){ // aprime to displaced muons
    HepMCNodeReader *hr = new HepMCNodeReader();
    hr->set_particle_filter_on(true);
    hr->insert_particle_filter_pid(13); // filter muons
    hr->insert_particle_filter_pid(13 * -1);
    hr->Verbosity(verbosity);
    se->registerSubsystem(hr);
  }
  else if(do_aprime_electron){ // aprime to displaced electrons
    HepMCNodeReader *hr = new HepMCNodeReader();
    hr->set_particle_filter_on(true);
    hr->insert_particle_filter_pid(11);
    hr->insert_particle_filter_pid(11*-1);
    hr->Verbosity(verbosity);
    se->registerSubsystem(hr);
  }
  else if(do_trimuon){
    HepMCNodeReader *hr = new HepMCNodeReader();
    //hr->set_particle_filter_on(true);
    //hr->insert_particle_filter_pid(13);
    //hr->insert_particle_filter_pid(13 * -1);
    hr->Verbosity(verbosity);
    se->registerSubsystem(hr);
  }
  else if(do_gun){ // particle gun
    PHG4SimpleEventGenerator *genp = new PHG4SimpleEventGenerator("MUP");
    genp->add_particles("mu+", 1);  // mu+
    //genp->add_particles("mu-", 1); // mu-
    //genp->add_particles("e+", 1); // positron
    //genp->add_particles("pi+", 1); // pions
    //genp->add_particles("kaon0L", 1); // k0long
    //genp->add_particles("proton", 1); // protons
    //genp->add_particles("Upsilon", 1); // upsilon
    genp->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
					   PHG4SimpleEventGenerator::Uniform,
					   PHG4SimpleEventGenerator::Uniform);
    if(is_displaced){
      genp->set_vertex_distribution_mean(0.0, 0.0, 500.);
    } else {
      genp->set_vertex_distribution_mean(0.0, 0.0, target_coil_pos_z);
    }

    genp->set_vertex_distribution_width(0.0, 0.0, 0.0);
    genp->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
    genp->set_vertex_size_parameters(0.0, 0.0);
    genp->set_pxpypz_range(-1., 1., -1., 1., 0., 100.);
    genp->Verbosity(verbosity);
    se->registerSubsystem(genp);
  } else if (do_dy or do_jpsi) {
    PHPythia8 *pythia8 = new PHPythia8();
    pythia8->Verbosity(verbosity);
    if (do_dy)
      pythia8->set_config_file("phpythia8_DY.cfg");
    else
      pythia8->set_config_file("phpythia8_Jpsi.cfg");
    if (is_displaced) {
      pythia8->set_vertex_distribution_mean(0.0, 0.0, 500., 0);
    } else {
      pythia8->set_vertex_distribution_mean(0, 0, target_coil_pos_z, 0);
    }
    pythia8->set_embedding_id(1);
    se->registerSubsystem(pythia8);

    pythia8->set_trigger_AND();
    PHPy8ParticleTrigger *trigger_mup = new PHPy8ParticleTrigger();
    trigger_mup->AddParticles("-13");
    // trigger_mup->SetPxHighLow(7, 0.5);
    // trigger_mup->SetPyHighLow(6, -6);
    trigger_mup->SetPzHighLow(120, 30);
    pythia8->register_trigger(trigger_mup);

    PHPy8ParticleTrigger *trigger_mum = new PHPy8ParticleTrigger();
    trigger_mum->AddParticles("13");
    // trigger_mum->SetPxHighLow(-0.5, 7);
    // trigger_mum->SetPyHighLow(6, -6);
    trigger_mum->SetPzHighLow(120, 30);
    pythia8->register_trigger(trigger_mum);

    HepMCNodeReader *hr = new HepMCNodeReader();
    hr->set_particle_filter_on(true);
    hr->insert_particle_filter_pid(13);
    hr->insert_particle_filter_pid(-13);
    se->registerSubsystem(hr);
  } else if (do_cosmic) {
    SQCosmicGen *cosmicGen = new SQCosmicGen();
    se->registerSubsystem(cosmicGen);
  } else {
    std::cout << " No input! " << std::endl;
    return 0;
  }

  // fun4All G4 module
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
  g4Reco->SetWorldMaterial("G4_AIR"); // G4_Galactic, G4_AIR
  // G4 Physics list to use
  g4Reco->SetPhysicsList("FTFP_BERT");

  // setup detectors
  SetupInsensitiveVolumes(g4Reco, do_shielding, do_fmag, do_kmag,
                          do_absorber); // insensitive volumes
  SetupBeamline(
      g4Reco, do_collimator,
      target_coil_pos_z -
          302.36); // collimator, targer and shielding between target and FMag
  if (do_target) {
    SetupTarget(g4Reco, target_coil_pos_z, target_l, target_z, use_g4steps);
  }
  // sensitive elements of the spectrometer
  SetupSensitiveDetectors(g4Reco, do_dphodo, do_station1DC, "SQ_ArCO2",
                          "SQ_Scintillator", 1);
  if (doEMCal) {
    SetupEMCal(g4Reco, "EMCal", 0., 0., 1930.);
  }
  se->registerSubsystem(g4Reco);

  // save G4 truth info to the Node Tree
  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  // digitizer
  SQDigitizer *digitizer = new SQDigitizer("DPDigitizer", 0);
  digitizer->Verbosity(verbosity);
  digitizer->set_enable_st1dc(do_station1DC); // these two lines need to be in
                                              // sync with the parameters used
  digitizer->set_enable_dphodo(
      do_dphodo); // in the SetupSensitiveVolumes() function call above
  if (doEMCal) {
    digitizer->registerEMCal("EMCal", 100);
  }
  se->registerSubsystem(digitizer);

  // tracking module
  std::cout << "*********** Start SQReco step now..." << std::endl;
  SQReco* reco = new SQReco();
  reco->Verbosity(verbosity);
  reco->set_legacy_rec_container(legacy_rec_container);
  //reco->set_geom_file_name("support/geom.root"); // not needed as it's created on the fly
  reco->set_enable_KF(true);                      // Kalman filter not needed for the track finding, disabling KF saves a lot of initialization time
  reco->setInputTy(SQReco::E1039);                // options are SQReco::E906 and SQReco::E1039
  reco->setFitterTy(SQReco::KFREF);               // not relevant for the track finding, options are SQReco::KFREF and SQReco::LEGACY
  reco->set_evt_reducer_opt("none");              // if not provided, event reducer will be using JobOptsSvc to intialize; to turn off, set it to "none", for normal tracking, set to something like "aoc"
  reco->set_enable_eval(false);                   // set to true to generate evaluation file which includes final track candidates 
  reco->set_eval_file_name("eval.root");          // evaluation filename
  reco->set_enable_eval_dst(false);               // set to true to include final track candidates in the DST tree
  //reco->add_eval_list(3);                         // include back partial tracks in eval tree for debuging
  //reco->add_eval_list(2);                         // include station-3+/- in eval tree for debuging
  //reco->add_eval_list(1);                         // include station-2 in eval tree for debugging
  se->registerSubsystem(reco);

  // truth node maker after tracking
  TruthNodeMaker* truthMaker = new TruthNodeMaker();
  truthMaker->set_legacy_rec_container(legacy_rec_container);
  if(do_aprime_muon or do_aprime_electron){
    truthMaker->set_m_process_type(3); // set process type to 3 (A' -> di lepton) since we only have a 3 particle process instead of 0+1->2+3
  }
  if(do_trimuon){
    truthMaker->set_m_process_type(3);
  }
  truthMaker->Verbosity(verbosity);
  se->registerSubsystem(truthMaker);

  // trigger emulator
  // needs TruthNodeMaker to associate the trigger to SQEvent 
  DPTriggerAnalyzer* dptrigger = new DPTriggerAnalyzer();                                                                                                                                                 
  dptrigger->set_road_set_file_name("$E1039_RESOURCE/trigger/trigger_67.txt");
  dptrigger->Verbosity(verbosity);
  se->registerSubsystem(dptrigger);  

  // event filter
  EvtFilter *evt_filter = new EvtFilter();
  evt_filter->Verbosity(verbosity);
  evt_filter->set_trigger_req(1<<5);
  //se->registerSubsystem(evt_filter);

  // truth vertexing
  SQTruthVertexing* truthVtx = new SQTruthVertexing();
  truthVtx->set_legacy_rec_container(legacy_rec_container);
  truthVtx->set_vtx_smearing(50.); // smear the truth z_vertex to mimic resolution effect, default is 0.
  //se->registerSubsystem(truthVtx);

  // vertexing for dimuon information
  if(legacy_rec_container){
    VertexFit* vertexing = new VertexFit();
    se->registerSubsystem(vertexing);
  }

  // analysis module
  gSystem->Load("libsim_ana.so");
  SimAna *sim_ana = new SimAna();  
  sim_ana->Verbosity(verbosity);
  sim_ana->set_out_name(out_file);
  sim_ana->set_legacy_rec_container(legacy_rec_container);
  if(do_analysis){
    se->registerSubsystem(sim_ana);    
  }

  // input 
  if(do_aprime_muon or do_aprime_electron){
    // use hepmc input
    Fun4AllHepMCInputManager *in = new Fun4AllHepMCInputManager("HEPMCIN");
    se->registerInputManager(in);
    stringstream ssin;
    ssin << "$DIR_CMANTILL/../../lhe/output/displaced_Aprime_Muons/" << ifile
         << ".txt";
    in->fileopen(gSystem->ExpandPathName(ssin.str().c_str()));
    in->Verbosity(verbosity);
    se->registerInputManager(in);
  } 
  else if(do_trimuon){
    Fun4AllHepMCInputManager *in = new Fun4AllHepMCInputManager("HEPMCIN");
    se->registerInputManager(in);
    stringstream ssin;
    ssin << "$DIR_CMANTILL/../../lhe/output/trimuon_0.5MS0gS1.hepmc";
    in->fileopen(gSystem->ExpandPathName(ssin.str().c_str()));
    in->Verbosity(verbosity);
    se->registerInputManager(in);
  }
  else {
    // need a dummy input to drive the event loop
    Fun4AllInputManager *in = new Fun4AllDummyInputManager("DUMMY");
    in->Verbosity(verbosity);
    se->registerInputManager(in);
  }

  se->run(nevent);

  // export the geometry
  // PHGeomUtility::ExportGeomtry(se->topNode(),"geom.root");

  // finish job - close and save output files
  se->End();
  se->PrintTimer();
  std::cout << "All done" << std::endl;

  // cleanup - delete the server and exit
  delete se;
  gSystem->Exit(0);
  return 0;
}
