
#define GRID

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <TSystem.h>

#include "G4_SensitiveDetectors.C"
#include "G4_Target.C"

#ifndef GRID
#include "EventDisplay.C"
#endif

R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libPHPythia8)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libdptrigger)
R__LOAD_LIBRARY(libembedding)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libmodule_example)
#endif


using namespace std;

int Fun4DB(
    const int nevent = 1
    )
{
  const double target_coil_pos_z = -300;
  const int nmu = 1;
  int embedding_opt = 0;

  const bool do_collimator = true;
  const bool do_target = true;
  const bool do_e1039_shielding = true;
  const double target_l = 7.9; //cm
  const double target_z = (7.9-target_l)/2.; //cm
  const int use_g4steps = 1;

  const double FMAGSTR = -1.044; //-1.054;
  const double KMAGSTR = -1.025; //-0.951;

  const bool gen_pythia8  = true;
  const bool gen_gun      = false;
  const bool gen_particle = false;
  const bool read_hepmc   = false;

  gSystem->Load("libfun4all");
  gSystem->Load("libg4detectors");
  gSystem->Load("libg4testbench");
  gSystem->Load("libg4eval");
  gSystem->Load("libktracker.so");

  recoConsts *rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  rc->Print();

  JobOptsSvc *jobopt_svc = JobOptsSvc::instance();
  jobopt_svc->init("run6_sim.opts");

  GeomSvc *geom_svc = GeomSvc::instance();
  //const double x0_shift = 0.0; //cm 
  //std::cout << "D2X::X0: " << geom_svc->getDetectorX0("D2X") << std::endl;
  //geom_svc->setDetectorX0("D2X", geom_svc->getDetectorX0("D2X")+x0_shift);
  //std::cout << "D2X::X0: " << geom_svc->getDetectorX0("D2X") << std::endl;

  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(10);

  // pythia8
  if(gen_pythia8) {
    gSystem->Load("libPHPythia8.so");

    PHPythia8 *pythia8 = new PHPythia8();
    pythia8->Verbosity(1);
    pythia8->set_config_file("phpythia8_DY.cfg");
    pythia8->set_vertex_distribution_mean(0, 0, target_coil_pos_z, 0);
    pythia8->set_embedding_id(1);
    se->registerSubsystem(pythia8);

    pythia8->set_trigger_AND();

    PHPy8ParticleTrigger* trigger_mup = new PHPy8ParticleTrigger();
    trigger_mup->AddParticles("-13");
    //trigger_mup->SetPxHighLow(7, 0.5);
    //trigger_mup->SetPyHighLow(6, -6);
    trigger_mup->SetPzHighLow(120, 10);
    pythia8->register_trigger(trigger_mup);

    PHPy8ParticleTrigger* trigger_mum = new PHPy8ParticleTrigger();
    trigger_mum->AddParticles("13");
    //trigger_mum->SetPxHighLow(-0.5, -7);
    //trigger_mum->SetPyHighLow(6, -6);
    trigger_mum->SetPzHighLow(120, 10);
    pythia8->register_trigger(trigger_mum);
  }
  
  if(gen_pythia8 || read_hepmc) {
    HepMCNodeReader *hr = new HepMCNodeReader();
    hr->set_particle_filter_on(true);
    hr->insert_particle_filter_pid(13);
    hr->insert_particle_filter_pid(-13);
    se->registerSubsystem(hr);
  }

  // single gun
  if(gen_gun) {
    // 15.529542 * 21.995780 * -0.628401 * 1.3662691 * 60.306243
    PHG4ParticleGun *gun = new PHG4ParticleGun("GUN");
    gun->set_name("mu-");
    gun->set_vtx(16, 22, 550);
    gun->set_mom(-0.6, 1.4, 60);
    se->registerSubsystem(gun);
  }

  // multi particle gun
  if(gen_particle) {
    PHG4SimpleEventGenerator *genp = new PHG4SimpleEventGenerator("MUP");
    //genp->set_seed(123);
    genp->add_particles("mu+", nmu);  // mu+,e+,proton,pi+,Upsilon
    genp->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
        PHG4SimpleEventGenerator::Uniform,
        PHG4SimpleEventGenerator::Uniform);
    genp->set_vertex_distribution_mean(0.0, 0.0, target_coil_pos_z);
    genp->set_vertex_distribution_width(0.0, 0.0, 0.0);
    genp->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
    genp->set_vertex_size_parameters(0.0, 0.0);

    if(FMAGSTR>0)
      //genp->set_pxpypz_range(0,6, -6,6, 10,100);
      genp->set_pxpypz_range(-3,6, -3,3, 10,100);
    else
      //genp->set_pxpypz_range(-6,0, -6,6, 10,100);
      genp->set_pxpypz_range(-6,3, -3,3, 10,100);


    genp->Verbosity(0);
    se->registerSubsystem(genp);
  }

  if(gen_particle) {
    PHG4SimpleEventGenerator *genm = new PHG4SimpleEventGenerator("MUP");
    //genm->set_seed(123);
    genm->add_particles("mu-", nmu);  // mu+,e+,proton,pi+,Upsilon
    genm->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
        PHG4SimpleEventGenerator::Uniform,
        PHG4SimpleEventGenerator::Uniform);
    genm->set_vertex_distribution_mean(0.0, 0.0, target_coil_pos_z);
    genm->set_vertex_distribution_width(0.0, 0.0, 0.0);
    genm->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
    genm->set_vertex_size_parameters(0.0, 0.0);

    if(FMAGSTR>0)
      //genm->set_pxpypz_range(-6,0, -6,6, 10,100);
      genm->set_pxpypz_range(-6,3, -3,3, 10,100);
    else
      //genm->set_pxpypz_range(0,6, -6,6, 10,100);
      genm->set_pxpypz_range(-3,6, -3,3, 10,100);

    genm->Verbosity(0);
    se->registerSubsystem(genm);
  }

  // Fun4All G4 module
  PHG4Reco *g4Reco = new PHG4Reco();
  //PHG4Reco::G4Seed(123);
  //g4Reco->set_field(5.);
  g4Reco->set_field_map(
      jobopt_svc->m_fMagFile+" "+
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
  SetupSensitiveDetectors(g4Reco);

  se->registerSubsystem(g4Reco);

  // save truth info to the Node Tree
  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  // digitizer
  DPDigitizer *digitizer = new DPDigitizer("DPDigitizer", 0);
  digitizer->Verbosity(0);
  se->registerSubsystem(digitizer);

  // PatternDBGen
  PatternDBGen *dbgen = new PatternDBGen();
  dbgen->Verbosity(0);
  se->registerSubsystem(dbgen);

  // trakcing module
  gSystem->Load("libktracker.so");
  KalmanFastTrackingWrapper *ktracker = new KalmanFastTrackingWrapper();
  ktracker->Verbosity(0);
  ktracker->set_enable_event_reducer(true);
  ktracker->set_DS_level(0);
  se->registerSubsystem(ktracker);

  // evaluation module
  gSystem->Load("libmodule_example.so");
  TrkEval *trk_eval = new TrkEval();
  trk_eval->Verbosity(0);
  trk_eval->set_hit_container_choice("Vector");
  trk_eval->set_out_name("trk_eval.root");
  se->registerSubsystem(trk_eval);

  ///////////////////////////////////////////
  // Output
  ///////////////////////////////////////////

  // input - we need a dummy to drive the event loop
  if(read_hepmc) {
    Fun4AllHepMCInputManager *in = new Fun4AllHepMCInputManager("HEPMCIN");
    in->Verbosity(10);
    in->set_vertex_distribution_mean(0,0,target_coil_pos_z,0);
    se->registerInputManager(in);
    in->fileopen("hepmcout.txt");
  } else {
    Fun4AllInputManager *in = new Fun4AllDummyInputManager("DUMMY");
    se->registerInputManager(in);
  }

  //Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", "DST.root");
  //se->registerOutputManager(out);

  if(gen_pythia8 && !read_hepmc) {
    Fun4AllHepMCOutputManager *out = new Fun4AllHepMCOutputManager("HEPMCOUT", "hepmcout.txt");
    out->set_embedding_id(1);
    se->registerOutputManager(out);
  }

  if (nevent >= 0)
  {
    se->run(nevent);

    PHGeomUtility::ExportGeomtry(se->topNode(),"geom.root");

    // finish job - close and save output files
    se->End();
    se->PrintTimer();
    std::cout << "All done" << std::endl;

    // cleanup - delete the server and exit
    delete se;
    gSystem->Exit(0);
  } else { // TEve event display
#ifndef GRID
    gROOT->LoadMacro("EventDisplay.C");
    EventDisplay(nevent);
#endif
  } 

  return 0;
}

PHG4ParticleGun *get_gun(const char *name = "PGUN")
{
  Fun4AllServer *se = Fun4AllServer::instance();
  PHG4ParticleGun *pgun = (PHG4ParticleGun *) se->getSubsysReco(name);
  return pgun;
}
