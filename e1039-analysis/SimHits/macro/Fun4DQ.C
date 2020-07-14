#include <TSystem.h>

#include <top/G4_Beamline.C>
#include <top/G4_Target.C>
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
#include "G4_EMCal.C"

R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libPHPythia8)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libdptrigger)
R__LOAD_LIBRARY(libevt_filter)
R__LOAD_LIBRARY(libanamodule)

using namespace std;

int Fun4DQ(const int nevent = 10,
	   std::string ifile = "Brem_1.04_z500_600_eps_-6.4")
{
  const bool simple = false;

  const bool do_collimator = true;
  const bool do_target     = true;
  const bool do_shielding  = true;
  const bool do_fmag       = true;
  const bool do_kmag       = true;
  const bool do_absorber   = true;

  const double collimator_pos_z = -602.36;
  const double target_coil_pos_z = -300.;
  const double target_l = 7.9; //cm
  const double target_z = (7.9-target_l)/2.; //cm

  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;

  recoConsts *rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  rc->Print();

  JobOptsSvc* jobopt_svc = JobOptsSvc::instance();
  jobopt_svc->init("run7_sim.opts");

  GeomSvc::UseDbSvc(true);
  GeomSvc* geom_svc = GeomSvc::instance();

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  if(simple){
    PHG4SimpleEventGenerator* genp = new PHG4SimpleEventGenerator("E+"); 
    genp->set_seed(123);
    genp->add_particles("e+", 1);  // mu+,e+,proton,pi+,Upsilon
    genp->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform, PHG4SimpleEventGenerator::Uniform, PHG4SimpleEventGenerator::Uniform);
    genp->set_vertex_distribution_mean(0.0, 0.0, 550.);
    genp->set_vertex_distribution_width(0.0, 0.0, 0.0);
    genp->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
    genp->set_vertex_size_parameters(0.0, 0.0);
    genp->set_pxpypz_range(-6., 6., -3. ,3., 10., 100.);
    se->registerSubsystem(genp);
  }
  else{
    int idLep=11;
    HepMCNodeReader *hr = new HepMCNodeReader();
    hr->set_particle_filter_on(true);
    hr->insert_particle_filter_pid(idLep);
    hr->insert_particle_filter_pid(idLep*-1);
    se->registerSubsystem(hr);
  }

  // Fun4All G4 module
  PHG4Reco* g4Reco = new PHG4Reco();
  //PHG4Reco::G4Seed(123);
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
  g4Reco->SetWorldShape("G4BOX"); // shape of our world - it is a box
  g4Reco->SetWorldMaterial("G4_AIR");  // // this is what our world is filled with G4_Galactic, G4_AIR
  g4Reco->SetPhysicsList("FTFP_BERT"); // Geant4 Physics list to use

  SetupBeamline(g4Reco, do_collimator, collimator_pos_z);
  SetupTarget(g4Reco, target_coil_pos_z, target_l, target_z, 1, 0);
  SetupInsensitiveVolumes(g4Reco, do_shielding, do_fmag, do_kmag, do_absorber);
  SetupSensitiveDetectors(g4Reco);
  SetupEMCal(g4Reco, "EMCal", 0., -110., 1930.);
  se->registerSubsystem(g4Reco);

  // save truth info to the Node Tree
  PHG4TruthSubsystem* truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  // digitizer
  SQDigitizer* digitizer = new SQDigitizer("Digitizer", 0);
  //digitizer->Verbosity(99);
  digitizer->registerEMCal("EMCal", 100);
  se->registerSubsystem(digitizer);

  // analyzer module
  AnaModule* ana = new AnaModule();
  ana->set_output_filename("ana.root");
  se->registerSubsystem(ana);

  // input - we need a dummy to drive the event loop
  if(simple){
    Fun4AllInputManager *in = new Fun4AllDummyInputManager("DUMMY");
    se->registerInputManager(in);
  }
  else{
    Fun4AllHepMCInputManager *in = new Fun4AllHepMCInputManager("HEPMCIN");
    se->registerInputManager(in);
    in->Verbosity(10);
    stringstream ssin; ssin << ifile << ".txt";
    in->fileopen(gSystem->ExpandPathName(ssin.str().c_str()));
  }

  // DST output manager, tunred off to save disk by default
  Fun4AllDstOutputManager* out = new Fun4AllDstOutputManager("DSTOUT", "DST.root");
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
