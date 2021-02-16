#include <TSystem.h>

#include "G4_EMCal.C"

R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libPHPythia8)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libdptrigger)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libevt_filter)
R__LOAD_LIBRARY(libanamodule)

using namespace std;

int Fun4EMCal(const int nevent = 1)
{
  GeomSvc::UseDbSvc(true);
  GeomSvc* p_geomSvc = GeomSvc::instance();

  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(0);

  PHG4ParticleGun* gun = new PHG4ParticleGun("GUN");
  gun->set_name("e-");
  gun->set_vtx(0., 0., -300.);
  gun->set_mom(0., 0., 2.);
  se->registerSubsystem(gun);

  // Fun4All G4 module
  PHG4Reco* g4Reco = new PHG4Reco();
  g4Reco->SetWorldSizeX(1000);      //unit is cm
  g4Reco->SetWorldSizeY(1000);
  g4Reco->SetWorldSizeZ(1000);
  g4Reco->SetWorldShape("G4BOX");   // shape of our world - it is a box
  g4Reco->SetWorldMaterial("G4_AIR"); // this is what our world is filled with, eg. G4_Galactic, G4_AIR

  // sub-systems
  SetupEMCal(g4Reco);
  se->registerSubsystem(g4Reco);

  // save truth info to the Node Tree
  PHG4TruthSubsystem* truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  SQDigitizer* digitizer = new SQDigitizer("Digitizer", 0);
  digitizer->Verbosity(0);
  digitizer->registerEMCal("EMCal", 100);
  se->registerSubsystem(digitizer);

  AnaModule* ana = new AnaModule();
  ana->set_output_filename("ana.root");
  se->registerSubsystem(ana);

  // input - we need a dummy to drive the event loop
  Fun4AllInputManager* in = new Fun4AllDummyInputManager("DUMMY");
  se->registerInputManager(in);

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
