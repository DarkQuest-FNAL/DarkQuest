#include <TSystem.h>

R__LOAD_LIBRARY(libinterface_main)
R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libanamodule)

int DSTToTree(const int nEvents = 18518)
{
  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(1000);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("SimDst");
  in->Verbosity(10);
  se->registerInputManager(in);
  in->fileopen("result.root");

  AnaModule* ana = new AnaModule();
  ana->set_output_filename("ana.root");
  se->registerSubsystem(ana);

  se->run(nEvents);
  se->End();
  delete se;
  return 0;
}
