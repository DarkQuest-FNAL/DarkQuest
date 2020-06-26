R__LOAD_LIBRARY(libsim_eval)
#include <iostream>
#include <sstream>
using namespace std;

int Fun4SimDst(std::string f_dst="Brem_0.36_z500_600_eps_-5.60")
{
  gSystem->Load("libsim_eval.so");

  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(1);
  Fun4AllInputManager *in = new Fun4AllDstInputManager("SimDst");
  se->registerInputManager(in);
  //stringstream ssin; ssin << f_dst << "_dst.root";
  stringstream ssin; ssin << "/pnfs/e906/scratch/users/cmantill/SimHits/Brem_0.36_z500_600_eps_-5.6-10000/1/out/Brem_0.36_z500_600_eps_-5.60_dst.root";
  in->fileopen(ssin.str().c_str());

  SimEval *sim_eval = new SimEval();
  //sim_eval->Verbosity(99);                                                                                                                                                        
  sim_eval->set_hit_container_choice("Vector");
  //stringstream ssout; ssout << "sim_eval_" << f_dst << ".root";
  stringstream ssout; ssout << "sim_eval_" << "Brem_0.36_z500_600_eps_-5.60" << ".root";
  sim_eval->set_out_name(ssout.str().c_str());
  se->registerSubsystem(sim_eval);

  se->run();
  se->End();
  se->PrintTimer();
  delete se;
  return 0;
}
