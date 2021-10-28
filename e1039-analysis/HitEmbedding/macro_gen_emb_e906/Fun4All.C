//R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libana_embedding)
using namespace std;

int Fun4All(const int run, const string fname, const int n_evt=0)
{
  recoConsts* rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", run);

  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllSRawEventInputManager* in = new Fun4AllSRawEventInputManager("SRawEventIM");
  in->Verbosity(0);
  in->enable_E1039_translation();
  se->registerInputManager(in);

  se->registerSubsystem(new FilterE906Nim3());
  se->registerSubsystem(new CalibXTv2());
  se->registerSubsystem(new GenEmbeddingData());

  in->fileopen(fname);
  se->run(n_evt);
  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
