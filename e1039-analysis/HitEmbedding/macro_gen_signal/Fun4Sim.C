#include <top/G4_Beamline.C>
#include <top/G4_Target.C>
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libSQPrimaryGen)
using namespace std;

int Fun4Sim(const int job_id=0, const int n_evt=0)
{
  recoConsts *rc = recoConsts::instance();
  Fun4AllServer *se = Fun4AllServer::instance();

  ///
  /// Global parameters
  ///
  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  rc->set_CharFlag("VTX_GEN_MATERIAL_MODE", "Target");

  ///
  /// Event generator
  ///
  const int gen_switch = 2;
  SQPrimaryParticleGen* sq_gen = new SQPrimaryParticleGen();
  switch (gen_switch) {
  case 1: // Drell-Yan: 5000 events = 12 hour
    sq_gen->set_massRange(1.5, 9.0);
    sq_gen->set_xfRange(0.00, 0.95);
    sq_gen->enableDrellYanGen();
    break;
  case 2: // J/psi: 5000 events = 20 hours
    sq_gen->set_xfRange(0.2, 1.0);
    sq_gen->enableJPsiGen();
    break;
  case 3: // psi'
    sq_gen->set_xfRange(0.2, 1.0); // Not tuned yet
    sq_gen->enablePsipGen();
    break;
  }
  se->registerSubsystem(sq_gen);

  ///
  /// Detector setting
  ///
  PHG4Reco *g4Reco = new PHG4Reco();
  g4Reco->set_field_map(
      rc->get_CharFlag("fMagFile")+" "+
      rc->get_CharFlag("kMagFile")+" "+
      Form("%f",FMAGSTR) + " " +
      Form("%f",KMAGSTR) + " " +
      "5.0",
      PHFieldConfig::RegionalConst);
  g4Reco->SetWorldSizeX(1000);
  g4Reco->SetWorldSizeY(1000);
  g4Reco->SetWorldSizeZ(5000);
  g4Reco->SetWorldShape("G4BOX");
  g4Reco->SetWorldMaterial("G4_AIR"); //G4_Galactic, G4_AIR
  g4Reco->SetPhysicsList("FTFP_BERT");

  SetupInsensitiveVolumes(g4Reco);
  SetupBeamline(g4Reco);
  SetupTarget(g4Reco);
  SetupSensitiveDetectors(g4Reco);

  se->registerSubsystem(g4Reco);

  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  /// digitizer
  SQDigitizer *digitizer = new SQDigitizer("DPDigitizer", 0);
  se->registerSubsystem(digitizer);

  /// Save only events that are in the geometric acceptance.
  SQGeomAcc* geom_acc = new SQGeomAcc();
  geom_acc->SetMuonMode(SQGeomAcc::PAIR);
  geom_acc->SetPlaneMode(SQGeomAcc::HODO_CHAM);
  geom_acc->SetNumOfH1EdgeElementsExcluded(4);
  se->registerSubsystem(geom_acc);

  // Make SQ nodes for truth info
  TruthNodeMaker* tnm = new TruthNodeMaker();
  tnm->SetJobID(job_id);
  se->registerSubsystem(tnm);

  ///
  /// Reconstruction
  ///
  SQReco* reco = new SQReco();
  reco->set_evt_reducer_opt("none");
  se->registerSubsystem(reco);

  VertexFit* vertexing = new VertexFit();
  vertexing->enable_fit_target_center();
  se->registerSubsystem(vertexing);

  ///
  /// Input, output and execution
  ///
  Fun4AllInputManager* man_in = new Fun4AllDummyInputManager("DUMMY");
  se->registerInputManager(man_in);

  Fun4AllDstOutputManager *man_out = new Fun4AllDstOutputManager("DSTOUT", "DST.root");
  se->registerOutputManager(man_out);
  man_out->AddNode("SQEvent");
  man_out->AddNode("SQHitVector");
  man_out->AddNode("SRecEvent");
  man_out->AddNode("SQMCEvent");
  man_out->AddNode("SQTruthTrackVector");
  man_out->AddNode("SQTruthDimuonVector");

  const bool count_only_good_events = true;
  se->run(n_evt, count_only_good_events);
  
  se->End();
  se->PrintTimer();
  rc->WriteToFile("recoConsts.tsv");
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
