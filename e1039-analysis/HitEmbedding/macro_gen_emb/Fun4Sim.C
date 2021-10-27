#include <top/G4_Beamline.C>
#include <top/G4_Target.C>
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libana_embedding)
using namespace std;

int Fun4Sim(const int nevent = 10)
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

  ///
  /// Event generator
  ///
  const double target_coil_pos_z = -300;
  PHG4SimpleEventGenerator *genp = new PHG4SimpleEventGenerator("MUP");
  genp->add_particles("mu+", 2);
  genp->set_vertex_distribution_mean(0.0, 0.0, target_coil_pos_z);
  if (FMAGSTR > 0) genp->set_pxpypz_range(+1, +3,  -0.5, 0.5,  40, 60);
  else             genp->set_pxpypz_range(-3, -1,  -0.5, 0.5,  40, 60);
  se->registerSubsystem(genp);

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

  // save truth info to the Node Tree
  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  // digitizer
  SQDigitizer *digitizer = new SQDigitizer("DPDigitizer", 0);
  se->registerSubsystem(digitizer);

  /// Save only events that are in the geometric acceptance.
  SQGeomAcc* geom_acc = new SQGeomAcc();
  geom_acc->SetMuonMode(SQGeomAcc::SINGLE);
  geom_acc->SetPlaneMode(SQGeomAcc::HODO_CHAM);
  geom_acc->SetNumOfH1EdgeElementsExcluded(4);
  se->registerSubsystem(geom_acc);

  // Make SQ nodes for truth info
  se->registerSubsystem(new TruthNodeMaker());

  se->registerSubsystem(new GenEmbeddingData());

  ///
  /// Input, output and run
  ///
  Fun4AllInputManager *in = new Fun4AllDummyInputManager("DUMMY");
  se->registerInputManager(in);

  se->run(nevent);
  
  PHGeomUtility::ExportGeomtry(se->topNode(),"geom.root");
  se->End();
  se->PrintTimer();
  rc->WriteToFile("recoConsts.tsv");
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
