#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <g4detectors/PHG4CollimatorSubsystem.h>
class SubsysReco;
R__LOAD_LIBRARY(libg4detectors)
#endif

void SetupEMCal(
  PHG4Reco* g4Reco,
  const std::string name = "EMCal",
  const double place_x = 0.,
  const double place_y = 0.,
  const double place_z = 0.,
  const int n_super_h = 6,
  const int n_super_v = 3,
  const int verbose   = 0.) 
{
  PHG4EMCalSubsystem* emcal = new PHG4EMCalSubsystem(name.c_str(), 0);
  emcal->SuperDetector(name.c_str());
  emcal->set_int_param("n_towers_x", n_super_h*6*2);   //1 super module = 6x6 module, 1 module = 2x2 tower
  emcal->set_int_param("n_towers_y", n_super_v*6*2);
  emcal->set_double_param("place_x", place_x);
  emcal->set_double_param("place_y", place_y);
  emcal->set_double_param("place_z", place_z);
  emcal->Verbosity(verbose);
  //emcal->set_int_param("absorberactive", 1);
  g4Reco->registerSubsystem(emcal);
 
  return;
}

