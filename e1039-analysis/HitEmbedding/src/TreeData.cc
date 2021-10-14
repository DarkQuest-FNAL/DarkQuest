#include "TreeData.h"
using namespace std;

EventData::EventData()
  : job_id(0)
  , event_id(0)
  , trig_bits(0)
  , rfp01(0)
  , rfp00(0)
  , rfm01(0)
  , weight(1.0)
  , rec_stat(0)
  , n_h1x(0)
  , n_h2x(0)
  , n_h3x(0)
  , n_h4x(0)
  , n_d1 (0)
  , n_d2 (0)
  , n_d3 (0)
{
  ;
}

TrackData::TrackData() 
  : charge(0)
{
  ;
}
  
DimuonData::DimuonData() 
  : mass(0)
  , pT(0)
  , x1(0)
  , x2(0)
  , xF(0)
  , costh_cs(0)
  , phi_cs(0)
{
  ;
}
