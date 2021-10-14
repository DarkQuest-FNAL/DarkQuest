#include <iomanip>
#include <TGraphErrors.h>
#include <interface_main/SQParamDeco.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/recoConsts.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
#include <geom_svc/CalibParamXT.h>
#include <geom_svc/CalibParamInTimeTaiwan.h>
#include "CalibXTv2.h"
using namespace std;

CalibXTv2::CalibXTv2(const std::string& name)
  : SubsysReco(name)
  , m_vec_hit(0)
  , m_cal_xt (0)
  , m_cal_int(0)
{
  ;
}

CalibXTv2::~CalibXTv2()
{
  if (m_cal_xt ) delete m_cal_xt ;
  if (m_cal_int) delete m_cal_int;
}

int CalibXTv2::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CalibXTv2::InitRun(PHCompositeNode* topNode)
{
  SQRun*       run_header = findNode::getClass<SQRun      >(topNode, "SQRun");
  SQParamDeco* param_deco = findNode::getClass<SQParamDeco>(topNode, "SQParamDeco");
  m_vec_hit               = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!run_header || !m_vec_hit) return Fun4AllReturnCodes::ABORTEVENT;

  if (! m_cal_xt) m_cal_xt = new CalibParamXT();
  m_cal_xt->SetMapIDbyDB(run_header->get_run_id());
  m_cal_xt->ReadFromDB();

  if (! m_cal_int) m_cal_int = new CalibParamInTimeTaiwan();
  m_cal_int->SetMapIDbyDB(run_header->get_run_id());
  m_cal_int->ReadFromDB();

  if (param_deco) {
    param_deco->set_variable(m_cal_xt ->GetParamID(), m_cal_xt ->GetMapID());
    param_deco->set_variable(m_cal_int->GetParamID(), m_cal_int->GetMapID());
  }

  recoConsts* rc = recoConsts::instance();
  rc->set_CharFlag(m_cal_xt ->GetParamID(), m_cal_xt ->GetMapID());
  rc->set_CharFlag(m_cal_int->GetParamID(), m_cal_int->GetMapID());

  if (Verbosity() > 0) {
    cout << Name() << ": " << m_cal_xt ->GetParamID() << " = " << m_cal_xt ->GetMapID() << "\n"
         << Name() << ": " << m_cal_int->GetParamID() << " = " << m_cal_int->GetMapID() << "\n";
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CalibXTv2::process_event(PHCompositeNode* topNode)
{
  GeomSvc* geom = GeomSvc::instance();
  for (SQHitVector::Iter it = m_vec_hit->begin(); it != m_vec_hit->end(); it++) {
    SQHit* hit = *it;
    int det = hit->get_detector_id();
    int ele = hit->get_element_id();
    hit->set_pos(geom->getMeasurement(det, ele));

    TGraphErrors* gr_t2x;
    TGraphErrors* gr_t2dx;
    if (m_cal_xt->Find(det, gr_t2x, gr_t2dx)) {
      double center, width;
      //if (! m_cal_int->Find(det, ele, center, width)) {
      if (! m_cal_int->Find(det, 0, center, width)) { // until e1039-core is updated
        cerr << "  WARNING:  Cannot find the in-time parameter for det=" << det << " ele=" << ele << " in CalibXTv2.\n";
        continue;
      }
      float t1 = center - width / 2;
      float t0 = center + width / 2;
      float   tdc_time = hit->get_tdc_time();
      float drift_time = t0 - tdc_time;
      float drift_dist = gr_t2x->Eval(drift_time);
      if (drift_dist < 0) drift_dist = 0;
      hit->set_in_time(t1 <= tdc_time && tdc_time <= t0);
      hit->set_drift_distance(drift_dist);
    } else { // Hodoscope
      hit->set_drift_distance(0);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CalibXTv2::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
