#ifndef __CALIB_XT_V2_H__
#define __CALIB_XT_V2_H__
#include <fun4all/SubsysReco.h>
class SQHitVector;
class CalibParamXT;
class CalibParamInTimeTaiwan;

/// SubsysReco module to calibrate (i.e. update) the drift distance and the in-time flag of chamber and prop tube hits.
/**
 * It is a temporary version made from `e1039-core/online/decoder_maindaq/CalibXT.h`.
 * It should be removed once the original class is updated.
 */
class CalibXTv2: public SubsysReco {
  SQHitVector* m_vec_hit;
  CalibParamXT* m_cal_xt;
  CalibParamInTimeTaiwan* m_cal_int;

 public:
  CalibXTv2(const std::string &name = "CalibXTv2");
  virtual ~CalibXTv2();
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
};

#endif // __CALIB_XT_V2_H__
