#ifndef _FILTER_E906_NIM3__H_
#define _FILTER_E906_NIM3__H_
#include <vector>
#include <fun4all/SubsysReco.h>
class SQEvent;

/// A SubsysReco module to filter E906 NIM3 events.
/** 
 * It accepts only good NIM3 (i.e. random RF) events.
 * A list of good spills is read from a list file.
 */
class FilterE906Nim3: public SubsysReco {
  std::string m_list_name;
  std::vector<int> m_list_spill_ok;
  SQEvent* mi_evt;

  int m_n_evt_all;
  int m_n_evt_spill;
  int m_n_evt_nim3;

 public:
  FilterE906Nim3(const std::string name="FilterE906Nim3");
  virtual ~FilterE906Nim3() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void        SetSpillListName(const std::string name) { m_list_name = name; }
  std::string GetSpillListName() const { return m_list_name; }
};

#endif // _FILTER_E906_NIM3__H_
