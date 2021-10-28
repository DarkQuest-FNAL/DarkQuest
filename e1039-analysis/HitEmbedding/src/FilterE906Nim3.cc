#include <fstream>
#include <algorithm>
#include <phool/getClass.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <interface_main/SQEvent.h>
#include "FilterE906Nim3.h"
using namespace std;

FilterE906Nim3::FilterE906Nim3(const std::string name)
  : SubsysReco(name)
  , m_list_name("list_spill_good.txt")
  , m_n_evt_all  (0)
  , m_n_evt_spill(0)
  , m_n_evt_nim3 (0)
{
  ;
}

int FilterE906Nim3::Init(PHCompositeNode* topNode)
{
  cout << "FilterE906Nim3::Init():  Read the good spill list from '" << m_list_name << "'.\n";
  ifstream ifs(m_list_name.c_str());
  if (! ifs.is_open()) {
    cout << "!!ERROR!!  Cannot open the file.  Abort." << endl;
    exit(1);
  }
  int sp;
  while (ifs >> sp) m_list_spill_ok.push_back(sp);
  ifs.close();
  if (m_list_spill_ok.size() == 0) {
    cout << "!!ERROR!!  No good spill was found.  Abort." << endl;
    exit(1);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int FilterE906Nim3::InitRun(PHCompositeNode* topNode)
{
  mi_evt = findNode::getClass<SQEvent>(topNode, "SQEvent");
  if (!mi_evt) return Fun4AllReturnCodes::ABORTEVENT;
  return Fun4AllReturnCodes::EVENT_OK;
}

int FilterE906Nim3::process_event(PHCompositeNode* topNode)
{
  m_n_evt_all++;
  int spill_id = mi_evt->get_spill_id();
  if (find(m_list_spill_ok.begin(), m_list_spill_ok.end(), spill_id) == m_list_spill_ok.end()) return Fun4AllReturnCodes::ABORTEVENT;
  m_n_evt_spill++;
  if (! mi_evt->get_trigger(SQEvent::NIM3)) return Fun4AllReturnCodes::ABORTEVENT;
  m_n_evt_nim3++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int FilterE906Nim3::End(PHCompositeNode* topNode)
{
  ofstream ofs("stat.txt");
  ofs << m_n_evt_all << "\n" << m_n_evt_spill << "\n" << m_n_evt_nim3 << "\n";
  ofs.close();
  return Fun4AllReturnCodes::EVENT_OK;
}
