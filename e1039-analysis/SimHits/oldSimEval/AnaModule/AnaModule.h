#ifndef _ANA_Module__H_
#define _ANA_Module__H_

#include <map>
#include <fun4all/SubsysReco.h>
#include <TString.h>

class TFile;
class TTree;
class SQHitVector;

class AnaModule: public SubsysReco 
{
public:
  AnaModule() {;}
  virtual ~AnaModule() {;}

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  void set_output_filename(const TString& n);

private:
  int GetNodes(PHCompositeNode* topNode);
  void MakeTree();

  // Input
  SQHitVector* hitVector;

  // Output
  TString saveName;
  TFile* saveFile;
  TTree* saveTree;

  int eventID;
  int nHits;
  double edep[2000];
  int hID[2000];
  int vID[2000];
};

#endif
