#ifndef _ANA_Module__H_
#define _ANA_Module__H_

#include <map>
#include <fun4all/SubsysReco.h>
#include <TString.h>
#include <TVector3.h>
#include <ktracker/SRawEvent.h>

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

  void set_output_filename(const TString& n) { saveName = n; }

private:
  int GetNodes(PHCompositeNode* topNode);
  void MakeTree();

  // Input
  SRawEvent* rawEvent;

  // Output
  TString saveName;
  TFile* saveFile;
  TTree* saveTree;

  Int_t runID, spillID, eventID;
  Int_t triggerBits;
  Int_t nHits;
  Int_t detectorID[15000], elementID[15000];
  Double_t tdcTime[15000], driftDistance[15000], pos[15000];

};

#endif
