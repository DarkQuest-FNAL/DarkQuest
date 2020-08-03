#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TH1I.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TH2D.h>

#include "SRawEvent.h"

using namespace std;

int main(int argc, char *argv[])
{

  SRawEvent* rawEvent = new SRawEvent();
  TFile* dataFile = new TFile(argv[1], "READ");
  TTree* dataTree = (TTree*)dataFile->Get("save");
  dataTree->SetBranchAddress("rawEvent", &rawEvent);

  Int_t runID, spillID, eventID;
  Int_t triggerBits;
  Int_t nHits;
  Int_t detectorID[15000], elementID[15000];
  Double_t tdcTime[15000], driftDistance[15000], pos[15000];
  

  TFile* saveFile = new TFile(argv[2], "recreate");
  TTree* saveTree = new TTree("Events", "Events");

  //saveTree->Branch("runID", &runID);
  //saveTree->Branch("spillID", &spillID);
  //saveTree->Branch("eventID", &eventID);
  //saveTree->Branch("triggerBits", &triggerBits);
  saveTree->Branch("nHits", &nHits);
  saveTree->Branch("detectorID", detectorID, "detectorID[nHits]/I");
  saveTree->Branch("elementID", elementID, "elementID[nHits]/I");
  saveTree->Branch("tdcTime", tdcTime, "tdcTime[nHits]/D");
  saveTree->Branch("driftDistance", driftDistance, "driftDistance[nHits]/D");
  saveTree->Branch("pos", pos, "pos[nHits]/D");

  for(Int_t i = 0; i < dataTree->GetEntries(); ++i) {
      dataTree->GetEntry(i);
      if(i % 1000 == 0) cout << i << endl;

      // The NIM1 and NIM3 bits are random triggers. 
      int nim1TriggerMask = 32; 
      int nim3TriggerMask = 128;  

      if (rawEvent->getTriggerBits()>0 && (rawEvent->getTriggerBits() & (nim1TriggerMask|nim3TriggerMask) != 0)){
	runID = rawEvent->getRunID();
	spillID = rawEvent->getSpillID();
	eventID = rawEvent->getEventID();
	triggerBits = rawEvent->getTriggerBits();
	nHits = rawEvent->getNHitsAll();	

	for(Int_t k = 0; k < rawEvent->getNHitsAll(); ++k)
	  {
	    Hit h = rawEvent->getHit(k);
	    detectorID[k] = h.detectorID; 
	    elementID[k] = h.elementID;
	    pos[k] = h.pos;
	   
	    // detector ID refers to the detector number as seen here:
	    // st1-drift chambers| D0: 1-6, D1: 7-12
	    // st2-drift chambers| D2: 13-18
	    // st3-drift chambers| D3p: 19-24, D3m: 25-30
	    // h1-hodoscope: H1B/T: 31/22 H1L/R: 33/34
	    // h2-hodoscope: H2L/R: 35/36 H2B/T: 37/38
	    // h3-hodoscope: H3B/T: 39/40
	    // h4-hodoscope: H4Y1L/R: 41/42 H4Y2L/R: 43/44 H4B/T: 45/46
	    // proto-tubes: 47-54
	    // dp-stations: DP1:55-58, DP2:59-6 

	    // elementID refers to the hodoscope bar or drift chamber channel hit, this is detector dependent but can give us an idea of position

	    tdcTime[k] = h.tdcTime;
	    driftDistance[k] = h.driftDistance;

	  }
        if(nHits > 0) saveTree->Fill();
      }

      rawEvent->clear();
  }

  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  return EXIT_SUCCESS;
}
