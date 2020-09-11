#include <iomanip>
#include <TFile.h>
#include <TTree.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include "AnaModule.h"

int AnaModule::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  eventID = 0;
  MakeTree();
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::process_event(PHCompositeNode* topNode)
{

  nHits = rawEvent->getNHitsAll();

  // The NIM1 and NIM3 bits are random triggers. 
  int nim1TriggerMask = 32; 
  int nim3TriggerMask = 128;  

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

  if (rawEvent->getTriggerBits()>0 && (rawEvent->getTriggerBits() & (nim1TriggerMask|nim3TriggerMask) != 0)){
    saveTree->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::End(PHCompositeNode* topNode)
{
  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::GetNodes(PHCompositeNode* topNode)
{
  rawEvent = findNode::getClass<SRawEvent>(topNode, "SRawEvent");
  if(!rawEvent)
    {
      rawEvent = nullptr;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

void AnaModule::MakeTree()
{
  saveFile = new TFile(saveName, "RECREATE");

  saveTree = new TTree("rawEvents", "Raw Tree Created by AnaModule");
  saveTree->Branch("nHits", &nHits);
  saveTree->Branch("detectorID", detectorID, "detectorID[nHits]/I");
  saveTree->Branch("elementID", elementID, "elementID[nHits]/I");
  saveTree->Branch("tdcTime", tdcTime, "tdcTime[nHits]/D");
  saveTree->Branch("driftDistance", driftDistance, "driftDistance[nHits]/D");
  saveTree->Branch("pos", pos, "pos[nHits]/D");
}
