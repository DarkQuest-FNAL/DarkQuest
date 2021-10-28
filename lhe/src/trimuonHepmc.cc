#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#include <TMath.h>
#include <TRandom3.h>
#include <TLorentzVector.h>

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/IO_GenEvent.h"

using namespace HepMC;

struct stdhep_entry {
  int isthep;      /* status code */
  int idhep;       /* The particle id */
  int mom1, mom2;  /* The position of the mother particle */
  int color1, color2; 
  double phep[5];  /* 4-Momentum, mass */
  double lifetime;
  double spin;
};

struct stdhep_event {
  stdhep_entry *particle;
  stdhep_entry *scalar;
  stdhep_entry *posparticle;
  stdhep_entry *negparticle;
};

inline void SetPxPyPzM(TLorentzVector& t, double x, double y, double z, double m){
  t.SetPxPyPzE(x,y,z, sqrt(x*x+y*y+z*z+m*m) );
}

int main(int argc,char** argv)
{

  double mass(-1), lifetime(-1);

  std::vector<stdhep_event> input_events;

  // read lhe
  ifstream myfile("/seaquest/users/cmantill/DarkQuest/lhe/data/trimuon/signalProcess0.5MS0gS1.lhe");
  if (myfile.is_open()) {
    std::string line;
    while (getline(myfile, line)){
      if (line.empty() || line[0] == '#')
	continue;
      string eventline = "<event";
      if (line.find(eventline) != string::npos) {
	struct stdhep_event temp_event = (struct stdhep_event) {0,0,0};                                                                                                                                                                                                       
	std::string newline;
	// ignore line right after event
	getline(myfile, newline);
	bool found_particle = false;
	while (getline(myfile, newline)){
	  struct stdhep_entry *temp = new struct stdhep_entry;                                                                                                                                                                                           
	  std::stringstream ss(newline);
	  if( ss >> temp->idhep >> temp->isthep >> temp->mom1 >> temp->mom2 >> temp->color1 >> temp->color2 >> temp->phep[0] >> temp->phep[1] >> temp->phep[2] >> temp->phep[3] >> temp->phep[4] >> temp->lifetime >> temp->spin ) {
	    found_particle = true;
	  }
	  else break;
	  if (temp->isthep==2 && temp->idhep==1001) {
	    temp_event.scalar = temp; 
	    if(mass==-1) mass = temp->phep[4];
	    if(lifetime==-1) lifetime = temp->lifetime;
	  }
	  else if (temp->isthep==1 && temp->idhep==13 && temp->mom1==3) {
	    temp_event.negparticle = temp;
	  }
	  else if (temp->isthep==1 && temp->idhep==-13 && temp->mom1==3) {
	    temp_event.posparticle = temp;
	  }
	  else if (temp->isthep==1 && temp->idhep==13 && temp->mom1<3 && temp->mom1>0) {
	    temp_event.particle = temp;
	  }
          else if (temp->isthep==-1 && temp->idhep==13){
	    std::cout << " momentum of incoming muon " << temp->phep[3] << std::endl;
	  }
	  else delete temp;
	}
	if (temp_event.scalar && temp_event.negparticle && temp_event.posparticle && temp_event.particle) {
	  input_events.push_back(temp_event);
        }
      }
    }
  }
  myfile.close();

  cout << "n input events " << input_events.size() << endl;
  
  // assuming they all come from the target or dump area?
  // production vertex: beamspot at y=2 cm; guess z=50 cm for mean interaction position (dump face at 25 cm, interaction length 16.77 cm)
  double vx_production[3] = {0.0, 2.0, 50.0}; 
  
  // write hepmc
  string outFile = "trimuon_0.5MS0gS1.hepmc";
  IO_GenEvent output_file(outFile);

  //int n_hepmc = 10e3;
  int event_number = 0;
  for (vector<stdhep_event>::iterator event = input_events.begin(); event!=input_events.end();++event) {

    // create HepMC evt                                                                                                                                                                                                                                                     
    GenEvent* evt = new GenEvent(Units::GEV, Units::CM);
    evt->set_event_number(event_number);

    TLorentzVector mu_vtx;
    SetPxPyPzM(mu_vtx, event->particle->phep[0]+vx_production[0],event->particle->phep[1]+vx_production[1],event->particle->phep[2]+vx_production[2], 0.105);

    std::cout << " event particle pz " << event->particle->phep[2] << " scalar pz " << event->scalar->phep[3] << " negtrack " << event->negparticle->phep[3] << " pos " << event->posparticle->phep[3] << std::endl;
    GenParticle* particle = new GenParticle(FourVector(event->particle->phep[0],event->particle->phep[1],event->particle->phep[2],event->particle->phep[3]), event->particle->idhep, 1);
    GenParticle* pscalar = new GenParticle(FourVector(event->scalar->phep[0],event->scalar->phep[1],event->scalar->phep[2],event->scalar->phep[3]), event->scalar->idhep, event->scalar->isthep);

    GenVertex* vpart = new GenVertex(FourVector(mu_vtx.Px(), mu_vtx.Py(), mu_vtx.Pz(), mu_vtx.E()) );
    evt->add_vertex( vpart );
    vpart->add_particle_out( pscalar );
    vpart->add_particle_out( particle );

    GenParticle* pposparticle = new GenParticle(FourVector(event->posparticle->phep[0],event->posparticle->phep[1],event->posparticle->phep[2],event->posparticle->phep[3]), event->posparticle->idhep, event->posparticle->isthep);
    pposparticle->set_status(1);
    GenParticle* pnegparticle =new GenParticle(FourVector(event->negparticle->phep[0],event->negparticle->phep[1],event->negparticle->phep[2],event->negparticle->phep[3]), event->negparticle->idhep, event->negparticle->isthep);
    pnegparticle->set_status(1);

    // create vertex
    TLorentzVector vtx;
    SetPxPyPzM(vtx, mu_vtx.Px(), mu_vtx.Py(), mu_vtx.Pz(), mass);

    GenVertex* vscalar = new GenVertex(FourVector(vtx.Px(), vtx.Py(), vtx.Pz(), vtx.E()) );
    evt->add_vertex( vscalar );
    vscalar->add_particle_out( pposparticle );
    vscalar->add_particle_out( pnegparticle );

    // write file
    output_file.write_event(evt);
    event_number++;
  }

}
