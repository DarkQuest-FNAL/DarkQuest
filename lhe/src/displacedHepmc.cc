#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <unistd.h>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

using namespace std;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <TSpline.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/IO_GenEvent.h"

using namespace HepMC;


struct stdhep_entry {
  int isthep;     /* status code */
  int idhep;      /* The particle id */
  int jmohep[2];  /* The position of the mother particle */
  double phep[5]; /* 4-Momentum, mass */
};

struct stdhep_event {
  stdhep_entry *aprime;
  stdhep_entry *postrack;
  stdhep_entry *negtrack;
};

// decay width calculations follow arXiv:1311.3870, Eq. 27 and 28
// all decay widths must be multiplied by epsilon^2
// decay width of A'->e+e- in units of GeV
double width_dielectron(double m_aprime) {
  double alpha = 1.0/137.036;
  double m_electron = 5.11e-3;
  double massratio_sq = m_electron*m_electron/(m_aprime*m_aprime);
  if (m_aprime>2*m_electron)
    return (alpha*m_aprime/3.0) * sqrt(1 - 4*massratio_sq) * (1 + 2*massratio_sq);
  else return 0.0;
}

// decay width of A'->mu+mu- in units of GeV
double width_dimuon(double m_aprime) {
  double alpha = 1.0/137.036;
  double m_muon = 105.66e-3;
  double massratio_sq = m_muon*m_muon/(m_aprime*m_aprime);
  if (m_aprime>2*m_muon)
    return (alpha*m_aprime/3.0) * sqrt(1 - 4*massratio_sq) * (1 + 2*massratio_sq);
  else return 0.0;
}

// decay width of A'->hadrons in units of GeV
double width_dihadron(TSpline * r_spline, double m_aprime) {
  double alpha = 1.0/137.036;
  double m_pion = 139.57e-3;
  if (m_aprime>2*m_pion)
    return (alpha*m_aprime/3.0) * r_spline->Eval(m_aprime);
  else return 0.0;
}

// total decay width in units of GeV
double width_total(TSpline * r_spline, double m_aprime) {
  return width_dielectron(m_aprime)+width_dimuon(m_aprime)+width_dihadron(r_spline,m_aprime);
}

// mean lifetime in units of cm
// must be divided by epsilon^2
double get_ctau(TSpline * r_spline, double m_aprime) {
  //1e9*elementary charge = conversion factor from GeV to joules
  return TMath::Ccgs()*TMath::Hbar()/(1e9*TMath::Qe()*width_total(r_spline,m_aprime));
}

// branching fraction to dimuon
double branching_dimuon(TSpline * r_spline, double m_aprime) {
  return width_dimuon(m_aprime)/width_total(r_spline,m_aprime);
}

// branching fraction to dielectron
double branching_dielectron(TSpline * r_spline, double m_aprime) {
  return width_dielectron(m_aprime)/width_total(r_spline,m_aprime);
}

int main(int argc,char** argv)
{
    int nevhep;             /* The event number */
    int rseed = 0;

    double eps = 1e-6;
    string mech = "";
    std::stringstream stream;    
    string epsStr = "";

    //int n_repeat = 2000; //number of times to sample the decay distribution for each input event
    int n_repeat = 2;
    double vx_production[3] = {0.0, 2.0, 50.0}; //beamspot at y=2 cm; guess z=50 cm for mean interaction position (dump face at 25 cm, interaction length 16.77 cm)
    float min_vz = 300.0;
    float max_vz = 800.0;

    float mass = -1.0;
    float ctau = 0.0;

    bool write_tree = true;
    bool ismuons = false;
    bool iselectrons = false;
    int c;

    int aflag = 0;
    int bflag = 0;
    char *cvalue = NULL;
    int index;

    if(argc<10) {
      printf("Wrong number of arguments, usage: <input stdhep filename> <output stdhep filename> mech muon/electron seed epsPow mass minVz maxVz");
    }

    mech = argv[3];
    if(strcmp(argv[4], "muon")==0) ismuons =true;
    else if(strcmp(argv[4],"electron")==0) iselectrons =true;
    else return(-1);
    rseed = atoi(argv[5]);
    eps = pow(10.0,atof(argv[6]));
    stream << atof(argv[6]);
    stream.precision(2);
    epsStr = stream.str();
    mass = atof(argv[7]);
    min_vz = atof(argv[8]);
    max_vz = atof(argv[9]);

    double e_cm[2000], r_ratio[2000];
    int n_rpoints = 0;
    FILE * r_file;
    r_file = fopen("data/rpp2017-hadronicrpp_page1001.dat","r");
    char line[1000];
    int numvals;
    while (fgets(line,1000,r_file)!=NULL) {
        numvals = sscanf(line," %lf %*f %*f %lf",&e_cm[n_rpoints], &r_ratio[n_rpoints]);
        if (numvals==2 && e_cm[n_rpoints]>e_cm[n_rpoints-1] && e_cm[n_rpoints]<10.0) { //reject repeated points with the same E, and only go up to 10 GeV
            n_rpoints++;
        }
    }
    //printf("loaded R data with %d points\n",n_rpoints);
    fclose(r_file);
    TSpline3 * r_spline = new TSpline3("R_ratio", e_cm, r_ratio, n_rpoints, "", 0.0, r_ratio[n_rpoints-1]);
    r_spline->SetName("R_ratio");
    if (mass>0) {
        ctau = get_ctau(r_spline,mass)/(eps*eps);
        //printf("mass=%f GeV, epsilon=%e, ctau=%f cm\n",mass, eps, ctau);
    }
    //if(ismuons){
    //  std::cout << "BR dimuon " << branching_dimuon(r_spline,mass) << std::endl;
    //}
    //if(iselectrons){
    //  std::cout << "BR dielectron "<< branching_dielectron(r_spline,mass) << std::endl;
    //}


    FILE * in_file;
    in_file = fopen(argv[optind],"r");

    //initialize the RNG
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();

    T = gsl_rng_mt19937;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,rseed);

    int ostream = 0;

    Double_t px0, py0, pz0;
    Double_t pz1, y1, ty1, x1_st1, tx1_st1, x1, tx1;
    Double_t pz2, y2, ty2, x2_st1, tx2_st1, x2, tx2;
    vector<stdhep_event> input_events;

    nevhep = 1;
    int n_accepted_events = 0;
    

    while (true) {

      
        char line[1000];
        bool found_event = false;
	
        while (fgets(line,1000,in_file)!=NULL) {
            if (strstr(line,"<event")!=NULL) {
                found_event = true;
                break;
            }
        }
        if (!found_event) {
            fclose(in_file);
            break;
        }

        int nup, idprup; //number of particles, process ID
        double xwgtup; //event weight
	
	int pdgID; //pdgID of decay
	if (ismuons){
	  pdgID = 13;
	}
	else if (iselectrons){
          pdgID= 11;
	}
	else{
	  break;
	}

	// read LHE
        struct stdhep_event temp_event = (struct stdhep_event) {0,0,0};
        while (true) {
            struct stdhep_entry *temp = new struct stdhep_entry;
            fgets(line,1000,in_file);
            char blah[1000];
            int n_tokens = sscanf(line,"%d %d %d %*d %lf %lf %lf %lf %lf",
				  &(temp->idhep),&(temp->isthep),&(temp->jmohep[0]),&(temp->phep[3]),&(temp->phep[0]),&(temp->phep[1]),&(temp->phep[2]),&(temp->phep[4]));
            if (n_tokens<8) break;
            switch (temp->isthep) {//translate between status conventions for HEPEUP (LHE) and HEPEVT (StdHep)
                case 1:
                case 2:
                    break;
                case -1:
                    temp->isthep = 3;
                    break;
                default:
                    temp->isthep = 0;
            }
            if (temp->isthep==2 && temp->idhep==666) {// intermediate particle, PDG ID 666
                if (temp_event.aprime) printf("WARNING: multiple A'\n");
                temp_event.aprime = temp;
                if (mass<0) { //only read mass on the first event (assume all events in a file have the same mass)
                    mass = temp->phep[4];
                    ctau = get_ctau(r_spline,mass)/(eps*eps);
                    printf("mass=%f GeV, epsilon=%e, ctau=%f cm\n",mass, eps, ctau);
                }
            }
            else if (temp->isthep==3 && temp->idhep==pdgID) {// final state particle. -
                if (temp_event.negtrack) printf("WARNING: multiple -\n");
                temp_event.negtrack = temp;
            }
            else if (temp->isthep==3 && temp->idhep==(pdgID*-1)) {// final state particle, +
                if (temp_event.postrack) printf("WARNING: multiple +\n");
                temp_event.postrack = temp;
            }
            else delete temp;
        }
        if (temp_event.aprime && temp_event.negtrack && temp_event.postrack) {
            input_events.push_back(temp_event);
        }
        else printf("WARNING: missing A' decays\n");
        //if (nevhep%1000==0) printf("%d\n",nevhep);
        nevhep++;
    }


    // writing HepMC file
    stream.str("");
    stream << std::fixed << std::setprecision(2) << mass;
    std::string massStr = stream.str();
    std::string lepStr = "Muons";
    if(iselectrons) lepStr = "Electrons";
    string outFile = "displaced_Aprime_"+lepStr+"/"+mech+"_"+massStr+"_z"+std::to_string((int)min_vz)+"_"+std::to_string((int)max_vz)+"_eps_"+epsStr+".txt";
    //WriterAscii output_file(outFile);
    IO_GenEvent output_file(outFile);
    int n_extra_repeats = 0;

    do {
      
        for (vector<stdhep_event>::iterator event = input_events.begin(); event!=input_events.end();++event) {
            double gamma, beta;
            gamma = event->aprime->phep[3]/event->aprime->phep[4];
            beta = sqrt(1.0-pow(gamma,-2.0));
            double decay_length = beta*gamma*ctau;
            double p = 0.0;
            for (int j=0;j<3;j++) p += event->aprime->phep[j]*event->aprime->phep[j];
            p = sqrt(p);

            px0 = event->aprime->phep[0];
            py0 = event->aprime->phep[1];
            pz0 = event->aprime->phep[2];
	    
	    double px1 = event->postrack->phep[0];
	    double py1 = event->postrack->phep[1];
	    double pz1 = event->postrack->phep[2];
	    double pt1 = event->postrack->phep[3];
	    
	    double px2 = event->negtrack->phep[0];
	    double py2 = event->negtrack->phep[1];
	    double pz2 = event->negtrack->phep[2];
	    double pt2 = event->negtrack->phep[3];
	    
            for (int i=0;i<n_repeat;i++) {
	      n_extra_repeats++;

	      double vx[4];
	      
	      double vtx_displacement = gsl_ran_exponential(r,decay_length);
	      for (int j=0;j<3;j++) vx[j] = vtx_displacement*event->aprime->phep[j]/p + vx_production[j];
	      if (vx[2]<min_vz || vx[2]>max_vz) continue;

	      // this does not include fmag or kmag kick
	      n_accepted_events++;
	      
	      // create HepMC evt
	      //GenEvent evt = GenEvent(Units::GEV, Units::CM);
	      GenEvent* evt = new GenEvent(Units::GEV, Units::CM);
	      //GenEvent evt = GenEvent(Units::GEV, Units::MM);
	      //evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);
	      evt->set_event_number(n_accepted_events);
	      // create A' particle 
	      // px      py        pz       e     pdgid status  
	      //GenParticlePtr paprime = std::make_shared<GenParticle>( FourVector(px0,py0,pz0,event->aprime->phep[3]), event->aprime->idhep, event->aprime->isthep);
	      GenParticle* paprime = new GenParticle(FourVector(px0,py0,pz0,event->aprime->phep[3]), event->aprime->idhep, event->aprime->isthep);
	      // create postrack particle
	      //GenParticlePtr ppostrack = std::make_shared<GenParticle>( FourVector(px1,py1,pz1,pt1), event->postrack->idhep, event->postrack->isthep);
	      GenParticle* ppostrack = new GenParticle(FourVector(px1,py1,pz1,pt1), event->postrack->idhep, event->postrack->isthep);
	      ppostrack->set_status(1);
	      // create negtrack particle
	      //GenParticlePtr pnegtrack = std::make_shared<GenParticle>( FourVector(px2,py2,pz2,pt2), event->negtrack->idhep, event->negtrack->isthep);
	      GenParticle* pnegtrack = new GenParticle(FourVector(px2,py2,pz2,pt2), event->negtrack->idhep, event->negtrack->isthep);
	      pnegtrack->set_status(1);
	      
	      // create A' vertex
	      // need to know where the vertex is (vx)
	      vx[3] = sqrt(vx[0]*vx[0] + vx[1]*vx[1] + vx[2]*vx[2] + event->aprime->phep[4]*event->aprime->phep[4]);
	      
	      //std::cout <<  vx[0] << vx[1] << vx[2] << vx[3] << std::endl;
	      //GenVertexPtr vaprime = std::make_shared<GenVertex>( FourVector(vx[0], vx[1], vx[2], vx[3]) );
	      GenVertex* vaprime = new GenVertex(FourVector(vx[0], vx[1], vx[2], vx[3]) );
	      evt->add_vertex( vaprime );
	      vaprime->add_particle_in( paprime );
	      vaprime->add_particle_out( ppostrack );
	      vaprime->add_particle_out( pnegtrack );
	      
	      // write file
	      output_file.write_event(evt);
	      //if (n_accepted_events%1000==0) printf("accepted %d\n",n_accepted_events);
            }
        }
    } while (n_accepted_events<10000 and n_extra_repeats <= 10000);//Go for 10k events
    //printf("%d events accepted by cuts after %d samples, %d samples of %d events\n",n_accepted_events,n_extra_repeats,n_repeat,int(input_events.size()));
    printf("%d %d %f %e %f %f %s \n",n_accepted_events,n_extra_repeats,mass,eps,min_vz,max_vz,mech.c_str());
}
