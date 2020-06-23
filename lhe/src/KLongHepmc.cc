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

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TF1.h>

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/IO_GenEvent.h"

using namespace HepMC;


inline void SetPxPyPzM(TLorentzVector& t, double x, double y, double z, double m){
    t.SetPxPyPzE(x,y,z, sqrt(x*x+y*y+z*z+m*m) );
}
void Decay2(TRandom* r, TLorentzVector m, TLorentzVector& d1, TLorentzVector& d2){
    double mdaught = 0.13957;
    double pdaught = sqrt(pow(m.M(),2)-pow(2*mdaught,2))/2.;
    double x,y,z;
    r->Sphere(x,y,z, pdaught);
    const TVector3 v(x,y,z);
    d1.SetVectM(v,mdaught);
    d2.SetVectM(-v,mdaught);
    d1.Boost( m.BoostVector() );
    d2.Boost( m.BoostVector() );
}

int main(int argc,char** argv)
{

    int n_hepmc = 10e3;

    double vx_prod[3] = {0.0, 0.0, 0.0};
    double vx_dump[3] = {0.0, 0.0, 500.0};

    // get kinematics
    TFile* f_kin = new TFile("../data/klong/kinematics.root","read");
    TH2D* h_kin = (TH2D*) f_kin->Get("klong");
    if (!h_kin){
        cout << "failed to retrieve kinematics file. did you generate it?" << endl;
        exit(0);
    }

    // writing HepMC file
    string outFile = "klong.hepmc";
    IO_GenEvent output_file(outFile);

    int ei=0;
    double pt,pz,phi, decay_length;
    double m = 0.497611; // pdg
    double ctau = 1534.; // lifetime is 5.116e-8 s
    TF1 *f_exp = new TF1("f_exp","exp(-x/[0])",500,600);
    TRandom3* rand = new TRandom3(2020);
    
    // kL decay vertex
    TLorentzVector kl, d1, d2, d3, vtx;
    double vz;

    while(ei < n_hepmc){

        h_kin->GetRandom2(pt,pz);
        phi = rand->Uniform(0,2*TMath::Pi());
        SetPxPyPzM(kl, pt*sin(phi), pt*cos(phi), pz, m);
        f_exp->SetParameter(0, kl.Beta()*kl.Gamma()*ctau);

        vz = f_exp->GetRandom();
        SetPxPyPzM(vtx, vz*kl.Px()/kl.Pz(), vz*kl.Px()/kl.Pz(), vz, m);
      	
        // create HepMC evt
        GenEvent* evt = new GenEvent(Units::GEV, Units::CM);
        evt->set_event_number(ei);
	
        // create KL particle and daughters
        // px      py        pz       e     pdgid status  
        GenParticle* gen_kl = new GenParticle(FourVector(kl.Px(), kl.Py(), kl.Pz(), kl.E()),130,2);
        // Decay2(rand,kl,d1,d2);
        // GenParticle* gen_d1 = new GenParticle(FourVector(d1.Px(), d1.Py(), d1.Pz(), d1.E()), 211,1);
        // GenParticle* gen_d2 = new GenParticle(FourVector(d2.Px(), d2.Py(), d2.Pz(), d2.E()),-211,1);

        // create KL vertex
        GenVertex* vaprime = new GenVertex(FourVector(vtx.Px(), vtx.Py(), vtx.Pz(), vtx.E()) );
        evt->add_vertex( vaprime );
        vaprime->add_particle_in( gen_kl );
        vaprime->add_particle_out( gen_kl );
        // vaprime->add_particle_out( gen_d1 );
        // vaprime->add_particle_out( gen_d2 );
	
        // write file
        output_file.write_event(evt);
        ei++;
        if (ei % 1000 ==0) cout << "processed " << ei << endl;
    }// generated all events      
}
