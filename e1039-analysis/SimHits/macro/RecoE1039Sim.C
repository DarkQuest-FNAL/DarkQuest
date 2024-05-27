#include "G4_EMCal.C"
#include <top/G4_Beamline.C>
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
#include <top/G4_Target.C>

R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libPHPythia8)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libdptrigger)
R__LOAD_LIBRARY(libevt_filter)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libSQPrimaryGen)
R__LOAD_LIBRARY(libsim_ana)

#include <iostream>
#include <sstream>
#include <string>
using namespace std;

/*
 * Macro used to analyze simulation in SpinQuest
 * isim = 1 to run on Aprime to dimuon signal
 * isim = 2 to run on Aprime to dielectron signal
 * isim = 3 to run on single particle gun:
 * - igun = 1: muon
 * - igun = 2: electron
 * - igun = 3: positron
 * - igun = 4: proton
 * - igun = 5: gamma
 * - igun = 6: pi+
 * - igun = 7: pi-
 * - igun = 8: klong
 * - igun = 9: pi0
 * isim = 4 to run on DY to dimuon sample generated with Pythia
 * isim = 5 to run on J/psi to dimuon sample generated with Pythia
 * isim = 6 to run on cosmic sample
 * isim = 7 to run on trimuon sample
 * isim = 8 to run on background-simulation only
 *
 * do_displaced_tracking: to run tracking - if particle is e+/e-/gamma then electron_tracking is set to true by default too
 * do_analysis: to produce the analysis ntuple
 * for Aprime signal, is_displaced to always set to True
 */

int RecoE1039Sim(const int nevents = 200,
                 const int isim = 1,
                 const int igun = 0,
                 const double zvertex = -300, // target_coil_pos_z
                 const bool do_displaced_tracking = true,
                 const bool do_analysis = true,
                 bool run_pileup = false,
                 std::string input_file = "Brem_0.550000_z500_600_eps_-5.4",
                 std::string input_path = "$HepMC_DIR/",
                 std::string out_file = "output.root",
                 std::string out_path = "./",
                 std::string pileup_file = "/pnfs/e1039/persistent/users/apun/bkg_study/e1039pythiaGen_26Oct21/10_bkge1039_pythia_wshielding_100M.root",
                 const int verbosity = 0)
{
    // input simulation
    bool do_aprime_muon{false}, do_aprime_electron{false};
    bool do_gun{false};
    bool do_dy{false}, do_jpsi{false}, do_cosmic{false}, do_trimuon{false};
    bool do_bkgOnly{false};

    // tracking options
    bool electron_tracking{false};

    // gun options
    std::string particle_name;

    switch (isim)
    {
    case 1:
        do_aprime_muon = true;
        std::cout << " DISPLACED A' TO MUONS " << std::endl;
        break;
    case 2:
        do_aprime_electron = true;
        electron_tracking = true;
        std::cout << " DISPLACED A' TO ELECTRONS " << std::endl;
        break;
    case 3:
        do_gun = true;
        switch (igun)
        {
        case 1: // muon gun
            particle_name = "mu-";
            break;
        case 2: // electron gun
            particle_name = "e-";
            electron_tracking = true;
            break;
        case 3: // positron gun
            particle_name = "e+";
            electron_tracking = true;
            break;
        case 4: // proton gun
            particle_name = "proton";
            break;
        case 5: // photon gun
            particle_name = "gamma";
            electron_tracking = true;
            break;
        case 6: // pi+ gun
            particle_name = "pi+";
            break;
        case 7: // pi- gun
            particle_name = "pi-";
            break;
        case 8: // klong gun
            particle_name = "kaon0L";
            break;
        case 9: // pi0 gun
            particle_name = "pi0";
            break;
        }
        std::cout << " " << particle_name << " GUN " << std::endl;
        break;
    case 4:
        do_dy = true;
        std::cout << " DO DY " << std::endl;
        break;
    case 5:
        do_jpsi = true;
        std::cout << " DO J/PSI " << std::endl;
        break;
    case 6:
        do_cosmic = true;
        std::cout << " DO COSMIC " << std::endl;
        break;
    case 7:
        do_trimuon = true;
        std::cout << " DO TRIMUON " << std::endl;
        break;
    case 8:
        do_bkgOnly = true;
        run_pileup = true;
        std::cout << " DO Background Only. Set the run_pileup to true " << std::endl;
        break;
    }

    /** Verbosity (https://github.com/E1039-Collaboration/e1039-core/blob/master/framework/fun4all/Fun4AllBase.h#L33-L55)
     *  the verbosity of different modules can also be modified separately for debugging
     */
    bool isDEBUG = false;
    if (verbosity > 0)
        isDEBUG = true;

    // legacy rec container
    const bool legacy_rec_container = true; // false is for e1039 format

    // save dst file
    const bool save_dst = false;

    // setup detectors in SpinQuest
    const bool do_collimator = true;
    const bool do_target = true;
    const bool do_shielding = true;
    const bool do_fmag = true;
    const bool do_kmag = true;
    const bool do_absorber = true;
    const bool do_dphodo = true;
    const bool do_station1DC = false; // station-1 drift chamber should be turned off by default
    const bool doEMCal = false;       // emcal turned off (for SpinQuest)

    // SpinQuest constants
    const double target_coil_pos_z = -300;
    const double target_l = 7.9;                   // cm
    const double target_z = (7.9 - target_l) / 2.; // cm
    const int use_g4steps = 1;
    // const double FMAGSTR = -1.054;
    // const double KMAGSTR = -0.951;
    const double FMAGSTR = -1.044;
    const double KMAGSTR = 1.025;

    // SpinQuest reco constants
    recoConsts *rc = recoConsts::instance();
    rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
    rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
    if (doEMCal)
    {
        rc->set_CharFlag(
            "AlignmentMille",
            "$DIR_TOP/data/alignment/align_mille.txt"); // alignment file needed for EMCAL
    }
    rc->set_CharFlag("fMagFile",
                     "$E1039_RESOURCE/geometry/magnetic_fields/tab.Fmag");
    rc->set_CharFlag("kMagFile",
                     "$E1039_RESOURCE/geometry/magnetic_fields/tab.Kmag");

    if (do_cosmic)
    {
        rc->init("cosmic");
        rc->set_BoolFlag("COARSE_MODE", true);
        rc->set_DoubleFlag("KMAGSTR", 0.);
        rc->set_DoubleFlag("FMAGSTR", 0.);
    }

    if (electron_tracking)
    {
        rc->set_BoolFlag(
            "TRACK_ELECTRONS",
            true); // track electrons by eliminating certain muon hit requirements
    }
    if (do_displaced_tracking)
    {
        rc->set_BoolFlag(
            "TRACK_DISPLACED",
            true); // track displaced particles by removing backwards extrapolation in st2+3 to st1 tracklet connection
    }

    if (isDEBUG)
    {
        rc->Print();
    }

    // geometry information
    GeomSvc::UseDbSvc(true);
    GeomSvc *geom_svc = GeomSvc::instance();
    if (isDEBUG)
    {
        std::cout << "print geometry information" << std::endl;
        geom_svc->printWirePosition();
        std::cout << " align printing " << std::endl;
        geom_svc->printAlignPar();
        std::cout << " table printing" << std::endl;
        geom_svc->printTable();
        std::cout << "done geometry printing" << std::endl;
    }

    // make the Server
    Fun4AllServer *se = Fun4AllServer::instance();
    se->Verbosity(verbosity);

    // input to the simulation
    if (do_aprime_muon)
    { // aprime to displaced muons
        HepMCNodeReader *hr = new HepMCNodeReader();
        hr->set_particle_filter_on(true);
        hr->insert_particle_filter_pid(13); // filter muons
        hr->insert_particle_filter_pid(13 * -1);
        hr->Verbosity(verbosity);
        se->registerSubsystem(hr);
    }
    else if (do_aprime_electron)
    { // aprime to displaced electrons
        HepMCNodeReader *hr = new HepMCNodeReader();
        hr->set_particle_filter_on(true);
        hr->insert_particle_filter_pid(11); // filter electrons
        hr->insert_particle_filter_pid(11 * -1);
        hr->Verbosity(verbosity);
        se->registerSubsystem(hr);
    }
    else if (do_trimuon)
    {
        HepMCNodeReader *hr = new HepMCNodeReader();
        hr->Verbosity(verbosity);
        se->registerSubsystem(hr);
    }
    else if (do_gun)
    { // single particle gun
        PHG4SimpleEventGenerator *genp = new PHG4SimpleEventGenerator("PARTICLEGUN");
        genp->add_particles(particle_name.c_str(), 1);

        genp->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                               PHG4SimpleEventGenerator::Uniform,
                                               PHG4SimpleEventGenerator::Uniform);
        genp->set_vertex_distribution_mean(10.0, 10.0, zvertex); // to set after FMAG: zvertex: 520
        genp->set_vertex_distribution_width(10.0, 10.0, 0.0);    // for protons set to 10.0 in z?
        genp->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
        genp->set_vertex_size_parameters(0.0, 0.0);

        genp->set_pxpypz_range(-.15, .15, -.15, .15, 10., 100.);

        genp->Verbosity(verbosity);
        se->registerSubsystem(genp);
    }
    else if (do_dy or do_jpsi)
    {
        PHPythia8 *pythia8 = new PHPythia8();
        pythia8->Verbosity(verbosity);
        if (do_dy)
            pythia8->set_config_file("$DIR_TOP/data/pythiaconfig/phpythia8_DY.cfg");
        else
            pythia8->set_config_file("$DIR_TOP/data/pythiaconfig/phpythia8_Jpsi.cfg");
        pythia8->set_vertex_distribution_mean(0.0, 0.0, zvertex, 0);
        pythia8->set_embedding_id(1);
        se->registerSubsystem(pythia8);

        pythia8->set_trigger_AND();
        PHPy8ParticleTrigger *trigger_mup = new PHPy8ParticleTrigger();
        trigger_mup->AddParticles("-13");
        trigger_mup->SetPzHighLow(120, 30);
        pythia8->register_trigger(trigger_mup);

        PHPy8ParticleTrigger *trigger_mum = new PHPy8ParticleTrigger();
        trigger_mum->AddParticles("13");
        trigger_mum->SetPzHighLow(120, 30);
        pythia8->register_trigger(trigger_mum);

        HepMCNodeReader *hr = new HepMCNodeReader();
        hr->set_particle_filter_on(true);
        hr->insert_particle_filter_pid(13);
        hr->insert_particle_filter_pid(-13);
        se->registerSubsystem(hr);
    }
    else if (do_cosmic)
    {
        SQCosmicGen *cosmicGen = new SQCosmicGen();
        se->registerSubsystem(cosmicGen);
    }
    else if (do_bkgOnly)
    {
        // no operation needed here
        ;
    }
    else
    {
        std::cout << " No input! " << std::endl;
        return 0;
    }

    // pileup
    if (run_pileup)
    {
        SQPileupGen *extgen = new SQPileupGen();
        // function
        TF1 *intensity_profile = new TF1("intensity_profile", "[0]*pow(x,[1])*exp(-[2]*x+exp(-[3]*x))+[4]", 0, 5000);
        intensity_profile->SetParameter(0, 6.35);
        intensity_profile->SetParameter(1, 1.38);
        intensity_profile->SetParameter(2, 4.9e-3);
        intensity_profile->SetParameter(3, 4.7e-3);
        intensity_profile->SetParameter(4, 178.8);
        extgen->set_beam_intensity_profile(intensity_profile);
        extgen->setExtInputFile(pileup_file);
        se->registerSubsystem(extgen);
    }

    // fun4All G4 module
    PHG4Reco *g4Reco = new PHG4Reco();
    g4Reco->set_field_map(
        rc->get_CharFlag("fMagFile") + " " + rc->get_CharFlag("kMagFile") + " " +
            Form("%f", FMAGSTR) + " " + Form("%f", KMAGSTR) + " " + "5.0",
        PHFieldConfig::RegionalConst);
    // size of the world - every detector has to fit in here
    g4Reco->SetWorldSizeX(1000);
    g4Reco->SetWorldSizeY(1000);
    g4Reco->SetWorldSizeZ(5000);
    // shape of our world - it is a tube
    g4Reco->SetWorldShape("G4BOX");
    // this is what our world is filled with
    g4Reco->SetWorldMaterial("G4_AIR"); // G4_Galactic, G4_AIR
    // G4 Physics list to use
    g4Reco->SetPhysicsList("FTFP_BERT");

    // setup detectors
    SetupInsensitiveVolumes(g4Reco, do_shielding, do_fmag, do_kmag,
                            do_absorber); // insensitive volumes

    SetupBeamline(
        g4Reco, do_collimator,
        target_coil_pos_z -
            302.36); // collimator, targer and shielding between target and FMag
    if (do_target)
    {
        SetupTarget(g4Reco, target_coil_pos_z, target_l, target_z, use_g4steps);
    }
    // sensitive elements of the spectrometer
    SetupSensitiveDetectors(g4Reco, do_dphodo, do_station1DC, "SQ_ArCO2",
                            "SQ_Scintillator", 0);
    if (doEMCal)
    {
        SetupEMCal(g4Reco, "EMCal", 0., 0., 1930.);
    }
    se->registerSubsystem(g4Reco);

    // save G4 truth info to the Node Tree
    PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
    g4Reco->registerSubsystem(truth);

    // digitizer
    SQDigitizer *digitizer = new SQDigitizer("DPDigitizer", 0);
    digitizer->Verbosity(verbosity);
    digitizer->set_enable_st1dc(do_station1DC); // these two lines need to be in
                                                // sync with the parameters used
    digitizer->set_enable_dphodo(
        do_dphodo); // in the SetupSensitiveVolumes() function call above
    if (doEMCal)
    {
        digitizer->registerEMCal("EMCal", 100);
    }
    se->registerSubsystem(digitizer);

    bool do_acceptance = false;
    if (do_acceptance)
    {
        // if turned on, only the events passing the geometric acceptance will be saved
        // https://github.com/E1039-Collaboration/e1039-core/blob/master/simulation/g4dst/SQGeomAcc.h
        SQGeomAcc *geom_acc = new SQGeomAcc();
        geom_acc->SetMuonMode(SQGeomAcc::PAIR);
        geom_acc->SetPlaneMode(SQGeomAcc::HODO_CHAM);
        // geom_acc->SetNumOfH1EdgeElementsExcluded(4);
        se->registerSubsystem(geom_acc);
    }

    // tracking module
    SQReco *reco = new SQReco();
    reco->Verbosity(verbosity);
    reco->set_legacy_rec_container(legacy_rec_container);
    // reco->set_geom_file_name("support/geom.root"); // not needed as it's created on the fly
    reco->set_enable_KF(true);         // Kalman filter not needed for the track finding, disabling KF saves a lot of initialization time
    reco->setInputTy(SQReco::E1039);   // options are SQReco::E906 and SQReco::E1039
    reco->setFitterTy(SQReco::KFREF);  // not relevant for the track finding, options are SQReco::KFREF and SQReco::LEGACY
    reco->set_evt_reducer_opt("none"); // if not provided, event reducer will be using JobOptsSvc to intialize; to turn off, set it to "none", for normal tracking, set to something like "aoc"
    // reco->set_evt_reducer_opt("r");              // if not provided, event reducer will be using JobOptsSvc to intialize; to turn off, set it to "none", for normal tracking, set to something like "aoc"
    reco->set_enable_eval(true);           // set to true to generate evaluation file which includes final track candidates
    reco->set_eval_file_name("eval.root"); // evaluation filename
    reco->set_enable_eval_dst(true);       // set to true to include final track candidates in the DST tree
    reco->add_eval_list(3);                // include back partial tracks in eval tree for debuging
    reco->add_eval_list(2);                // include station-3+/- in eval tree for debuging
    reco->add_eval_list(1);                // include station-2 in eval tree for debugging
    se->registerSubsystem(reco);

    // truth node maker after tracking
    TruthNodeMaker *truthMaker = new TruthNodeMaker();
    truthMaker->set_legacy_rec_container(legacy_rec_container);
    if (do_aprime_muon or do_aprime_electron)
    {
        truthMaker->set_m_process_type(3); // set process type to 3 (A' -> di lepton) since we only have a 3 particle process instead of 0+1->2+3
    }
    if (do_trimuon)
    {
        truthMaker->set_m_process_type(5);
    }
    truthMaker->Verbosity(verbosity);
    se->registerSubsystem(truthMaker);

    // trigger emulator
    // needs TruthNodeMaker to associate the trigger to SQEvent
    DPTriggerAnalyzer *dptrigger = new DPTriggerAnalyzer();
    dptrigger->set_road_set_file_name("$E1039_RESOURCE/trigger/trigger_67.txt");
    dptrigger->set_dproad_set_file_name("$DIR_TOP/data/trigger/DPTrigger_road16.txt");
    dptrigger->Verbosity(verbosity);
    se->registerSubsystem(dptrigger);

    // event filter
    EvtFilter *evt_filter = new EvtFilter();
    evt_filter->Verbosity(verbosity);
    evt_filter->set_trigger_req(1 << 5);
    // se->registerSubsystem(evt_filter);

    // truth vertexing
    SQTruthVertexing *truthVtx = new SQTruthVertexing();
    truthVtx->set_legacy_rec_container(legacy_rec_container);
    truthVtx->set_vtx_smearing(50.); // smear the truth z_vertex to mimic resolution effect, default is 0.
    // se->registerSubsystem(truthVtx);

    // vertexing for dimuon information
    if (legacy_rec_container)
    {
        VertexFit *vertexing = new VertexFit();
        vertexing->Verbosity(verbosity);
        se->registerSubsystem(vertexing);
    }

    // analysis module
    // gSystem->Load("libsim_ana.so");
    SimAna *sim_ana = new SimAna();
    sim_ana->Verbosity(verbosity);
    std::string ofile = out_path + out_file;
    sim_ana->set_out_name(ofile);
    sim_ana->set_legacy_rec_container(legacy_rec_container);
    sim_ana->save_secondaries(false); // set to true to save secondaries
    if (do_analysis)
    {
        se->registerSubsystem(sim_ana);
    }

    // input
    if (do_aprime_muon or do_aprime_electron)
    {
        // use hepmc input
        Fun4AllHepMCInputManager *in = new Fun4AllHepMCInputManager("HEPMCIN");
        se->registerInputManager(in);
        stringstream ssin;
        ssin << input_path << input_file << ".txt";
        std::cout << "Aprime Input path " << ssin.str().c_str() << std::endl;
        in->fileopen(gSystem->ExpandPathName(ssin.str().c_str()));
        in->Verbosity(verbosity);
        se->registerInputManager(in);
    }
    else if (do_trimuon)
    {
        Fun4AllHepMCInputManager *in = new Fun4AllHepMCInputManager("HEPMCIN");
        se->registerInputManager(in);
        stringstream ssin;
        // ssin << input_path << input_file << ".hepmc";
        ssin << input_path << input_file; // I changed the run_trimuon script to include the .hepmc postscript
        std::cout << "Trimuon Input path " << ssin.str().c_str() << std::endl;
        in->fileopen(gSystem->ExpandPathName(ssin.str().c_str()));
        in->Verbosity(verbosity);
        se->registerInputManager(in);
    }
    else
    {
        // need a dummy input to drive the event loop
        Fun4AllInputManager *in = new Fun4AllDummyInputManager("DUMMY");
        in->Verbosity(verbosity);
        se->registerInputManager(in);
    }

    // output (DST file)

    std::string dstfile = ofile;
    dstfile.resize(dstfile.size() - 5); // remove root from ending
    dstfile.append("_DST.root");

    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", dstfile.c_str());
    // these classes are needed for hit embedding
    out->AddNode("SQEvent");
    out->AddNode("SQHitVector");
    out->AddNode("SRecEvent");
    out->AddNode("SQMCEvent");
    out->AddNode("SQTruthTrackVector");
    out->AddNode("SQTruthDimuonVector");
    // add these to save same content of analysis ntuples
    // out->AddNode("SQRecTrackVector");
    // out->AddNode("SQRecDimuonVector");
    // out->AddNode("SQRecSt3TrackletVector");
    out->AddNode("G4TruthInfo");
    out->AddNode("PHG4HitContainer");
    out->AddNode("PHG4TruthInfoContainer");
    if (save_dst)
    {
        se->registerOutputManager(out);
    }

    se->run(nevents);

    // export the geometry
    // PHGeomUtility::ExportGeomtry(se->topNode(),"geom.root");

    // finish job - close and save output files
    se->End();
    se->PrintTimer();

    // cleanup - delete the server and exit
    delete se;
    gSystem->Exit(0);
    return 0;
}
