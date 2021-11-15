import argparse
import os

sim_to_isim = {'aprime-muon': 1,
               'aprime-electron': 2,
               'gun': 3,
               'DY': 4,
               'J-psi': 5,
               'cosmic': 6,
               'trimuon': 7,
               }
gun_to_igun = {'muon': 1,
               'electron': 2,
               'positron': 3,
               'proton': 4,
               'gamma': 5,
               'pi+': 6,
               'pi-': 7,
               'klong': 8,
               }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run Reco E1039 Sim ROOT script with options')
    parser.add_argument('-n', dest="n", type=int, help="Number of events", required=True)
    parser.add_argument('--sim', dest="sim", type=str, help="Simulation type", choices=sim_to_isim.keys(), required=True)
    parser.add_argument('--gun', dest="gun", type=str, help="Particle type for gun", choices=gun_to_igun.keys())
    parser.add_argument('--zvertex', dest="zvertex", type=float, help="Position of z vertex", default="300.")
    parser.add_argument('--displaced', action='store_true', default=False, help="Apply displaced tracking in reconstruction")
    parser.add_argument('--pileup', action='store_true', default=False, help="Mix pileup particles")
    parser.add_argument('--inputfile', default="", help="Input HEP MC file")
    parser.add_argument('--inputpath', default="", help="input directory")
    parser.add_argument('--outfile', dest="outfile", default="output.root", help="Name of output file")
    parser.add_argument('--outpath', dest="outpath", default="./", help="Output directory")
    parser.add_argument('--pudir', dest="pudir", default="/pnfs/e1039/persistent/users/apun/bkg_study/e1039pythiaGen_26Oct21/", help="Name of directory with pileup files")
    parser.add_argument('-v', '--verbose', default = False, action='store_true', dest="verbose", help="Verbose output, print all data to screen")

    args = parser.parse_args()
    
    # set options from args
    gun = gun_to_igun[args.gun] if args.gun >0 else 0
    
    cmd_options = "%i"%args.n # number of events
    cmd_options += ",%i"%sim_to_isim[args.sim] # isim
    cmd_options += ",%i"%gun # igun
    cmd_options += ",%.2f"%args.zvertex # z vertex position
    if args.displaced:
        cmd_options += ",true" # do displaced tracking
    else:
        cmd_options += ",false" # not do displaced tracking
    cmd_options += ",true" # do analysis ntuple
    if args.pileup:
        cmd_options += ",true" # mix pileup
    else:
        cmd_options += ",false"
    cmd_options += ',\\"%s\\"'%args.inputfile # input file (for strings ROOT has a weird formatting - requires strings to be within \)
    cmd_options += ',\\"%s\\"'%args.inputpath # input path
    cmd_options += ',\\"%s\\"'%args.outfile # output file
    cmd_options += ',\\"%s\\"'%args.outpath # output path

    # pileup: pick file from directory
    if args.pileup:
        import random
        pufile = args.pudir + "/" + random.choice( os.listdir(args.pudir) )
    else:
        pufile = "/pnfs/e1039/persistent/users/apun/bkg_study/e1039pythiaGen_26Oct21/10_bkge1039_pythia_wshielding_100M.root"
    cmd_options += ',\\"%s\\"'%pufile # pileup file

    if args.verbose:
        cmd_options += ",1" # verbosity
    else:
        cmd_options += ",0" # no verbosity

    # command
    cmd = 'root -b -q RecoE1039Sim.C\(%s\)'%cmd_options
    # where to execute command
    cwd = os.getcwd()
    
    print(cmd)
    #os.system(cmd)
