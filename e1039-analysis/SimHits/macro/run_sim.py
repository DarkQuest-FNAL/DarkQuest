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
    parser.add_argument('-n', dest="n", type=int, help="Number of events")
    parser.add_argument('--sim', dest="sim", type=str, help="Simulation type", choices=sim_to_isim.keys(), required=True)
    parser.add_argument('--gun', dest="gun", type=str, help="Particle type for gun", choices=gun_to_igun.keys())
    parser.add_argument('--zvertex', dest="zvertex", type=float, help="Position of z vertex", default="300.")
    parser.add_argument('--displaced', action=store_true, help="Apply displaced tracking in reconstruction")
    parser.add_argument('--inputfile', dest="inputfile", help="Input HEP MC file")
    parser.add_argument('--inputpath', dest="inputpath", help="input directory", default="/seaquest/users/cmantill/DarkQuest/lhe/output/displaced_Aprime_Muons_z500-600/")
    parser.add_argument('--outfile', dest="outfile", help="Name of output file")
    parser.add_argument('--outpath', dest="outpath", help="Output directory")
    parser.add_argument('-v', '--verbose', default = False, action='store_true', dest="verbose", help="Verbose output, print all data to screen")

    args = parser.parse_args()
    
    # set options from args
    gun = gun_to_igun[args.gun] if args.gun >0 else 0
    
    cmd_options = "nevent=%i"%args.n
    cmd_options += ",isim=%i"%sim_to_isim[args.sim]
    cmd_options += ",igun=%i"%gun
    cmd_options += ",zvertex=%.2f"%args.zvertex
    if args.displaced:
        cmd_options += ",do_displaced_tracking=true"
    else:
        cmd_options += ",do_displaced_tracking=false"
    cmd_options += ",do_analysis=true"
    cmd_options += ",input_file=%s"%args.inputfile
    cmd_options += ",input_path=%s"%args.inputpath
    cmd_options += ",out_file=%s"%args.outfile
    cmd_options += ",out_path=%s"%args.outpath
    if args.verbose:
        cmd_options += ",verbosity=1"
    else:
        cmd_options += ",verbosity=0"

    # start process
    env = os.environ.copy()
    from subprocess import PIPE, Popen
    # command
    cmd = 'root -b -q RecoE1039Sim.C\(%s\)'%cmd_options
    # where to execute command
    cwd = './'

    proc = Popen(cmd, cwd=cwd, stdout=PIPE, universal_newlines=True, env=env)

