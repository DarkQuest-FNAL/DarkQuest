import glob
import subprocess
import os
import multiprocessing
import time

cores = 12
nevent = 10000 #input(f"How many events of input file to process? (1000 is a good number to start) ")


#files = glob.glob("/cms/data/seaquest/users/cmantill/DarkQuest/lhe/output/displaced_Aprime_Muons_z500-600/Eta*.txt")
#files = glob.glob("/cms/data/seaquest/users/cmantill/DarkQuest/lhe/output/displaced_Aprime_Muons_z500-600/Brem*.txt")
files = glob.glob("/cms/data/seaquest/users/dsperka/DarkQuest/output/displaced_Aprime_Muons/Brem*.txt")

files = [f.split("/")[-1][:-4] for f in files]
for f in files:
    mech = f.split("_")[0]
    mass = f.split("_")[1]
    logeps = f.split("_")[5]
    if (float(mass)>1.5):
        files.remove(f)
        continue
    if (float(logeps)>-5):
        files.remove(f)
        continue
    if (mech=="Brem" and float(mass)<0.6 and float(logeps)<-7):
        files.remove(f)
        continue
    if (mech=="Eta" and float(logeps)<-7):
        files.remove(f)
        continue    
    
print(len(files))


def process(f):
        input_file_size = os.path.getsize(f'/cms/data/seaquest/users/dsperka/DarkQuest/output/displaced_Aprime_Muons/{f}.txt')
        if input_file_size == 0:
                print(f'\t\t\t{f}\t	Input file is empty',flush=True)
                os.system(f'touch /home/dsperka/data/{f}.no_input_placeholder')
                return 'error'
        command = f'mkdir /tmp/processing/{f} && cd /tmp/processing/{f} && cp /home/dsperka/DarkQuest-FNAL/DarkQuest/e1039-analysis/SimHits/macro/DPTrigger_road16.txt . && root -b -q /home/dsperka/DarkQuest-FNAL/DarkQuest/e1039-analysis/SimHits/macro/RecoE1039Sim.C\({nevent},1,0,500.0,true,true,false,\\"{f}\\",\\"/cms/data/seaquest/users/dsperka/DarkQuest/output/displaced_Aprime_Muons/\\"\) > command_output.txt 2>command_errors.txt ; mv output.root /cms/data/seaquest/users/dsperka/DarkQuest/16May2024/{f}.root'
        #the last chained command is added with ; because currently root crashes right at the end
        start_time = time.time()
        return_code = os.system(command)
        runtime = time.time()-start_time
        print(f'Time = {round(runtime,2)}\terror={return_code}\t{f}',flush=True)
        return runtime
#if input("About to delete /tmp/processing/* /home/caspian/data/* is that ok? (y/n) ") != 'y':
#	quit()
subprocess.call("rm -r /tmp/processing".split(),shell=False)
subprocess.call("mkdir /tmp/processing".split(),shell=False)
pool = multiprocessing.Pool(processes=cores)
print('Starting', flush=True)
times = pool.map(process,files)
print(times,flush=True)

quit()

