"""
script to generate and submit condor jobs
"""
import os
import time

def GenerateExecutable(ijob):
    idir = "/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_10_6_0/src/PostCorrNTuple/forCombine_default/Fiducial/test0/commands/gof/"
    os.system("mkdir -p " + idir + "/tmp_" + str(ijob))
    jobname = f"{idir}/tmp_{ijob}/run_gof_{ijob}"
    with open(jobname + ".sh", "w") as bashscript:
        bashscript.write("#!/bin/bash\n")
        bashscript.write("\n")
        bashscript.write("ulimit -s unlimited\n")
        bashscript.write("set -e\n")
        bashscript.write("cd /afs/cern.ch/work/y/yofeng/public/WpT/higgsCombine/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit\n")
        bashscript.write("eval $(scramv1 runtime -sh);\n")
        bashscript.write("\n")
        bashscript.write(f"cd {idir}/\n")
        
        bashscript.write(f"combine -M GoodnessOfFit datacard.txt --algo=saturated -t 100 -s {ijob} --toysFreq -n _tmp_{ijob}\n")

    # change it to executable
    os.system("chmod +x " + jobname + ".sh")

    job_desc = """Universe = vanilla
Executable = {jobname}.sh
Log        = {jobname}.log
Output     = {jobname}.out
Error      = {jobname}.error
+JobFlavour = "longlunch"
queue 1\n
""".format(jobname = jobname)

    with open(jobname + ".condor", 'w') as outfile:
        outfile.write(job_desc)
        outfile.close()

    os.system("condor_submit " + jobname + ".condor")
    time.sleep(0.1)
    
    
if __name__ == "__main__":
   for ijob in range(100):
       GenerateExecutable(ijob) 