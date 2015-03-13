#!/bin/bash
cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubUpgrade/Working_Dir/620SLHC23/src/FastTiming/Generic_Analizer/macros
eval `scramv1 runtime -sh`
root -b -q .x VertexDeterminator_minuit.C
