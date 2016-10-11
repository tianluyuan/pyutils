#!/bin/bash

# Taken from a script output from MCConfigClass.py (by robj)
# Modified for genie script 11/8/13 T.Y.

echo "cd into " ${TMPDIR}
cd ${TMPDIR}

echo "Running ND280"
runND280 -t ${TMPDIR} -c %CFG_PATH%

/usr/local/adm/bin/CUStageOut *genie_evtrate* %EVTRATE%
/usr/local/adm/bin/CUStageOut *.log %LOG%
