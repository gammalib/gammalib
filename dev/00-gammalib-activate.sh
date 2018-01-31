#!/bin/sh

export GAMMALIB=${CONDA_PREFIX}
source ${GAMMALIB}/bin/gammalib-init.sh
unset PYTHONPATH
unset LD_LIBRARY_PATH
unset DYLD_LIBRARY_PATH
