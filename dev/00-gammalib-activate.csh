#!/bin/csh

setenv GAMMALIB ${CONDA_PREFIX}
source ${GAMMALIB}/bin/gammalib-init.csh
unsetenv PYTHONPATH
unsetenv LD_LIBRARY_PATH
unsetenv DYLD_LIBRARY_PATH
