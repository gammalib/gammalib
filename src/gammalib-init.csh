# Filename: gammalib-init.csh
# Description: C-shell flavor initialization for GammaLib
#              runs gammalib-setup to generate a temporary csh script
#              tailored specifically to this user and GammaLib
#              installation, then source that.
# Author/Date: Juergen Knoedlseder, IRAP, February 1, 2011
#
if(${?GAMMALIB} == 0) then 
  echo "gammalib-init.csh: ERROR -- set GAMMALIB before sourcing gammalib-init.csh"
else if(-x "$GAMMALIB/bin/gammalib-setup") then 
  set gammalib_init=`$GAMMALIB/bin/gammalib-setup csh`
  if($status == 0 && "x$gammalib_init" != x) then
    if(-f "$gammalib_init") then
      source $gammalib_init
    endif
    \rm -f $gammalib_init
  endif
  unset gammalib_init
else
  echo "gammalib-init.csh: ERROR -- cannot execute $GAMMALIB/bin/gammalib-setup"
endif
