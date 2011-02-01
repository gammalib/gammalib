# Filename: gammalib-init.sh
# Description: Bourne-shell flavor initialization for GammaLib.
#              Runs gammalib-setup to generate a sh script tailored
#              specifically to this user and GammaLib installation,
#              then source that.
# Author/Date: Juergen Knoedlseder, IRAP, February 1, 2011
#
if [ "x$GAMMALIB" = x ]; then 
  echo "gammalib-init.sh: ERROR -- set GAMMALIB before sourcing gammalib-init.sh"
elif [ -x "$GAMMALIB/bin/gammalib-setup" ]; then 
  gammalib_init=`$GAMMALIB/bin/gammalib-setup sh`
  if [ $? -eq 0 -a "x$gammalib_init" != x ]; then
    if [ -f "$gammalib_init" ]; then
      . $gammalib_init
    fi
    rm -f $gammalib_init
  fi
  unset gammalib_init
else
  echo "gammalib-init.sh: ERROR -- cannot execute $GAMMALIB/bin/gammalib-setup"
fi
