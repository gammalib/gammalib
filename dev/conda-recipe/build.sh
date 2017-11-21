#!/bin/bash

# Set Mac OS X minimum deployment target
if [ `uname` == Darwin ]; then
  export CXXFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
  export LDFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
fi

# Configure GammaLib
./configure --prefix="$PREFIX"

# Build
make

# Install
make install
