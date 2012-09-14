#!/bin/sh

if command -v libtoolize >/dev/null 2>&1; then
  libtoolize --copy
else
  if command -v glibtoolize >/dev/null 2>&1; then
    glibtoolize --copy
  fi
fi
aclocal -I m4
autoconf
autoheader
automake --add-missing --copy
