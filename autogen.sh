#!/bin/sh

#aclocal -I m4
if which libtoolize >/dev/null; then
  libtoolize --copy
else
  if which glibtoolize >/dev/null; then
    glibtoolize --copy
  fi
fi
aclocal -I m4
autoconf
autoheader
automake --add-missing --copy
