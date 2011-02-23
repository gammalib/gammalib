#!/bin/sh

aclocal -I m4
if which libtoolize >/dev/null; then
  libtoolize --copy
else
  if which glibtoolize >/dev/null; then
    glibtoolize --copy
  fi
fi
autoconf
autoheader
automake --add-missing --copy
