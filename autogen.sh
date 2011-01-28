#!/bin/sh

#autoreconf --force --install
#autoreconf --force --install --symlink \
#&& aclocal -I m4 \
#&& automake --add-missing \
#&& autoconf
aclocal -I m4
libtoolize --copy
autoconf
autoheader
automake --add-missing --copy

