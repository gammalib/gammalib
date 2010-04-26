#!/bin/sh

#autoreconf --force --install -I config -I m4
libtoolize --force \
&& aclocal \
&& automake --add-missing \
&& autoconf
