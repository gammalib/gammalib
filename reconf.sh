#!/bin/sh
rm -f config.cache
rm -f acconfig.h
echo "- autoreconf."
autoreconf --force --install --symlink
echo "- aclocal."
aclocal -I m4
echo "- autoconf."
autoconf
echo "- automake."
automake -a -c
exit
