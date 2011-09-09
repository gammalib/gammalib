#! /usr/bin/env python

#from distutils.core import setup, Extension
from distutils import sysconfig
#import glob
#import os
#import sys
import commands


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
	"""
	"""
	# Get configuration variables
	cvars = sysconfig.get_config_vars()
	
	# Dump configuration variables
	#print cvars
	for key in cvars:
		print key, "=", cvars[key]

	# Get library information
	print commands.getoutput("file build/temp.*/gammalib_wrap.o")
	print commands.getoutput("file ../src/.libs/libgamma_python.a")
	print commands.getoutput("file ../src/.libs/libgamma.a")
	print commands.getoutput("file ../src/.libs/libgamma.dylib")
	print commands.getoutput("file build/_gammalib.so")
	print commands.getoutput("file build/lib.*/_gammalib.so")
	print commands.getoutput("file /usr/local/gamma/lib/libcfitsio.dylib")
	print commands.getoutput("file /opt/local/lib/libcfitsio.dylib")
	print commands.getoutput("file /usr/lib/libreadline.dylib")