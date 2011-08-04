#! /usr/bin/env python

from gammalib import *

def test_one(a):
	a[0]=1
	a[1]=2
	a[2]=3
	print a
	return a

def test_two(a):
	a.remove(0,0)
	a.remove(0,1)
	a.insert(0,2)
	a.remove(1,2)
	a.remove(1,1)
	a.remove(0,1)
	print a
	return a

def test_three(a):
	a.insert(0,10)
	print a
	a[7]=99
	a.remove(8,2)
	a[0]=99
	a.remove(1,5)
	print a
	return a


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Perform testing.
	"""
	# Perform tests
	#a = GFitsTableBitCol("Test",3)
	#a = GFitsTableBoolCol("Test",3)
	#a = GFitsTableByteCol("Test",3)
	#a = GFitsTableDoubleCol("Test",3)
	#a = GFitsTableFloatCol("Test",3)
	#a = GFitsTableLongCol("Test",3)
	#a = GFitsTableLongLongCol("Test",3)
	#a = GFitsTableShortCol("Test",3)
	#a = GFitsTableStringCol("Test",3,10)
	#a = GFitsTableULongCol("Test",3)
	a = GFitsTableUShortCol("Test",3)
	a = test_one(a)
	a = test_two(a)
	a = test_three(a)
