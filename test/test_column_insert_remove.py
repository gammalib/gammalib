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

def set(a):
	#a[0] = 1.0
	a[1] = 65536
	a[2] = True
	#a[3] = "Hallo"
	print a[1]


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Perform testing.
	"""
	# Perform tests
	#a = GFitsTableBitCol("Test",4)
	#a = GFitsTableBoolCol("Test",4)
	#a = GFitsTableByteCol("Test",4)
	#a = GFitsTableDoubleCol("Test",4)
	#a = GFitsTableFloatCol("Test",4)
	#a = GFitsTableLongCol("Test",4)
	#a = GFitsTableLongLongCol("Test",4)
	#a = GFitsTableShortCol("Test",4)
	#a = GFitsTableStringCol("Test",4,10)
	#a = GFitsTableULongCol("Test",4)
	a = GFitsTableUShortCol("Test",4)
	a = test_one(a)
	a = test_two(a)
	a = test_three(a)
	#set(a)
	
