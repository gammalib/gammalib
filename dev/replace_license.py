#! /usr/bin/env python

import os
import glob


# =========================== #
# Read license text from file #
# =========================== #
def read_license(filename):
	"""
	Read license text from file.
	"""
	# Initialise lines
	lines = []
	
	# Loop over file
	for line in open(filename, 'r'):
		if len(line.rstrip("\n").strip(" ")) > 0:
			lines.append(line)
	
	# Return lines
	return lines


# ==================== #
# Replace license text #
# ==================== #
def replace_license(filename, license_old, license_new):
	"""
	Replace license text in file.
	"""
	# Set temporary name
	tmpname = "test.txt"
	
	# Open temporary file
	file = open(tmpname, "w")
	
	# Initialise old text line finder
	ifind = 0
	
	# Loop over file
	for line in open(filename, 'r'):
		
		# Search line for license text
		if line == license_old[ifind]:
			ifind += 1
			if ifind >= len(license_old):
				for replace in license_new:
					file.write(replace)
				print "License text upated in file "+filename+"."
				ifind = 0
		
		# Handle partial license text
		elif ifind > 0:
			print "Partial old license text found in file "+filename+". Ignored."
			ifind = 0
		
		# ... otherwise write out line
		else:
			file.write(line)
	
	# Close temporary file
	file.close()
	
	# Replace input file
	os.rename(tmpname, filename)
	
	# Return
	return


# ======================================================== #
# Replace license text in files with name given by pattern #
# ======================================================== #
def update_files(pattern, license_old, license_new):
	"""
	Replace license text in files with name given by pattern.
	"""
	# Get sorted list of files
	files = glob.glob(pattern)
	files.sort()

	# Loop over all files
	for file in files:
		replace_license(file, license_old, license_new)
	
	# Return
	return
	

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
	"""
	Replace old by new license text.
	"""
	# Load license files
	license_old = read_license("license_old.txt")
	license_new = read_license("license_new.txt")
	
	# Replace license text
	update_files("include/*.hpp", license_old, license_new)
	update_files("pyext/*.i", license_old, license_new)
	update_files("src/*/*.*pp", license_old, license_new)
	update_files("inst/*/include/*.hpp", license_old, license_new)
	update_files("inst/*/pyext/*.i", license_old, license_new)
	update_files("inst/*/src/*.cpp", license_old, license_new)
	update_files("inst/*/test/*.*pp", license_old, license_new)
	update_files("test/*.*pp", license_old, license_new)

