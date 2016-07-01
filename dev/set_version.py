#! /usr/bin/env python
# ==========================================================================
# Set the version number of the software
#
# Copyright (C) 2016 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import sys


# ========================== #
# Get current version number #
# ========================== #
def get_current():
    """
    Get current version number

    Returns
    -------
    version : str
        Version number
    """
    # Initialse version number of empty string
    version = ''

    # Open configure.ac file
    f = open('configure.ac', 'r')

    # Read all lines
    for line in f:

        # Search for pattern
        if line.startswith('AC_INIT('):

            # Extract version number
            start   = line.find(',')+1
            stop    = line.find(',', start)
            version = line[start:stop].strip().lstrip('[').rstrip(']')

    # Close file
    f.close()

    # Return version
    return version


# ========================== #
# Get libtool version number #
# ========================== #
def get_libtool():
    """
    Get libtools version number

    Returns
    -------
    version : str
        Version number
    """
    # Initialse version number of empty string
    version = ''

    # Open configure.ac file
    f = open('configure.ac', 'r')

    # Read all lines
    for line in f:

        # Search for pattern
        if line.startswith('GAMMALIB_LT_VERSION='):

            # Extract version number
            start   = line.find('"')+1
            stop    = line.find('"', start)
            version = line[start:stop].strip()

    # Close file
    f.close()

    # Return version
    return version


# ============================== #
# Set version number in one file #
# ============================== #
def set_version(filename, version, current):
    """
    Set version number

    Parameters
    ----------
    filename : str
        File name
    version : str
        Version number
    current : str
        Current version number
    """
    # Open file in read mode
    f = open(filename, 'r')

    # Read file content
    content = f.read()

    # Close file
    f.close()

    # Replace version number
    new_content = content.replace(current, version)

    # Open file in write mode
    f = open(filename, 'w')

    # Write file content
    f.write(new_content)

    # Close file
    f.close()

    # Return
    return


# =============================== #
# Set version number in all files #
# =============================== #
def set_versions(version, current):
    """
    Set version numbers in all files

    Parameters
    ----------
    version : str
        Version number
    current : str
        Current version number
    """
    # List of all relevant files, including their access paths from the
    # source directory
    filenames = ['configure.ac',
                 'README.md',
                 'sonar-project.properties',
                 'doc/Doxyfile',
                 'doc/source/conf.py']

    # Set versions by loop over all files
    for filename in filenames:
        set_version(filename, version, current)

    # Return
    return


# =============================================== #
# Set libtool version number in configure.ac file #
# =============================================== #
def set_libtool_version(library=False, added=False, removed=False, changed=False):
    """
    Set special version numbers in configure.ac file

    Parameters
    ----------
    library : bool, optional
        Library source code changed
    added : bool, optional
        Interfaces added
    removed : bool, optional
        Interfaces removed
    changed : bool, optional
        Interfaces changed
    """
    # Get libtool version number
    old   = get_libtool()
    split = old.split(':')

    # Decompose libtool version number
    current  = int(split[0])
    revision = int(split[1])
    age      = int(split[2])

    # If the library source code has changed at all since the last update, then
    # increment revision
    if library:
        revision += 1

    # If any interfaces have been added, removed, or changed since the last
    # update, increment current, and set revision to 0
    if added or removed or changed:
        current  += 1
        revision  = 0

    # If any interfaces have been added since the last public release, then
    # increment age
    if added:
        age += 1

    # If any interfaces have been removed or changed since the last public
    # release, then set age to 0
    if removed or changed:
        age = 0

    # Build new libtool version number
    version = '%d:%d:%d' % (current, revision, age)

    # Set libtool version number
    set_version('configure.ac', version, old)

    # Return
    return


# ===================================================== #
# Set special version number in doc/source/conf.py file #
# ===================================================== #
def set_special_conf(version, current):
    """
    Set special version numbers in doc/source/conf.py file

    Parameters
    ----------
    version : str
        Version number
    current : str
        Current version number
    """
    # Generate current version number
    split = current.split('.')
    cur   = '%s.%s' % (split[0], split[1])

    # Generate new version number
    split = version.split('.')
    ver   = '%s.%s' % (split[0], split[1])

    # Set version number
    set_version('doc/source/conf.py', ver, cur)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Check for version number
    if len(sys.argv) < 3:
        print('Usage: set_version.py version [larc]')
        print('       l : soucre code changed')
        print('       a : interface(s) added')
        print('       r : interface(s) removed')
        print('       c : interface(s) changed')
        sys.exit()

    # Extract version number
    version = sys.argv[1]

    # Extract code change flags
    if len(sys.argv) == 3:
        flags = sys.argv[2]

    # Set change flags
    library = False
    added   = False
    removed = False
    changed = False
    if flags.find('l') != -1:
        library = True
    if flags.find('a') != -1:
        added = True
    if flags.find('r') != -1:
        removed = True
    if flags.find('c') != -1:
        changed = True

    # Get the current version
    current = get_current()

    # Set versions
    set_versions(version, current)

    # Set libtool version
    set_libtool_version(library=library, added=added, removed=removed,
                        changed=changed)

    # Set special version number in doc/source/conf.py file
    set_special_conf(version, current)

