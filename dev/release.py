#! /usr/bin/env python
# ==========================================================================
# GammaLib release manager
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
import os
import sys


# ========================== #
# Get current version number #
# ========================== #
def get_current_version():
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
def get_libtool_version():
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
                 'gammalib.pc.in',
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
    old   = get_libtool_version()
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
def set_special_doc_source_conf(version, current):
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


# ======================= #
# Set the package version #
# ======================= #
def set_package_version(version, library=False, added=False, removed=False,
                        changed=False):
    """
    Set the package version

    Parameters
    ----------
    version : str
        Version number
    """
    # Get the current package version
    current = get_current_version()

    # Set versions
    set_versions(version, current)

    # Set libtool version
    set_libtool_version(library=library, added=added, removed=removed,
                        changed=changed)

    # Set special version number in doc/source/conf.py file
    set_special_doc_source_conf(version, current)

    # Return
    return


# ======================================= #
# Check whether version string is correct #
# ======================================= #
def check_version(version):
    """
    Check whether version string is correct

    Parameters
    ----------
    version : str
        Version number

    Returns
    -------
    check : bool
        True if version string is correct
    """
    # Initialise check with False
    check = False

    # Split version string
    split = version.split('.')

    # If there are three elements or four elements we may have a valid version
    # string
    if len(split) == 3 or len(split) == 4:

        # If the three first elements are integers we may have a valid version
        # string
        if split[0].isdigit() and split[1].isdigit() and split[2].isdigit():

            # If we have three elements the version is valid and we are done
            if len(split) == 3:
                check = True

            # Otherwise check if the 4th element is of format devn, where n
            # is an integer number
            elif split[3][0:3] == 'dev':
                if split[3][3:].isdigit():
                    check = True

    # Return check
    return check


# ================ #
# Manage main menu #
# ================ #
def main_menu():
    """
    Manage main menu
    """
    # Print main menu
    print('[1] Make a new package release')
    print('[2] Set the package version')
    print('[q] Quit')

    # Wait for the input
    waiting = True
    while waiting:
        choice = str(raw_input('Enter your choice: '))
        if choice == '1' or choice == '2' or choice == 'q':
            waiting = False

    # Return choice
    return choice


# ====================== #
# Manage package release #
# ====================== #
def release_menu():
    """
    Manage package release
    """
    # Get current package version
    current = get_current_version()

    # Create a new release branch
    
    # Set the package version
    print(' Current GammaLib version is "%s"' % current)
    if confirm(' Do you want to change the package version?'):
        version_menu()

    # Commit changes

    # Create tarball
    
    # Return
    return


# ============================== #
# Manage package version setting #
# ============================== #
def version_menu():
    """
    Manage package version setting
    """
    # Get the current package version
    current = get_current_version()

    # Print current package version
    print(' Current GammaLib version is "%s"' % current)

    # Looping until the version change is accepted
    looping = True
    while looping:

        # Get package version from input
        waiting = True
        while waiting:
            version = str(raw_input('  Please enter new GammaLib version: '))
            if version == 'q':
                sys.exit()
            if check_version(version) == False:
                print('  *** Invalid GammaLib version. '
                      'Please enter the version in the format x.y.z[.devn]')
            else:
                waiting = False

        # If the version number change is confirmed then do the change
        if confirm('  Change version to "%s"?' % version):

            # Determine package changes
            library = confirm('  Has the source code changed since last release?')
            added   = confirm('  Have interfaces been added since last release?')
            removed = confirm('  Have interfaces been removed since last release?')
            changed = confirm('  Have interfaces changed since last release?')

            # Set package version
            set_package_version(version, library=library, added=added,
                                removed=removed, changed=changed)

            # Get new package and libtool versions
            version = get_current_version()
            libtool = get_libtool_version()

            # Inform about version change
            print(' GammaLib version changed to "%s", libtool version set to "%s"' %
                  (version, libtool))

            # Ask whether changes should be committed
            if confirm(' Should the changes be committed?'):
                commit('GammaLib package version changed to %s' % version)

            # Signal that we can exit looping
            looping = False

    # Return
    return


# ================= #
# Confirm something #
# ================= #
def confirm(text):
    """
    Confirm something

    Parameters
    ----------
    text : str
        Thing to confirm

    Returns
    -------
    flag : bool
        True if something is confirmed
    """
    # Confirmation loop
    waiting = True
    while waiting:
        confirmation = str(raw_input(text+' (y/n): '))
        if confirmation == 'q':
            sys.exit()
        elif confirmation == 'y':
            flag    = True
            waiting = False
        elif confirmation == 'n':
            flag    = False
            waiting = False

    # Return confirmation flag
    return flag


# ============== #
# Commit changes #
# ============== #
def commit(message):
    """
    Commit changes into Git repo

    Parameters
    ----------
    message : str
        Commit message
    """
    # Show changes
    print(' The following files have been changed:')
    os.system('git status -s')

    # Ask again whether the changes should be committed
    if confirm(' Commit all changes?'):
    
        # Stage all files
        os.system('git add -u')

        # Commit all changes
        os.system('git commit -m "%s"' % message)
    
    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Clear console
    os.system('clear')

    # Print header
    print('GammaLib release manager')
    print('========================')
    print('')
    
    # Enter endless loop
    while True:

        # Show main menu
        choice = main_menu()

        # Dispatch according to choice
        if choice == '1':
            release_menu()
        elif choice == '2':
            version_menu()
        elif choice == 'q':
            break

