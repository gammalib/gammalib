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
import commands


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
    Commit changes into git repo

    Parameters
    ----------
    message : str
        Commit message
    """
    # Show changes
    print("The following files have been changed:")
    os.system('git status -s')

    # Ask again whether the changes should be committed
    if confirm("Commit all changes?"):

        # Optionally enter a commit message
        if len(message) == 0:
            message = str(raw_input("Please enter a commit message: "))

        # Stage all files
        os.system('git add -u')

        # Commit all changes
        os.system('git commit -m "%s"' % message)
    
    # Return
    return


# ========================== #
# Check if git branch exists #
# ========================== #
def git_branch_exists(branch):
    """
    Check if git branch exists

    Parameters
    ----------
    branch : str
        Branch name

    Returns
    -------
    exists : bool
        True if git branch exists
    """
    # Check if git branch exists
    rc = os.system('git rev-parse --verify %s >/dev/null 2>&1' % branch)

    # Set exists flag
    if rc == 0:
        exists = True
    else:
        exists = False

    # Return flag
    return exists


# ================= #
# Create git branch #
# ================= #
def create_git_branch(branch):
    """
    Create a git branch

    Parameters
    ----------
    branch : str
        Branch name
    """
    # If branch exists already then ask whether it should be deleted
    if git_branch_exists(branch):
        if confirm("Git branch '%s' exists already, do you want to delete "
                   "it?" % branch):
            os.system('git branch -D %s' % branch)

    # If branch exists already then switch to the new branch, otherwise create
    # new branch and switch to it
    if git_branch_exists(branch):
        os.system('git checkout %s' % branch)
    else:
        print("Created new branch '%s'" % branch)
        os.system('git checkout -b %s' % branch)

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
    set_version_in_file('configure.ac', version, old)

    # Return
    return


# ============================== #
# Set version number in one file #
# ============================== #
def set_version_in_file(filename, version, current):
    """
    Set version number in one file

    The function replaces all current version numbers in a file by a new
    version number.

    Parameters
    ----------
    filename : str
        File name
    version : str
        New version number
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


# ======================= #
# Set the package version #
# ======================= #
def set_package_version(version):
    """
    Set the package version

    Parameters
    ----------
    version : str
        Version number
    """
    # List of all files in which the version numbers should be changed. The
    # file name should include their access paths with respect to the source
    # directory
    filenames = ['configure.ac',
                 'gammalib.pc.in',
                 'README.md',
                 'sonar-project.properties',
                 'doc/Doxyfile',
                 'doc/source/conf.py']


    # Get the current package version
    current = get_current_version()

    # Set versions by loop over all files
    for filename in filenames:
        set_version_in_file(filename, version, current)

    # Generate current version number
    split = current.split('.')
    cur   = '%s.%s' % (split[0], split[1])

    # Generate new version number
    split = version.split('.')
    ver   = '%s.%s' % (split[0], split[1])

    # Set version number in doc/source/conf.py file
    set_version_in_file('doc/source/conf.py', ver, cur)

    # Return
    return


# ============= #
# Build tarball #
# ============= #
def build_tarball(branch, folder):
    """
    Build tarball using 'make dist'

    Parameters
    ----------
    branch : str
        Git branch
    folder : str
        Folder name
    """
    # Set base path
    base_path = os.getcwd()

    # Enable/disable logging
    if confirm("Log actions in logfile?"):
        verbosity = ' >> %s/%s.log 2>&1' % (base_path, folder)
        os.system('rm -rf %s.log >/dev/null 2>&1' % folder)
    else:
        verbosity = ' >/dev/null 2>&1'

    # Check package
    print("Create tarball")
    rc = os.system('make dist%s' % verbosity)
    if rc == 0:
        print("Tarball creation successful")
    else:
        print("*** Failure in creation of tarball")

    # Return
    return


# ============= #
# Check tarball #
# ============= #
def check_tarball(filename):
    """
    Check tarball using 'make distcheck'

    Parameters
    ----------
    filename : str
        Tarball file name
    """
    # Set base path
    base_path = os.getcwd()

    # Set base folder
    folder = filename.rstrip('.tar.gz')

    # Enable/disable logging
    if confirm("Log check in logfile?"):
        name      = filename.rstrip('.tar.gz')
        verbosity = ' >> %s/%s.check 2>&1' % (base_path, folder)
        os.system('rm -rf %s.check >/dev/null 2>&1' % folder)
    else:
        verbosity = ' >/dev/null 2>&1'

    # Check tarball
    print("Check tarball")
    rc = os.system('make distcheck%s' % verbosity)
    if rc == 0:
        print("Tarball checking successful")
    else:
        print("*** Failure in checking of tarball")

    # Return
    return


# ================ #
# Manage main menu #
# ================ #
def main_menu():
    """
    Manage main menu
    """
    # Print main menu
    print('[1] Make a new release')
    print('[2] Set the package version in current branch')
    print('[3] Set the libtool version in current branch')
    print('[4] Commit changes in current branch')
    print('[5] Create tarball from current branch')
    print('[6] Check tarball')
    print('[q] Quit')

    # Wait for the input
    waiting = True
    while waiting:
        choice = str(raw_input('Enter your choice: '))
        if choice == '1' or choice == '2' or choice == '3' or choice == '4' or \
           choice == '5' or choice == '6' or choice == 'q':
            waiting = False

    # Return choice
    return choice


# ====================== #
# Manage package release #
# ====================== #
def release_menu(branch='release'):
    """
    Manage package release

    Parameters
    ----------
    branch : str, optional
        Release branch name
    """
    # Annonce actions
    print("")
    print("Make a new release")
    print("------------------")

    # Step 1: Create a release branch
    if confirm("Step 1: Do you want to create a '%s' branch?" % branch):
        create_git_branch(branch)
    else:
        if git_branch_exists(branch):
            os.system('git checkout %s' % branch)
        else:
            requested = branch
            branch    = commands.getoutput('git rev-parse --abbrev-ref HEAD')
            print("WARNING: No '%s' branch found, continue working on "
                  "current branch '%s'" % (requested, branch))

    # Get current package version
    current = get_current_version()

    # Step 2: Set package version
    print("")
    if confirm("Step 2: Current GammaLib version is '%s'. Do you want to "
               "change the package version?" % current):
        package_version_menu()

    # Get current libtool version
    current = get_libtool_version()

    # Step 3: Set libtool version
    print("")
    if confirm("Step 3: Current Libtool version is '%s'. Do you want to "
               "change the libtool version?" % current):
        libtool_version_menu()

    # Step 4: Commit changes
    print("")
    if confirm("Step 4: Commit changes?"):
        package = get_current_version()
        libtool = get_libtool_version()
        message = "GammaLib package version set to '%s' and libtool version "\
                  "set to '%s'" % (package, libtool)
        commit(message)

    # Step 5: Push branch into git
    print("")
    if confirm("Step 5: Push changes?"):
        os.system('git push origin %s' % branch)

    # Step 6: Build tarball
    print("")
    if confirm("Step 6: Build tarball?"):
        folder = "gammalib-%s" % get_current_version()
        build_tarball(branch, folder)

    # Step 7: Check tarball
    print("")
    if confirm("Step 7: Check tarball?"):
        filename = "gammalib-%s.tar.gz" % get_current_version()
        check_tarball(filename)

    # Print separator
    print("")

    # Return
    return


# ============================== #
# Manage package version setting #
# ============================== #
def package_version_menu():
    """
    Manage package version setting
    """
    # Get the current package version
    current = get_current_version()

    # Looping until the version change is accepted
    looping = True
    while looping:

        # Get package version from input
        waiting = True
        while waiting:
            version = str(raw_input("Current GammaLib version is '%s'. Please "
                                    "enter new GammaLib version: " % current))
            if version == "q":
                sys.exit()
            if check_version(version) == False:
                print("*** Invalid GammaLib version. "
                      "Please enter the version in the format x.y.z[.devn]")
            else:
                waiting = False

        # If the version number change is confirmed then do the change
        if confirm("Change version to '%s'?" % version):

            # Set package version
            set_package_version(version)

            # Get new package and libtool versions
            version = get_current_version()

            # Inform about version change
            print("GammaLib version changed to '%s'" % (version))

            # Signal that we can exit looping
            looping = False

    # Return
    return


# ============================== #
# Manage libtool version setting #
# ============================== #
def libtool_version_menu():
    """
    Manage libtool version setting
    """
    # Get current libtool version
    current = get_libtool_version()

    # Continue only if libtool version should be changed
    if confirm("Current libtool version is '%s'. Do you want to change the "
               "version?" % current):

        # Determine package changes
        library = confirm(" Has the source code changed since last release?")
        added   = confirm(" Have interfaces been added since last release?")
        removed = confirm(" Have interfaces been removed since last release?")
        changed = confirm(" Have interfaces changed since last release?")

        # Set libtool version
        set_libtool_version(library=library, added=added, removed=removed,
                            changed=changed)

        # Get new libtool version
        version = get_libtool_version()

        # Print new libtool version
        print("Libtool version changed to '%s'" % version)

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
            print('')
        elif choice == '2':
            package_version_menu()
            print('')
        elif choice == '3':
            libtool_version_menu()
            print('')
        elif choice == '4':
            commit('')
            print('')
        elif choice == '5':
            branch = commands.getoutput('git rev-parse --abbrev-ref HEAD')
            folder = "gammalib-%s" % get_current_version()
            build_tarball(branch, folder)
            print('')
        elif choice == '6':
            filename = "gammalib-%s.tar.gz" % get_current_version()
            check_tarball(filename)
            print('')
        elif choice == 'q':
            break

