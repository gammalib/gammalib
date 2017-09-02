#! /usr/bin/env python
# ==========================================================================
# GammaLib code generator
#
# Copyright (C) 2017 Juergen Knoedlseder
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
import glob
import commands
import fileinput
from datetime import date
import gammalib


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


# ============ #
# Get response #
# ============ #
def response(text, confirm_response=False):
    """
    Get response

    Parameters
    ----------
    text : str
        Thing we want to get a response to
    confirm_response : bool
        True if response should be confirmed

    Returns
    -------
    response : str
        Reponse
    """
    # Looping until response is accepted
    looping = True
    while looping:

        # Get response from input
        waiting = True
        while waiting:
            response = str(raw_input('%s: ' % text))
            if response == 'q':
                sys.exit()
            else:
                waiting = False

        # If the version number change is confirmed then do the change
        if confirm_response:
            if confirm('Is "%s" ok?' % response):
                looping = False
        else:
            looping = False

    # Return response
    return response


# ======== #
# Add file #
# ======== #
def add_file(infile, outfile, tokens):
    """
    Add a class file, replacing all tokens

    Parameters
    ----------
    infile : str
        Input filename
    outfile : str
        Output filename
    tokens : list of dict
        Token list
    """
    # Create output file
    file = open(outfile, "w")

    # Initialise number of header lines
    n_header = 0

    # Loop over input file
    for line in open(infile, 'r'):

        # Signal if this is a header line
        is_header = False
        if len(line) > 75:
            is_header = line[1] == '*' and line[75] == '*'
            n_header += 1

        # Replace tokens in line
        for token in tokens:
            line = line.replace(token['pattern'], token['string'])

        # If we have a header line than format the line
        if is_header:
            text = gammalib.strip_whitespace(line[2:69])
            if n_header == 2:
                line = ' *'+gammalib.centre(text,73)+'*\n'
            elif n_header == 4:
                line = ' *'+gammalib.left('  '+text,73)+'*\n'

        # Write out line
        file.write(line)

    # Close file
    file.close()

    # Return
    return


# ==================== #
# Determine base class #
# ==================== #
def get_base_class(filename):
    """
    Determine base class of a file

    Parameters
    ----------
    filename : str
        Filename

    Returns
    -------
    baseclass : str
        Base class name
    """
    # Initialise empty baseclass
    baseclass = ''

    # Loop over file
    for line in open(filename, 'r'):

        # Search for base class name
        if 'class' in line:
            i = line.find('public')
            if i != -1:
                baseclass = line[i+7:].rstrip('{ \n')
                break

    # Return
    return baseclass


# ========== #
# Set tokens #
# ========== #
def set_tokens(name, instrument, author):
    """
    Set replacement tokens
    """
    # Get current year
    year = str(date.today().year)

    # Set tokens
    tokens = [{'pattern': 'XXX', 'string': name.upper()},
              {'pattern': '[INSTRUMENT]', 'string': instrument},
              {'pattern': '[AUTHOR]', 'string': author},
              {'pattern': '[YEAR]', 'string': year}]

    # Return tokens
    return tokens


# ======================== #
# Manage instrument module #
# ======================== #
def module_menu():
    """
    Manage instrument module
    """
    # Annonce actions
    print("")
    print("Add instrument module")
    print("---------------------")

    # Stay in loop until there is a final confirmation
    while True:

        # Enter a module name
        while True:
            name = response('Please enter a 3 digit-name for the module').lower()
            if len(name) == 3:
                break
            else:
                print('*** Error: Module name "%s" has not 3 digits.' % name)

        # If module exists already then ask if classes should be added, and
        # quit if this is not the case
        if os.path.isdir('inst/%s' % name):
            if not confirm('Module "%s" exists already. Do you want to add '
                           'classes to this module?' % name):
                return

        # Enter other information
        instrument = response('Please enter the instrument name (e.g. "Fermi/LAT")')
        author     = response('Please enter your name (e.g. "Joe Public")')

        # Ask to confirm module summary
        print('All right. Have now:')
        print(' Module name ....: "%s"' % name)
        print(' Instrument name : "%s"' % instrument)
        print(' Your name ......: "%s"' % author)
        if confirm('Is this correct?'):
            break

    # Set tokens
    tokens = set_tokens(name, instrument, author)

    # If the module exists already then determine all base classes that are
    # not yet implemented
    if os.path.isdir('inst/%s' % name):

        # Determine which of the baseclasses are not yet implemented
        baseclasses = ['GObservation', 'GResponse',
                       'GEventAtom', 'GEventBin', 'GEventCube', 'GEventList',
                       'GInstDir', 'GRoi']
        files = glob.glob('inst/%s/include/*.hpp' % name)
        for file in files:
            baseclass = get_base_class(file)
            if baseclass in baseclasses:
                baseclasses.remove(baseclass)

        # If all baseclasses are already there then exist now
        if len(baseclasses) == 0:
            print('All right. The module exists already and all base classes '
                  'are already implemented. Return now to main menu.')
            return
        else:
            print("All right. Now let's add some classes:")

        # Loop over all missing baseclasses
        for baseclass in baseclasses:
            classname = 'G'+name.upper()+baseclass[1:]
            if confirm('Add "%s" class?' % classname):
                add_instrument_class(name, tokens, classname)

    # Return
    return


# ======================= #
# Add an instrument class #
# ======================= #
def add_instrument_class(name, tokens, classname):
    """
    Add an instrument class

    Parameters
    ----------
    name : str
        Module name
    tokens : str
        Tokens for replacement
    classname : str
        Class name
    """
    # Set destination file names
    incfile = 'inst/%s/include/%s.hpp' % (name, classname)
    srcfile = 'inst/%s/src/%s.cpp' % (name, classname)
    pyfile  = 'inst/%s/pyext/%s.i' % (name, classname)
    
    # Add files
    add_file('src/template/GXXXRoi.hpp', incfile, tokens)
    add_file('src/template/GXXXRoi.cpp', srcfile, tokens)
    add_file('src/template/GXXXRoi.i',   pyfile, tokens)

    # Update module header
    filename   = 'inst/%s/include/G%sLib.hpp' % (name, name.upper())
    afterline  = '#include "G%sObservation.hpp"' % (name.upper())
    insertline = '#include "%s.hpp"\n' % (classname)
    for line in fileinput.FileInput(filename,inplace=1):
        if afterline in line:
            line = line.replace(line, line+insertline)
        print line,

    # Update module Makefile.am
    filename    = 'inst/%s/Makefile.am' % (name)
    insertline1 = '          src/%s.cpp ' % (classname)
    insertline2 = '                     include/%s.hpp ' % (classname)
    for line in fileinput.FileInput(filename,inplace=1):
        if 'sources=' in line.replace(' ', ''):
            last = line[-2]
            line = line.replace(line, line+insertline1+last+'\n')
        if 'pkginclude_HEADERS=' in line.replace(' ', ''):
            last = line[-2]
            line = line.replace(line, line+insertline2+last+'\n')
        print line,

    # Update Python module
    filename = 'inst/%s/pyext/%s.i' % (name, name)
    with open(filename, 'a') as f:
        f.write('%%include "%s.i"' % classname)

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
    print('[1] Add instrument module')
    print('[2] Add spectral model')
    print('[3] Add spatial model')
    print('[4] Add temporal model')
    print('[4] Add generic class')
    print('[q] Quit')

    # Wait for the input
    waiting = True
    while waiting:
        choice = str(raw_input('Enter your choice: '))
        if choice == '1' or choice == '2' or choice == '3' or choice == '4' or \
           choice == 'q':
            waiting = False

    # Return choice
    return choice


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Check if script has been started from GammaLib root directory
    if not os.path.isfile('gammalib.pc.in'):
        print('"codgen.py" script needs to be started from GammaLib source code '
              'root directory, you are somewhere else. Quit now.')
        sys.exit()

    # Clear console
    os.system('clear')

    # Print header
    print('GammaLib code generator')
    print('=======================')
    print('')

    # Enter endless loop
    while True:

        # Show main menu
        choice = main_menu()

        # Dispatch according to choice
        if choice == '1':
            module_menu()
            print('')
        elif choice == '2':
            #package_version_menu()
            print('')
        elif choice == '3':
            #libtool_version_menu()
            print('')
        elif choice == '4':
            #commit('')
            print('')
        elif choice == 'q':
            break
