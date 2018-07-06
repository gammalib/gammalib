#! /usr/bin/env python
# ==========================================================================
# GammaLib code generator
#
# Copyright (C) 2017-2018 Juergen Knoedlseder
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
            text = line.strip(' *\n')
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


# ===================== #
# Set base class tokens #
# ===================== #
def set_base_tokens(classname, author, what, object):
    """
    Set replacement tokens for a generic base class
    """
    # Get current year
    year = str(date.today().year)

    # Set tokens
    tokens = [{'pattern': 'GTPLBase', 'string': classname},
              {'pattern': 'GTPLBASE', 'string': classname.upper()},
              {'pattern': '[WHAT]', 'string': what},
              {'pattern': '[AUTHOR]', 'string': author},
              {'pattern': 'TPL_OBJECT', 'string': object.lower()},
              {'pattern': '[YEAR]', 'string': year}]

    # Return tokens
    return tokens


# ========================== #
# Set container class tokens #
# ========================== #
def set_container_tokens(classname, author, what, object):
    """
    Set replacement tokens for a generic container class
    """
    # Get current year
    year = str(date.today().year)

    # Set tokens
    tokens = [{'pattern': 'GTPLContainer', 'string': classname},
              {'pattern': 'GTPLContainer', 'string': classname.upper()},
              {'pattern': '[WHAT]', 'string': what},
              {'pattern': '[AUTHOR]', 'string': author},
              {'pattern': 'TPL_CONTAINER', 'string': object.lower()},
              {'pattern': '[YEAR]', 'string': year}]

    # Return tokens
    return tokens


# ===================== #
# Set instrument tokens #
# ===================== #
def set_inst_tokens(name, instrument, author):
    """
    Set replacement tokens for an instrument
    """
    # Get current year
    year = str(date.today().year)

    # Set tokens
    tokens = [{'pattern': 'xxx', 'string': name.lower()},
              {'pattern': 'XXX', 'string': name.upper()},
              {'pattern': '[INSTRUMENT]', 'string': instrument},
              {'pattern': '[AUTHOR]', 'string': author},
              {'pattern': '[YEAR]', 'string': year}]

    # Return tokens
    return tokens


# ======================= #
# Add an instrument class #
# ======================= #
def add_instrument_class(name, tokens, classname, tmpname=None):
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
    tmpname : str, optional
        Template name
    """
    # Set template classname
    if tmpname == None:
        tmpname = classname.replace(name.upper(), 'XXX')

    # Set destination file names
    incfile = 'inst/%s/include/%s.hpp' % (name, classname)
    srcfile = 'inst/%s/src/%s.cpp'     % (name, classname)
    pyfile  = 'inst/%s/pyext/%s.i'     % (name, classname)
    
    # Add files
    add_file('src/template/%s.hpp' % tmpname, incfile, tokens)
    add_file('src/template/%s.cpp' % tmpname, srcfile, tokens)
    add_file('src/template/%s.i'   % tmpname, pyfile, tokens)

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
        f.write('%%include "%s.i"\n' % classname)

    # Return
    return


# ============================================== #
# Adjust configuration for new instrument module #
# ============================================== #
def adjust_config(name):
    """
    Parameters
    ----------
    name : str
        Module name
    """
    # Adjust inst/Makefile.am
    config = 'if WITH_INST_%s\nINST_%s = %s\nendif\n' % (name.upper(), name.upper(), name)
    for line in fileinput.FileInput('inst/Makefile.am',inplace=1):
        if 'Detect configuration' in line:
            line = line.replace(line, line+config)
        elif 'SUBDIRS=' in line.replace(' ', ''):
            line = line.replace(line, line[:-1]+' $(INST_%s)\n' % name.upper())
        print line,

    # Adjust configure.ac
    checks = '# Checks for %s interface\n' \
             'AC_ARG_WITH([%s],\n' \
             '            [AS_HELP_STRING([--with-%s],\n' \
             '                            [compile in %s specific interface [default=yes]])],\n' \
             '            [],\n' \
             '            [with_%s=yes])\n' \
             'if test "x$with_%s" = "xyes"; then\n' \
             '  AC_CHECK_FILE([$srcdir/inst/%s/README.md],\n' \
             '                [with_%s=yes],\n' \
             '                [with_%s=no])\n' \
             'fi\n' \
             'AM_CONDITIONAL(WITH_INST_%s, test "x$with_%s" = "xyes")\n' % \
             (name.upper(), name, name, name.upper(), name, name, name,
              name, name, name.upper(), name)
    wrap   = '  if test "x$with_%s" = "xyes"; then\n' \
             '    AC_CHECK_FILES([$srcdir/pyext/gammalib/%s_wrap.cpp\n' \
             '                    $srcdir/pyext/gammalib/%s.py],,\n' \
             '                   [has_wrappers="no"])\n' \
             '  fi' % (name, name, name)
    files  = 'if test "x$with_%s" = "xyes"; then\n' \
             '  AC_CONFIG_FILES([inst/%s/Makefile])\n' \
             '  AC_CONFIG_FILES([inst/%s/test/Makefile])\n' \
             'fi' % (name, name, name)
    info   = 'if test "x$with_%s" = "xyes"; then\n' \
             '  echo "  * %s interface                (yes)"\n' \
             'else\n' \
             '  echo "  - %s interface                (no)"\n' \
             'fi' % (name, name.upper(), name.upper())
    for line in fileinput.FileInput('configure.ac',inplace=1):
        if '# Insert new checks here' in line:
            print(checks)
        elif '# Insert new wrappers here' in line:
            print(wrap)
        elif '# Inset new Makefiles here' in line:
            print(files)
        elif '# Insert new interface informations here' in line:
            print(info)
        print line,

    # Return
    return


# ===================================================== #
# Adjust Python configuration for new instrument module #
# ===================================================== #
def adjust_python_config(name, instrument):
    """
    Parameters
    ----------
    name : str
        Module name
    instrument : str
        Instrument
    """
    # Adjust pyext/setup.py.in
    module = 'if \'@WITH_INST_%s_TRUE@\' != \'#\':\n' \
             '    print(\'Include %s instrument interface in gammalib Python module.\')\n' \
             '    inst_modules.append(\'%s\')' % \
             (name.upper(), instrument, name)
    for line in fileinput.FileInput('pyext/setup.py.in',inplace=1):
        if '# Insert new module here' in line:
            print(module)
        print line,

    # Adjust pyext/Makefile.am
    target  = 'if WITH_INST_%s\n' \
              '  %s_SWIG_TARGETS = $(wrapperdir)/gammalib/%s_wrap.cpp \\\n' \
              '                     $(wrapperdir)/gammalib/%s.py\n' \
              '  %s_TESTS        = $(top_srcdir)/inst/%s/test/test_%s.py\n' \
              'endif' % (name.upper(), name.upper(), name, name,
                         name.upper(), name, name.upper())
    targets = '                    $(%s_SWIG_TARGETS) ' % name.upper()
    rule    = '# Rule for %s module\n' \
              'if WITH_INST_%s\n' \
              '$(wrapperdir)/gammalib/%s.py: $(wrapperdir)/gammalib/%s_wrap.cpp\n' \
              '$(wrapperdir)/gammalib/%s_wrap.cpp: $(top_srcdir)/inst/%s/pyext/%s.i\n' \
              '	if $(SWIGCOMPILE) -MMD -MF "gammalib/%s.Tpi" -I$(top_srcdir)/pyext -o gammalib/%s_wrap.cpp -outdir gammalib $<; \\\n' \
              '	then mv -f "gammalib/%s.Tpi" "gammalib/%s.Pi"; else rm -f "gammalib/%s.Tpi"; exit 1; fi\n' \
              'endif\n' % (name.upper(), name.upper(), name, name, name,
                           name, name, name, name, name, name, name)
    for line in fileinput.FileInput('pyext/Makefile.am',inplace=1):
        if '# Insert new target here' in line:
            print(target)
        elif 'INST_SWIG_TARGETS=' in line.replace(' ', ''):
            last = line[-2]
            line = line.replace(line, line+targets+last+'\n')
        elif '$(MWL_TESTS)' in line:
            line = line.replace(line, line[:-1]+' $(%s_TESTS)\n' % name.upper())
        if '# Insert new rule here' in line:
            print(rule)
        print line,

    # Adjust pyext/gammalib/__init__.py.in
    module = ('from gammalib.%s import *' % name)
    for line in fileinput.FileInput('pyext/gammalib/__init__.py.in',inplace=1):
        if '# Insert new module here' in line:
            print(module)
        print line,

    # Return
    return


# =================================================== #
# Adjust test configuration for new instrument module #
# =================================================== #
def adjust_test_config(name, instrument):
    """
    Parameters
    ----------
    name : str
        Module name
    instrument : str
        Instrument
    """
    # Adjust test/Makefile.am
    test    = 'if WITH_INST_%s\n' \
              '  INST_%s     = $(top_builddir)/inst/%s/test/test_%s\n' \
              '  TEST_%s     = :$(top_srcdir)/inst/%s/test\n' \
              '  TEST_%s_ENV = TEST_%s_DATA=$(top_srcdir)/inst/%s/test/data\n' \
              'endif' % (name.upper(), name.upper(), name, name.upper(),
                         name.upper(), name, name.upper(), name.upper(),
                         name)
    testenv = '                    $(TEST_%s_ENV) ' % name.upper()
    tests   = '        $(INST_%s) \\\n' % name.upper()
    for line in fileinput.FileInput('test/Makefile.am',inplace=1):
        if '# Insert new test here' in line:
            print(test)
        elif 'TEST_PYTHON_ENV=' in line.replace(' ', ''):
            line = line.replace(line, line[:-1]+'$(TEST_%s)\n' % name.upper())
        elif 'TEST_DATA=$(top_srcdir)/test/data' in line:
            last = line[-2]
            line = line.replace(line, line+testenv+last+'\n')
        elif '$(TEST_PYTHON_SCRIPT)' in line:
            line = line.replace(line, tests+line)
        print line,

    # Adjust test/test_python.py
    test  = '# Try importing %s tests\n' \
            'try:\n' \
            '    import test_%s\n' \
            '    has_%s = True\n' \
            'except:\n' \
            '    has_%s = False\n' % \
            (instrument, name.upper(), name, name)
    data  = '        os.environ[\'TEST_%s_DATA\'] = \'%s/data\'' % \
            (name.upper(), name)
    suite = '    # Optionally handle %s suite\n' \
            '    if has_%s:\n' \
            '        suite_%s = test_%s.Test()\n' \
            '        suite_%s.set()\n' \
            '        suites.append(suite_%s)\n' % \
            (instrument, name, name, name.upper(), name, name)
    cpy   = '        os.system(\'cp -r %%s %%s\' %% (head+\'/%s\',  \'%s\'))' % \
            (name, name)
    for line in fileinput.FileInput('test/test_python.py',inplace=1):
        if '# Insert new test here' in line:
            print(test)
        if '# Insert new environment here' in line:
            print(data)
        if '# Insert new suite here' in line:
            print(suite)
        if '# Copy over test data' in line:
            print(cpy)
        print line,

    # Return
    return


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
        print('\nAll right. Have now:')
        print('Module name .....: "%s"' % name)
        print('Instrument name .: "%s"' % instrument)
        print('Your name .......: "%s"' % author)
        if confirm('Is this correct?'):
            break

    # Set tokens
    tokens = set_inst_tokens(name, instrument, author)

    # If the module exists already then determine all base classes that are
    # not yet implemented and propose to add them one by one ...
    if os.path.isdir('inst/%s' % name):

        # Determine which of the baseclasses are not yet implemented
        baseclasses = ['GObservation', 'GResponse',
                       'GEventList', 'GEventAtom', 'GEventCube', 'GEventBin',
                       'GInstDir', 'GRoi']
        files = glob.glob('inst/%s/include/*.hpp' % name)
        for file in files:
            baseclass = get_base_class(file)
            if baseclass in baseclasses:
                baseclasses.remove(baseclass)

        # If all baseclasses are already there then exist now
        if len(baseclasses) == 0:
            print('\nAll right. The module exists already and all base classes '
                  'are already implemented. Return now to main menu.')
            return
        else:
            print('\nAll right. Now let\'s add some classes:')

        # Loop over all missing baseclasses
        for baseclass in baseclasses:
            classname = 'G'+name.upper()+baseclass[1:]
            if confirm('Add "%s" class?' % classname):
                add_instrument_class(name, tokens, classname)

    # ... otherwise add a new module
    else:

        # Setup base classes
        baseclasses = ['GObservation', 'GResponse']

        # Ask further questions
        print('\nAll right. You want a new "%s" instrument module.' % name)
        while True:
            add_list = confirm('Do you want event list support?')
            add_bin  = confirm('Do you want binned event data support?')
            if not add_list and not add_bin:
                print('*** Warning: You need no support for event lists and '
                      'binned event data. This means that you have no data :o')
                if confirm('Is this true? You want a module that cannot handle data?'):
                    break
            else:
                break

        # Add base classes according to choice
        if add_list or add_bin:
            baseclasses.append('GInstDir')
        if add_list:
            baseclasses.extend(['GEventList', 'GEventAtom', 'GRoi'])
        if add_bin:
            baseclasses.extend(['GEventCube', 'GEventBin'])

        # Create directory structure
        os.mkdir('inst/%s' % name)
        os.mkdir('inst/%s/caldb' % name)
        os.mkdir('inst/%s/include' % name)
        os.mkdir('inst/%s/pyext' % name)
        os.mkdir('inst/%s/src' % name)
        os.mkdir('inst/%s/test' % name)
        os.mkdir('inst/%s/test/data' % name)

        # Add some files
        add_file('src/template/inst_Makefile.am',
                 'inst/%s/Makefile.am' % name, tokens)
        add_file('src/template/inst_README.md',
                 'inst/%s/README.md' % name, tokens)
        add_file('src/template/GXXXLib.hpp',
                 'inst/%s/include/G%sLib.hpp' % (name, name.upper()), tokens)
        add_file('src/template/xxx.i',
                 'inst/%s/pyext/%s.i' % (name, name), tokens)
        add_file('src/template/inst_test_Makefile.am',
                 'inst/%s/test/Makefile.am' % name, tokens)
        add_file('src/template/inst_test_XXX.hpp',
                 'inst/%s/test/test_%s.hpp' % (name, name.upper()), tokens)
        add_file('src/template/inst_test_XXX.cpp',
                 'inst/%s/test/test_%s.cpp' % (name, name.upper()), tokens)
        add_file('src/template/inst_test_XXX.py',
                 'inst/%s/test/test_%s.py' % (name, name.upper()), tokens)

        # Adjust configuration
        adjust_config(name)

        # Adjust Python configuration
        adjust_python_config(name, instrument)

        # Adjust test configuration
        adjust_test_config(name, instrument)

        # Loop over all baseclasses
        for baseclass in baseclasses:
            classname = 'G'+name.upper()+baseclass[1:]
            add_instrument_class(name, tokens, classname)

    # Return
    return


# ==================== #
# Manage generic class #
# ==================== #
def generic_class_menu():
    """
    Manage generic class
    """
    # Annonce actions
    print("")
    print("Add generic class")
    print("-----------------")

    # Stay in loop until there is a final confirmation
    while True:

        # Get class directory
        while True:
            dir = response('Please enter directory where class should reside '
                           '(e.g. "inst/cta", "src/obs")').lower()
            if os.path.isdir(dir):
                break
            else:
                print('*** Error: Directory not found.')

        # Enter other information
        classname = response('Please enter a class name (e.g. "GEnergy")')
        what      = response('Please say what the class is for (e.g. "Energy")')
        object    = response('Please say how an instance of the class should be named (e.g. "energy")')
        author    = response('Please enter your name (e.g. "Joe Public")')

        # Ask to confirm module summary
        print('\nAll right. Have now:')
        print('Class name ......: "%s"' % classname)
        print('Class descriptor : "%s"' % what)
        print('Class instance ..: "%s"' % object)
        print('Your name .......: "%s"' % author)
        if confirm('Is this correct?'):
            break

    # Set tokens
    tokens = set_base_tokens(classname, author, what, object)

    # If we have an instrument model then add a class to the instrument
    if 'inst/' in dir:
        name = dir[5:8]
        add_instrument_class(name, tokens, classname, 'GTPLBase')

    # Otherwise we have a core module
    else:
        print('Core modules not supported yet, sorry ;)')

    # Return
    return


# ====================== #
# Manage container class #
# ====================== #
def container_class_menu():
    """
    Manage container class
    """
    # Annonce actions
    print("")
    print("Add container class")
    print("-------------------")

    # Stay in loop until there is a final confirmation
    while True:

        # Get class directory
        while True:
            dir = response('Please enter directory where class should reside '
                           '(e.g. "inst/cta", "src/obs")').lower()
            if os.path.isdir(dir):
                break
            else:
                print('*** Error: Directory not found.')

        # Enter other information
        basename  = response('Please enter class name for the objects in the '
                             'container (e.g. "GEnergy")')
        baseinst  = response('Please specify how an instance of the class '
                             'objects should be named (e.g. "energy")')
        classname = response('Please enter class name for the container '
                             '(e.g. "GEnergies")')
        classinst = response('Please specify how an instance of the container '
                             'class should be named (e.g. "energies")')
        what      = response('Please say what the class contains (e.g. "Energy")')
        author    = response('Please enter your name (e.g. "Joe Public")')

        # Ask to confirm module summary
        print('\nAll right. Have now:')
        print('Container class name .........: "%s"' % classname)
        print('Container class instance .....: "%s"' % classinst)
        print('Object class in container ....: "%s"' % basename)
        print('Object instance in container .: "%s"' % baseinst)
        print('Object descriptor ............: "%s"' % what)
        print('Your name ....................: "%s"' % author)
        if confirm('Is this correct?'):
            break

    # Set tokens
    base_tokens      = set_base_tokens(basename, author, what, baseinst)
    container_tokens = set_container_tokens(classname, author, what, classinst)
    container_tokens.extend(base_tokens)

    # If we have an instrument model then add a class to the instrument
    if 'inst/' in dir:
        name = dir[5:8]
        if confirm('Do you want to add the base class "%s"?' % basename):
            add_instrument_class(name, base_tokens, basename, 'GTPLBase')
        add_instrument_class(name, container_tokens, classname, 'GTPLContainer')

    # Otherwise we have a core module
    else:
        print('Core modules not supported yet, sorry ;)')

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
    print('[1] Add generic class')
    print('[2] Add container class')
    print('[3] Add instrument module')
    print('[q] Quit')

    # Wait for the input
    waiting = True
    while waiting:
        choice = str(raw_input('Enter your choice: '))
        if choice == '1' or choice == '2' or choice == '3' or choice == 'q':
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
            generic_class_menu()
            print('')
        elif choice == '2':
            container_class_menu()
            print('')
        elif choice == '3':
            module_menu()
            print('')
        elif choice == 'q':
            break
