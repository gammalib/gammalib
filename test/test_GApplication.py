# ==========================================================================
# This module performs unit tests for the GammaLib application module
#
# Copyright (C) 2012-2017 Juergen Knoedlseder
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
import gammalib
import test_support


# ========================================== #
# Test class for GammaLib application module #
# ========================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib application module
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Set members
        self._logfile = 'test_logger.log'

        # Return
        return

    # Check log file versus a reference
    def _check_log_file(self, ref, calls):
        """
        Check log file versus a reference

        Parameters
        ----------
        ref : list of str
            Reference lines
        calls : list of str
            Methods that are tested
        """
        # Open log file, read all line and close file
        file  = open(self._logfile, 'r')
        lines = file.readlines()
        file.close()

        # Check all lines versus reference. There is a special kluge in case
        # that the reference is "true\n" as in some Python/SWIG version this
        # is displayed by "1\n". Both variants are tested to avoid version
        # dependent test results.
        for i, line in enumerate(lines):
            if ref[i] == 'true\n':
                self.test_assert(line == ref[i] or line == '1\n',
                                 'Check '+calls[i]+': '+line)
            else:
                self.test_value(line, ref[i], 'Check '+calls[i]+': '+line)

        # Return
        return

    # Setup GApplicationPars container
    def _setup_pars(self):
        """
        Setup GApplicationPars container

        Returns
        -------
        pars : `~gammalib.GApplicationPars`
            Parameter container
        """
        # Setup application parameter container
        pars = gammalib.GApplicationPars()
        for i in range(10):
            name = '%s' % i
            par  = gammalib.GApplicationPar(name,'r','a','1.0','0.0','2.0','Dummy')
            pars.append(par)

        # Return application parameter container
        return pars

    # Test GLog class
    def _test_log(self):
        """
        Test GLog class
        """
        # Allocate logger
        log = gammalib.GLog(self._logfile, True)

        # Test print methods
        log('Test __call__(std::string)\n')
        log(True)
        log('\n')
        log(41)
        log('\n')
        log(3.1415)
        log('\n')
        log.parformat('This is a parameter')
        log('\n')
        log.toupper('upper case')
        log('\n')
        log.tolower('LOWER CASE')
        log('\n')
        log.fill('Higgs', 3)
        log('\n')
        log.left('Left', 10)
        log('\n')
        log.right('Right', 10)
        log('\n')
        log.centre('Centre', 10)
        log('\n')
        log.header0('Header 0')
        log('\n')
        log.header1('Header 1')
        log.header2('Header 2')
        log.header3('Header 3')

        # Close logger
        log.close()

        # Set reference
        ref = []
        ref.append('Test __call__(std::string)\n')
        ref.append('true\n') # New swig version (3.x.y)
        ref.append('41\n')
        ref.append('3.1415\n')
        ref.append(' This is a parameter .......: \n')
        ref.append('UPPER CASE\n')
        ref.append('lower case\n')
        ref.append('HiggsHiggsHiggs\n')
        ref.append('Left      \n')
        ref.append('     Right\n')
        ref.append('  Centre  \n')
        ref.append('\n')
        ref.append('+==========+\n')
        ref.append('| Header 1 |\n')
        ref.append('+==========+\n')
        ref.append('+----------+\n')
        ref.append('| Header 2 |\n')
        ref.append('+----------+\n')
        ref.append('=== Header 3 ===\n')

        # Set calls for error messages
        calls = []
        calls.append('__call__(std::string)')
        calls.append('__call__(bool)')
        calls.append('__call__(int)')
        calls.append('__call__(double)')
        calls.append('parformat()')
        calls.append('toupper()')
        calls.append('tolower()')
        calls.append('fill()')
        calls.append('left()')
        calls.append('right()')
        calls.append('centre()')
        calls.append('header0()')
        calls.append('header1()')
        calls.append('header1()')
        calls.append('header1()')
        calls.append('header2()')
        calls.append('header2()')
        calls.append('header2()')
        calls.append('header3()')

        # Check log file
        self._check_log_file(ref, calls)

        # Create copy of logger
        cpy_log = log.copy()

        # Print one more line into copy
        cpy_log.header3('Header 3bis')

        # Append reference and call for error message
        ref.append('=== Header 3bis ===\n')
        calls.append('header3()')

        # Check log file again
        self._check_log_file(ref, calls)

        # Return
        return

    # Test GApplication class
    def _test_app(self):
        """
        Test GApplication class
        """
        # Set PFILES environment variable
        os.environ['PFILES'] = os.environ['TEST_DATA']

        # Allocate test application
        app = gammalib.GApplication('test_GApplication', '1.1.0')

        # Check proper allocation of test application
        self.test_value(app._name(), 'test_GApplication', 'Check application name')
        self.test_value(app._version(), '1.1.0', 'Check application version')
        self.test_value(app.pars().size(), 10, 'Check number of parameters')

        # Check that we properly can loop over all parameters. The npars
        # variable counts the number of iterations. If the parameter "chatter"
        # is encounter then set the chattiness to 4
        npars = 0
        for par in app:
            npars += 1
            if par.name() == 'chatter':
                par.integer(4)

        # Check outcome of loop
        self.test_value(npars, 10, 'Check number of parameters in loop')
        self.test_value(app['chatter'].integer(), 4, 'Check that "chatter" '
                        'parameter could be stored correctly')

        # Check access by index
        self.test_value(app[1].integer(), 1)
        self.test_value(app[-5].integer(), 4)

        # Check setting by index
        par = app[1]
        par.integer(5)
        app[1] = par
        par = app[-5]
        par.integer(4)
        app[-5] = par
        self.test_value(app[1].integer(), 5)
        self.test_value(app[-5].integer(), 4)

        # Check setting and getting of parameters
        app['real']     = 100.0
        app['integer']  = 10
        app['string']   = 'GAL'
        app['filename'] = 'newfile.fits'
        app['time']     = gammalib.GTime('2000-01-01T12:00:00')
        self.test_value(app['real'].real(), 100.0)
        self.test_value(app['integer'].integer(), 10)
        self.test_value(app['string'].string(), 'GAL')
        self.test_value(app['filename'].filename().url(), 'newfile.fits')
        self.test_value(app['time'].time().utc(), '2000-01-01T12:00:00')

        # Check boolean exception (does not work on older Linux systems)
        #self.test_try('Test GApplication boolean parameter exception')
        #try:
        #    app['real'] = False
        #    self.test_try_failure('Exception not thrown')
        #except:
        #    self.test_try_success()
        #else:
        #    self.test_try_failure('This should never happen')

        # Check integer exception
        self.test_try('Test GApplication integer parameter exception')
        try:
            app['string'] = 1
            self.test_try_failure('Exception not thrown')
        except:
            self.test_try_success()
        else:
            self.test_try_failure('This should never happen')

        # Check string exception
        self.test_try('Test GApplication string parameter exception')
        try:
            app['real'] = 'Unknown'
            self.test_try_failure('Exception not thrown')
        except:
            self.test_try_success()
        else:
            self.test_try_failure('This should never happen')
        self.test_try('Test GApplication string parameter exception')
        try:
            app['time'] = 'Unknown'
            self.test_try_failure('Exception not thrown')
        except:
            self.test_try_success()
        else:
            self.test_try_failure('This should never happen')
        self.test_try('Test GApplication string parameter exception')
        try:
            app['clobber'] = 'Unknown'
            self.test_try_failure('Exception not thrown')
        except:
            self.test_try_success()
        else:
            self.test_try_failure('This should never happen')

        # Check time exception
        self.test_try('Test GApplication time parameter exception')
        try:
            app['real'] = gammalib.GTime()
            self.test_try_failure('Exception not thrown')
        except:
            self.test_try_success()
        else:
            self.test_try_failure('This should never happen')

        # Check handling of undefined parameters
        app['real']    = 'INDEF'
        app['integer'] = 'INDEF'
        app['time']    = 'INDEF'
        self.test_assert(app['real'].is_undefined(), 'INDEF real parameter')
        self.test_assert(app['integer'].is_undefined(), 'INDEF integer parameter')
        self.test_assert(app['time'].is_undefined(), 'INDEF time parameter')
        app['real']    = 'NONE'
        app['integer'] = 'NONE'
        app['time']    = 'NONE'
        self.test_assert(app['real'].is_undefined(), 'NONE real parameter')
        self.test_assert(app['integer'].is_undefined(), 'NONE integer parameter')
        self.test_assert(app['time'].is_undefined(), 'NONE time parameter')
        app['real']    = 'UNDEF'
        app['integer'] = 'UNDEF'
        app['time']    = 'UNDEF'
        self.test_assert(app['real'].is_undefined(), 'UNDEF real parameter')
        self.test_assert(app['integer'].is_undefined(), 'UNDEF integer parameter')
        self.test_assert(app['time'].is_undefined(), 'UNDEF time parameter')
        app['real']    = 'UNDEFINED'
        app['integer'] = 'UNDEFINED'
        app['time']    = 'UNDEFINED'
        self.test_assert(app['real'].is_undefined(), 'UNDEFINED real parameter')
        self.test_assert(app['integer'].is_undefined(), 'UNDEFINED integer parameter')
        self.test_assert(app['time'].is_undefined(), 'UNDEFINED time parameter')

        # Check handling of NaN parameters
        app['real']    = 'INF'
        app['integer'] = 'INF'
        self.test_assert(app['real'].is_notanumber(), 'INF real parameter')
        self.test_assert(app['integer'].is_valid(), 'INF integer parameter')
        app['real']    = 'INFINITY'
        app['integer'] = 'INFINITY'
        self.test_assert(app['real'].is_notanumber(), 'INFINITY real parameter')
        self.test_assert(app['integer'].is_valid(), 'INFINITY integer parameter')
        app['real']    = 'NAN'
        app['integer'] = 'NAN'
        self.test_assert(app['real'].is_notanumber(), 'NAN real parameter')
        self.test_assert(app['integer'].is_valid(), 'NAN integer parameter')

        # Return
        return

    # Test GApplicationPars class
    def _test_pars(self):
        """
        Test GApplicationPars class
        """
        # Test GApplicationPars constructor with bad filename
        self.test_try('Test GApplicationPars constructor with bad filename')
        try:
            pars = gammalib.GApplicationPars('testme.par')
            self.test_try_failure('Exception not thrown')
        except:
            self.test_try_success()
        else:
            self.test_try_failure('This should never happen')

        # Return
        return

    # Test GApplicationPars class access operators
    def _test_pars_access(self):
        """
        Test GApplicationPars class model access
        """
        # Setup application parameter container
        pars = self._setup_pars()

        # Perform application parameter access tests
        test_support._container_access_index(self, pars)

        # Check application parameter setting by index from start
        par     = gammalib.GApplicationPar('98','r','a','1.0','0.0','2.0','Dummy')
        pars[3] = par
        self.test_value(pars[3].name(),  '98')

        # Check application parameter setting by index from end
        par      = gammalib.GApplicationPar('99','r','a','1.0','0.0','2.0','Dummy')
        pars[-2] = par
        self.test_value(pars[-2].name(), '99')

        # Return
        return

    # Test GApplicationPars class slicing
    def _test_pars_slicing(self):
        """
        Test GApplicationPars class slicing
        """
        # Setup application parameter container
        pars = self._setup_pars()

        # Perform slicing tests
        test_support._container_slicing(self, pars)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('app')

        # Append tests
        self.append(self._test_log, 'Test GLog')
        self.append(self._test_app, 'Test GApplication')
        self.append(self._test_pars, 'Test GApplicationPars')
        self.append(self._test_pars_access, 'Test GApplicationPars parameter access')
        self.append(self._test_pars_access, 'Test GApplicationPars slicing')

        # Return
        return
