# ==========================================================================
# This module performs unit tests for the GammaLib application module.
#
# Copyright (C) 2012-2015 Juergen Knoedlseder
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
import gammalib


# ========================================== #
# Test class for GammaLib application module #
# ========================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib application module.
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Set members
        self.logfile = "test_logger.log"

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("app")

        # Append tests
        self.append(self.test_log, "Test GLog")
        self.append(self.test_pars, "Test GApplicationPars")

        # Return
        return

    # Test GLog
    def test_log(self):
        """
        Test GLog.
        """
        # Allocate logger
        log = gammalib.GLog(self.logfile, True)

        # Test print methods
        log("Test __call__(std::string)")
        log("\n")
        log(True)
        log("\n")
        log(41)
        log("\n")
        log(3.1415)
        log("\n")
        log.parformat("This is a parameter")
        log("\n")
        log.toupper("upper case")
        log("\n")
        log.tolower("LOWER CASE")
        log("\n")
        log.fill("Higgs", 3)
        log("\n")
        log.left("Left", 10)
        log("\n")
        log.right("Right", 10)
        log("\n")
        log.centre("Centre", 10)
        log("\n")
        log.header0("Header 0")
        log("\n")
        log.header1("Header 1")
        log.header2("Header 2")
        log.header3("Header 3")

        # Close logger
        log.close()

        # Set reference
        ref = []
        ref.append("Test __call__(std::string)\n")
        #ref.append("1\n")    # Old swig version
        ref.append("true\n") # New swig version (3.x.y)
        ref.append("41\n")
        ref.append("3.1415\n")
        ref.append(" This is a parameter .......: \n")
        ref.append("UPPER CASE\n")
        ref.append("lower case\n")
        ref.append("HiggsHiggsHiggs\n")
        ref.append("Left      \n")
        ref.append("     Right\n")
        ref.append("  Centre  \n")
        ref.append("\n")
        ref.append("+==========+\n")
        ref.append("| Header 1 |\n")
        ref.append("+==========+\n")
        ref.append("+----------+\n")
        ref.append("| Header 2 |\n")
        ref.append("+----------+\n")
        ref.append("=== Header 3 ===\n")

        # Set calls for error messages
        calls = []
        calls.append("__call__(std::string)")
        calls.append("__call__(bool)")
        calls.append("__call__(int)")
        calls.append("__call__(double)")
        calls.append("parformat()")
        calls.append("toupper()")
        calls.append("tolower()")
        calls.append("fill()")
        calls.append("left()")
        calls.append("right()")
        calls.append("centre()")
        calls.append("header0()")
        calls.append("header1()")
        calls.append("header1()")
        calls.append("header1()")
        calls.append("header2()")
        calls.append("header2()")
        calls.append("header2()")
        calls.append("header3()")

        # Check file
        file = open(self.logfile, "r")
        lines = file.readlines()
        file.close()
        for i, line in enumerate(lines):
            if ref[i] == "true\n":
                self.test_assert(line == ref[i] or line == "1",
                                 "Test "+calls[i]+": "+line)
            else:
                self.test_assert(line == ref[i], "Test "+calls[i]+": "+line)

        # Return
        return

    # Test GApplicationPars
    def test_pars(self):
        """
        Test GApplicationPars.
        """
        # Test GApplicationPars constructor with bad filename.
        self.test_try("Test GApplicationPars constructor with bad filename")
        try:
            pars = gammalib.GApplicationPars("testme.par")
            self.test_try_failure("Exception not thrown.")
        except:
            self.test_try_success()
        else:
            self.test_try_failure("This should never happen.")

        # Return
        return
