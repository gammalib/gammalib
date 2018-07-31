# ==========================================================================
# This module performs unit tests for the GammaLib COMPTEL module.
#
# Copyright (C) 2012-2018 Juergen Knoedlseder
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


# ====================================== #
# Test class for GammaLib COMPTEL module #
# ====================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib COMPTEL module
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Test class pickeling
    def _test_pickeling(self):
        """
        Test class pickeling
        """
        # Set CALDB
        os.environ['CALDB'] = os.environ['TEST_COM_DATA']+'/../caldb'

        # Perform pickeling tests of empty classes
        test_support.pickeling(self, gammalib.GCOMD1Response())
        test_support.pickeling(self, gammalib.GCOMD2Response())
        test_support.pickeling(self, gammalib.GCOMDri())
        test_support.pickeling(self, gammalib.GCOMEventAtom())
        test_support.pickeling(self, gammalib.GCOMEventBin())
        test_support.pickeling(self, gammalib.GCOMEventCube())
        test_support.pickeling(self, gammalib.GCOMEventList())
        #test_support.pickeling(self, gammalib.GCOMIaq())
        #test_support.pickeling(self, gammalib.GCOMInstChars())
        test_support.pickeling(self, gammalib.GCOMInstDir())
        test_support.pickeling(self, gammalib.GCOMModelDRBFitting())
        test_support.pickeling(self, gammalib.GCOMOad())
        test_support.pickeling(self, gammalib.GCOMOads())
        test_support.pickeling(self, gammalib.GCOMObservation())
        test_support.pickeling(self, gammalib.GCOMResponse())
        test_support.pickeling(self, gammalib.GCOMRoi())
        #test_support.pickeling(self, gammalib.GCOMSelection())
        #test_support.pickeling(self, gammalib.GCOMStatus())
        test_support.pickeling(self, gammalib.GCOMTim())

        # Setup test
        caldb   = gammalib.GCaldb('cgro','comptel')
        dir     = gammalib.GSkyDir()
        instdir = gammalib.GCOMInstDir(dir, 10.0)
        energy  = gammalib.GEnergy(1.0,'MeV')
        time    = gammalib.GTime(1.0,'sec')
        map     = gammalib.GSkyMap('CAR','GAL',0.0,0.0,1.0,1.0,10,10,10)
        dri     = gammalib.GCOMDri(map, 0.0, 2.0, 10)
        atom    = gammalib.GCOMEventAtom()
        ebds    = gammalib.GEbounds(1,gammalib.GEnergy(0.75,'MeV'),gammalib.GEnergy(30.0,'MeV'))
        gti     = gammalib.GGti(gammalib.GTime(1.0,'secs'),gammalib.GTime(3.0,'secs'))
        dre     = gammalib.GCOMDri(map, 0.0, 2.0, 10)
        obs     = gammalib.GObservations(os.environ['TEST_COM_DATA']+'/obs.xml')
        model   = gammalib.GModels(os.environ['TEST_COM_DATA']+'/crab.xml')
        dre.ebounds(ebds)
        dre.gti(gti)
        atom.dir(instdir)
        atom.energy(energy)
        atom.time(time)
        atom.e1(100.0)
        atom.e2(200.0)
        atom.phibar(11.1)
        atom.theta(22.2)
        atom.phi(33.3)
        atom.eha(44.4)
        atom.psd(55.5)
        atom.tof(66.6)
        atom.modcom(1)
        atom.reflag(2)
        atom.veto(3)
        list = gammalib.GCOMEventList()
        list.append(atom)
        bin = gammalib.GCOMEventBin()
        bin.dir(instdir)
        bin.energy(energy)
        bin.time(time)
        bin.counts(11.0)
        bin.solidangle(11.0)
        bin.ewidth(energy)
        bin.ontime(111.0)
        oad = gammalib.GCOMOad()
        oad.tstart(gammalib.GTime(1.0,'secs'))
        oad.tstop(gammalib.GTime(3.0,'secs'))
        oad.zaxis(dir)
        oad.xaxis(dir)
        oad.tjd(1000)
        oad.tics(100)
        oad.gcaz(10.0)
        oad.gcel(20.0)
        oad.georad(30.0)
        oads = gammalib.GCOMOads()
        oads.append(oad)

        # Take provision for rounding errors
        fits = gammalib.GFits()
        dre.write(fits)
        dre.read(fits[0])

        # Perform pickeling tests of filled classes
        test_support.pickeling(self, gammalib.GCOMD1Response(caldb,'DEFAULT'))
        test_support.pickeling(self, gammalib.GCOMD2Response(caldb,'DEFAULT'))
        test_support.pickeling(self, gammalib.GCOMDri(dri))
        test_support.pickeling(self, gammalib.GCOMEventAtom(atom))
        test_support.pickeling(self, gammalib.GCOMEventBin(bin))
        test_support.pickeling(self, gammalib.GCOMEventCube(dre))
        test_support.pickeling(self, gammalib.GCOMEventList(list))
        #test_support.pickeling(self, gammalib.GCOMIaq())
        #test_support.pickeling(self, gammalib.GCOMInstChars(caldb,'DEFAULT'))
        test_support.pickeling(self, gammalib.GCOMInstDir(instdir))
        test_support.pickeling(self, gammalib.GCOMModelDRBFitting(model['Background']))
        test_support.pickeling(self, gammalib.GCOMOad(oad))
        test_support.pickeling(self, gammalib.GCOMOads(oads))
        test_support.pickeling(self, gammalib.GCOMObservation(obs[0]))
        test_support.pickeling(self, gammalib.GCOMResponse(caldb, 'UNH(1.0-3.0)MeV'))
        test_support.pickeling(self, gammalib.GCOMRoi(instdir, 10.0, 1.0, 5.0))
        #test_support.pickeling(self, gammalib.GCOMSelection())
        #test_support.pickeling(self, gammalib.GCOMStatus())
        test_support.pickeling(self, gammalib.GCOMTim(gti))

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('COMPTEL')

        # Append tests
        self.append(self._test_pickeling, 'Test COMPTEL class pickeling')

        # Return
        return
