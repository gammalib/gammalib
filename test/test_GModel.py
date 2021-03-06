# ==========================================================================
# This module performs unit tests for the GammaLib model module
#
# Copyright (C) 2012-2020 Juergen Knoedlseder
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


# ==================================== #
# Test class for GammaLib model module #
# ==================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib model module
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

    # Setup GModels container
    def _setup_models(self):
        """
        Setup GModels container

        Returns
        -------
        models : `~gammalib.GModels`
            Models container
        """
        # Setup model container
        models = gammalib.GModels()
        model  = gammalib.GModelSky()
        for i in range(10):
            model.name('%s' % i)
            models.append(model)

        # Return model container
        return models

    # Test GModels class access operators
    def _test_models_access(self):
        """
        Test GModels class model access
        """
        # Setup model container and sky model
        models = self._setup_models()
        model  = gammalib.GModelSky()

        # Perform model access tests
        test_support.container_access_index(self, models)

        # Check model getting by name
        self.test_value(models['7'].name(), '7')

        # Check model setting by name
        model.name('Test name')
        models['3'] = model
        self.test_value(models['Test name'].name(), 'Test name')

        # Check model setting by index from start
        model.name('Test index from start')
        models[1] = model
        self.test_value(models[1].name(), 'Test index from start')

        # Check model setting by index from end
        model.name('Test index from end')
        models[-2] = model
        self.test_value(models[-2].name(), 'Test index from end')

        # Return
        return

    # Test GModels class slicing
    def _test_models_slicing(self):
        """
        Test GModels class slicing
        """
        # Setup models container
        models = self._setup_models()

        # Perform slicing tests
        test_support.container_slicing(self, models)

        # Return
        return

    # Test class pickeling
    def _test_pickeling(self):
        """
        Test class pickeling
        """
        # Perform pickeling tests of empty classes
        test_support.pickeling(self, gammalib.GModelPar())
        test_support.pickeling(self, gammalib.GModelAssociation())
        test_support.pickeling(self, gammalib.GModelAssociations())
        test_support.pickeling(self, gammalib.GModels())
        test_support.pickeling(self, gammalib.GModelSky())
        test_support.pickeling(self, gammalib.GModelSpatialComposite())
        test_support.pickeling(self, gammalib.GModelSpatialDiffuseConst())
        test_support.pickeling(self, gammalib.GModelSpatialDiffuseCube())
        test_support.pickeling(self, gammalib.GModelSpatialDiffuseMap())
        test_support.pickeling(self, gammalib.GModelSpatialEllipticalDisk())
        test_support.pickeling(self, gammalib.GModelSpatialEllipticalGauss())
        test_support.pickeling(self, gammalib.GModelSpatialPointSource())
        test_support.pickeling(self, gammalib.GModelSpatialRadialDisk())
        test_support.pickeling(self, gammalib.GModelSpatialRadialRing())
        test_support.pickeling(self, gammalib.GModelSpatialRadialGauss())
        test_support.pickeling(self, gammalib.GModelSpatialRadialProfileDMBurkert())
        test_support.pickeling(self, gammalib.GModelSpatialRadialProfileDMEinasto())
        test_support.pickeling(self, gammalib.GModelSpatialRadialProfileDMZhao())
        test_support.pickeling(self, gammalib.GModelSpatialRadialProfileGauss())
        test_support.pickeling(self, gammalib.GModelSpatialRadialShell())
        test_support.pickeling(self, gammalib.GModelSpectralBrokenPlaw())
        test_support.pickeling(self, gammalib.GModelSpectralComposite())
        test_support.pickeling(self, gammalib.GModelSpectralConst())
        test_support.pickeling(self, gammalib.GModelSpectralExpInvPlaw())
        test_support.pickeling(self, gammalib.GModelSpectralExponential())
        test_support.pickeling(self, gammalib.GModelSpectralExpPlaw())
        test_support.pickeling(self, gammalib.GModelSpectralFunc())
        test_support.pickeling(self, gammalib.GModelSpectralGauss())
        test_support.pickeling(self, gammalib.GModelSpectralLogParabola())
        test_support.pickeling(self, gammalib.GModelSpectralMultiplicative())
        test_support.pickeling(self, gammalib.GModelSpectralBins())
        test_support.pickeling(self, gammalib.GModelSpectralNodes())
        test_support.pickeling(self, gammalib.GModelSpectralPlaw())
        test_support.pickeling(self, gammalib.GModelSpectralPlawEnergyFlux())
        test_support.pickeling(self, gammalib.GModelSpectralPlawPhotonFlux())
        test_support.pickeling(self, gammalib.GModelSpectralSmoothBrokenPlaw())
        test_support.pickeling(self, gammalib.GModelSpectralSuperExpPlaw())
        test_support.pickeling(self, gammalib.GModelSpectralTable())
        test_support.pickeling(self, gammalib.GModelTemporalConst())
        test_support.pickeling(self, gammalib.GModelTemporalLightCurve())
        test_support.pickeling(self, gammalib.GModelTemporalPhaseCurve())

        # Setup for tests
        fct     = gammalib.GFilename(os.environ['TEST_DATA']+'/filefunction.txt')
        lcrv    = gammalib.GFilename(os.environ['TEST_DATA']+'/model_temporal_lightcurve.fits')
        pcrv    = gammalib.GFilename(os.environ['TEST_DATA']+'/model_temporal_phasecurve.fits')
        cube    = gammalib.GFilename(os.environ['TEST_DATA']+'/test_cube.fits')
        table   = gammalib.GFilename(os.environ['TEST_DATA']+'/model_point_table.fits')
        dir     = gammalib.GSkyDir()
        dir.radec_deg(83.6331, 22.0145)
        time    = gammalib.GTime(1.0,'sec')
        pivot   = gammalib.GEnergy(1.0,'TeV')
        ptsrc   = gammalib.GModelSpatialPointSource(1.0,2.0)
        plaw    = gammalib.GModelSpectralPlaw(1.0,-2.0,pivot)
        cnst    = gammalib.GModelTemporalConst(3.0)
        sky     = gammalib.GModelSky(ptsrc, plaw, cnst)
        engs    = gammalib.GEnergies(10,pivot,gammalib.GEnergy(10.0,'TeV'))
        ebds    = gammalib.GEbounds(10,gammalib.GEnergy(1.0,'TeV'),gammalib.GEnergy(10.0,'TeV'))
        circle  = gammalib.GSkyRegionCircle(83.6331, 22.0145, 3.0)
        spacom  = gammalib.GModelSpatialComposite()
        spacom.append(ptsrc, 'Src1')
        spacom.append(ptsrc, 'Src2')
        speexp  = gammalib.GModelSpectralExponential()
        speexp.exponent(plaw)
        specom  = gammalib.GModelSpectralComposite()
        specom.append(plaw, 'Src1')
        specom.append(plaw, 'Src2')
        spemul  = gammalib.GModelSpectralMultiplicative()
        spemul.append(plaw, 'Src1')
        spemul.append(plaw, 'Src2')
        spebin  = gammalib.GModelSpectralBins(plaw, ebds)
        spenod  = gammalib.GModelSpectralNodes(plaw, engs)
        spamap  = gammalib.GModelSpatialDiffuseMap(cube, 2.0)
        spamap.mc_cone(circle)
        spacube = gammalib.GModelSpatialDiffuseCube(cube, 3.0)
        spacube.mc_cone(circle)
        xml1    = gammalib.GXmlElement('association name="Crab"')
        xml1.append(gammalib.GXmlElement('property name="RA" value="83.0"'))
        xml1.append(gammalib.GXmlElement('property name="DEC" value="22.0" error="1.0"'))
        xml2    = gammalib.GXmlElement('association name="Vela"')
        xml2.append(gammalib.GXmlElement('property name="RA" value="128.50" error="2.1"'))
        xml2.append(gammalib.GXmlElement('property name="DEC" value="-45.83"'))
        xml3    = gammalib.GXmlElement('associations')
        xml3.append(xml1)
        xml3.append(xml2)
        filefct = gammalib.GModelSpectralFunc()
        filefct.append(gammalib.GEnergy(1.0, 'MeV'),  1.0)
        filefct.append(gammalib.GEnergy(10.0, 'MeV'), 10.0)

        # Perform pickeling tests of filled classes
        test_support.pickeling(self, gammalib.GModelPar('Par',2.0,5.0))
        test_support.pickeling(self, gammalib.GModelAssociation(xml1))
        test_support.pickeling(self, gammalib.GModelAssociations(xml3))
        test_support.pickeling(self, gammalib.GModels(self._setup_models()))
        test_support.pickeling(self, gammalib.GModelSky(sky))
        test_support.pickeling(self, gammalib.GModelSpatialComposite(spacom))
        test_support.pickeling(self, gammalib.GModelSpatialDiffuseConst(2.0))
        test_support.pickeling(self, spacube)
        test_support.pickeling(self, spamap)
        test_support.pickeling(self, gammalib.GModelSpatialEllipticalDisk(dir,0.2,0.3,45.0))
        test_support.pickeling(self, gammalib.GModelSpatialEllipticalGauss(dir,0.2,0.3,45.0))
        test_support.pickeling(self, gammalib.GModelSpatialPointSource(ptsrc))
        test_support.pickeling(self, gammalib.GModelSpatialRadialDisk(dir,0.2))
        test_support.pickeling(self, gammalib.GModelSpatialRadialRing(dir,0.3,0.1))
        test_support.pickeling(self, gammalib.GModelSpatialRadialGauss(dir,0.2))
        #test_support.pickeling(self, gammalib.GModelSpatialRadialProfileDMBurkert())
        #test_support.pickeling(self, gammalib.GModelSpatialRadialProfileDMEinasto())
        #test_support.pickeling(self, gammalib.GModelSpatialRadialProfileDMZhao())
        test_support.pickeling(self, gammalib.GModelSpatialRadialProfileGauss(dir,0.2))
        test_support.pickeling(self, gammalib.GModelSpatialRadialShell(dir,0.3,0.1))
        test_support.pickeling(self, gammalib.GModelSpectralBrokenPlaw(1.0,-2.0,pivot,-2.1))
        test_support.pickeling(self, gammalib.GModelSpectralComposite(specom))
        test_support.pickeling(self, gammalib.GModelSpectralConst(3.0))
        test_support.pickeling(self, gammalib.GModelSpectralExpInvPlaw(1.0,-2.0,pivot,0.9))
        test_support.pickeling(self, gammalib.GModelSpectralExponential(speexp))
        test_support.pickeling(self, gammalib.GModelSpectralExpPlaw(1.0,-2.0,pivot,pivot))
        test_support.pickeling(self, gammalib.GModelSpectralFunc(fct, 2.0))
        test_support.pickeling(self, gammalib.GModelSpectralFunc(filefct))
        test_support.pickeling(self, gammalib.GModelSpectralGauss(1.0,pivot,pivot))
        test_support.pickeling(self, gammalib.GModelSpectralLogParabola(1.0,-2.0,pivot,0.9))
        test_support.pickeling(self, gammalib.GModelSpectralMultiplicative(spemul))
        test_support.pickeling(self, gammalib.GModelSpectralBins(spebin))
        test_support.pickeling(self, gammalib.GModelSpectralNodes(spenod))
        test_support.pickeling(self, gammalib.GModelSpectralPlaw(plaw))
        test_support.pickeling(self, gammalib.GModelSpectralPlawEnergyFlux(1.0,-2.0,pivot,pivot))
        test_support.pickeling(self, gammalib.GModelSpectralPlawPhotonFlux(1.0,-2.0,pivot,pivot))
        test_support.pickeling(self, gammalib.GModelSpectralSmoothBrokenPlaw(1.0,-2.0,pivot,-2.1,pivot,0.9))
        test_support.pickeling(self, gammalib.GModelSpectralSuperExpPlaw(1.0,-2.0,pivot,pivot,0.9))
        test_support.pickeling(self, gammalib.GModelSpectralTable(table, 2.0))
        test_support.pickeling(self, gammalib.GModelTemporalConst(cnst))
        test_support.pickeling(self, gammalib.GModelTemporalLightCurve(lcrv, 2.0))
        test_support.pickeling(self, gammalib.GModelTemporalPhaseCurve(pcrv, time, 0.5, 1.0, 1.0, 1.0, 2.0))

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('model')

        # Append tests
        self.append(self._test_models_access, 'Test GModels model access')
        self.append(self._test_models_slicing, 'Test GModels slicing')
        self.append(self._test_pickeling, 'Test model class pickeling')

        # Return
        return
