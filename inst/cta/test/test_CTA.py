# ==========================================================================
# This module performs unit tests for the GammaLib CTA module.
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
import math
import test_support


# ================================== #
# Test class for GammaLib CTA module #
# ================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib CTA module
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Set test directories
        self._data  = os.environ['TEST_CTA_DATA']
        self._caldb = self._data + '/../caldb'

        # Return
        return

    # Setup GCTAEventList container
    def _setup_eventlist(self):
        """
        Setup GCTAEventList container

        Returns
        -------
        list : `~gammalib.GCTAEventList`
            GCTAEventList container
        """
        # Setup event list container
        list = gammalib.GCTAEventList()
        for i in range(10):
            dir    = gammalib.GCTAInstDir(gammalib.GSkyDir(),float(i),float(i))
            energy = gammalib.GEnergy(float(i),'MeV')
            time   = gammalib.GTime(float(i),'sec')
            atom   = gammalib.GCTAEventAtom(dir, energy, time)
            list.append(atom)

        # Return event list container
        return list

    # Setup GCTAEventCube
    def _setup_eventcube(self):
        """
        Setup GCTAEventCube

        Returns
        -------
        cube : `~gammalib.GCTAEventCube`
            GCTAEventCube
        """
        # Setup event cube
        map    = gammalib.GSkyMap('CAR','CEL',0.0,0.0,1.0,1.0,5,5,2)
        ebds   = gammalib.GEbounds(2, gammalib.GEnergy(1.0, 'TeV'),
                                      gammalib.GEnergy(10.0, 'TeV'))
        gti    = gammalib.GGti(gammalib.GTime(0.0,'sec'),
                               gammalib.GTime(1.0,'sec'))
        cube   = gammalib.GCTAEventCube(map, ebds, gti)
        counts = 0.0
        for bin in cube:
            counts += 1.0
            bin.counts(1.0)

        # Return event cube
        return cube

    # Test GCTAEventList class access operators
    def _test_eventlist_access(self):
        """
        Test GCTAEventList class observation access
        """
        # Setup event list container
        list = self._setup_eventlist()

        # Perform event list access tests
        test_support.energy_container_access_index(self, list)

        # Return
        return

    # Test GCTAEventList class slicing
    def _test_eventlist_slicing(self):
        """
        Test GCTAEventList class slicing
        """
        # Setup event list container
        list = self._setup_eventlist()

        # Perform slicing tests
        test_support.energy_container_slicing(self, list)

        # Return
        return

    # Test effective area response
    def _test_aeff(self):
        """
        Test GCTAAeff classes
        """
        # Test GCTAAeff2D file constructor
        filename = self._caldb + '/prod1_gauss.fits'
        aeff     = gammalib.GCTAAeff2D(filename)

        # Test Aeff values
        self.test_value(aeff(0.0, 0.0), 5535774176.75, 0.1,
                        'Test reference effective area value')
        self.test_value(aeff(1.0, 0.0), 20732069462.7, 0.1,
                        'Test reference effective area value')
        self.test_value(aeff(0.0, 0.01745), 5682897797.76, 0.1,
                        'Test reference effective area value')
        self.test_value(aeff(1.0, 0.01745), 18446656815.1, 0.1,
                        'Test reference effective area value')

        # Test that Aeff values outside boundaries are zero
        self.test_value(aeff(-1.80001, 0.0), 0.0, 1.0e-6,
                        'Test that effective area is zero for energy below'
                        ' minimum energy')
        self.test_value(aeff(+2.20001, 0.0), 0.0, 1.0e-6,
                        'Test that effective area is zero for energy above'
                        ' maximum energy')
        self.test_value(aeff(0.0, -0.00001), 0.0, 1.0e-6,
                        'Test that effective area is zero for offset angle'
                        ' below minimum offset angle')
        self.test_value(aeff(0.0,  0.13963), 0.0, 1.0e-6,
                        'Test that effective area is zero for offset angle'
                        ' above maximum offset angle')

        # Test GCTAAeffPerfTable file constructor
        filename = self._caldb + '/cta_dummy_irf.dat'
        aeff     = gammalib.GCTAAeffPerfTable(filename)

        # Test Aeff values
        self.test_value(aeff(0.0, 0.0), 2738898000.0, 0.1)
        self.test_value(aeff(1.0, 0.0), 16742420500.0, 0.1)
        self.test_value(aeff(0.0, 0.01745), 2590995083.29, 0.1)
        self.test_value(aeff(1.0, 0.01745), 15838314971.2, 0.1)

        # Test GCTAAeffArf file constructor
        filename = self._caldb + '/dc1/arf.fits'
        aeff     = gammalib.GCTAAeffArf(filename)

        # Test Aeff values
        self.test_value(aeff(0.0, 0.0), 1607246236.98, 0.1)
        self.test_value(aeff(1.0, 0.0), 4582282342.98, 0.1)
        self.test_value(aeff(0.0, 0.01745), 1607246236.98, 0.1)
        self.test_value(aeff(1.0, 0.01745), 4582282342.98, 0.1)

        # Return
        return

    # Test point spread function response
    def _test_psf(self):
        """
        Test GCTAPsf classes
        """
        # Test GCTAPsf2D file constructor
        filename = self._caldb + '/prod1_gauss.fits'
        psf      = gammalib.GCTAPsf2D(filename)

        # Test PSF values
        self.test_value(psf(0.0, 0.0, 0.0), 163782.469465, 1.0e-6)
        self.test_value(psf(0.001, 0.0, 0.0), 97904.9307797, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.0), 616076.98558, 1.0e-6)
        self.test_value(psf(0.001, 1.0, 0.0), 88932.681708, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.01745), 433247.309504, 1.0e-6)
        self.test_value(psf(0.001, 1.0, 0.01745), 111075.0692681, 1.0e-6)

        # Test GCTAPsfKing file constructor
        filename = self._caldb + '/prod1_king.fits'
        psf      = gammalib.GCTAPsfKing(filename)

        # Test PSF values
        self.test_value(psf(0.0, 0.0, 0.0), 213616.312600672, 1.0e-6)
        self.test_value(psf(0.001, 0.0, 0.0), 90918.3030269623, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.0), 1126804.99931516, 1.0e-5)
        self.test_value(psf(0.001, 1.0, 0.0), 54873.6646449112, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.01745), 660972.636049452, 1.0e-6)
        self.test_value(psf(0.001, 1.0, 0.01745), 80272.4048345619, 1.0e-6)

        # Test GCTAPsfPerfTable file constructor
        filename = self._caldb + '/cta_dummy_irf.dat'
        psf      = gammalib.GCTAPsfPerfTable(filename)

        # Test PSF values
        self.test_value(psf(0.0, 0.0, 0.0), 537853.354917, 1.0e-6)
        self.test_value(psf(0.001, 0.0, 0.0), 99270.360144, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.0), 1292604.7473727, 1.0e-6)
        self.test_value(psf(0.001, 1.0, 0.0), 22272.4258111, 1.0e-6)
        self.test_value(psf(0.0, 1.0, 0.01745), 1292604.7473727, 1.0e-6)
        self.test_value(psf(0.001, 1.0, 0.01745), 22272.4258111, 1.0e-6)

        # Test GCTAPsfVector file constructor
        filename = self._caldb + '/dc1/psf_magic.fits'
        psf      = gammalib.GCTAPsfVector(filename)

        # Print
        self.test_value(psf(0.0, -1.0, 0.0), 42263.9572394, 1.0e-6)
        self.test_value(psf(0.001, -1.0, 0.0), 37008.8652966, 1.0e-6)
        self.test_value(psf(0.0, 0.0, 0.0), 208989.164294, 1.0e-6)
        self.test_value(psf(0.001, 0.0, 0.0), 108388.031915, 1.0e-6)
        self.test_value(psf(0.0, 0.0, 0.01745), 208989.164294, 1.0e-6)
        self.test_value(psf(0.001, 0.0, 0.01745), 108388.031915, 1.0e-6)

        # Return
        return

    # Test energy dispersion
    def _test_edisp(self):
        """
        Test GCTAEdisp classes
        """
        # Set some energies
        eng0    = gammalib.GEnergy(1.0, 'TeV')
        eng01   = gammalib.GEnergy(1.2589254, 'TeV')
        eng001  = gammalib.GEnergy(1.023293, 'TeV')
        eng0001 = gammalib.GEnergy(1.0023052, 'TeV')
        eng1001 = gammalib.GEnergy(10.023052, 'TeV')
        eng1    = gammalib.GEnergy(10.0, 'TeV')
        eng30   = gammalib.GEnergy(30.0, 'TeV')

        # Test GCTAEdispRmf file constructor
        filename = self._caldb + '/dc1/rmf.fits'
        edisp    = gammalib.GCTAEdispRmf(filename)

        # Test energy dispersion values
        self.test_value(edisp(eng30, eng0), 0.0, 1.0e-9)
        self.test_value(edisp(eng0, eng30), 0.0, 1.0e-9)

        # Test GCTAEdispPerfTable file constructor
        filename = self._caldb + '/cta_dummy_irf.dat'
        edisp    = gammalib.GCTAEdispPerfTable(filename)

        # Test energy dispersion values
        self.test_value(edisp(eng0,    eng0),      4.3386871e-06, 1.0e-6)
        self.test_value(edisp(eng0001, eng0),      4.3273514e-06, 1.0e-6)
        self.test_value(edisp(eng001,  eng0),      4.1090489e-06, 1.0e-6)
        self.test_value(edisp(eng01,   eng0),      1.4984969e-07, 1.0e-6)
        self.test_value(edisp(eng1,    eng1),      7.8454726e-07, 1.0e-6)
        self.test_value(edisp(eng1001, eng1, 0.0), 7.8194077e-07, 1.0e-6)

        # Test GCTAResponseIrf file constructor
        db  = gammalib.GCaldb(self._caldb)
        irf = 'cta_dummy_irf'
        rsp = gammalib.GCTAResponseIrf(irf, db)

        # Test nedisp computations
        #dir  = gammalib.GSkyDir()
        #pnt  = gammalib.GCTAPointing()
        #time = gammalib.GTime()
        #self.test_value(rsp.nedisp(dir, GEnergy(3.7, "TeV"), time, pnt, 
        #                           GEbounds(GEnergy(0.1, "TeV"),
        #                                    GEnergy(10.0, "TeV"))),
        #                1.0, 0.005)
        #self.test_value(rsp.nedisp(dir, GEnergy(3.7, "TeV"), time, pnt, 
        #                           GEbounds(GEnergy(2.72345,  "TeV"),
        #                                    GEnergy(5.026615, "TeV"))),
        #                1.0, 0.005)
        #self.test_value(rsp.nedisp(dir, GEnergy(3.7, "TeV"), time, pnt, 
        #                           GEbounds(GEnergy(3.7, "TeV"),
        #                                    GEnergy(10.0, "TeV"))),
        #                0.5, 0.005)

    # Test response
    def _test_response(self):
        """
        Test response classes
        """
        # Load 1DC CTA observation (ARF, PSF, RMF)
        filename = self._data + '/irf_1dc.xml'
        obs      = gammalib.GObservations(filename)

        # Return
        return

    # Test On/Off analysis
    def _test_onoff(self):
        """
        Test On/Off analysis
        """
        # Load model container
        models = gammalib.GModels(self._data + '/crab_irf.xml')

        # Create On region
        ondir = gammalib.GSkyDir()
        ondir.radec_deg(83.6331, 22.0145)
        on = gammalib.GSkyRegions()
        on.append(gammalib.GSkyRegionCircle(ondir, 0.2))

        # Create Off region
        offdir = gammalib.GSkyDir()
        offdir.radec_deg(83.6331, 23.5145)
        off = gammalib.GSkyRegions()
        off.append(gammalib.GSkyRegionCircle(offdir, 0.5))

        # Set energy binning
        etrue = gammalib.GEbounds(40, gammalib.GEnergy(0.1,  'TeV'),
                                      gammalib.GEnergy(10.0, 'TeV'))
        ereco = gammalib.GEbounds(20, gammalib.GEnergy(0.1,  'TeV'),
                                      gammalib.GEnergy(10.0, 'TeV'))

        # Create On/Off observations from CTA observations
        filename = self._data + '/irf_unbinned.xml'
        inobs    = gammalib.GObservations(filename)
        outobs   = gammalib.GObservations()
        for run in inobs:
            onoff = gammalib.GCTAOnOffObservation(run, models, 'Crab',
                                                  etrue, ereco, on, off)
            outobs.append(onoff)

        # Load On/Off models and attach them to the observations
        models = gammalib.GModels(self._data + '/onoff_model.xml')
        outobs.models(models)

        # Perform maximum likelihood fit
        lm = gammalib.GOptimizerLM()
        outobs.optimize(lm)
        outobs.errors(lm)

        # Test On/Off model fitting results
        sky = outobs.models()['Crab']
        bgd = outobs.models()['Background']
        self.test_value(sky['Prefactor'].value(), 6.456877e-16, 1.0e-18,
                        'Check sky model prefactor value')
        self.test_value(sky['Prefactor'].error(), 2.176260e-17, 1.0e-20,
                        'Check sky model prefactor error')
        self.test_value(sky['Index'].value(), -2.575639, 1.0e-4,
                        'Check sky model index value')
        self.test_value(sky['Index'].error(), 0.030702, 1.0e-4,
                        'Check sky model index error')
        self.test_value(bgd['Prefactor'].value(), 1.182291, 0.01,
                        'Check background model prefactor value')
        self.test_value(bgd['Prefactor'].error(), 0.152625, 0.01,
                        'Check background model prefactor error')
        self.test_value(bgd['Index'].value(), 0.520937, 0.01,
                        'Check background model index value')
        self.test_value(bgd['Index'].error(), 0.086309, 0.01,
                        'Check background model index error')

        # Save PHA, ARF and RMFs
        for run in outobs:
            run.on_spec().save('test_cta_onoff_pha_on.fits', True)
            run.off_spec().save('test_cta_onoff_pha_off.fits', True)
            run.arf().save('test_cta_onoff_arf.fits', True)
            run.rmf().save('test_cta_onoff_rmf.fits', True)

        # Save On/Off observations
        outobs.save('test_cta_onoff.xml')

        # Return
        return

    # Test class pickeling
    def _test_pickeling(self):
        """
        Test class pickeling
        """
        # Set CALDB
        os.environ['CALDB'] = os.environ['TEST_DATA']+'/caldb'

        # Perform pickeling tests of empty classes
        test_support.pickeling(self, gammalib.GCTAAeff2D())
        #test_support.pickeling(self, gammalib.GCTAAeffArf())
        #test_support.pickeling(self, gammalib.GCTAAeffPerfTable())
        test_support.pickeling(self, gammalib.GCTABackground3D())
        #test_support.pickeling(self, gammalib.GCTABackgroundPerfTable())
        test_support.pickeling(self, gammalib.GCTACubeBackground())
        test_support.pickeling(self, gammalib.GCTACubeEdisp())
        test_support.pickeling(self, gammalib.GCTACubeExposure())
        test_support.pickeling(self, gammalib.GCTACubePsf())
        #test_support.pickeling(self, gammalib.GCTACubeSourceDiffuse())
        #test_support.pickeling(self, gammalib.GCTACubeSourcePoint())
        test_support.pickeling(self, gammalib.GCTAEdisp2D())
        #test_support.pickeling(self, gammalib.GCTAEdispPerfTable())
        #test_support.pickeling(self, gammalib.GCTAEdispRmf())
        test_support.pickeling(self, gammalib.GCTAEventAtom())
        test_support.pickeling(self, gammalib.GCTAEventBin())
        test_support.pickeling(self, gammalib.GCTAEventCube())
        test_support.pickeling(self, gammalib.GCTAEventList())
        test_support.pickeling(self, gammalib.GCTAInstDir())
        test_support.pickeling(self, gammalib.GCTAModelSkyCube())
        test_support.pickeling(self, gammalib.GCTAModelBackground())
        test_support.pickeling(self, gammalib.GCTAModelAeffBackground())
        test_support.pickeling(self, gammalib.GCTAModelCubeBackground())
        test_support.pickeling(self, gammalib.GCTAModelIrfBackground())
        test_support.pickeling(self, gammalib.GCTAModelRadialAcceptance())
        test_support.pickeling(self, gammalib.GCTAModelRadialGauss())
        test_support.pickeling(self, gammalib.GCTAModelRadialPolynom())
        test_support.pickeling(self, gammalib.GCTAModelRadialProfile())
        test_support.pickeling(self, gammalib.GCTAModelSpatialGradient())
        test_support.pickeling(self, gammalib.GCTAModelSpatialMultiplicative())
        test_support.pickeling(self, gammalib.GCTAObservation())
        test_support.pickeling(self, gammalib.GCTAOnOffObservation())
        test_support.pickeling(self, gammalib.GCTAPointing())
        test_support.pickeling(self, gammalib.GCTAPsf2D())
        test_support.pickeling(self, gammalib.GCTAPsfKing())
        #test_support.pickeling(self, gammalib.GCTAPsfPerfTable())
        test_support.pickeling(self, gammalib.GCTAPsfTable())
        #test_support.pickeling(self, gammalib.GCTAPsfVector())
        test_support.pickeling(self, gammalib.GCTAResponseCube())
        test_support.pickeling(self, gammalib.GCTAResponseIrf())
        test_support.pickeling(self, gammalib.GCTAResponseTable())
        test_support.pickeling(self, gammalib.GCTARoi())

        # Setup test
        list    = self._setup_eventlist()
        cube    = self._setup_eventcube()
        atom    = list[0]
        dir     = gammalib.GSkyDir()
        instdir = gammalib.GCTAInstDir(dir, 2.0, -3.0)
        pivot   = gammalib.GEnergy(1.0,'TeV')
        plaw    = gammalib.GModelSpectralPlaw(1.0,-2.0,pivot)
        irfname = os.environ['TEST_DATA']+'/caldb/data/cta/prod2/bcf/North_0.5h/irf_file.fits.gz'
        bin     = gammalib.GCTAEventBin()
        bin.dir(instdir)
        bin.energy(pivot)
        bin.time(gammalib.GTime(1.0,'sec'))
        bin.counts(1.0)
        bin.solidangle(0.1)
        bin.ewidth(pivot)
        bin.ontime(100.0)
        bin.weight(1.0)
        emin      = gammalib.GEnergy(1.0,'TeV')
        emax      = gammalib.GEnergy(10.0,'TeV')
        engs      = gammalib.GEnergies(10,emin,emax)
        ebds      = gammalib.GEbounds(2,emin,emax)
        region    = gammalib.GSkyRegionCircle(dir, 0.2)
        regs      = gammalib.GSkyRegions()
        regs.append(region)
        models    = gammalib.GModels(self._data + '/crab_irf.xml')
        expcube   = gammalib.GCTACubeExposure('CAR','CEL',0.,0.,0.1,0.1,10,10,engs)
        psfcube   = gammalib.GCTACubePsf('CAR','CEL',0.,0.,0.1,0.1,10,10,engs,1.0,10)
        bgdcube   = gammalib.GCTACubeBackground('CAR','CEL',0.,0.,0.1,0.1,10,10,ebds)
        edispcube = gammalib.GCTACubeEdisp('CAR','CEL',0.,0.,0.1,0.1,10,10,engs,1.0,10)
        caldb     = gammalib.GCaldb('cta','prod2')
        rspirf    = gammalib.GCTAResponseIrf('North_0.5h', caldb)
        rspcube1  = gammalib.GCTAResponseCube(expcube, psfcube, bgdcube)
        rspcube2  = gammalib.GCTAResponseCube(expcube, psfcube, edispcube, bgdcube)
        bin.energy(pivot)
        obs1      = gammalib.GCTAObservation()
        obs1.events(list)
        obs1.response(rspirf)
        obs2      = gammalib.GCTAObservation()
        obs2.events(cube)
        obs2.response(rspirf)
        obs3      = gammalib.GCTAObservation()
        obs3.events(cube)
        obs3.response(rspcube1)
        obs4      = gammalib.GCTAObservation()
        obs4.events(cube)
        obs4.response(rspcube2)
        obs5      = gammalib.GCTAOnOffObservation(obs1, models, 'Crab', ebds, ebds, regs, regs)
        radgauss  = gammalib.GCTAModelRadialGauss(1.0)
        radacc    = gammalib.GCTAModelRadialAcceptance(radgauss, plaw)
        multi     = gammalib.GCTAModelSpatialMultiplicative()
        multi.append(radgauss, 'Src1')
        multi.append(radgauss, 'Src2')
        kingname  = self._caldb + '/prod1_king.fits'
        hessaeff  = self._data  + '/irf_hess_aeff.fits.gz'
        hesspsf   = self._data  + '/irf_hess_psf.fits.gz'
        hessedisp = self._data  + '/irf_hess_edisp.fits.gz'
        hessbkg   = self._data  + '/irf_hess_bkg.fits.gz'
        hessirf   = gammalib.GCTAResponseIrf()
        hessirf.load_aeff(hessaeff)
        hessirf.load_psf(hesspsf)
        hessirf.load_edisp(hessedisp)
        hessirf.load_background(hessbkg)
        skycube = self._data + '/crab_modcube.fits.gz'

        # Perform pickeling tests of filled classes
        test_support.pickeling(self, gammalib.GCTAAeff2D(irfname))
        #test_support.pickeling(self, gammalib.GCTAAeffArf())
        #test_support.pickeling(self, gammalib.GCTAAeffPerfTable())
        test_support.pickeling(self, gammalib.GCTABackground3D(irfname))
        #test_support.pickeling(self, gammalib.GCTABackgroundPerfTable())
        test_support.pickeling(self, gammalib.GCTACubeBackground(bgdcube))
        test_support.pickeling(self, gammalib.GCTACubeEdisp(edispcube))
        test_support.pickeling(self, gammalib.GCTACubeExposure(expcube))
        test_support.pickeling(self, gammalib.GCTACubePsf(psfcube))
        #test_support.pickeling(self, gammalib.GCTACubeSourceDiffuse())
        #test_support.pickeling(self, gammalib.GCTACubeSourcePoint())
        test_support.pickeling(self, gammalib.GCTAEdisp2D(irfname))
        #test_support.pickeling(self, gammalib.GCTAEdispPerfTable())
        #test_support.pickeling(self, gammalib.GCTAEdispRmf())
        test_support.pickeling(self, gammalib.GCTAEventAtom(atom))
        test_support.pickeling(self, gammalib.GCTAEventBin(bin))
        test_support.pickeling(self, gammalib.GCTAEventCube(cube))
        test_support.pickeling(self, gammalib.GCTAEventList(list))
        test_support.pickeling(self, gammalib.GCTAInstDir(instdir))
        test_support.pickeling(self, gammalib.GCTAModelSkyCube(skycube,plaw))
        test_support.pickeling(self, gammalib.GCTAModelBackground(radgauss,plaw))
        test_support.pickeling(self, gammalib.GCTAModelAeffBackground(plaw))
        test_support.pickeling(self, gammalib.GCTAModelCubeBackground(plaw))
        test_support.pickeling(self, gammalib.GCTAModelIrfBackground(plaw))
        test_support.pickeling(self, gammalib.GCTAModelRadialAcceptance(radacc))
        test_support.pickeling(self, gammalib.GCTAModelRadialGauss(1.0))
        test_support.pickeling(self, gammalib.GCTAModelRadialPolynom([1.0,2.0]))
        test_support.pickeling(self, gammalib.GCTAModelRadialProfile(1.0,2.0,3.0))
        test_support.pickeling(self, gammalib.GCTAModelSpatialGradient(1.0,2.0))
        test_support.pickeling(self, gammalib.GCTAModelSpatialMultiplicative(multi))
        test_support.pickeling(self, gammalib.GCTAObservation(obs1))
        test_support.pickeling(self, gammalib.GCTAObservation(obs2))
        test_support.pickeling(self, gammalib.GCTAObservation(obs3))
        test_support.pickeling(self, gammalib.GCTAObservation(obs4))
        test_support.pickeling(self, gammalib.GCTAOnOffObservation(obs5))
        test_support.pickeling(self, gammalib.GCTAPointing(dir))
        test_support.pickeling(self, gammalib.GCTAPsf2D(irfname))
        test_support.pickeling(self, gammalib.GCTAPsfKing(kingname))
        #test_support.pickeling(self, gammalib.GCTAPsfPerfTable())
        test_support.pickeling(self, gammalib.GCTAPsfTable(hesspsf))
        #test_support.pickeling(self, gammalib.GCTAPsfVector())
        test_support.pickeling(self, gammalib.GCTAResponseCube(rspcube1))
        test_support.pickeling(self, gammalib.GCTAResponseCube(rspcube2))
        test_support.pickeling(self, gammalib.GCTAResponseIrf(rspirf))
        test_support.pickeling(self, gammalib.GCTAResponseIrf(hessirf))
        #test_support.pickeling(self, gammalib.GCTAResponseTable()) # No constructor
        test_support.pickeling(self, gammalib.GCTARoi(instdir,2.0))

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('CTA')

        # Append tests
        self.append(self._test_eventlist_access, 'Test GCTAEventList event access')
        self.append(self._test_eventlist_slicing, 'Test GCTAEventList slicing')
        self.append(self._test_aeff, 'Test CTA effective area classes')
        self.append(self._test_psf, 'Test CTA PSF classes')
        self.append(self._test_edisp, 'Test CTA energy dispersion classes')
        self.append(self._test_response, 'Test CTA response classes')
        self.append(self._test_onoff, 'Test CTA On/Off analysis')
        self.append(self._test_pickeling, 'Test CTA class pickeling')

        # Return
        return
