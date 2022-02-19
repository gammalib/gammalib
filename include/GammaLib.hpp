/***************************************************************************
 *          GammaLib.hpp - Gamma-Ray Astronomy Library Header file         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2022 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GammaLib.hpp
 * @brief GammaLib definitions
 * @author Juergen Knoedlseder
 */

#ifndef GAMMALIB_HPP
#define GAMMALIB_HPP


/***************************************************************************
 *                              Core services                              *
 ***************************************************************************/

/* __ Typemaps ___________________________________________________________ */
#include "GTypemaps.hpp"

/* __ Interface classes __________________________________________________ */
#include "GBase.hpp"
#include "GContainer.hpp"
#include "GRegistry.hpp"

/* __ Common tools _______________________________________________________ */
#include "GException.hpp"
#include "GNodeArray.hpp"
#include "GBilinear.hpp"
#include "GCsv.hpp"
#include "GRan.hpp"
#include "GUrl.hpp"
#include "GUrlFile.hpp"
#include "GUrlString.hpp"
#include "GFilename.hpp"
#include "GDaemon.hpp"

/* __ Linear algebra module ______________________________________________ */
#include "GVector.hpp"
#include "GMatrixBase.hpp"
#include "GMatrix.hpp"
#include "GMatrixSparse.hpp"
#include "GMatrixSymmetric.hpp"

/* __ Numerics module ____________________________________________________ */
#include "GIntegral.hpp"
#include "GIntegrals.hpp"
#include "GDerivative.hpp"
#include "GFunction.hpp"
#include "GFunctions.hpp"
#include "GMath.hpp"
#include "GNdarray.hpp"
#include "GFft.hpp"
#include "GFftWavetable.hpp"

/* __ FITS module ________________________________________________________ */
#include "GFits.hpp"
#include "GFitsHDU.hpp"
#include "GFitsHeader.hpp"
#include "GFitsHeaderCard.hpp"
#include "GFitsImage.hpp"
#include "GFitsImageByte.hpp"
#include "GFitsImageSByte.hpp"
#include "GFitsImageUShort.hpp"
#include "GFitsImageShort.hpp"
#include "GFitsImageULong.hpp"
#include "GFitsImageLong.hpp"
#include "GFitsImageLongLong.hpp"
#include "GFitsImageFloat.hpp"
#include "GFitsImageDouble.hpp"
#include "GFitsAsciiTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableCol.hpp"
#include "GFitsTableBitCol.hpp"
#include "GFitsTableByteCol.hpp"
#include "GFitsTableBoolCol.hpp"
#include "GFitsTableStringCol.hpp"
#include "GFitsTableUShortCol.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableULongCol.hpp"
#include "GFitsTableLongCol.hpp"
#include "GFitsTableLongLongCol.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GFitsTableCFloatCol.hpp"
#include "GFitsTableCDoubleCol.hpp"

/* __ XML module _________________________________________________________ */
#include "GXml.hpp"
#include "GXmlNode.hpp"
#include "GXmlDocument.hpp"
#include "GXmlText.hpp"
#include "GXmlElement.hpp"
#include "GXmlComment.hpp"
#include "GXmlAttribute.hpp"
#include "GXmlPI.hpp"

/* __ VO module __________________________________________________________ */
#include "GVOClient.hpp"
#include "GVOHub.hpp"
#include "GVOTable.hpp"

/* __ XSPEC module _______________________________________________________ */
#include "GArf.hpp"
#include "GPha.hpp"
#include "GRmf.hpp"

/* __ Application module _________________________________________________ */
#include "GApplication.hpp"
#include "GLog.hpp"
#include "GApplicationPars.hpp"
#include "GApplicationPar.hpp"

/* __ Optimizer module ___________________________________________________ */
#include "GOptimizer.hpp"
#include "GOptimizerLM.hpp"
#include "GOptimizerPar.hpp"
#include "GOptimizerPars.hpp"
#include "GOptimizerFunction.hpp"

/* __ Unit Test class___________________________________________________ */
#include "GTestCase.hpp"
#include "GTestSuite.hpp"
#include "GTestSuites.hpp"

/***************************************************************************
 *                        Analysis support services                        *
 ***************************************************************************/

/* __ Sky handling _______________________________________________________ */
#include "GSkyDir.hpp"
#include "GSkyDirs.hpp"
#include "GHorizDir.hpp"
#include "GSkyPixel.hpp"
#include "GSkyMap.hpp"
#include "GSkyRegions.hpp"
#include "GSkyRegion.hpp"
#include "GSkyRegionCircle.hpp"
#include "GSkyRegionRectangle.hpp"
#include "GSkyRegionMap.hpp"
#include "GSkyProjection.hpp"
#include "GHealpix.hpp"
#include "GWcsRegistry.hpp"
#include "GWcs.hpp"
#include "GWcsAIT.hpp"
#include "GWcsARC.hpp"
#include "GWcsAZP.hpp"
#include "GWcsCAR.hpp"
#include "GWcsGLS.hpp"
#include "GWcsMER.hpp"
#include "GWcsMOL.hpp"
#include "GWcsSIN.hpp"
#include "GWcsSFL.hpp"
#include "GWcsSTG.hpp"
#include "GWcsTAN.hpp"

/* __ Observation handling _______________________________________________ */
#include "GEnergy.hpp"
#include "GEnergies.hpp"
#include "GTime.hpp"
#include "GTimes.hpp"
#include "GTimeReference.hpp"
#include "GCaldb.hpp"
#include "GObservations.hpp"
#include "GObservation.hpp"
#include "GObservationRegistry.hpp"
#include "GEvents.hpp"
#include "GEventList.hpp"
#include "GEventCube.hpp"
#include "GEvent.hpp"
#include "GEventAtom.hpp"
#include "GEventBin.hpp"
#include "GGti.hpp"
#include "GRoi.hpp"
#include "GEbounds.hpp"
#include "GPhases.hpp"
#include "GInstDir.hpp"
#include "GResponse.hpp"
#include "GResponseCache.hpp"
#include "GResponseVectorCache.hpp"
#include "GPhotons.hpp"
#include "GPhoton.hpp"
#include "GSource.hpp"
#include "GPulsar.hpp"
#include "GPulsarEphemeris.hpp"
#include "GEphemerides.hpp"

/* __ Model handling _____________________________________________________ */
#include "GModelPar.hpp"
#include "GModelAssociation.hpp"
#include "GModelAssociations.hpp"
#include "GModels.hpp"
#include "GModel.hpp"
#include "GModelRegistry.hpp"
#include "GModelSky.hpp"
#include "GModelData.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelSpatialRadialDisk.hpp"
#include "GModelSpatialRadialGauss.hpp"
#include "GModelSpatialRadialGeneralGauss.hpp"
#include "GModelSpatialRadialShell.hpp"
#include "GModelSpatialRadialRing.hpp"
#include "GModelSpatialRadialProfile.hpp"
#include "GModelSpatialRadialProfileGauss.hpp"
#include "GModelSpatialRadialProfileDMBurkert.hpp"
#include "GModelSpatialRadialProfileDMEinasto.hpp"
#include "GModelSpatialRadialProfileDMZhao.hpp"
#include "GModelSpatialElliptical.hpp"
#include "GModelSpatialEllipticalDisk.hpp"
#include "GModelSpatialEllipticalGauss.hpp"
#include "GModelSpatialEllipticalGeneralGauss.hpp"
#include "GModelSpatialDiffuse.hpp"
#include "GModelSpatialDiffuseConst.hpp"
#include "GModelSpatialDiffuseCube.hpp"
#include "GModelSpatialDiffuseMap.hpp"
#include "GModelSpatialComposite.hpp"
#include "GModelSpectral.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelSpectralBrokenPlaw.hpp"
#include "GModelSpectralComposite.hpp"
#include "GModelSpectralConst.hpp"
#include "GModelSpectralExpPlaw.hpp"
#include "GModelSpectralExpInvPlaw.hpp"
#include "GModelSpectralExponential.hpp"
#include "GModelSpectralSuperExpPlaw.hpp"
#include "GModelSpectralFunc.hpp"
#include "GModelSpectralGauss.hpp"
#include "GModelSpectralLogParabola.hpp"
#include "GModelSpectralMultiplicative.hpp"
#include "GModelSpectralBins.hpp"
#include "GModelSpectralNodes.hpp"
#include "GModelSpectralPlaw.hpp"
#include "GModelSpectralPlawPhotonFlux.hpp"
#include "GModelSpectralPlawEnergyFlux.hpp"
#include "GModelSpectralSmoothBrokenPlaw.hpp"
#include "GModelSpectralTable.hpp"
#include "GModelSpectralTablePar.hpp"
#include "GModelSpectralTablePars.hpp"
#include "GModelTemporal.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GModelTemporalConst.hpp"
#include "GModelTemporalLightCurve.hpp"
#include "GModelTemporalPhaseCurve.hpp"

#endif /* GAMMALIB_HPP */
