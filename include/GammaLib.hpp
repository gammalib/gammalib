/***************************************************************************
 *         GammaLib.hpp  -  Gamma-Ray Astronomy Library Header file        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2013 by Juergen Knoedlseder                         *
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

/* __ Interface classes __________________________________________________ */
#include "GBase.hpp"
#include "GRegistry.hpp"

/* __ Common tools _______________________________________________________ */
#include "GException.hpp"
#include "GNodeArray.hpp"
#include "GCsv.hpp"
#include "GRan.hpp"

/* __ Numerics module ____________________________________________________ */
#include "GVector.hpp"
#include "GMatrixBase.hpp"
#include "GMatrix.hpp"
#include "GSymMatrix.hpp"
#include "GSparseMatrix.hpp"
#include "GIntegral.hpp"
#include "GDerivative.hpp"
#include "GFunction.hpp"
#include "GNumerics.hpp"

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

/* __ Application module _________________________________________________ */
#include "GApplication.hpp"
#include "GLog.hpp"
#include "GPars.hpp"
#include "GPar.hpp"

/* __ Optimizer module ___________________________________________________ */
#include "GOptimizer.hpp"
#include "GOptimizerLM.hpp"
#include "GOptimizerPars.hpp"
#include "GOptimizerFunction.hpp"

/* __ Unit Test class___________________________________________________ */
#include "GTestCase.hpp"
#include "GTestSuite.hpp"
#include "GTestSuites.hpp"

/***************************************************************************
 *                        Analysis support services                        *
 ***************************************************************************/

/* __ Skymap handling ____________________________________________________ */
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"
#include "GSkymap.hpp"
#include "GWcs.hpp"
#include "GWcsRegistry.hpp"
#include "GWcslib.hpp"
#include "GWcsAZP.hpp"
#include "GWcsCAR.hpp"
#include "GWcsMER.hpp"
#include "GWcsTAN.hpp"
#include "GWcsSTG.hpp"
#include "GWcsHPX.hpp"

/* __ Observation handling _______________________________________________ */
#include "GEnergy.hpp"
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
#include "GInstDir.hpp"
#include "GPointing.hpp"
#include "GResponse.hpp"
#include "GPhotons.hpp"
#include "GPhoton.hpp"
#include "GSource.hpp"

/* __ Model handling _____________________________________________________ */
#include "GModelPar.hpp"
#include "GModels.hpp"
#include "GModel.hpp"
#include "GModelRegistry.hpp"
#include "GModelSky.hpp"
#include "GModelData.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelRadialDisk.hpp"
#include "GModelRadialGauss.hpp"
#include "GModelRadialShell.hpp"
#include "GModelSpatialConst.hpp"
#include "GModelSpatialCube.hpp"
#include "GModelSpatialMap.hpp"
#include "GModelSpectral.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelSpectralConst.hpp"
#include "GModelSpectralExpPlaw.hpp"
#include "GModelSpectralFunc.hpp"
#include "GModelSpectralNodes.hpp"
#include "GModelSpectralPlaw.hpp"
#include "GModelSpectralPlaw2.hpp"
#include "GModelSpectralLogParabola.hpp"
#include "GModelTemporal.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GModelTemporalConst.hpp"

#endif /* GAMMALIB_HPP */
