/***************************************************************************
 *         GammaLib.hpp  -  Gamma-Ray Astronomy Library Header file        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2011 by Jurgen Knodlseder                           *
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

#ifndef GAMMALIB_HPP
#define GAMMALIB_HPP

/* __ Includes ___________________________________________________________ */


/***************************************************************************
 *                              Core services                              *
 ***************************************************************************/

/* __ Common tools _______________________________________________________ */
#include "GException.hpp"
#include "GNodeArray.hpp"
#include "GCsv.hpp"
#include "GRan.hpp"

/* __ Numerics module ____________________________________________________ */
#include "GVector.hpp"
#include "GMatrixBase.hpp"
#include "GMatrixTools.hpp"
#include "GMatrix.hpp"
#include "GSymMatrix.hpp"
#include "GSparseMatrix.hpp"
#include "GIntegral.hpp"
#include "GIntegrand.hpp"
#include "GDerivative.hpp"
#include "GFunction.hpp"

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


/***************************************************************************
 *                        Analysis support services                        *
 ***************************************************************************/

/* __ Skymap handling ____________________________________________________ */
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"
#include "GSkymap.hpp"
#include "GWcs.hpp"
#include "GWcsCAR.hpp"
#include "GWcsHPX.hpp"

/* __ Observation handling _______________________________________________ */
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GObservations.hpp"
#include "GObservation.hpp"
#include "GEvents.hpp"
#include "GEventList.hpp"
#include "GEventCube.hpp"
#include "GEvent.hpp"
#include "GGti.hpp"
#include "GRoi.hpp"
#include "GEbounds.hpp"
#include "GInstDir.hpp"
#include "GPointing.hpp"
#include "GResponse.hpp"
#include "GPhoton.hpp"

/* __ Model handling _____________________________________________________ */
#include "GModelPar.hpp"
#include "GModels.hpp"
#include "GModel.hpp"
#include "GModelRegistry.hpp"
#include "GModelSky.hpp"
#include "GModelData.hpp"
#include "GModelPointSource.hpp"
#include "GModelExtendedSource.hpp"
#include "GModelDiffuseSource.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GModelSpatialConst.hpp"
#include "GModelSpatialCube.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GModelRadial.hpp"
#include "GModelRadialRegistry.hpp"
#include "GModelRadialDisk.hpp"
#include "GModelRadialGauss.hpp"
#include "GModelRadialShell.hpp"
#include "GModelSpectral.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelSpectralConst.hpp"
#include "GModelSpectralExpPlaw.hpp"
#include "GModelSpectralFunc.hpp"
#include "GModelSpectralPlaw.hpp"
#include "GModelSpectralPlaw2.hpp"
#include "GModelTemporal.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GModelTemporalConst.hpp"

#endif /* GAMMALIB_HPP */
