/***************************************************************************
 *                         gammalib - SWIG file                            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
 * ----------------------------------------------------------------------- *
 * Usage:                                                                  *
 * swig -c++ -python -Wall gammalib.i                                      *
 ***************************************************************************/
%module gammalib

/* __ Support module _____________________________________________________ */
%include "GException.i"
%include "GNodeArray.i"
%include "GCsv.i"
%include "GRan.i"

/* __ Numerics module ____________________________________________________ */
%include "GVector.i"
%include "GMatrixBase.i"
//%include "GMatrixTools.i"
%include "GMatrix.i"
%include "GSymMatrix.i"
%include "GSparseMatrix.i"
//%include "GIntegral.i"
//%include "GIntegrand.i"

/* __ FITS module ________________________________________________________ */
%include "GFits.i"
%include "GFitsHDU.i"
%include "GFitsHeader.i"
%include "GFitsHeaderCard.i"
%include "GFitsImage.i"
%include "GFitsImageByte.i"
%include "GFitsImageSByte.i"
%include "GFitsImageUShort.i"
%include "GFitsImageShort.i"
%include "GFitsImageULong.i"
%include "GFitsImageLong.i"
%include "GFitsImageLongLong.i"
%include "GFitsImageFloat.i"
%include "GFitsImageDouble.i"
%include "GFitsTable.i"
%include "GFitsAsciiTable.i"
%include "GFitsBinTable.i"
%include "GFitsTableCol.i"
%include "GFitsTableBitCol.i"
%include "GFitsTableByteCol.i"
%include "GFitsTableBoolCol.i"
%include "GFitsTableStringCol.i"
%include "GFitsTableUShortCol.i"
%include "GFitsTableShortCol.i"
%include "GFitsTableULongCol.i"
%include "GFitsTableLongCol.i"
%include "GFitsTableLongLongCol.i"
%include "GFitsTableFloatCol.i"
%include "GFitsTableDoubleCol.i"
%include "GFitsTableCFloatCol.i"
%include "GFitsTableCDoubleCol.i"

/* __ XML module _________________________________________________________ */
%include "GXml.i"
%include "GXmlNode.i"
%include "GXmlDocument.i"
%include "GXmlText.i"
%include "GXmlElement.i"
%include "GXmlComment.i"
%include "GXmlAttribute.i"
%include "GXmlPI.i"

/* __ Application module _________________________________________________ */
%include "GApplication.i"
%include "GLog.i"
%include "GPars.i"
%include "GPar.i"

/* __ Optimizer module ___________________________________________________ */
%include "GOptimizer.i"
%include "GOptimizerLM.i"
%include "GOptimizerPars.i"
//%include "GOptimizerFunction.i"


/***************************************************************************
 *                        Analysis support services                        *
 ***************************************************************************/

/* __ Skymap handling ____________________________________________________ */
%include "GSkyDir.i"
%include "GSkyPixel.i"
%include "GSkymap.i"
%include "GWcs.i"
%include "GWcsRegistry.i"
%include "GWcslib.i"
%include "GWcsCAR.i"
%include "GWcsHPX.i"

/* __ Observation handling _______________________________________________ */
%include "GObservations.i"
%include "GObservation.i"
%include "GObservationRegistry.i"
%include "GEvents.i"
%include "GEventList.i"
%include "GEventCube.i"
%include "GEvent.i"
%include "GEventAtom.i"
%include "GEventBin.i"
%include "GInstDir.i"
%include "GEnergy.i"
%include "GTime.i"
%include "GRoi.i"
%include "GEbounds.i"
%include "GGti.i"
%include "GPointing.i"
%include "GResponse.i"
%include "GPhoton.i"

/* __ Model handling _____________________________________________________ */
%include "GModelPar.i"
%include "GModels.i"
%include "GModel.i"
%include "GModelRegistry.i"
%include "GModelSky.i"
%include "GModelData.i"
%include "GModelPointSource.i"
%include "GModelExtendedSource.i"
%include "GModelDiffuseSource.i"
%include "GModelSpatial.i"
%include "GModelSpatialRegistry.i"
%include "GModelSpatialConst.i"
%include "GModelSpatialCube.i"
%include "GModelSpatialPtsrc.i"
%include "GModelRadial.i"
%include "GModelRadialRegistry.i"
%include "GModelRadialDisk.i"
%include "GModelRadialGauss.i"
%include "GModelRadialShell.i"
%include "GModelSpectral.i"
%include "GModelSpectralRegistry.i"
%include "GModelSpectralConst.i"
%include "GModelSpectralExpPlaw.i"
%include "GModelSpectralFunc.i"
%include "GModelSpectralPlaw.i"
%include "GModelSpectralPlaw2.i"
%include "GModelTemporal.i"
%include "GModelTemporalRegistry.i"
%include "GModelTemporalConst.i"

/* __ Instrument specific ________________________________________________ */
#ifdef WITH_INST_MWL
%include "mwllib.i"
#endif
#ifdef WITH_INST_CTA
%include "ctalib.i"
#endif
#ifdef WITH_INST_LAT
%include "latlib.i"
#endif
