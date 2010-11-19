/***************************************************************************
 *                         gammalib - SWIG file                            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 * Usage:                                                                  *
 * swig -c++ -python -Wall -includeall gammalib.i                          *
 ***************************************************************************/
%module gammalib

/* __ Common tools _______________________________________________________ */
#include "GException.i"
#include "GNodeArray.i"

/* __ Numerics module ____________________________________________________ */
#include "GVector.i"
#include "GMatrixBase.i"
//#include "GMatrixTools.hpp"
#include "GMatrix.i"
#include "GSymMatrix.i"
#include "GSparseMatrix.i"
//#include "GIntegral.hpp"
//#include "GIntegrand.hpp"

/* __ FITS module ________________________________________________________ */
#include "GFits.i"
#include "GFitsHDU.i"
#include "GFitsHeader.i"
#include "GFitsHeaderCard.i"
#include "GFitsImage.i"
//#include "GFitsImageByte.i"
#include "GFitsTable.i"
#include "GFitsAsciiTable.i"
#include "GFitsBinTable.i"
#include "GFitsTableCol.i"
#include "GFitsTableBitCol.i"
#include "GFitsTableBoolCol.i"
#include "GFitsTableDoubleCol.i"
#include "GFitsTableFloatCol.i"
#include "GFitsTableULongCol.i"
#include "GFitsTableLongCol.i"
#include "GFitsTableUShortCol.i"
#include "GFitsTableLongLongCol.i"
#include "GFitsTableShortCol.i"
#include "GFitsTableStringCol.i"

/* __ Application module _________________________________________________ */
#include "GApplication.i"
//#include "GLog.i"
#include "GPars.i"
#include "GPar.i"

/* __ Optimizer module ___________________________________________________ */
#include "GOptimizer.i"
#include "GOptimizerLM.i"
#include "GOptimizerPars.i"
//#include "GOptimizerFunction.i"


/***************************************************************************
 *                        Analysis support services                        *
 ***************************************************************************/

/* __ Skymap handling ____________________________________________________ */
#include "GSkyDir.i"
#include "GSkyPixel.i"
#include "GSkymap.i"
#include "GWcs.i"
#include "GWcsCAR.i"
#include "GWcsHPX.i"

/* __ Observation handling _______________________________________________ */
//#include "GEnergy.i"
//#include "GTime.i"
#include "GObservations.i"
#include "GObservation.i"
//#include "GEvents.i"
//#include "GEventList.i"
//#include "GEventCube.i"
//#include "GEvent.i"
//#include "GGti.i"
//#include "GRoi.i"
//#include "GEbounds.i"
//#include "GInstDir.i"
//#include "GPointing.i"
//#include "GResponse.i"

/* __ Model handling _____________________________________________________ */
#include "GModels.i"
#include "GModel.i"
#include "GModelPar.i"
#include "GModelSpatial.i"
#include "GModelSpatialPtsrc.i"
#include "GModelSpectral.i"
#include "GModelSpectralPlaw.i"
#include "GModelTemporal.i"
#include "GModelTemporalConst.i"

/* __ CTA ________________________________________________________________ */
#include "GCTAObservation.i"



