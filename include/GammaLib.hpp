/***************************************************************************
 *         GammaLib.hpp  -  Gamma-Ray Astronomy Library Header file        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
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

/* __ Numerics module ____________________________________________________ */
#include "GVector.hpp"
#include "GMatrixBase.hpp"
#include "GMatrixTools.hpp"
#include "GMatrix.hpp"
#include "GSymMatrix.hpp"
#include "GSparseMatrix.hpp"

/* __ FITS module ________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFits.hpp"
#include "GFitsHDU.hpp"
#include "GFitsHeader.hpp"
#include "GFitsHeaderCard.hpp"
#include "GFitsData.hpp"
#include "GFitsImage.hpp"
#include "GFitsImageFlt.hpp"
#include "GFitsImageDbl.hpp"
#include "GFitsAsciiTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableCol.hpp"
#include "GFitsTableBitCol.hpp"
#include "GFitsTableLogCol.hpp"
#include "GFitsTableStrCol.hpp"
#include "GFitsTableShtCol.hpp"
#include "GFitsTableLngCol.hpp"
#include "GFitsTableLlgCol.hpp"
#include "GFitsTableFltCol.hpp"
#include "GFitsTableDblCol.hpp"

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
#include "GEbounds.hpp"
#include "GResponse.hpp"

/* __ Model handling _____________________________________________________ */
#include "GModels.hpp"
#include "GModel.hpp"
#include "GModelPar.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GModelSpectral.hpp"
#include "GModelSpectralPlaw.hpp"
#include "GModelTemporal.hpp"
#include "GModelTemporalConst.hpp"


/***************************************************************************
 *                      Instrument specific services                       *
 ***************************************************************************/

/* __ CTA specific code __________________________________________________ */
#include "GCTAObservation.hpp"
#include "GCTAEventList.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTAResponse.hpp"

/* __ Fermi/LAT specific code ____________________________________________ */
#include "GLATObservation.hpp"
#include "GLATEventList.hpp"
#include "GLATEventCube.hpp"
#include "GLATEventAtom.hpp"
#include "GLATEventBin.hpp"
#include "GLATResponse.hpp"

#endif /* GAMMALIB_HPP */
