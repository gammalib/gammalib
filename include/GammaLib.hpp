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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


/***************************************************************************
 *                              Core services                              *
 ***************************************************************************/

/* __ Common tools _______________________________________________________ */
#include "GException.hpp"
#include "GNodeArray.hpp"

/* __ Numerics module ____________________________________________________ */
#include "GMatrixTools.hpp"
#include "GVector.hpp"
#include "GMatrixBase.hpp"
#include "GMatrix.hpp"
#include "GSymMatrix.hpp"
#include "GSparseMatrix.hpp"

/* __ FITS module ________________________________________________________ */
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
#include "GFitsTableLogCol.hpp"
#include "GFitsTableStrCol.hpp"
#include "GFitsTableShtCol.hpp"
#include "GFitsTableLngCol.hpp"
#include "GFitsTableFltCol.hpp"
#include "GFitsTableDblCol.hpp"

/* __ Parameter interface module _________________________________________ */


/* __ Reporting module ___________________________________________________ */


/* __ Application module _________________________________________________ */


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
#include "GData.hpp"
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

/* __ Fermi/LAT specific code ____________________________________________ */
#include "GLATObservation.hpp"
#include "GLATEventList.hpp"
#include "GLATEventCube.hpp"
#include "GLATEventAtom.hpp"
#include "GLATEventBin.hpp"
#include "GLATResponse.hpp"

#endif /* GAMMALIB_HPP */
