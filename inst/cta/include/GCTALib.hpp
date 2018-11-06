/***************************************************************************
 *                 GCTALib.hpp - CTA Support Header files                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @file GCTALib.hpp
 * @brief Collection of CTA support header files
 * @author Juergen Knoedlseder
 */

#ifndef GCTALIB_HPP
#define GCTALIB_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"

/* __ CTA typemaps _______________________________________________________ */
#include "GCTATypemaps.hpp"

/* __ CTA specific headers _______________________________________________ */
#include "GCTAException.hpp"
#include "GCTAObservation.hpp"
#include "GCTAOnOffObservation.hpp"
#include "GCTAEventList.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTAEventCube.hpp"
#include "GCTAEventBin.hpp"
#include "GCTAInstDir.hpp"
#include "GCTARoi.hpp"
#include "GCTAPointing.hpp"
#include "GCTAResponse.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTAResponseCube.hpp"
#include "GCTAResponseTable.hpp"
#include "GCTAResponseCache.hpp"
#include "GCTAAeff.hpp"
#include "GCTAAeffPerfTable.hpp"
#include "GCTAAeffArf.hpp"
#include "GCTAAeff2D.hpp"
#include "GCTAPsf.hpp"
#include "GCTAPsfPerfTable.hpp"
#include "GCTAPsfVector.hpp"
#include "GCTAPsf2D.hpp"
#include "GCTAPsfKing.hpp"
#include "GCTAPsfTable.hpp"
#include "GCTAEdisp.hpp"
#include "GCTAEdispPerfTable.hpp"
#include "GCTAEdispRmf.hpp"
#include "GCTAEdisp2D.hpp"
#include "GCTABackground.hpp"
#include "GCTABackgroundPerfTable.hpp"
#include "GCTABackground3D.hpp"
#include "GCTACubeExposure.hpp"
#include "GCTACubeBackground.hpp"
#include "GCTACubePsf.hpp"
#include "GCTACubeEdisp.hpp"
#include "GCTACubeSource.hpp"
#include "GCTACubeSourcePoint.hpp"
#include "GCTACubeSourceDiffuse.hpp"
#include "GCTAModelBackground.hpp"
#include "GCTAModelSpatial.hpp"
#include "GCTAModelSpatialGradient.hpp"
#include "GCTAModelSpatialMultiplicative.hpp"
#include "GCTAModelSpatialRegistry.hpp"
#include "GCTAModelCubeBackground.hpp"
#include "GCTAModelIrfBackground.hpp"
#include "GCTAModelAeffBackground.hpp"
#include "GCTAModelRadial.hpp"
#include "GCTAModelRadialRegistry.hpp"
#include "GCTAModelRadialGauss.hpp"
#include "GCTAModelRadialPolynom.hpp"
#include "GCTAModelRadialProfile.hpp"
#include "GCTAModelRadialAcceptance.hpp"

/* __ CTA specific definitions ___________________________________________ */
#define G_CTA_MJDREF 51544.5                 //!< Reference of CTA time frame

#endif /* GCTALIB_HPP */
