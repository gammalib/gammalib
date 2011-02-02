/***************************************************************************
 *                 GCTASupport.hpp  -  CTA support functions               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTASupport.hpp
 * @brief Definition of support function used by CTA classes
 * @author J. Knodlseder
 */

#ifndef GCTASUPPORT_HPP
#define GCTASUPPORT_HPP

/* __ Includes ___________________________________________________________ */

/* __ Namespaces _________________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Prototypes _________________________________________________________ */
double cta_roi_arclength(const double& rad,     const double& dist,
                         const double& cosdist, const double& sindist,
                         const double& roi,     const double& cosroi);

#endif /* GCTASUPPORT_HPP */
