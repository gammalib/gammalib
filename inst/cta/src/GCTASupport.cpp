/***************************************************************************
 *                 GCTASupport.cpp  -  CTA support functions               *
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
 * @file GCTASupport.cpp
 * @brief Implementation of support function used by CTA classes
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include <iostream>
#include "GCTASupport.hpp"
#include "GTools.hpp"

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_CHECK_FOR_NAN 0


/***********************************************************************//**
 * @brief Returns length of circular arc within circular ROI
 *
 * @param[in] rad Circle radius in radians (<pi).
 * @param[in] dist Circle centre distance to ROI centre (<pi).
 * @param[in] cosdist Cosine of circle centre distance to ROI centre.
 * @param[in] sindist Sinus of circle centre distance to ROI centre.
 * @param[in] roi Radius of ROI in radians.
 * @param[in] cosroi Cosine of ROI radius.
 *
 * This method returns the arclength in radians of a circle of radius 'rad'
 * with a centre that is offset by 'dist' from the ROI centre, where the ROI
 * radius is given by 'roi'. To speed-up computations, the cosines and sinus
 * of 'roi' and 'psf' should be calculated by the client and be passed to
 * the method.
 ***************************************************************************/
double cta_roi_arclength(const double& rad,     const double& dist,
                         const double& cosdist, const double& sindist,
                         const double& roi,     const double& cosroi)
{
    // Declare arclength
    double arclength;

    // Handle special case that circle centre matches ROI centre
    if (dist == 0.0) {
        if (rad > roi) arclength = 0.0;   // Circle outside ROI
        else           arclength = twopi; // Circle inside ROI
    }

    // ... otherwise circle and ROI centres are not identical
    else {

        // Handle special case that we evaluate exactly at the circle
        // centre. In this case we have in fact a point, and if this point
        // falls within the ROI it has a formal arclength of 2pi.
        //
        if (rad == 0.0) {
            if (dist > roi) arclength = 0.0;   // Circle centre outside ROI
            else            arclength = twopi; // Circle centre inside ROI
        }

        // ... otherwise we have to handle the general case
        else {
            double d = roi - dist;
            if (-rad >= d)
                arclength = 0.0;
            else if (rad <= d)
                arclength = twopi;
            else {
                double cosrad = std::cos(rad);
                double sinrad = std::sin(rad);
                double cosang = (cosroi - cosdist*cosrad) / (sindist*sinrad);
                arclength     = 2.0 * arccos(cosang);
                #if G_CHECK_FOR_NAN
                if (std::isnan(arclength)) {
                    std::cout << "cta_roi_arclength: NaN occured";
                    std::cout << " rad=" << rad;
                    std::cout << " sinrad=" << sinrad;
                    std::cout << " cosrad=" << cosrad;
                    std::cout << " sindist=" << sindist;
                    std::cout << " cosdist=" << cosdist;
                    std::cout << " cosang=" << cosang;
                    std::cout << std::endl;
                }
                #endif
            }
        }

    } // endelse: Circle and ROI centres were not identical

    // Return arclength
    return arclength;
}
