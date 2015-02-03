/***************************************************************************
 *      GCTACubeSourcePoint.cpp - CTA cube analysis point source class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Juergen Knoedlseder                         *
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
 * @file GCTACubeSourcePoint.cpp
 * @brief CTA cube analysis point source class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTACubeSourcePoint.hpp"
#include "GTools.hpp"
#include "GModelSpatial.hpp"
#include "GObservation.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GCTAEventCube.hpp"
#include "GCTAResponseCube.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_SET       "GCTACubeSourcePoint::set(GModelSpatial&, GObservation&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTACubeSourcePoint::GCTACubeSourcePoint(void) : GCTACubeSource()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] source Point source cube.
 ***************************************************************************/
GCTACubeSourcePoint::GCTACubeSourcePoint(const GCTACubeSourcePoint& source) :
                     GCTACubeSource(source)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(source);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTACubeSourcePoint::~GCTACubeSourcePoint(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] source Point source cube.
 * @return Point source cube.
 ***************************************************************************/
GCTACubeSourcePoint& GCTACubeSourcePoint::operator=(const GCTACubeSourcePoint& source)
{
    // Execute only if object is not identical
    if (this != &source) {

        // Copy base class members
        this->GCTACubeSource::operator=(source);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(source);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GCTACubeSourcePoint::clear(void)
{
    // Free class members
    free_members();
    this->GCTACubeSource::free_members();

    // Initialise members
    this->GCTACubeSource::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Deep copy of point source cube.
 ***************************************************************************/
GCTACubeSourcePoint* GCTACubeSourcePoint::clone(void) const
{
    // Return deep copy
    return new GCTACubeSourcePoint(*this);
}


/***********************************************************************//**
 * @brief Set point source cube for a given observation
 *
 * @param[in] name Model name.
 * @param[in] model Spatial model.
 * @param[in] obs Observation.
 ***************************************************************************/
void GCTACubeSourcePoint::set(const std::string&   name,
                              const GModelSpatial& model,
                              const GObservation&  obs)
{
    // Get pointer to model source model
    const GModelSpatialPointSource* ptsrc = dynamic_cast<const GModelSpatialPointSource*>(&model);
    if (ptsrc == NULL) {
        std::string msg = "Model is not a spatial point source model.";
        throw GException::invalid_value(G_SET, msg);
    }

    // Get pointer on CTA event cube
    const GCTAEventCube* cube = dynamic_cast<const GCTAEventCube*>(obs.events());
    if (cube == NULL) {
        std::string msg = "Observation does not contain a CTA event cube.";
        throw GException::invalid_value(G_SET, msg);
    }

    // Get pointer on CTA response cube
    const GCTAResponseCube* rsp = dynamic_cast<const GCTAResponseCube*>(obs.response());
    if (rsp == NULL) {
        std::string msg = "Observation does not contain a CTA response cube.";
        throw GException::invalid_value(G_SET, msg);
    }

    // Get point source attributes
    m_name        = name;
    m_dir         = ptsrc->dir();
    GTime srcTime = cube->time();

    // Set PSF deltas in radians
    m_deltas = rsp->psf().deltas();
    for (int i = 0; i < m_deltas.size(); ++i) {
        m_deltas[i] *= gammalib::deg2rad;
    }

    // Initialise data members
    m_aeff.assign(cube->ebins(), 0.0);
    m_delta_map.assign(cube->npix(), 0.0);
    m_psf.assign(cube->ebins()*m_deltas.size(), 0.0);

    // Compute deadtime corrected effective area
    for (int i = 0; i < cube->ebins(); ++i) {

        // Get source energy
        const GEnergy& srcEng = cube->energy(i);

        // Get exposure
        double aeff = rsp->exposure()(m_dir, srcEng);

        // Divide by ontime as the binned likelihood function is later
        // multiplying by ontime
        aeff /= obs.ontime();

        // Apply deadtime correction
        aeff *= obs.deadc(srcTime);

        // Store deadtime corrected effective area in vector
        m_aeff[i] = aeff;

    } // endfor: looped over energies

    // Compute distance map
    for (int i = 0; i < cube->npix(); ++i) {

        // Get cube pixel sky direction
        GSkyDir obsDir = cube->map().inx2dir(i);

        // Determine angular distance between point source direction and
        // cube pixel sky direction (radians)
        double delta = m_dir.dist(obsDir);

        // Store angular distance in vector
        m_delta_map[i] = delta;
        
    } // endfor: looped over spatial cube pixels

    // Compute point spread function
    for (int i = 0, inx = 0; i < cube->ebins(); ++i) {

        // Get source energy
        const GEnergy& srcEng = cube->energy(i);

        // Loop over delta bins
        for (int k = 0; k < m_deltas.size(); ++k, ++inx) {

            // Get PSF
            m_psf[inx] = rsp->psf()(m_dir, m_deltas[k], srcEng);

        } // endfor: looped over delta bins

    } // endfor: looped over PSF

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set point source cube for a given observation
 *
 * @param[in] ieng Energy layer.
 * @param[in] delta Distance between true and measured photon direction
 *                  (radians).
 ***************************************************************************/
double GCTACubeSourcePoint::psf(const int& ieng, const double& delta) const
{
    // Set node array interpolation values
    m_deltas.set_value(delta);

    // Compute offset
    int offset = ieng * m_deltas.size();

    // Compute PSF by bi-linear interpolation
    double psf = m_deltas.wgt_left()  * m_psf[offset+m_deltas.inx_left()] +
                 m_deltas.wgt_right() * m_psf[offset+m_deltas.inx_right()];

    // Reset negative PSF values
    if (psf < 0.0) {
        psf = 0.0;
    }

    // Return
    return psf;
}


/***********************************************************************//**
 * @brief Print point source cube information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing point source cube information.
 ***************************************************************************/
std::string GCTACubeSourcePoint::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTACubeSourcePoint ===");
        result.append("\n"+gammalib::parformat("Source name") + name());

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTACubeSourcePoint::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_aeff.clear();
    m_delta_map.clear();
    m_psf.clear();
    m_deltas.clear();
   
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] source Point source cube.
 ***************************************************************************/
void GCTACubeSourcePoint::copy_members(const GCTACubeSourcePoint& source)
{
    // Copy members
    m_dir       = source.m_dir;
    m_aeff      = source.m_aeff;
    m_delta_map = source.m_delta_map;
    m_psf       = source.m_psf;
    m_deltas    = source.m_deltas;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTACubeSourcePoint::free_members(void)
{
    // Return
    return;
}
