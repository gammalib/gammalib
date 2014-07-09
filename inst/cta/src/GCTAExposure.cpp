/***************************************************************************
 *            GCTAExposure.cpp - CTA mean point spread function class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Chia-Chun Lu                                     *
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
 * @file GCTAExposure.cpp
 * @brief CTA mean point spread function class implementation
 * @author Chia-Chun Lu
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTAExposure.hpp"
#include "GCTAObservation.hpp"
#include "GMath.hpp"
using namespace gammalib;
/* __ Method name definitions ____________________________________________ */

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
GCTAExposure::GCTAExposure(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] exp Point spread function.
 ***************************************************************************/
GCTAExposure::GCTAExposure(const GCTAExposure& expcube)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(expcube);

    // Return
    return;
}

/***********************************************************************//**
 * @brief constructor taking definitions of 3D bins.
 *
 * @param[in] obs observation container.
 * @param[in] x   sky coordinate.
 * @param[in] y   sky coordinate.
 * @param[in] dx  pixel size.
 * @param[in] dy  pixel size.
 * @param[in] nx  number of x pixels.
 * @param[in] ny  number of y pixels.
 * @param[in] emin minimum energy.
 * @param[in] emax maximum energy.
 * @param[in] nebins number of energy bins.
 *
 * This constructor resample an effective area talbe and saves it as a 3D cube.
 * The index of the map is 
 * The bin size is decided evenly in delta^2 space.
 ***************************************************************************/
GCTAExposure::GCTAExposure(const GObservations& obs, 
			 const double& x, const double& y, 
			 const double& dx, const double& dy,
			 const int& nx, const int& ny,
			 const double& emin, const double& emax, const int& nebins)
{
    // Initialise class members
    init_members();
    m_obs = obs;
    m_nebins = nebins;
    m_cube = GSkymap("CAR", "CEL", x, y, dx, dy,
		     nx, ny, nebins);
    // Initialise energy bounds
    GEnergy gemin = GEnergy(emin, "TeV");
    GEnergy gemax = GEnergy(emax, "TeV");
    m_ebounds = GEbounds(nebins, gemin, gemax, true);

    //
    set_expcube();
    // Return
    return;
}

/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAExposure::~GCTAExposure(void)
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
 * @param[in] exp exposure cube
 * @return Exposure cube.
 ***************************************************************************/
GCTAExposure& GCTAExposure::operator= (const GCTAExposure& expcube)
{
    // Execute only if object is not identical
    if (this != &expcube) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(expcube);

    } // endif: object was not identical

    // Return this object
    return *this;
}

/***********************************************************************//**
 * @brief Set the Exposure cube
 * This function takes the effective area table given in the observation xml
 * to produce a 4D Exposure cube.
 ***************************************************************************/
void GCTAExposure::set_expcube(void)
{
  for (int i = 0 ; i< m_obs.size(); i++){
    GCTAObservation *obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
    GCTAResponse rsp = obs->response();
    GSkyDir pnt = obs->pointing().dir();
    for (int ie = 0 ; ie < m_nebins ; ie++){
      for (int pix = 0 ; pix < m_cube.npix() ; pix++){
	  GSkyDir dir = m_cube.inx2dir(pix);
	  double logeng = m_ebounds.emean(ie).log10TeV();
	  double theta = pnt.dist(dir); // radian
	  m_cube(pix,ie) = rsp.aeff(theta, 
				    0.0, 0.0, 0.0, 
				    logeng);
	}
      }
    save("expcube.fits",true);
  }
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
void GCTAExposure::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    //this->GCTAExposure::free_members();

    // Initialise members
    // this->GCTAExposure::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA event cube into FITS file.
 *
 * @param[in] fits FITS file.
 ***************************************************************************/
void GCTAExposure::write(GFits& fits) const
{
    // Write cube
    m_cube.write(fits);

    // Write energy boundaries
    m_ebounds.write(fits);

    // Return
    return;
}




/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Deep copy of point spread function instance.
 ***************************************************************************/
GCTAExposure* GCTAExposure::clone(void) const
{
    return new GCTAExposure(*this);
}

/***********************************************************************//**
 * @brief Load point spread function from performance table
 *
 * @param[in] filename Performance table file name.
 *
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method loads the point spread function information from an effective
 * area response table.
 ***************************************************************************/
void GCTAExposure::load(const std::string& filename)
{
    return;
}

void GCTAExposure::save(const std::string& filename,
		       const bool& clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write event cube
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);
    return;
}


/***********************************************************************//**
 * @brief Print point spread function information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing point spread function information.
 ***************************************************************************/
std::string GCTAExposure::print(const GChatter& chatter) const
{
    std::string result;
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
void GCTAExposure::init_members(void)
{
    // Initialise members
    m_obs.clear();
    m_cube.clear();
    m_ebounds.clear();
    m_nebins = 0;
   
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] expcube Exposure cube
 ***************************************************************************/
void GCTAExposure::copy_members(const GCTAExposure& cube)
{
    // Initialise members
    m_obs = cube.m_obs;
    m_cube = cube.m_cube;
    m_ebounds = cube.m_ebounds;
    m_nebins = m_nebins;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAExposure::free_members(void)
{
    // Return
    return;
}
