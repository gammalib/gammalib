/***************************************************************************
 *            GCTAMeanPsf.cpp - CTA mean point spread function class       *
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
 * @file GCTAMeanPsf.cpp
 * @brief CTA mean point spread function class implementation
 * @author Chia-Chun Lu
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTAMeanPsf.hpp"
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
GCTAMeanPsf::GCTAMeanPsf(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
GCTAMeanPsf::GCTAMeanPsf(const GCTAMeanPsf& psfcube)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(psfcube);

    // Return
    return;
}

/***********************************************************************//**
 * @brief constructor taking definitions of 4D bins.
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
 * @param[in] min minimum delta.
 * @param[in] max maximum delta.
 * @param[in] nbins number of delta bins.
 *
 * This constructor resample a Psf function and saves it as a 4D skymap.
 * The index of the map is 
 * The bin size is decided evenly in delta^2 space.
 ***************************************************************************/
GCTAMeanPsf::GCTAMeanPsf(const GObservations& obs, 
			 const double& x, const double& y, 
			 const double& dx, const double& dy,
			 const int& nx, const int& ny,
			 const double& emin, const double& emax, const int& nebins,
			 const double& min, const double& max, const int& nbins)
{
    // Initialise class members
    init_members();
    m_obs = obs;
    m_nbins = nbins;
    m_nebins = nebins;
    m_cube = GSkymap("CAR", "CEL", x, y, dx, dy,
		     nx, ny, nebins*nbins);
    // Initialise energy bounds
    GEnergy gemin = GEnergy(emin, "TeV");
    GEnergy gemax = GEnergy(emax, "TeV");
    m_ebounds = GEbounds(nebins, gemin, gemax, true);

    // Initialise delta node array
    std::vector<double> deltas;
    for (int i = 0 ; i < nbins; i++){
      double binsize = (max*max - min*min)/nbins;
      double delta = std::sqrt(binsize*0.5*i + min*min);
      deltas.push_back(delta);
    }
    m_deltas = GNodeArray(deltas);

    //
    set_psfcube();
    // Return
    return;
}

/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAMeanPsf::~GCTAMeanPsf(void)
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
 * @param[in] psf Point spread function.
 * @return Point spread function.
 ***************************************************************************/
GCTAMeanPsf& GCTAMeanPsf::operator= (const GCTAMeanPsf& psf)
{
    // Execute only if object is not identical
    if (this != &psf) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(psf);

    } // endif: object was not identical

    // Return this object
    return *this;
}

/***********************************************************************//**
 * @brief Set the psf skymap
 * This function takes the psf response function given in the observation xml
 * to produce a 4D psf cube.
 ***************************************************************************/
void GCTAMeanPsf::set_psfcube(void)
{
  for (int i = 0 ; i< m_obs.size(); i++){
    GCTAObservation *obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
    GCTAResponse rsp = obs->response();
    GSkyDir pnt = obs->pointing().dir();
    for (int ie = 0 ; ie < m_nebins ; ie++){
      for (int pix = 0 ; pix < m_cube.npix() ; pix++){
	for (int idelta = 0 ; idelta < m_nbins ; idelta++){
	  GSkyDir dir = m_cube.inx2dir(pix);
	  double logeng = m_ebounds.emean(ie).log10TeV();
	  double theta = pnt.dist(dir); // radian
	  double delta = m_deltas[idelta]*deg2rad; // radian
	  int inx = idelta + ie * m_nbins;
	  m_cube(pix,inx) = rsp.psf(delta, theta, 
				    0.0, 0.0, 0.0, 
				    logeng);
	}
      }
    }
    save("test.fits",true);
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
void GCTAMeanPsf::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    //this->GCTAMeanPsf::free_members();

    // Initialise members
    // this->GCTAMeanPsf::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA event cube into FITS file.
 *
 * @param[in] fits FITS file.
 ***************************************************************************/
void GCTAMeanPsf::write(GFits& fits) const
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
GCTAMeanPsf* GCTAMeanPsf::clone(void) const
{
    return new GCTAMeanPsf(*this);
}

/***********************************************************************//**
 * @brief Load point spread function from performance table
 *
 * @param[in] filename Performance table file name.
 *
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method loads the point spread function information from a PSF
 * response table.
 ***************************************************************************/
void GCTAMeanPsf::load(const std::string& filename)
{
    return;
}

void GCTAMeanPsf::save(const std::string& filename,
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
std::string GCTAMeanPsf::print(const GChatter& chatter) const
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
void GCTAMeanPsf::init_members(void)
{
    // Initialise members
    m_obs.clear();
    m_cube.clear();
    m_ebounds.clear();
    m_deltas.clear();
    m_nbins = 0;
    m_nebins = 0;
   
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] psfcube psf cube
 ***************************************************************************/
void GCTAMeanPsf::copy_members(const GCTAMeanPsf& cube)
{
    // Initialise members
    m_obs = cube.m_obs;
    m_cube = cube.m_cube;
    m_ebounds = cube.m_ebounds;
    m_deltas = cube.m_deltas;
    m_nbins = m_nbins;
    m_nebins = m_nebins;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAMeanPsf::free_members(void)
{
    // Return
    return;
}
