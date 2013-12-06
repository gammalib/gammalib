/***************************************************************************
 *             GLATLtCubeMap.cpp - Fermi/LAT livetime cube map             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GLATLtCubeMap.cpp
 * @brief Fermi/LAT livetime cube map class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATLtCubeMap.hpp"
#include "GTools.hpp"
#include "GMath.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_COSTHETA                            "GLATLtCubeMap::costheta(int&)"
#define G_PHI                                      "GLATLtCubeMap::phi(int&)"

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
GLATLtCubeMap::GLATLtCubeMap(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param map Livetime cube map.
 ***************************************************************************/
GLATLtCubeMap::GLATLtCubeMap(const GLATLtCubeMap& map)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(map);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATLtCubeMap::~GLATLtCubeMap(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Operators                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] map Livetime cube map.
 ***************************************************************************/
GLATLtCubeMap& GLATLtCubeMap::operator= (const GLATLtCubeMap& map)
{
    // Execute only if object is not identical
    if (this != &map) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(map);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Sum function multiplied by livetime over zenith angle
 *
 * @param[in] dir Sky direction.
 * @param[in] fct Function to evaluate.
 *
 * Computes
 * \f[\sum_{\cos \theta} T_{\rm live}(\cos \theta) f(\cos \theta)\f]
 * where
 * \f$T_{\rm live}(\cos \theta)\f$ is the livetime as a function of the
 * cosine of the zenith angle, and
 * \f$f(\cos \theta)\f$ is a function that depends on the cosine of the
 * zenith angle.
 * This method assumes that \f$T_{\rm live}(\cos \theta)\f$ is stored as a
 * set of maps.
 ***************************************************************************/
double GLATLtCubeMap::operator() (const GSkyDir& dir, _ltcube_ctheta fct)
{
    // Get map index
    int pixel = m_map.dir2pix(dir);

    // Initialise sum
    double sum = 0.0;

    // Loop over zenith angles
    for (int i = 0; i < m_num_ctheta; ++i)
        sum += m_map(pixel, i) * (*fct)(costheta(i));

    // Return sum
    return sum;
}


/***********************************************************************//**
 * @brief Sum function multiplied by livetime over zenith and azimuth angles
 *
 * @param[in] dir Sky direction.
 * @param[in] fct Function to evaluate.
 *
 * Computes
 * \f[\sum_{\cos \theta, \phi} T_{\rm live}(\cos \theta, \phi) 
 *    f(\cos \theta, \phi)\f]
 * where
 * \f$T_{\rm live}(\cos \theta, \phi)\f$ is the livetime as a function of
 * the cosine of the zenith and of the azimuth angle, and
 * \f$f(\cos \theta, \phi)\f$ is a function that depends on the cosine of
 * the zenith angle and of the azimuth angle.
 * This method assumes that \f$T_{\rm live}(\cos \theta, \phi)\f$ is
 * stored as a set of maps in a 2D array with \f$\cos \theta\f$ being the
 * most rapidely varying parameter and with the first map starting at
 * index m_num_ctheta (the first m_num_ctheta maps are the livetime cube
 * maps without any \f$\phi\f$ dependence).
 ***************************************************************************/
double GLATLtCubeMap::operator() (const GSkyDir& dir, _ltcube_ctheta_phi fct)
{
    // Get map index
    int pixel = m_map.dir2pix(dir);

    // Initialise sum
    double sum = 0.0;

    // Loop over azimuth and zenith angles. Note that the map index starts
    // with m_num_ctheta as the first m_num_ctheta maps correspond to an
    // evaluation without any phi-dependence.
    for (int iphi = 0, i = m_num_ctheta; iphi < m_num_phi; ++iphi) {
        double p = phi(iphi);
        for (int itheta = 0; itheta < m_num_ctheta; ++itheta, ++i) {
            sum += m_map(pixel, i) * (*fct)(costheta(itheta), p);
        }
    }

    // Return sum
    return sum;
}


/***********************************************************************//**
 * @brief Sum effective area multiplied by livetime over zenith and
 *        (optionally) azimuth angles
 *
 * @param[in] dir True sky direction.
 * @param[in] energy True photon energy.
 * @param[in] aeff Effective area.
 *
 * Computes
 * \f[\sum_{\cos \theta, \phi} T_{\rm live}(\cos \theta, \phi) 
 *    A_{\rm eff}(\log E, \cos \theta, \phi)\f]
 * where
 * \f$T_{\rm live}(\cos \theta, \phi)\f$ is the livetime as a function of
 * the cosine of the zenith and the azimuth angle, and
 * \f$A_{\rm eff}(\log E, \cos \theta, \phi)\f$ is the effective area that
 * depends on
 * the log10 of the energy (in MeV),
 * the cosine of the zenith angle, and
 * the azimuth angle.
 * This method assumes that \f$T_{\rm live}(\cos \theta, \phi)\f$ is
 * stored as a set of maps in a 2D array with \f$\cos \theta\f$ being the
 * most rapidely varying parameter and with the first map starting at
 * index m_num_ctheta (the first m_num_ctheta maps are the livetime cube
 * maps without any \f$\phi\f$ dependence).
 ***************************************************************************/
double GLATLtCubeMap::operator() (const GSkyDir& dir, const GEnergy& energy,
                                  const GLATAeff& aeff)
{
    // Get map index
    int pixel = m_map.dir2pix(dir);

    // Initialise sum
    double sum = 0.0;

    // Circumvent const correctness
    GLATAeff* fct = ((GLATAeff*)&aeff);

    // If livetime cube and response have phi dependence then sum over
    // zenith and azimuth. Note that the map index starts with m_num_ctheta
    // as the first m_num_ctheta maps correspond to an evaluation without
    // any phi-dependence.
    if (hasphi() && aeff.hasphi()) {
        for (int iphi = 0, i = m_num_ctheta; iphi < m_num_phi; ++iphi) {
            double p = phi(iphi);
            for (int itheta = 0; itheta < m_num_ctheta; ++itheta, ++i)
                sum += m_map(pixel, i) * (*fct)(energy.log10MeV(), costheta(i), p);
        }
    }

    // ... otherwise sum only over zenith angle
    else {
        for (int i = 0; i < m_num_ctheta; ++i)
            sum += m_map(pixel, i) * (*fct)(energy.log10MeV(), costheta(i));
    }

    // Return sum
    return sum;
}


/***********************************************************************//**
 * @brief Sum effective area multiplied by livetime over zenith angles
 *
 * @param[in] dir True sky direction.
 * @param[in] energy True photon energy.
 * @param[in] offset Offset from true direction (deg).
 * @param[in] psf Point spread function.
 * @param[in] aeff Effective area.
 *
 * Computes
 * \f[\sum_{\cos \theta} T_{\rm live}(\cos \theta) 
 *    PSF(\log E, \delta, \cos \theta) A_{\rm eff}(\cos \theta, \phi)\f]
 * where
 * \f$T_{\rm live}(\cos \theta)\f$ is the livetime as a function of the
 * cosine of the zenith angle, and
 * \f$PSF(\log E, \delta, \cos \theta)\f$ is the point spread function that
 * depends on
 * the log10 of the energy (in MeV),
 * the offset angle from the true direction (in degrees), and
 * the cosine of the zenith angle, and
 * \f$A_{\rm eff}(\cos \theta, \phi)\f$ is the effective area that depends
 * on the cosine of the zenith angle and (optionally) of the azimuth angle.
 ***************************************************************************/
double GLATLtCubeMap::operator() (const GSkyDir& dir, const GEnergy& energy,
                                  const double& offset, const GLATPsf& psf,
                                  const GLATAeff& aeff)
{
    // Get map index
    int pixel = m_map.dir2pix(dir);

    // Initialise sum
    double sum = 0.0;

    // Circumvent const correctness
    GLATPsf*  fpsf  = ((GLATPsf*)&psf);
    GLATAeff* faeff = ((GLATAeff*)&aeff);

    // Get log10 of energy
    double logE = energy.log10MeV();

    // If livetime cube and response have phi dependence then sum over
    // zenith and azimuth. Note that the map index starts with m_num_ctheta
    // as the first m_num_ctheta maps correspond to an evaluation without
    // any phi-dependence.
    if (hasphi() && aeff.hasphi()) {
        for (int iphi = 0, i = m_num_ctheta; iphi < m_num_phi; ++iphi) {
            double p = phi(iphi);
            for (int itheta = 0; itheta < m_num_ctheta; ++itheta, ++i)
                sum += m_map(pixel, i) * (*faeff)(logE, costheta(i), p) *
                       (*fpsf)(offset, logE, costheta(i));
        }
    }

    // ... otherwise sum only over zenith angle
    else {
        for (int i = 0; i < m_num_ctheta; ++i)
            sum += m_map(pixel, i) * (*faeff)(logE, costheta(i)) *
                   (*fpsf)(offset, logE, costheta(i));
    }

    // Return sum
    return sum;
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
void GLATLtCubeMap::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GLATLtCubeMap* GLATLtCubeMap::clone(void) const
{
    return new GLATLtCubeMap(*this);
}


/***********************************************************************//**
 * @brief Load livetime cube from FITS file
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
void GLATLtCubeMap::read(const GFitsTable* hdu)
{
    // Clear object
    clear();

    // Load skymap
    m_map.read(*hdu);

    // Set costheta binning scheme
    std::string scheme = gammalib::strip_whitespace(gammalib::toupper(hdu->string("THETABIN")));
    m_sqrt_bin = (scheme == "SQRT(1-COSTHETA)");

    // Read attributes
    m_num_ctheta = hdu->integer("NBRBINS");
    m_num_phi    = hdu->integer("PHIBINS");
    m_min_ctheta = hdu->real("COSMIN");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save livetime cube map into FITS file
 *
 * @param[in] file FITS file.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
void GLATLtCubeMap::write(GFits* file) const
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return cos theta value for an index
 *
 * @param[in] index Bin index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Bin index outside value range.
 *
 * The cos theta value is computed using
 * \f[\cos \theta = 1 - \left( \frac{{\rm index} + 0.5}{\rm ncostheta}
 *    \right)^N (1 - {\rm costhetamin})\f]
 * where
 * \f${\rm index}\f$ is the bin index,
 * \f${\rm ncostheta}\f$ is the number of cos theta bins,
 * \f$N\f$ is either 1 or 2, and
 * \f${\rm costhetamin}\f$ is the minimum cos theta value.
 * Default values for LAT are \f$N=2\f$ and \f${\rm costhetamin}=0\f$.
 ***************************************************************************/
double GLATLtCubeMap::costheta(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_num_ctheta)
        throw GException::out_of_range(G_COSTHETA, index, 0, m_num_ctheta-1);
    #endif

    // Set cos theta scale
    double f = (index+0.5)/m_num_ctheta;
    if (m_sqrt_bin)
        f = f*f;

    // Set cos theta value
    double costheta = 1.0 - f * (1.0 - m_min_ctheta); 

    // Return costheta
    return costheta;
}


/***********************************************************************//**
 * @brief Return phi value (in radians) for an index
 *
 * @param[in] index Bin index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Bin index outside value range.
 *
 * The phi value is computed using
 * \f[\phi = \frac{{\rm index} + 0.5}{\rm nphi} \frac{\pi}{4}\f]
 * where
 * \f${\rm index}\f$ is the bin index, and
 * \f${\rm nphi}\f$ is the number of phi bins.
 ***************************************************************************/
double GLATLtCubeMap::phi(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_num_phi)
        throw GException::out_of_range(G_PHI, index, 0, m_num_phi-1);
    #endif

    // Set phi value
    double phi = (index+0.5) / m_num_phi * gammalib::pi / 4.0;

    // Return phi
    return phi;
}


/***********************************************************************//**
 * @brief Return cos theta binning scheme
 *
 * Returns either "SQRT(1-COSTHETA)" or "COSTHETA".
 ***************************************************************************/
std::string GLATLtCubeMap::costhetabin(void) const
{
    // Set string
    std::string result = (m_sqrt_bin) ? "SQRT(1-COSTHETA)" : "COSTHETA";

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Print lifetime cube map information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing lifetime cube map information.
 ***************************************************************************/
std::string GLATLtCubeMap::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GLATLtCubeMap ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of cos theta bins") +
                      gammalib::str(ncostheta()));
        result.append("\n"+gammalib::parformat("Number of phi bins") +
                      gammalib::str(nphi()));
        result.append("\n"+gammalib::parformat("Cos theta binning"));
        if (m_sqrt_bin) {
            result.append("sqrt");
        }
        else {
            result.append("linear");
        }
        result.append("\n"+gammalib::parformat("Minimum cos theta") +
                      gammalib::str(costhetamin()));
        result.append("\n"+m_map.print(chatter));

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATLtCubeMap::init_members(void)
{
    // Initialise members
    m_map.clear();
    m_num_ctheta = 0;
    m_num_phi    = 0;
    m_min_ctheta = 0.0;
    m_sqrt_bin   = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param map Livetime cube map.
 ***************************************************************************/
void GLATLtCubeMap::copy_members(const GLATLtCubeMap& map)
{
    // Copy members
    m_map        = map.m_map;
    m_num_ctheta = map.m_num_ctheta;
    m_num_phi    = map.m_num_phi;
    m_min_ctheta = map.m_min_ctheta;
    m_sqrt_bin   = map.m_sqrt_bin;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATLtCubeMap::free_members(void)
{
    // Return
    return;
}
