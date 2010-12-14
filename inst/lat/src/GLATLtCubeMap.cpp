/***************************************************************************
 *            GLATLtCubeMap.cpp  -  Fermi LAT lifetime cube map            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATLtCubeMap.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS_OP1                    "GLATLtCubeMap::operator() (double&)"
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
 * @param map Lifetime cube map.
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
 * @param[in] map Lifetime cube map.
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
 * most rapidely varying parameter.
 ***************************************************************************/
double GLATLtCubeMap::operator() (const GSkyDir& dir, _ltcube_ctheta_phi fct)
{
    // Get map index
    int pixel = m_map.dir2pix(dir);

    // Initialise sum
    double sum = 0.0;

    // Loop over azimuth and zenith angles
    for (int iphi = 0, i = 0; iphi < m_num_phi; ++iphi) {
        double p = phi(iphi);
        for (int itheta = 0; itheta < m_num_ctheta; ++itheta, ++i)
            sum += m_map(pixel, i) * (*fct)(costheta(itheta), p);
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
 * @brief Load lifetime cube from FITS file
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
void GLATLtCubeMap::read(const GFitsTable* hdu)
{
    // Clear object
    clear();

    // Load skymap
    m_map.read(hdu);

    // Set costheta binning scheme
    std::string scheme = strip_whitespace(toupper(hdu->string("THETABIN")));
    m_sqrt_bin = (scheme == "SQRT(1-COSTHETA)");

    // Read attributes
    m_num_ctheta = hdu->integer("NBRBINS");
    m_num_phi    = hdu->integer("PHIBINS");
    m_min_ctheta = hdu->real("COSMIN");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save lifetime cube map into FITS file
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
 ***************************************************************************/
double GLATLtCubeMap::phi(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_num_phi)
        throw GException::out_of_range(G_PHI, index, 0, m_num_phi-1);
    #endif

    // Set phi value
    double phi = (index+0.5) / m_num_phi * pi / 4.0;

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
 ***************************************************************************/
std::string GLATLtCubeMap::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GLATLtCubeMap ===");
    result.append("\n"+parformat("Number of cos theta bins")+str(ncostheta()));
    result.append("\n"+parformat("Number of phi bins")+str(nphi()));
    result.append("\n"+parformat("Cos theta binning"));
    if (m_sqrt_bin)
        result.append("sqrt");
    else
        result.append("linear");
    result.append("\n"+parformat("Minimum cos theta")+str(costhetamin()));
    result.append("\n"+m_map.print());

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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param map Lifetime cube map.
 ***************************************************************************/
void GLATLtCubeMap::copy_members(const GLATLtCubeMap& map)
{
    // Copy members
    m_map        = map.m_map;
    m_num_ctheta = map.m_num_ctheta;
    m_num_phi    = map.m_num_phi;

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


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] map Lifetime cube map.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GLATLtCubeMap& map)
{
     // Write lifetime cube map in output stream
    os << map.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] map Lifetime cube map.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GLATLtCubeMap& map)
{
    // Write lifetime cube map into logger
    log << map.print();

    // Return logger
    return log;
}
