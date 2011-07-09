/***************************************************************************
 *               GException_sky.cpp  -  sky exception handlers             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GException_sky.cpp
 * @brief Sky module exception classes implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief General GSkymap error
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Error message.
 ***************************************************************************/
GException::skymap::skymap(std::string origin,
                           std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Bad sky map parameter
 *
 * @param[in] origin Method that throws the error.
 * @param[in] par Map parameter.
 * @param[in] message Error message.
 ***************************************************************************/
GException::skymap_bad_par::skymap_bad_par(std::string origin,
                                           int         par,
                                           std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid map parameter ("+str(par)+"). " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Incompatible data size
 *
 * @param[in] origin Method that throws the error.
 * @param[in] size Encountered data size.
 * @param[in] expected Expected data size.
 * @param[in] message Error message.
 ***************************************************************************/
GException::skymap_bad_size::skymap_bad_size(std::string origin,
                                             int         size,
                                             int         expected,
                                             std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Data size ("+str(size)+") is incompatible with "
                "expected size ("+str(expected)+"). " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Incompatible ctype keywords
 *
 * @param[in] origin Method that throws the error.
 * @param[in] ctype1 Encountered ctype1 keyword.
 * @param[in] ctype2 Encountered ctype2 keyword.
 * @param[in] message Error message.
 ***************************************************************************/
GException::skymap_bad_ctype::skymap_bad_ctype(std::string origin,
                                               std::string ctype1,
                                               std::string ctype2,
                                               std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "CTYPE1 "+ctype1+") and CTYPE2 ("+ctype2+") header keywords"
                " are not compatible. " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Bad image dimension
 *
 * @param[in] origin Method that throws the error.
 * @param[in] naxis Encountered number of axes.
 * @param[in] message Error message.
 ***************************************************************************/
GException::skymap_bad_image_dim::skymap_bad_image_dim(std::string origin,
                                                       int         naxis,
                                                       std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid number of image dimensions (naxis="+str(naxis)+")."
                " " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief General GWcs error
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Error message.
 ***************************************************************************/
GException::wcs::wcs(std::string origin,
                     std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid WCS code
 *
 * @param[in] origin Method that throws the error.
 * @param[in] wcs WCS.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::wcs_invalid::wcs_invalid(std::string origin,
                                     std::string wcs,
                                     std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid World Coordinate System code ("+wcs+"). " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Coordinate system invalid
 *
 * @param[in] origin Method that throws the error.
 * @param[in] coordsys Specified coordinate system.
 ***************************************************************************/
GException::wcs_bad_coords::wcs_bad_coords(std::string origin,
                                           std::string coordsys)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid WCS coordinate system ("+coordsys+"). "
                "Should be one of EQU/CEL/GAL/ECL/HEL/SGL.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief No projection function declared
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::wcs_no_proj_fct::wcs_no_proj_fct(std::string origin,
                                             std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "No projection function declared. " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Healpix resolution invalid
 *
 * @param[in] origin Method that throws the error.
 * @param[in] nside Specified nside parameter.
 ***************************************************************************/
GException::wcs_hpx_bad_nside::wcs_hpx_bad_nside(std::string origin,
                                                 int         nside)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid Healpix nside value ("+str(nside)+"). "
                "Should be one of 1,2,4,8,16,32,64,128,256,512,1024,2048,"
                "4096,8192.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Healpix pixel ordering invalid
 *
 * @param[in] origin Method that throws the error.
 * @param[in] ordering Specified pixel ordering.
 ***************************************************************************/
GException::wcs_hpx_bad_ordering::wcs_hpx_bad_ordering(std::string origin,
                                                       std::string ordering)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid pixel ordering ("+ordering+"). "
                "Should be either 'RING' or 'NESTED'.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Singular matrix encountered
 *
 * @param[in] origin Method that throws the error.
 * @param[in] naxis Number of matrix axes.
 * @param[in] mat Matrix (stored in vector).
 ***************************************************************************/
GException::wcs_singular_matrix::wcs_singular_matrix(std::string                origin,
                                                     int                        naxis,
                                                     const std::vector<double>& mat)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Singular matrix encountered: mat=[";
    for (int i = 0, ij = 0; i < naxis; ++i) {
        if (i > 0)
            m_message += ",";
        m_message += "[";
        for (int j = 0; j < naxis; ++j, ++ij) {
            if (j > 0)
                m_message += ",";
            m_message += str(mat[ij]);
        }
        m_message += "]";
    }
    m_message += "].";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid WCS parameter
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::wcs_invalid_parameter::wcs_invalid_parameter(std::string origin,
                                                         std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid WCS parameter encountered.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid (x,y) coordinate(s)
 *
 * @param[in] origin Method that throws the error.
 * @param[in] num Number of invalid coordinates.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::wcs_invalid_x_y::wcs_invalid_x_y(std::string origin,
                                             int         num,
                                             std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = str(num)+" (x,y) coordinates were invalid.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid (phi,theta) coordinate(s)
 *
 * @param[in] origin Method that throws the error.
 * @param[in] num Number of invalid coordinates.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::wcs_invalid_phi_theta::wcs_invalid_phi_theta(std::string origin,
                                                         int         num,
                                                         std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = str(num)+" (phi,theta) coordinates were invalid.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}

