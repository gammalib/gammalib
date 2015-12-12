/***************************************************************************
 *                      GFilename.cpp - Filename class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Juergen Knoedlseder                              *
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
 * @file GFilename.cpp
 * @brief Filename class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GFilename.hpp"
#include "GTools.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_SET_FILENAME                "GFilename::set_filename(std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GFilename::GFilename(void)
{
    // Initialise members
    init_members();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Filename constructor
 *
 * @param[in] filename File name.
 *
 * Constructs GFilename object by assigning a @p filename.
 ***************************************************************************/
GFilename::GFilename(const std::string& filename)
{
    // Initialise members
    init_members();

    // Set filename
    set_filename(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Filename constructor
 *
 * @param[in] filename File name.
 *
 * Constructs GFilename object by assigning a @p filename.
 ***************************************************************************/
GFilename::GFilename(const char* filename)
{
    // Initialise members
    init_members();

    // Set filename
    set_filename(std::string(filename));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] filename File name.
 ***************************************************************************/
GFilename::GFilename(const GFilename& filename)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFilename::~GFilename(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] filename File name.
 * @return File name.
 ***************************************************************************/
GFilename& GFilename::operator=(const GFilename& filename)
{
    // Execute only if object is not identical
    if (this != &filename) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(filename);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear file name
 ***************************************************************************/
void GFilename::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone file name
 *
 * @return Pointer to deep copy of file name.
 ***************************************************************************/
GFilename* GFilename::clone(void) const
{
    // Clone object
    return new GFilename(*this);
}


/***********************************************************************//**
 * @brief Return extension name
 *
 * @param[in] defaultname Default extension name (default: "").
 * @return String containing extension name.
 *
 * Returns the extension name. If no extension name is given the name
 * provided by the @p defaultname argument will be used. By default, the
 * @p defaultname parameter is an empty string.
 ***************************************************************************/
std::string GFilename::extname(const std::string& defaultname) const
{
    // Set the extension name
    std::string extname = m_extname;

    // If no extension name is set, the use the default extension name
    if (!has_extname()) {
        extname = defaultname;
    }

    // Return extension name
    return (extname);
}


/***********************************************************************//**
 * @brief Return extension number
 *
 * @param[in] defaultextno Default extension number (default: 1).
 * @return Integer containing extension number.
 *
 * Returns the extension number. If no extension number is given the number
 * provided by the @p defaultextno argument will be used. By default, the
 * @p defaultextno parameter is set to 1.
 ***************************************************************************/
int GFilename::extno(const int& defaultextno) const
{
    // Set the extension number
    int extno = m_extno;

    // If no extension number is set, the use the default extension number
    if (!has_extno()) {
        extno = defaultextno;
    }

    // Return extension number
    return (extno);
}


/***********************************************************************//**
 * @brief Return extension version number
 *
 * @param[in] defaultextver Default extension version number (default: 1).
 * @return Integer containing extension version number.
 *
 * Returns the extension version number. If no extension version number is given the number
 * provided by the @p defaultextver argument will be used. By default, the
 * @p defaultextver parameter is set to 1.
 ***************************************************************************/
int GFilename::extver(const int& defaultextver) const
{
    // Set the extension version number
    int extver = m_extver;

    // If no extension number is set, the use the default extension version number
    if (!has_extver()) {
        extver = defaultextver;
    }

    // Return extension number
    return (extver);
}



/***********************************************************************//**
 * @brief Print file name information
 *
 * @param[in] chatter Chattiness (default: NORMAL).
 * @return String containing file name information.
 ***************************************************************************/
std::string GFilename::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GFilename ===");

        // Append file name
        result.append("\n"+gammalib::parformat("File name"));
        result.append(m_filename);

        // Append extension name
        result.append("\n"+gammalib::parformat("Extension name"));
        if (has_extname()) {
            result.append(m_extname);
        }
        else {
            result.append("Not provided");
        }

        // Append extension number
        result.append("\n"+gammalib::parformat("Extension number"));
        if (has_extno()) {
            result.append(gammalib::str(m_extno));
        }
        else {
            result.append("Not provided");
        }

        // Append extension version
        result.append("\n"+gammalib::parformat("Extension version"));
        if (has_extver()) {
            result.append(gammalib::str(m_extver));
        }
        else {
            result.append("Not provided");
        }

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
void GFilename::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_extname.clear();
    m_extno  = -1;
    m_extver = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] filename File name.
 ***************************************************************************/
void GFilename::copy_members(const GFilename& filename)
{
    // Copy members
    m_filename = filename.m_filename;
    m_extname  = filename.m_extname;
    m_extno    = filename.m_extno;
    m_extver   = filename.m_extver;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFilename::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set file name
 *
 * @param[in] filename File name.
 *
 * @exception GException::invalid_argument
 *            Invalid file name specified.
 *
 * Sets all attributes of the object by parsing the file name. The file name
 * can have one of the following formats:
 *
 *      file.fits            No extension
 *      file.fits[1]         Extension number 1 (starting from 0)
 *      file.fits[EVENTS]    EVENTS extension
 *      file.fits[EVENTS,2]  Version 2 of the EVENTS extension
 ***************************************************************************/
void GFilename::set_filename(const std::string& filename)
{
    // Clear file name
    clear();

    // Strip any whitespace
    std::string fname = gammalib::strip_whitespace(filename);

    // Continue only if filename is not empty
    if (!fname.empty()) {

        // Check for [ symbol in filename
        size_t start = fname.find_first_of("[");

        // If we have an extension then separate filename into filename
        // part and extension part and fill the information of th object
        if (start != std::string::npos) {

            // Check for ] symbol in filename
            size_t stop = fname.find_first_of("]", start);
            
            // If there is no ] symbol then throw an exception
            if (stop == std::string::npos) {
                std::string msg = "Missing ] symbol in filename \""+
                                  fname+"\" extension. Please correct the "
                                  "filename.";
                throw GException::invalid_argument(G_SET_FILENAME, msg);
            }

            // If ] is not the last character then throw an exception
            if (stop < fname.length()-1) {
                std::string msg = "Characters beyond ] symbol in filename \""+
                                  fname+"\". Please correct the filename.";
                throw GException::invalid_argument(G_SET_FILENAME, msg);
            }

            // Extract extension name
            std::string extname = gammalib::strip_whitespace(fname.substr(start+1, stop-start-1));

            // Check if there is a , symbol that separates extension name
            // and extension version
            size_t sep = extname.find_first_of(",");

            // If we have a separator, then extract the extension version
            if (sep != std::string::npos) {

                // Get extension version string
                std::string extver = gammalib::strip_whitespace(extname.substr(sep+1, std::string::npos));

                // If string is empty then throw an exception
                if (extver.empty()) {
                    std::string msg = "No extension version found after , "
                                      "symbol in filename \""+fname+
                                      "\". Please correct the filename.";
                    throw GException::invalid_argument(G_SET_FILENAME, msg);
                }

                // Extract extension version
                m_extver = gammalib::toint(extver);

                // Throw an exception if the extension version is not
                // positive
                if (m_extver <= 0) {
                    std::string msg = "Non-positive extension version "
                                      "encountered in filename \""+fname+
                                      "\". Please correct the filename.";
                    throw GException::invalid_argument(G_SET_FILENAME, msg);
                }

                // Update the extension name
                extname = gammalib::strip_whitespace(extname.substr(0,sep));

            } // endif: extension version provided

            // If we have a purely numerical extension then convert the
            // extension name into an extension number
            if (extname.find_first_of("0123456789") != std::string::npos) {

                // Extract extension number
                m_extno = gammalib::toint(extname);

                // Throw an exception if the extension version is not
                // positive
                if (m_extno < 0) {
                    std::string msg = "Negative extension number encountered "
                                      "in filename \""+fname+"\". Please "
                                      "correct the filename.";
                    throw GException::invalid_argument(G_SET_FILENAME, msg);
                }

            } // endif: extension number provided

            // ... otherwise store the extension name
            else {

                // Store extension name
                m_extname = gammalib::toupper(extname);

                // Throw an exception if extension name is empty
                if (m_extname.empty()) {
                    std::string msg = "Empty extension name encountered "
                                      "in filename \""+fname+"\". Please "
                                      "correct the filename.";
                    throw GException::invalid_argument(G_SET_FILENAME, msg);
                }

            } // endelse: extension name provided

            // Store the file name
            m_filename = gammalib::strip_whitespace(fname.substr(0,start));

        } // endif: had extension

        // ... otherwise simply set the filename
        else {
            m_filename = fname;
        }

    } // endif: filename was not empty

    // Return
    return;
}
