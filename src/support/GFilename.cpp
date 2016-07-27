/***************************************************************************
 *                      GFilename.cpp - Filename class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2016 by Juergen Knoedlseder                         *
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
#if defined(HAVE_LIBCFITSIO)
#if defined(HAVE_CFITSIO_FITSIO_H)
#include <cfitsio/fitsio.h>
#elif defined(HAVE_LIBCFITSIO0_FITSIO_H)
#include <libcfitsio0/fitsio.h>
#else
#include <fitsio.h>
#endif
#endif
#include <sys/stat.h>             // for stat structure and S_ISREG
#include <cstdio>                 // for std::remove()
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
 * @brief Checks whether file exists
 *
 * @return True if file exists.
 *
 * Checks whether a file exists on disk. In case that the file is a FITS
 * file, the method also checks whether a compressed version of the file
 * (with a .gz, .Z, .z, or .zip extension) exists on disk (see is_fits()
 * method).
 ***************************************************************************/
bool GFilename::exists(void) const
{
    // First check if file is a FITS file
    bool exists = is_fits();

    // If the file is not a FITS file, then check whether a file with the
    // given name exists
    if (!exists) {

        // Allocate file information structure
        struct stat info;

        // Get file information structure
        int ret = stat(url().c_str(), &info);

        // If the file is a regular file then signal that it exists
        if (ret == 0 && S_ISREG(info.st_mode)) {
            exists = true;
        }

    } // endif: file was not a FITS file

    // Return result
    return (exists);
}


/***********************************************************************//**
 * @brief Checks whether file is a FITS file
 *
 * @return True if file is a FITS file.
 *
 * Test if the file or a compressed version of the file (with a .gz, .Z, .z,
 * or .zip extension) is a FITS file. This method is thread safe.
 ***************************************************************************/
bool GFilename::is_fits(void) const
{
    // Initialise result
    bool is_fits = false;

    // Check now for a FITS file. This works only if cfitsio is available.
    // Put the code into a critical zone as it might be called from within
    // a parallelized thread.
    #if defined(HAVE_LIBCFITSIO)
    #pragma omp critical(GFilename_is_fits)
    {
        int       status = 0;
        fitsfile* fptr   = NULL;
        status           = ffopen(&fptr, url().c_str(), 0, &status);
        if (status == 0) {
            is_fits = true;
        }
        ffclos(fptr, &status);
    }
    #endif

    // Return result
    return (is_fits);
}


/***********************************************************************//**
 * @brief Remove file from disk
 *
 * Removes file or a compressed version of the file (with a .gz, .Z, .z,
 * or .zip extension) if it exists on disk.
 ***************************************************************************/
void GFilename::remove(void) const
{
    // Set compressed file name variants
    GFilename name_gz  = url()+".gz";
    GFilename name_Z   = url()+".Z";
    GFilename name_z   = url()+".z";
    GFilename name_zip = url()+".zip";

    // Remove file if one of the variants exists on disk
    if (name_gz.exists()) {
        std::remove(name_gz.url().c_str());
    }
    else if (name_Z.exists()) {
        std::remove(name_Z.url().c_str());
    }
    else if (name_z.exists()) {
        std::remove(name_z.url().c_str());
    }
    else if (name_zip.exists()) {
        std::remove(name_zip.url().c_str());
    }
    else if (exists()) {
        std::remove(url().c_str());
    }

    // Return
    return;
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
 * @param[in] defaultno Default extension number (default: -11).
 * @return Integer containing extension number.
 *
 * Returns the extension number. If no extension number is given the number
 * provided by the @p defaultno argument will be used. By default, the
 * @p defaultno parameter is set to -1.
 ***************************************************************************/
int GFilename::extno(const int& defaultno) const
{
    // Set the extension number
    int extno = m_extno;

    // If no extension number is set, the use the default extension number
    if (!has_extno()) {
        extno = defaultno;
    }

    // Return extension number
    return (extno);
}


/***********************************************************************//**
 * @brief Return extension version number
 *
 * @param[in] defaultver Default extension version number (default: 0).
 * @return Integer containing extension version number.
 *
 * Returns the extension version number. If no extension version number is
 * given the number provided by the @p defaultver argument will be used. By
 * default, the @p defaultextver parameter is set to 0.
 ***************************************************************************/
int GFilename::extver(const int& defaultver) const
{
    // Set the extension version
    int extver = m_extver;

    // If no extension version is set, the use the default extension version
    if (!has_extver()) {
        extver = defaultver;
    }

    // Return extension version
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

        // Append filename
        result.append("\n"+gammalib::parformat("File name"));
        result.append(m_filename);

        // Append URL
        result.append("\n"+gammalib::parformat("URL"));
        result.append(m_url);

        // Append protocol
        result.append("\n"+gammalib::parformat("Protocol"));
        result.append(m_protocol);

        // Append access path
        result.append("\n"+gammalib::parformat("Access path"));
        result.append(m_path);

        // Append file name
        result.append("\n"+gammalib::parformat("File"));
        result.append(m_file);

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

        // Append expression
        result.append("\n"+gammalib::parformat("Selection expression"));
        if (has_expression()) {
            result.append(m_expression);
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
    m_url.clear();
    m_protocol.clear();
    m_path.clear();
    m_file.clear();
    m_extname.clear();
    m_extno  = -1;
    m_extver = 0;
    m_expression.clear();

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
    m_filename   = filename.m_filename;
    m_url        = filename.m_url;
    m_protocol   = filename.m_protocol;
    m_path       = filename.m_path;
    m_file       = filename.m_file;
    m_extname    = filename.m_extname;
    m_extno      = filename.m_extno;
    m_extver     = filename.m_extver;
    m_expression = filename.m_expression;

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

    // Strip any whitespace and store the full name
    m_filename = gammalib::strip_whitespace(filename);

    // Continue only if filename is not empty
    if (!m_filename.empty()) {

        // Check for [ symbol in filename
        size_t start = m_filename.find_first_of("[");

        // If we have an extension then separate filename into filename
        // part and extension part and fill the information of th object
        if (start != std::string::npos) {

            // Check for ] symbol in filename
            size_t stop = m_filename.find_first_of("]", start);
            
            // If there is no ] symbol then throw an exception
            if (stop == std::string::npos) {
                std::string msg = "Missing ] symbol in file name \""+
                                  m_filename+"\" extension. Please correct "
                                  "the file name.";
                throw GException::invalid_argument(G_SET_FILENAME, msg);
            }

            // If ] is not the last character then check if an opening bracket is next
            if (stop < m_filename.length()-1) {

                // Get character
                std::string character = gammalib::strip_whitespace(m_filename.substr(stop+1,1));

                // If no bracket is coming, throw an exception
                if (character != "[") {
                    std::string msg = "Non-bracket character \""+character+
                                      "\" beyond ] symbol in file name \""+
                                      m_filename+"\" expression. Please "
                                      "correct the file name.";
                    throw GException::invalid_argument(G_SET_FILENAME, msg);
                }
             }

            // Extract extension name
            std::string extname = gammalib::strip_whitespace(m_filename.substr(start+1, stop-start-1));

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
                                      "symbol in file name \""+m_filename+
                                      "\". Please correct the file name.";
                    throw GException::invalid_argument(G_SET_FILENAME, msg);
                }

                // Extract extension version
                m_extver = gammalib::toint(extver);

                // Throw an exception if the extension version is not
                // positive
                if (m_extver <= 0) {
                    std::string msg = "Non-positive extension version "
                                      "encountered in file name \""+
                                      m_filename+ "\". Please correct the "
                                      "file name.";
                    throw GException::invalid_argument(G_SET_FILENAME, msg);
                }

                // Update the extension name
                extname = gammalib::strip_whitespace(extname.substr(0,sep));

            } // endif: extension version provided

            // If we have an empty extension name then throw an exception
            if (extname.empty()) {
                std::string msg = "An empty extension has been specified in "
                                  "file name \""+m_filename+"\". Please "
                                  "correct the file name.";
                throw GException::invalid_argument(G_SET_FILENAME, msg);
            }

            // If we have a purely numerical extension then convert the
            // extension name into an extension number
            if (extname.find_first_not_of("+-0123456789") == std::string::npos) {

                // Extract extension number
                m_extno = gammalib::toint(extname);

                // Throw an exception if the extension version is not
                // positive
                if (m_extno < 0) {
                    std::string msg = "Negative extension number encountered "
                                      "in file name \""+m_filename+"\". "
                                      "Please correct the file name.";
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
                                      "in file name \""+m_filename+"\". "
                                      "Please correct the file name.";
                    throw GException::invalid_argument(G_SET_FILENAME, msg);
                }

            } // endelse: extension name provided

            // Store the file name
            m_url = gammalib::strip_whitespace(m_filename.substr(0,start));

            // Check if there is an expression behind the extension name
            size_t expr_start = m_filename.find_first_of("[", stop+1);

            if (expr_start != std::string::npos) {

               // Check for ] symbol in expression
               size_t expr_stop = m_filename.find_first_of("]", expr_start);

               // If there is no ] symbol then throw an exception
               if (expr_stop == std::string::npos) {
                   std::string msg = "Missing ] symbol in file name \""+
                                     m_filename+"\" expression. Please "
                                     "correct the file name.";
                   throw GException::invalid_argument(G_SET_FILENAME, msg);
               }

               // If ] is not the last character then throw an exception
               if (expr_stop < m_filename.length()-1) {
                   std::string msg = "Characters beyond ] symbol in file "
                                     "name \""+m_filename+"\" expression. "
                                     "Please correct the file name.";
                   throw GException::invalid_argument(G_SET_FILENAME, msg);
                }

               // Store expression
               m_expression = gammalib::strip_whitespace(m_filename.substr(expr_start+1, expr_stop - expr_start-1));
            }

        } // endif: had extension

        // ... otherwise simply set the filename
        else {
            m_url = m_filename;
        }

        // Extract access protocol
        size_t proto_end = m_filename.find_first_of(":");
        if (proto_end != std::string::npos) {
            m_protocol = m_filename.substr(0, proto_end);
            proto_end++;
        }
        else {
            proto_end = 0;
        }

        // Set start and end of file name
        size_t file_start = m_filename.find_last_of("/");
        size_t file_end   = (start == std::string::npos) ? m_filename.length() : start;

        // If a "/" symbol was found then extract the path and file
        if (file_start != std::string::npos) {
        
            // Extract path
            size_t path_start = m_filename.find_first_not_of("/", proto_end);
            if (path_start != std::string::npos) {
                if (path_start > 0) {
                    path_start--;
                }
                int length = file_start - path_start + 1;
                if (length > 0) {
                    m_path = m_filename.substr(path_start, length);
                }
            }

            // Extract file name
            int length = file_end - file_start - 1;
            if (length > 0) {
                m_file = m_filename.substr(file_start+1, length);
            }
        }

        // ... otherwise there is no path and we only extract the file
        else {
            int length = file_end - proto_end;
            if (length > 0) {
                m_file = m_filename.substr(proto_end, length);
            }
        }

    } // endif: filename was not empty

    // Return
    return;
}
