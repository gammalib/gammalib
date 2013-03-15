/***************************************************************************
 *                       GUrlFile.hpp - File URL class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GUrlFile.cpp
 * @brief File URL class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdarg>      // std::va_list
#include <cstdio>       // std::fopen, std::fgets, std::fclose, etc...
#include "GUrlFile.hpp"
#include "GTools.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPEN                   "GUrlFile::open(std::string&, std::string&)"

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
GUrlFile::GUrlFile(void) : GUrl()
{
    // Initialise members
    init_members();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Opening constructor
 *
 * @param[in] url File name.
 * @param[in] mode File mode.
 *
 * Constructs GUrlFile object by opening a file @p url in the specified
 * @p mode. Any environment variable present in the file name will be
 * automatically expanded.
 ***************************************************************************/
GUrlFile::GUrlFile(const std::string& url, const std::string& mode) : GUrl()
{
    // Initialise members
    init_members();

    // Open file
    open(url, mode);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] url URL.
 ***************************************************************************/
GUrlFile::GUrlFile(const GUrlFile& url) : GUrl(url)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(url);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GUrlFile::~GUrlFile(void)
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
 * @param[in] url URL.
 * @return URL.
 ***************************************************************************/
GUrlFile& GUrlFile::operator=(const GUrlFile& url)
{
    // Execute only if object is not identical
    if (this != &url) {

        // Copy base class members
        this->GUrl::operator=(url);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(url);

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
 * @brief Clear instance
 ***************************************************************************/
void GUrlFile::clear(void)
{
    // Free class members
    free_members();
    this->GUrl::free_members();

    // Initialise members
    this->GUrl::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of file URL.
 ***************************************************************************/
GUrlFile* GUrlFile::clone(void) const
{
    // Clone object
    return new GUrlFile(*this);
}


/***********************************************************************//**
 * @brief Open file
 *
 * @param[in] url File name.
 * @param[in] mode File mode.
 *
 * @exception GException::file_not_found
 *            File not found.
 * @exception GException::file_open_error
 *            Unable to open file.
 *
 * Opens a file @p url in the specified @p mode. Any environment variable
 * present in the filename will be automatically expanded.
 *
 * @todo Strip any file:// prefix
 ***************************************************************************/
void GUrlFile::open(const std::string& url, const std::string& mode)
{
    // First close any existing file
    close();

    // Expand environment variables
    std::string filename = expand_env(url);

    // Check if file exists
    if (!file_exists(filename)) {
        throw GException::file_not_found(G_OPEN, filename);
    }

    // Try opening file. Throw an exception if opening failed.
    m_fptr = std::fopen(filename.c_str(), mode.c_str());
    if (m_fptr == NULL) {
        throw GException::file_open_error(G_OPEN, filename);
    }

    // Store URL and mode
    m_url  = url;
    m_mode = mode;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close file
 ***************************************************************************/
void GUrlFile::close(void)
{
    // Close file
    if (m_fptr != NULL) {
        std::fclose(m_fptr);
    }

    // Reset members
    m_url.clear();
    m_mode.clear();
    m_fptr = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read block of data from file in buffer
 *
 * @param[in] buffer Data buffer.
 * @param[in] nbyte Number of Bytes to be read.
 * @return Number of Bytes that were effectively read.
 *
 * Reads @p nbyte Bytes from the file into a @p buffer. The position
 * indicator of the file is advanced by the total amount of bytes read.
 *
 * The total number of Bytes successfully read is returned. If this number
 * differs from the @p nbyte parameter, either a reading error occurred or
 * the end-of-file was reached while reading. In both cases, the proper
 * indicator is set.
 *
 * If either @p buffer is NULL or @p nbyte is zero, the method returns zero
 * and both the file state and the content pointed by @p buffer remain
 * unchanged.
 *
 * If no file has been opened, the method returns zero and both the file
 * state and the content pointed by @p buffer remain unchanged.
 ***************************************************************************/
int GUrlFile::read(void* buffer, const int& nbyte)
{
    // Initialise number of Bytes read
    int nread = 0;

    // Continue only if file is opened
    if (m_fptr != NULL) {

        // Continue only if buffer is valid and nbyte is positive
        if (buffer != NULL && nbyte > 0) {
            nread = std::fread(buffer, 1, nbyte, m_fptr);
        }

    } // endif: File was opened

    // Return number of Bytes read
    return nread;
}


/***********************************************************************//**
 * @brief Write block of data buffer into file
 *
 * @param[in] buffer Data buffer.
 * @param[in] nbyte Number of Bytes to be written.
 * @return Number of Bytes that were effectively written.
 *
 * Writes @p nbyte Bytes from a @p buffer into the file. The position
 * indicator of the file is advanced by the total amount of bytes written.
 *
 * The total number of Bytes successfully written is returned. If this number
 * differs from the @p nbyte parameter, a writing error prevented the
 * function from completing. 
 *
 * If either @p buffer is NULL or @p nbyte is zero, the method returns zero
 * and both the file state and the content pointed by @p buffer remain
 * unchanged.
 *
 * If no file has been opened, the method returns zero and both the file
 * state and the content pointed by @p buffer remain unchanged.
 ***************************************************************************/
int GUrlFile::write(const void* buffer, const int& nbyte)
{
    // Initialise number of Bytes written
    int nwritten = 0;

    // Continue only if file is opened
    if (m_fptr != NULL) {

        // Continue only if buffer is valid and nbyte is positive
        if (buffer != NULL && nbyte > 0) {
            nwritten = std::fwrite(buffer, 1, nbyte, m_fptr);
        }

    } // endif: File was opened

    // Return number of Bytes written
    return nwritten;
}


/***********************************************************************//**
 * @brief Return next character from file
 *
 * @return Next character in file.
 *
 * Returns the character currently pointed by the internal file position
 * indicator of the file. The internal file position indicator is then
 * advanced to the next character.
 *
 * If the stream is at the end-of-file when called, the function returns EOF
 * and sets the end-of-file indicator for the file.
 *
 * If a read error occurs, the method returns EOF.
 *
 * If no file has been opened, the method returns EOF.
 ***************************************************************************/
int GUrlFile::getchar(void)
{
    // Initialise character to EOF
    int character = EOF;

    // Continue only if file is opened
    if (m_fptr != NULL) {

        // Get next character
        character = std::fgetc(m_fptr);

    } // endif: File was opened

    // Return character
    return character;
}


/***********************************************************************//**
 * @brief Write character into file
 *
 * @param[in] character Character.
 *
 * Writes a character to the file and advances the position indicator.
 *
 * If no file has been opened, the method does nothing.
 ***************************************************************************/
void GUrlFile::putchar(const int& character)
{
    // Continue only if file is opened
    if (m_fptr != NULL) {

        // Write character
        std::fputc(character, m_fptr);

    } // endif: File was opened

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read formatted data from file
 *
 * @param[in] format Format.
 * @param[in] ... Optional parameters.
 *
 * Reads data from a file and stores them according to the parameter format
 * into the locations pointed by the additional arguments. The additional
 * arguments should point to already allocated objects of the type specified
 * by their corresponding format specifier within the format string.
 *
 * If no file has been opened, the method does nothing.
 ***************************************************************************/
void GUrlFile::scanf(const char* format, ...)
{
    // Continue only if file is opened
    if (m_fptr != NULL) {

        // Declare argument pointer
        std::va_list arg_ptr;

        // Set start argument
        va_start(arg_ptr, format);

        // Read data from file
        vfscanf(m_fptr, format, arg_ptr);

        // Stop argument reading
        va_end(arg_ptr);

    } // endif: File was opened

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write formatted data into file
 *
 * @param[in] format Format.
 * @param[in] ... Optional parameters.
 *
 * Writes the C string pointed by format to the file. If format includes
 * format specifiers (subsequences beginning with %), the additional
 * arguments following format are formatted and inserted in the resulting
 * string replacing their respective specifiers.
 *
 * After the format parameter, the function expects at least as many
 * additional arguments as specified by format.
 *
 * If no file has been opened, the method does nothing.
 ***************************************************************************/
void GUrlFile::printf(const char* format, ...)
{
    // Continue only if file is opened
    if (m_fptr != NULL) {

        // Declare argument pointer
        std::va_list arg_ptr;

        // Set start argument
        va_start(arg_ptr, format);

        // Write data into file
        std::vfprintf(m_fptr, format, arg_ptr);

        // Stop argument reading
        va_end(arg_ptr);

    } // endif: File was opened

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print URL information
 *
 * @return String containing URL information.
 ***************************************************************************/
std::string GUrlFile::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GUrlFile ===");

    // Append information
    if (m_fptr != NULL) {
    
        // Get information
        int pos = std::ftell(m_fptr);
        int err = std::ferror(m_fptr);
        
        // Append it
        result.append("\n"+parformat("File URL")+m_url);
        result.append("\n"+parformat("File mode")+m_mode);
        result.append("\n"+parformat("File position indicator")+str(pos));
        result.append("\n"+parformat("File error"));
        if (err == 0) {
            result.append("none");
        }
        else {
            result.append(str(err));
        }
    }
    else {
        result.append("\n"+parformat("URL")+"none");
    }

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
void GUrlFile::init_members(void)
{
    // Initialise members
    m_url.clear();
    m_mode.clear();
    m_fptr = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] url URL.
 *
 * @todo The simple copying of the file pointer is not clean as we have no
 * control anymore over the object. The file pointer should be cloned
 * somehow.
 ***************************************************************************/
void GUrlFile::copy_members(const GUrlFile& url)
{
    // Copy members
    m_url  = url.m_url;
    m_mode = url.m_mode;
    m_fptr = url.m_fptr;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GUrlFile::free_members(void)
{
    // Close file
    close();

    // Return
    return;
}
