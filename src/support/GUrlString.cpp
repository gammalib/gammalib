/***************************************************************************
 *                     GUrlString.cpp - String URL class                   *
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
 * @file GUrlString.cpp
 * @brief String URL class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdarg>         // std::va_list
#include <cstring>         // std::memset() function
#include <cstdio>          // std::fopen, std::fgets, std::fclose, etc...
#include "GUrlString.hpp"
#include "GTools.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPEN                 "GUrlString::open(std::string&, std::string&)"

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
GUrlString::GUrlString(void) : GUrl()
{
    // Initialise members
    init_members();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Text string constructor
 *
 * @param[in] string Text string.
 *
 * Constructs GUrlString object by assigning a text @p string to the
 * object's buffer.
 ***************************************************************************/
GUrlString::GUrlString(const std::string& string) : GUrl()
{
    // Initialise members
    init_members();

    // Open URL
    open(string);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] url Unified Resource Locator.
 ***************************************************************************/
GUrlString::GUrlString(const GUrlString& url) : GUrl(url)
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
GUrlString::~GUrlString(void)
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
 * @param[in] url Unified Resource Locator.
 * @return Unified Resource Locator.
 ***************************************************************************/
GUrlString& GUrlString::operator=(const GUrlString& url)
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
 * @brief Clear string URL
 ***************************************************************************/
void GUrlString::clear(void)
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
 * @brief Clone string URL
 *
 * @return Pointer to deep copy of string URL.
 ***************************************************************************/
GUrlString* GUrlString::clone(void) const
{
    // Clone object
    return new GUrlString(*this);
}


/***********************************************************************//**
 * @brief Open string URL
 *
 * @param[in] string Text string.
 * @param[in] mode Mode parameter (ignored).
 ***************************************************************************/
void GUrlString::open(const std::string& string, const std::string& mode)
{
    // First close any existing URL
    close();

    // Store text string and initialise index
    m_buffer = string;
    m_index  = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close file
 ***************************************************************************/
void GUrlString::close(void)
{
    // Clear buffer
    m_buffer.clear();

    // Initialise index
    m_index = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read block of data from string buffer
 *
 * @param[in] buffer Data buffer.
 * @param[in] nbyte Number of Bytes to be read.
 * @return Number of Bytes that were effectively read.
 *
 * Reads @p nbyte Bytes from the string into a @p buffer. The position
 * indicator of the string is advanced by the total amount of bytes read.
 *
 * The total number of Bytes successfully read is returned. If this number
 * differs from the @p nbyte parameter, the end of the string was reached.
 * The proper indicator is set.
 *
 * If either @p buffer is NULL or @p nbyte is zero, the method returns zero
 * and both the string state and the content pointed by @p buffer remain
 * unchanged.
 *
 * If no string exists, the method returns zero and both the string state
 * and the content pointed by @p buffer remain unchanged.
 ***************************************************************************/
int GUrlString::read(void* buffer, const int& nbyte)
{
    // Initialise number of Bytes read
    int nread = 0;

    // Continue only if buffer is valid and nbyte is positive
    if (buffer != NULL && nbyte > 0) {

        // Continue only if we have a string
        if (!m_buffer.empty()) {

            // Determine remaining number of Bytes in string
            int nremain = m_buffer.length() - m_index;

            // If Bytes are remaining then copy them to buffer
            if (nremain > 0) {

                // Set number of Bytes to extract
                nread = (nbyte < nremain) ? nbyte : nremain;

                // Extract Bytes
                const char* src = m_buffer.c_str() + m_index;
                char*       dst = reinterpret_cast<char*>(buffer);
                for (int i = 0; i < nread; ++i) {
                    *dst++ = *src++;
                }

                // Forward position indicator
                m_index += nread;

            } // endif: Bytes were remaining

        } // endif: there was a string

    } // endif: input parameters were valid

    // Return number of Bytes read
    return nread;
}


/***********************************************************************//**
 * @brief Write block of data buffer into string
 *
 * @param[in] buffer Data buffer.
 * @param[in] nbyte Number of Bytes to be written.
 * @return Number of Bytes that were effectively written.
 *
 * Writes @p nbyte Bytes from a @p buffer into the string. The position
 * indicator of the string is advanced by the total amount of bytes written.
 *
 * The total number of Bytes successfully written is returned. 
 *
 * If either @p buffer is NULL or @p nbyte is zero, the method returns zero
 * and both the string state and the content pointed by @p buffer remain
 * unchanged.
 ***************************************************************************/
int GUrlString::write(const void* buffer, const int& nbyte)
{
    // Initialise number of Bytes written
    int nwritten = 0;

    // Continue only if buffer is valid and nbyte is positive
    if (buffer != NULL && nbyte > 0) {

        // Create string to append
        std::string sbuffer;
        sbuffer.reserve(nbyte);
        const char* src = reinterpret_cast<const char*>(buffer);
        for (int i = 0; i < nbyte; ++i) {
            sbuffer.push_back(*src++);
        }

        // Append string to buffer
        m_buffer.append(sbuffer);

        // Set number of Bytes appended
        nwritten = nbyte;

        // Forward position indicator
        m_index += nwritten;

    } // endif: input parameters were valid

    // Return number of Bytes written
    return nwritten;
}


/***********************************************************************//**
 * @brief Return next character from string
 *
 * @return Next character in string.
 *
 * Returns the character currently pointed by the internal position
 * indicator. The internal position indicator is then advanced to the next
 * character.
 *
 * If the indicator is at the end of the string when called, the function
 * returns EOF.
 *
 * If no string exists, the method returns EOF.
 ***************************************************************************/
int GUrlString::getchar(void)
{
    // Initialise character to EOF
    int character = EOF;

    // Continue only if we have a string and if the index has not yet reached
    // the end of the buffer
    if (!m_buffer.empty() && m_index < m_buffer.length()) {

        // Get next character
        character = m_buffer[m_index];

        // Forward position indicator
        m_index++;

    } // endif: we had a string

    // Return character
    return character;
}


/***********************************************************************//**
 * @brief Write character into string
 *
 * @param[in] character Character.
 *
 * Writes a character to the string and advances the position indicator.
 ***************************************************************************/
void GUrlString::putchar(const int& character)
{
    // Append character to buffer
    m_buffer.push_back(character);

    // Forward position indicator
    m_index++;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read formatted data from string
 *
 * @param[in] format Format.
 * @param[in] ... Optional parameters.
 *
 * Reads data from a string and stores them according to the parameter format
 * into the locations pointed by the additional arguments. The additional
 * arguments should point to already allocated objects of the type specified
 * by their corresponding format specifier within the format string.
 *
 * If no string exists, the method does nothing.
 *
 * @todo The position indicator is not forwared as I have no idea how to do
 * this in fact!!!
 ***************************************************************************/
void GUrlString::scanf(const char* format, ...)
{
    // Continue only if we have a string and if the index has not yet reached
    // the end of the buffer
    if (!m_buffer.empty() && m_index < m_buffer.length()) {

        // Declare argument pointer
        std::va_list arg_ptr;

        // Set start argument
        va_start(arg_ptr, format);

        // Get pointer to source
        const char* src = m_buffer.c_str() + m_index;

        // Read data from string
        vsscanf(src, format, arg_ptr);

        // Stop argument reading
        va_end(arg_ptr);

    } // endif: we had a string

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write formatted data into string
 *
 * @param[in] format Format.
 * @param[in] ... Optional parameters.
 *
 * Writes the C string pointed by format to the string. If format includes
 * format specifiers (subsequences beginning with %), the additional
 * arguments following format are formatted and inserted in the resulting
 * string replacing their respective specifiers.
 *
 * After the format parameter, the function expects at least as many
 * additional arguments as specified by format.
 ***************************************************************************/
void GUrlString::printf(const char* format, ...)
{
    // Allocate buffer for printing. The length here is fixed, which is
    // not a good thing as we cannot really know how many Bytes arrive.
    // But how else can we deal with this???
    char buffer[8192];
    std::memset(&buffer, 0, sizeof(buffer));

    // Declare argument pointer
    std::va_list arg_ptr;

    // Set start argument
    va_start(arg_ptr, format);

    // Write data into buffer
    std::vsprintf(buffer, format, arg_ptr);

    // Stop argument reading
    va_end(arg_ptr);

    // Convert buffer to C++ string
    std::string sbuffer = std::string(buffer);

    // If not empty then append buffer
    int nwritten = sbuffer.length();
    if (nwritten > 0) {
        m_buffer.append(buffer);
        m_index += nwritten;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print URL information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing URL information.
 ***************************************************************************/
std::string GUrlString::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GUrlString ===");

        // Append URL information
        result.append("\n"+parformat("String size")+str(m_buffer.length()));
        result.append("\n"+parformat("String position indicator")+str(m_index));

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
void GUrlString::init_members(void)
{
    // Initialise members
    m_index = 0;
    m_buffer.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] url URL.
 ***************************************************************************/
void GUrlString::copy_members(const GUrlString& url)
{
    // Copy members
    m_index  = url.m_index;
    m_buffer = url.m_buffer;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GUrlString::free_members(void)
{
    // Close file
    close();

    // Return
    return;
}
