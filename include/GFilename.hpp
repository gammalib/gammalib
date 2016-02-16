/***************************************************************************
 *                      GFilename.hpp - Filename class                     *
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
 * @file GFilename.hpp
 * @brief Filename class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GFILENAME_HPP
#define GFILENAME_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @class GFilename
 *
 * @brief Filename class
 *
 * This class handles filenames. A filename is a string composed of an
 * optional protocol (http:, ftp:, file:), an absolute or relative access
 * path, a file, and optionally a FITS extension. Examples of valid file
 * names are
 *
 *     myfits.fits
 *     myfile.fits[EVENTS]
 *     ./data/myfile.fits
 *     ~/data/myfile.fits
 *     /home/myuser/data/myfile.fits
 *     http://www.irap.omp.eu/data/myfile.fits
 *     ftp://www.irap.omp.eu/data/myfile.fits
 *     file:///home/myuser/data/myfile.fits
 *
 * A filename without the optional FITS extension is called a Uniform
 * Resource Locator (URL) an is accessed using the url() method. The URL
 * can be decomposed into the protocol, access path and the filename using
 * the protocol(), path(), and file().
 *
 * The FITS extension is implemente using the GFitsExtension class.
 ***************************************************************************/
class GFilename : public GBase {

    // Friend functions
    friend std::string operator+(const GFilename& filename, const std::string& string);
    friend std::string operator+(const std::string& string, const GFilename& filename);
    friend bool        operator==(const GFilename &a, const GFilename &b);
    friend bool        operator!=(const GFilename &a, const GFilename &b);
    
public:
    // Constructors and destructors
    GFilename(void);
    GFilename(const std::string& filename);
    GFilename(const char* filename);
    GFilename(const GFilename& filename);
    virtual ~GFilename(void);

    // Operators
    GFilename&         operator=(const GFilename& filename);
                       operator std::string(void) const;

    // Methods
    void               clear(void);
    GFilename*         clone(void) const;
    std::string        classname(void) const;
    bool               is_empty(void) const;
    int                length(void) const;
    std::string        url(void) const;
    //std::string        protocol(void) const;
    //std::string&       path(void) const;
    //std::string&       file(void) const;
    bool               exists(void) const;
    bool               is_fits(void) const;
    void               remove(void) const;
    std::string        extname(const std::string& defaultname = "") const;
    const std::string& expression(void) const;
    int                extno(const int& defaultno = -1) const;
    int                extver(const int& defaultver = 0) const;
    bool               has_extname(void) const;
    bool               has_extno(void) const;
    bool               has_extver(void) const;
    bool               has_expression(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GFilename& filename);
    void free_members(void);
    void set_filename(const std::string& filename);

    // Protected members
    std::string m_filename;   //!< Full file name
    //std::string m_protocol;   //!< Access protocol
    std::string m_url;        //!< File name (with stripped extension info)
    //std::string m_path;       //!< Path
    std::string m_extname;    //!< Extension name ("": not set)
    int         m_extno;      //!< Extension number  (-1: not set)
    int         m_extver;     //!< Extension version (0: not set)
    std::string m_expression; //!< Selection expression ("": not set)
};


/***********************************************************************//**
 * @brief Implicit filename std::string convertor
 *
 * @return Full filename as std::string.
 *
 * Returns the full filename including any FITS extension as std::string.
 ***************************************************************************/
inline
GFilename::operator std::string(void) const
{
    return (gammalib::expand_env(m_filename));
}


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFilename").
 ***************************************************************************/
inline
std::string GFilename::classname(void) const
{
    return ("GFilename");
}


/***********************************************************************//**
 * @brief Signal if filename is empty
 *
 * @return True if filename is empty, false otherwise.
 ***************************************************************************/
inline
bool GFilename::is_empty(void) const
{
    return (m_filename.empty());
}


/***********************************************************************//**
 * @brief Return length of filename
 *
 * @return Length of filename.
 *
 * Returns the length of the filename, excluding any FITS extension.
 ***************************************************************************/
inline
int GFilename::length(void) const
{
    return (m_url.length());
}


/***********************************************************************//**
 * @brief Return Uniform Resource Locator (URL)
 *
 * @return Uniform Resource Locator without FITS extension.
 *
 * Returns the Uniform Resource Locator without FITS extension. Any
 * environment variable in the URL string will be expanded.
 ***************************************************************************/
inline
std::string GFilename::url(void) const
{
    return (gammalib::expand_env(m_url));
}


/***********************************************************************//**
 * @brief Return expression name
 *
 * @return String containing file expression.
 ***************************************************************************/
inline
const std::string& GFilename::expression(void) const
{
    return (m_expression);
}


/***********************************************************************//**
 * @brief Signal if filename has an extension name
 *
 * @return True if filename has an extension name, false otherwise.
 ***************************************************************************/
inline
bool GFilename::has_extname(void) const
{
    return (!m_extname.empty());
}


/***********************************************************************//**
 * @brief Signal if filename has an extension number
 *
 * @return True if filename has an extension number, false otherwise.
 ***************************************************************************/
inline
bool GFilename::has_extno(void) const
{
    return (m_extno >= 0);
}


/***********************************************************************//**
 * @brief Signal if filename has an extension version
 *
 * @return True if filename has an extension version, false otherwise.
 ***************************************************************************/
inline
bool GFilename::has_extver(void) const
{
    return (m_extver > 0);
}


/***********************************************************************//**
 * @brief Signal if filename has an expression
 *
 * @return True if filename has an expression, false otherwise.
 ***************************************************************************/
inline
bool GFilename::has_expression(void) const
{
    return (!m_expression.empty());
}


/***********************************************************************//**
 * @brief String addition operator
 *
 * @param[in] filename Filename.
 * @param[in] string String.
 * @return String with filename + string.
 ***************************************************************************/
inline
std::string operator+(const GFilename& filename, const std::string& string)
{
    return (std::string(filename)+string);
}


/***********************************************************************//**
 * @brief String addition operator
 *
 * @param[in] string String.
 * @param[in] filename Filename.
 * @return String with string + filename.
 ***************************************************************************/
inline
std::string operator+(const std::string& string, const GFilename& filename)
{
    return (string+std::string(filename));
}


/***********************************************************************//**
 * @brief Filename equality operator
 *
 * @param[in] a First filename.
 * @param[in] b Second filename.
 * @return True if filenames are equal.
 ***************************************************************************/
inline
bool operator==(const GFilename &a, const GFilename &b)
{
    return (std::string(a) == std::string(b));
}


/***********************************************************************//**
 * @brief Filename inequality operator
 *
 * @param[in] a First filename.
 * @param[in] b Second filename.
 * @return True if filenames are not equal.
 ***************************************************************************/
inline
bool operator!=(const GFilename &a, const GFilename &b)
{
    return (std::string(a) != std::string(b));
}

#endif /* GFILENAME_HPP */
