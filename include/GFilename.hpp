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

public:
    // Constructors and destructors
    GFilename(void);
    GFilename(const std::string& filename);
    GFilename(const char* filename);
    GFilename(const GFilename& filename);
    virtual ~GFilename(void);

    // Operators
    GFilename&         operator=(const GFilename& filename);
    const std::string& operator()(void) const;

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
    std::string m_fullname;   //!< Full file name
    std::string m_filename;   //!< File name (with stripped extension info)
    //std::string m_protocol;   //!< Access protocol
    //std::string m_url;        //!< File name (with stripped extension info)
    //std::string m_path;       //!< Path
    std::string m_extname;    //!< Extension name ("": not set)
    int         m_extno;      //!< Extension number  (-1: not set)
    int         m_extver;     //!< Extension version (0: not set)
    std::string m_expression; //!< Selection expression ("": not set)
};


/***********************************************************************//**
 * @brief Filename operator
 *
 * @return Full filename.
 *
 * Returns the full filename including any FITS extension.
 ***************************************************************************/
inline
const std::string& GFilename::operator()(void) const
{
    return (m_fullname);
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
    return (m_fullname.empty());
}


/***********************************************************************//**
 * @brief Return length of filename
 *
 * @return Length of filename.
 ***************************************************************************/
inline
int GFilename::length(void) const
{
    return (m_filename.length());
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
    return (gammalib::expand_env(m_filename));
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

#endif /* GFILENAME_HPP */
