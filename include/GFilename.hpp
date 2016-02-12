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


/***********************************************************************//**
 * @class GFilename
 *
 * @brief Filename class
 *
 * This class handles filenames.
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
    GFilename& operator=(const GFilename& filename);

    // Methods
    void               clear(void);
    GFilename*         clone(void) const;
    std::string        classname(void) const;
    bool               empty(void) const;
    int                size(void) const;
    int                length(void) const;
    const std::string& fullname(void) const;
    const std::string& filename(void) const;
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
    std::string m_extname;    //!< Extension name ("": not set)
    int         m_extno;      //!< Extension number  (-1: not set)
    int         m_extver;     //!< Extension version (0: not set)
    std::string m_expression; //!< Selection expression ("": not set)
};


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
bool GFilename::empty(void) const
{
    return (m_filename.empty());
}


/***********************************************************************//**
 * @brief Return size of filename
 *
 * @return Size of filename.
 *
 * The size of the file name is equal to its length.
 ***************************************************************************/
inline
int GFilename::size(void) const
{
    return (m_filename.size());
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
 * @brief Return full filename
 *
 * @return String containing filename with extension.
 *
 * Returns the full filename, including any extension.
 ***************************************************************************/
inline
const std::string& GFilename::fullname(void) const
{
    return (m_fullname);
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return String containing filename without extension.
 *
 * Returns the file name without any extension.
 ***************************************************************************/
inline
const std::string& GFilename::filename(void) const
{
    return (m_filename);
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
