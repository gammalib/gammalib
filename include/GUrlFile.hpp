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
 * @file GUrlFile.hpp
 * @brief File URL class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GURLFILE_HPP
#define GURLFILE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <cstdio>     // FILE
#include "GUrl.hpp"


/***********************************************************************//**
 * @class GUrlFile
 *
 * @brief File URL class
 ***************************************************************************/
class GUrlFile : public GUrl {

public:
    // Constructors and destructors
    GUrlFile(void);
    GUrlFile(const std::string& url, const std::string& mode);
    GUrlFile(const GUrlFile& url);
    virtual ~GUrlFile(void);

    // Operators
    GUrlFile& operator=(const GUrlFile& url);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GUrlFile*   clone(void) const;
    virtual void        open(const std::string& url, const std::string& mode);
    virtual void        close(void);
    virtual int         read(void* buffer, const int& nbyte);
    virtual int         write(const void* buffer, const int& nbyte);
    virtual int         getchar(void);
    virtual void        putchar(const int& character);
    virtual void        scanf(const char* format, ...);
    virtual void        printf(const char* format, ...);
    virtual std::string print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GUrlFile& url);
    void free_members(void);

    // Protected members
    std::string m_url;   //!< File URL
    std::string m_mode;  //!< File mode
    FILE*       m_fptr;  //!< File pointer
};

#endif /* GURLFILE_HPP */
