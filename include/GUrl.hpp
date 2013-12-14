/***************************************************************************
 *                    GUrl.hpp - Abstract URL base class                   *
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
 * @file GUrl.hpp
 * @brief Abstract URL base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GURL_HPP
#define GURL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GUrl
 *
 * @brief Abstract URL base class
 *
 * This class defines the abstract interface for all URL classes. The URL
 * classes implement handling of various URL types through a standard
 * interface. This allows to develop URL independent code.
 ***************************************************************************/
class GUrl : public GBase {

public:
    // Constructors and destructors
    GUrl(void);
    GUrl(const GUrl& url);
    virtual ~GUrl(void);

    // Operators
    GUrl& operator=(const GUrl& url);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GUrl*       clone(void) const = 0;
    virtual void        open(const std::string& url, const std::string& mode) = 0;
    virtual void        close(void) = 0;
    virtual int         read(void* buffer, const int& nbyte) = 0;
    virtual int         write(const void* buffer, const int& nbyte) = 0;
    virtual int         getchar(void) const = 0;
    virtual void        putchar(const int& character) = 0;
    virtual void        scanf(const char* format, ...) = 0;
    virtual void        printf(const char* format, ...) = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GUrl& url);
    void free_members(void);

};

#endif /* GURL_HPP */
