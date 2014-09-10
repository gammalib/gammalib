/***************************************************************************
 *                     GUrlString.hpp - String URL class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Juergen Knoedlseder                         *
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
 * @file GUrlString.hpp
 * @brief String URL class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GURLSTRING_HPP
#define GURLSTRING_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GUrl.hpp"


/***********************************************************************//**
 * @class GUrlString
 *
 * @brief String URL class
 *
 * This class implements the URL interface for a string.
 ***************************************************************************/
class GUrlString : public GUrl {

public:
    // Constructors and destructors
    GUrlString(void);
    GUrlString(const std::string& string);
    GUrlString(const GUrlString& url);
    virtual ~GUrlString(void);

    // Operators
    GUrlString& operator=(const GUrlString& url);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GUrlString* clone(void) const;
    virtual std::string classname(void) const;
    virtual void        open(const std::string& string, const std::string& mode = "");
    virtual void        close(void);
    virtual int         read(void* buffer, const int& nbyte);
    virtual int         write(const void* buffer, const int& nbyte);
    virtual int         get_char(void) const;
    virtual void        put_char(const int& character);
    virtual void        scanf(const char* format, ...);
    virtual void        printf(const char* format, ...);
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const std::string&  string(void) const { return m_buffer; }
    void                rewind(void) { m_index=0; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GUrlString& url);
    void free_members(void);

    // Protected members
    mutable int m_index;   //!< String position indicator
    std::string m_buffer;  //!< Text string
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GUrlString").
 ***************************************************************************/
inline
std::string GUrlString::classname(void) const
{
    return ("GUrlString");
}

#endif /* GURLSTRING_HPP */
