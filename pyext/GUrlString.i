/***************************************************************************
 *                      GUrlString.i - String URL class                    *
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
 * @file GUrlString.i
 * @brief String URL class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GUrlString.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GUrlString
 *
 * @brief String URL class
 ***************************************************************************/
class GUrlString : public GUrl {
public:
    // Constructors and destructors
    GUrlString(void);
    GUrlString(const std::string& string);
    GUrlString(const GUrlString& url);
    virtual ~GUrlString(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GUrlString* clone(void) const;
    virtual void        open(const std::string& string, const std::string& mode = "");
    virtual void        close(void);
    virtual int         getchar(void) const;
    virtual void        putchar(const int& character);

    // Other methods
    const std::string&  string(void) const;
    void                rewind(void);
};


/***********************************************************************//**
 * @brief GUrlString class extension
 ***************************************************************************/
%extend GUrlString {
    GUrlString copy() {
        return (*self);
    }
};
