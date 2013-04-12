/***************************************************************************
 *                     GUrl.i - Abstract URL base class                    *
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
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GUrl.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GUrl
 *
 * @brief Abstract URL base class
 ***************************************************************************/
class GUrl : public GBase {
public:
    // Constructors and destructors
    GUrl(void);
    GUrl(const GUrl& url);
    virtual ~GUrl(void);

    // Pure virtual methods
    virtual void  clear(void) = 0;
    virtual GUrl* clone(void) const = 0;
    virtual void  open(const std::string& url, const std::string& mode) = 0;
    virtual void  close(void) = 0;
    virtual int   getchar(void) = 0;
    virtual void  putchar(const int& character) = 0;
};


/***********************************************************************//**
 * @brief GUrl class extension
 ***************************************************************************/
%extend GUrl {
    char* read(const int& nbyte) {
        char* buffer  = new char[nbyte+1];
        int   nread   = self->read(buffer, nbyte);
        buffer[nread] = '\0';
        return buffer;
    }
    int write(const std::string& buffer, const int& nbyte) {
        return self->write(buffer.c_str(), nbyte);
    }
};
