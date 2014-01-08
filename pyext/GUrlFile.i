/***************************************************************************
 *                        GUrlFile.i - File URL class                      *
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
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GUrlFile.hpp"
#include "GTools.hpp"
%}


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

    // Implemented pure virtual base class methods
    virtual void      clear(void);
    virtual GUrlFile* clone(void) const;
    virtual void      open(const std::string& url, const std::string& mode);
    virtual void      close(void);
    virtual int       get_char(void) const;
    virtual void      put_char(const int& character);
};


/***********************************************************************//**
 * @brief GUrlFile class extension
 ***************************************************************************/
%extend GUrlFile {
    GUrlFile copy() {
        return (*self);
    }
};
