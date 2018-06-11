/***************************************************************************
 *                GFitsHeaderCard.i - FITS header card class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2018 by Juergen Knoedlseder                         *
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
 * @file GFitsHeaderCard.i
 * @brief FITS header card class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsHeaderCard.hpp"
#include "GTools.hpp"
%}
%include stl.i


/***********************************************************************//**
 * @class GFitsHeaderCard
 *
 * @brief Implements FITS header card class SWIG interface
 *
 * This class implements FITS header card. Each card consists of a
 * keyname (string), a value (string, floating pointer, integer or logical)
 * and a comment (string). COMMENT or HISTORY cards do not have any value.
 ***************************************************************************/
class GFitsHeaderCard : public GBase {

public:
    // Constructors & Destructors
    GFitsHeaderCard(void);
    GFitsHeaderCard(const std::string& keyname, const std::string& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const double& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const int& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const bool& value,
                    const std::string& comment);
    GFitsHeaderCard(const GFitsHeaderCard& card);
    virtual ~GFitsHeaderCard(void);

    // Methods
    void               clear(void);
    GFitsHeaderCard*   clone(void) const;
    std::string        classname(void) const;
    void               keyname(const std::string& keyname);
    const std::string& keyname(void) const;
    void               value(const std::string& value);
    void               value(const bool& value);
    void               value(const float& value);
    void               value(const double& value);
    void               value(const unsigned short& value);
    void               value(const short& value);
    void               value(const unsigned int& value);
    void               value(const int& value);
    void               value(const long& value);
    void               value(const unsigned long& value);
    void               value(const long long& value);
    const std::string& value(void) const;
    const int&         decimals(void) const;
    void               unit(const std::string& unit);
    const std::string& unit(void) const;
    void               comment(const std::string& comment);
    const std::string& comment(void) const;
    std::string        string(void) const;
    double             real(void) const;
    int                integer(void) const;
};


/***********************************************************************//**
 * @brief GFitsHeaderCard class SWIG extension
 ***************************************************************************/
%extend GFitsHeaderCard {
    GFitsHeaderCard copy() {
        return (*self);
    }
}
