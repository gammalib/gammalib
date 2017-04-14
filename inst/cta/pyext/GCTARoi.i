/***************************************************************************
 *                 GCTARoi.i - CTA region of interest class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2017 by Jurgen Knodlseder                           *
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
 * @file GCTARoi.i
 * @brief CTA region of interest class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTARoi.hpp"
%}


/***********************************************************************//**
 * @class GCTARoi
 *
 * @brief CTA region of interest class
 ***************************************************************************/
class GCTARoi : public GRoi {

public:
    // Constructors and destructors
    GCTARoi(void);
    explicit GCTARoi(const GXmlElement& xml);
    GCTARoi(const GCTAInstDir& centre, const double& radius);
    GCTARoi(const GCTARoi& roi);
    virtual ~GCTARoi(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCTARoi*    clone(void) const;
    virtual std::string classname(void) const;
    virtual bool        contains(const GEvent& event) const;

    // Other methods
    bool                is_valid(void) const;
    bool                contains(const GCTAInstDir& dir) const;
    const GCTAInstDir&  centre(void) const;
    const double&       radius(void) const;
    void                centre(const GCTAInstDir& centre);
    void                radius(const double& radius);
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GCTARoi class extension
 ***************************************************************************/
%extend GCTARoi {
    GCTARoi copy() {
        return (*self);
    }
};
