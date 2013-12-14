/***************************************************************************
 *              GCTAInstDir.i - CTA instrument direction class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAInstDir.i
 * @brief CTA instrument direction class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAInstDir.hpp"
%}


/***********************************************************************//**
 * @class GCTAInstDir
 *
 * @brief CTA instrument direction class
 ***************************************************************************/
class GCTAInstDir : public GInstDir {

public:
    // Constructors and destructors
    GCTAInstDir(void);
    explicit GCTAInstDir(const GSkyDir& dir);
    GCTAInstDir(const GCTAInstDir& dir);
    virtual ~GCTAInstDir(void);

    // Implemented pure virtual base class methods
    void         clear(void);
    GCTAInstDir* clone(void) const;

    // Other methods
    void     dir(const GSkyDir& dir);
    GSkyDir& dir(void);
};


/***********************************************************************//**
 * @brief GCTAInstDir class extension
 ***************************************************************************/
%extend GCTAInstDir {
    GCTAInstDir copy() {
        return (*self);
    }
    GCTAInstDir(GInstDir* dir) {
        GCTAInstDir* ptr = dynamic_cast<GCTAInstDir*>(dir);
        if (ptr != NULL) {
            return (ptr->clone());
        }
        else {
            throw GException::bad_type("GCTAInstDir(GInstDir*)",
                                       "GInstDir not of type GCTAInstDir");
        }
    }
};
