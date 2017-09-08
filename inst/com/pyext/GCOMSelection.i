/***************************************************************************
 *              GCOMSelection.i - COMPTEL selection set class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMSelection.i
 * @brief COMPTEL selection set class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMSelection.hpp"
%}


/***********************************************************************//**
 * @class GCOMSelection
 *
 * @brief COMPTEL selection set class
 ***************************************************************************/
class GCOMSelection : public GBase {

public:
    // Constructors and destructors
    GCOMSelection(void);
    GCOMSelection(const GCOMSelection& select);
    virtual ~GCOMSelection(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCOMSelection* clone(void) const;
    virtual std::string    classname(void) const;

    // Other methods
    void init_statistics(void) const;
    bool use_event(const GCOMEventAtom& event) const;
};


/***********************************************************************//**
 * @brief GCOMSelection class extension
 ***************************************************************************/
%extend GCOMSelection {
    GCOMSelection copy() {
        return (*self);
    }
};
