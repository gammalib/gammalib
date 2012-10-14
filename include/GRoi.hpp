/***************************************************************************
 *             GRoi.hpp  -  Region of interest interface class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GRoi.hpp
 * @brief Region of interest interface class definition
 * @author Juergen Knoedlseder
 */

#ifndef GROI_HPP
#define GROI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GRoi
 *
 * @brief Interface for the region of interest classes
 *
 * The region of interest class holds instrument specific information about
 * the spatial region in detector or telescopes coordinates that is used
 * for an analysis. In particular, the definition of a region of interest
 * is required for an unbinned analysis.
 ***************************************************************************/
class GRoi : public GBase {

    // Friend classes
    friend class GObservation;

public:
    // Constructors and destructors
    GRoi(void);
    GRoi(const GRoi& roi);
    virtual ~GRoi(void);

    // Operators
    virtual GRoi& operator= (const GRoi& roi);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GRoi*       clone(void) const = 0;
    virtual std::string print(void) const = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GRoi& roi);
    void free_members(void);
};

#endif /* GROI_HPP */
