/***************************************************************************
 *             GXXXRoi.hpp - [INSTRUMENT] region of interest class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file GXXXRoi.hpp
 * @brief [INSTRUMENT] region of interest class definition
 * @author [AUTHOR]
 */

#ifndef HEADER_PROTECT
#define HEADER_PROTECT

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRoi.hpp"


/***********************************************************************//**
 * @class GXXXRoi
 *
 * @brief [INSTRUMENT] region of interest class
 *
 * The [INSTRUMENT] region of interest class defines the event direction
 * region that is used for unbinned data analysis.
 *
 * @todo Complete the class description.
 ***************************************************************************/
class GXXXRoi : public GRoi {

public:
    // Constructors and destructors
    GXXXRoi(void);
    GXXXRoi(const GXXXRoi& roi);
    virtual ~GXXXRoi(void);

    // Operators
    GXXXRoi& operator=(const GXXXRoi& roi);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GXXXRoi*    clone(void) const;
    virtual std::string classname(void) const;
    virtual bool        contains(const GEvent& event) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    // TODO: Add any further methods that are needed

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXXXRoi& roi);
    void free_members(void);
    
    // Protected members
    // TODO: Add any data members that are necessary
    // Example:
    double m_radius; //!< Region of interest radius
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GXXXRoi").
 ***************************************************************************/
inline
std::string GXXXRoi::classname(void) const
{
    return ("GXXXRoi");
}

#endif /* HEADER_PROTECT */
