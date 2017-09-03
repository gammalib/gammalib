/***************************************************************************
 *              GCOMRoi.hpp - COMPTEL region of interest class             *
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
 * @file GCOMRoi.hpp
 * @brief COMPTEL region of interest class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMROI_HPP
#define GCOMROI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRoi.hpp"


/***********************************************************************//**
 * @class GCOMRoi
 *
 * @brief COMPTEL region of interest class
 *
 * The COMPTEL region of interest class defines the event direction
 * region that is used for unbinned data analysis.
 *
 * @todo Complete the class description.
 ***************************************************************************/
class GCOMRoi : public GRoi {

public:
    // Constructors and destructors
    GCOMRoi(void);
    GCOMRoi(const GCOMRoi& roi);
    virtual ~GCOMRoi(void);

    // Operators
    GCOMRoi& operator=(const GCOMRoi& roi);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMRoi*    clone(void) const;
    virtual std::string classname(void) const;
    virtual bool        contains(const GEvent& event) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    // TODO: Add any further methods that are needed

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMRoi& roi);
    void free_members(void);
    
    // Protected members
    // TODO: Add any data members that are necessary
    // Example:
    double m_radius; //!< Region of interest radius
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMRoi").
 ***************************************************************************/
inline
std::string GCOMRoi::classname(void) const
{
    return ("GCOMRoi");
}

#endif /* GCOMROI_HPP */
