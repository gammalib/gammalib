/***************************************************************************
 *              GLATEventList.hpp - Fermi/LAT event list class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
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
 * @file GLATEventList.hpp
 * @brief Fermi/LAT event list class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATEVENTLIST_HPP
#define GLATEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventList.hpp"
#include "GLATEventAtom.hpp"
#include "GLATRoi.hpp"
#include "GFits.hpp"
#include "GFitsHDU.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GLATEventList
 *
 * @brief Fermi/LAT event list class
 ***************************************************************************/
class GLATEventList : public GEventList {

public:
    // Constructors and destructors
    GLATEventList(void);
    GLATEventList(const GLATEventList& list);
    virtual ~GLATEventList(void);

    // Operators
    virtual GLATEventList&       operator= (const GLATEventList& list);
    virtual GLATEventAtom*       operator[](const int& index);
    virtual const GLATEventAtom* operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GLATEventList* clone(void) const;
    virtual std::string    classname(void) const;
    virtual int            size(void) const;
    virtual void           load(const std::string& filename);
    virtual void           save(const std::string& filename,
                                const bool& clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual void           roi(const GRoi& roi);
    virtual const GLATRoi& roi(void) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GLATEventList& list);
    void         free_members(void);
    virtual void set_energies(void) { return; }
    virtual void set_times(void) { return; }
    void         read_events(const GFitsTable& hdu);
    void         read_ds_keys(const GFitsHDU& hdu);

    // Protected members
    GLATRoi                    m_roi;            //!< Region of interest
    std::vector<GLATEventAtom> m_events;         //!< Events
    std::vector<std::string>   m_difrsp_label;   //!< Diffuse response model labels
    std::vector<std::string>   m_ds_type;        //!< Data selection types
    std::vector<std::string>   m_ds_unit;        //!< Data selection units
    std::vector<std::string>   m_ds_value;       //!< Data selection values
    std::vector<std::string>   m_ds_reference;   //!< Data selection references
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GLATEventList").
 ***************************************************************************/
inline
std::string GLATEventList::classname(void) const
{
    return ("GLATEventList");
}


/***********************************************************************//**
 * @brief Return number of events in list
 *
 * @return Number of events in list.
 ***************************************************************************/
inline
int GLATEventList::size(void) const
{
    return m_events.size();
}


/***********************************************************************//**
 * @brief Return number of events in list
 *
 * @return Number of events in list.
 ***************************************************************************/
inline
int GLATEventList::number(void) const
{
    return m_events.size();
}


/***********************************************************************//**
 * @brief Return Region of Interest
 *
 * @return Region of Interest.
 ***************************************************************************/
inline
const GLATRoi& GLATEventList::roi(void) const
{
    return m_roi;
}

#endif /* GLATEVENTLIST_HPP */
