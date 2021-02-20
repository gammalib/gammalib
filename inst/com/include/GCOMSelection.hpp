/***************************************************************************
 *             GCOMSelection.hpp - COMPTEL selection set class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2021 by Juergen Knoedlseder                         *
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
 * @file GCOMSelection.hpp
 * @brief COMPTEL selection set class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMSELECTION_HPP
#define GCOMSELECTION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"

/* __ Forward declarations _______________________________________________ */
class GFitsHDU;
class GCOMEventAtom;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMSelection
 *
 * @brief COMPTEL selection set class
 *
 * This class implements a COMPTEL selection set.
 ***************************************************************************/
class GCOMSelection : public GBase {

public:
    // Constructors and destructors
    GCOMSelection(void);
    GCOMSelection(const GCOMSelection& select);
    virtual ~GCOMSelection(void);

    // Operators
    GCOMSelection& operator=(const GCOMSelection& select);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCOMSelection* clone(void) const;
    virtual std::string    classname(void) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void init_statistics(void) const;
    bool use_event(const GCOMEventAtom& event) const;
    void read(const GFitsHDU& hdu);
    void write(GFitsHDU& hdu) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMSelection& select);
    void free_members(void);

    // Protected members
    double m_e1_min;       //!< Minimum D1 energy deposit (MeV)
    double m_e1_max;       //!< Maximum D1 energy deposit (MeV)
    double m_e2_min;       //!< Minimum D2 energy deposit (MeV)
    double m_e2_max;       //!< Maximum D2 energy deposit (MeV)
    double m_tof_min;      //!< Minimum TOF window
    double m_tof_max;      //!< Maximum TOF window
    double m_psd_min;      //!< Minimum PSD window
    double m_psd_max;      //!< Maximum PSD window
    double m_zeta_min;     //!< Minimum Earth horizon angle - Phibar window
    double m_zeta_max;     //!< Maximum Earth horizon angle - Phibar window
    int    m_reflag_min;   //!< Minimum rejection flag
    int    m_reflag_max;   //!< Maximum rejection flag
    int    m_vetoflag_min; //!< Minimum veto flag
    int    m_vetoflag_max; //!< Maximum veto flag

    // Selection statistics
    mutable int m_num_events_checked;  //!< Number of checked events
    mutable int m_num_events_used;     //!< Number of used events
    mutable int m_num_events_rejected; //!< Number of rejected events
    mutable int m_num_e1_min;          //!< Number of events below E1 threshold
    mutable int m_num_e1_max;          //!< Number of events above E1 threshold
    mutable int m_num_e2_min;          //!< Number of events below E2 threshold
    mutable int m_num_e2_max;          //!< Number of events above E2 threshold
    mutable int m_num_tof_min;         //!< Number of events below TOF threshold
    mutable int m_num_tof_max;         //!< Number of events above TOF threshold
    mutable int m_num_psd_min;         //!< Number of events below PSD threshold
    mutable int m_num_psd_max;         //!< Number of events above PSD threshold
    mutable int m_num_zeta_min;        //!< Number of events below Zeta threshold
    mutable int m_num_zeta_max;        //!< Number of events above Zeta threshold
    mutable int m_num_reflag_min;      //!< Number of events below rejection flag threshold
    mutable int m_num_reflag_max;      //!< Number of events above rejection flag threshold
    mutable int m_num_vetoflag_min;    //!< Number of events below veto flag threshold
    mutable int m_num_vetoflag_max;    //!< Number of events above veto flag threshold
    mutable int m_num_no_scatter;      //!< Number of events without scatter angle
    mutable int m_num_invalid_modcom;  //!< Number of events with invalid minitelescopes
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMSelection").
 ***************************************************************************/
inline
std::string GCOMSelection::classname(void) const
{
    return ("GCOMSelection");
}

#endif /* GCOMSELECTION_HPP */
