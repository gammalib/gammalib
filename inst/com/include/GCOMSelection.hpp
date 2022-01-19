/***************************************************************************
 *             GCOMSelection.hpp - COMPTEL selection set class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2022 by Juergen Knoedlseder                         *
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
#include "GPhases.hpp"
#include "GPulsar.hpp"
#include "GModelTemporalPhaseCurve.hpp"

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
    void           init_statistics(void) const;
    bool           use_event(const GCOMEventAtom& event) const;
    const double&  e1_min(void) const;
    void           e1_min(const double& e1_min);
    const double&  e1_max(void) const;
    void           e1_max(const double& e1_max);
    const double&  e2_min(void) const;
    void           e2_min(const double& e2_min);
    const double&  e2_max(void) const;
    void           e2_max(const double& e2_max);
    const int&     tof_min(void) const;
    void           tof_min(const int& tof_min);
    const int&     tof_max(void) const;
    void           tof_max(const int& tof_max);
    const int&     psd_min(void) const;
    void           psd_min(const int& psd_min);
    const int&     psd_max(void) const;
    void           psd_max(const int& psd_max);
    const int&     reflag_min(void) const;
    void           reflag_min(const int& reflag_min);
    const int&     reflag_max(void) const;
    void           reflag_max(const int& reflag_max);
    const int&     vetoflag_min(void) const;
    void           vetoflag_min(const int& vetoflag_min);
    const int&     vetoflag_max(void) const;
    void           vetoflag_max(const int& vetoflag_max);
    const int&     fpmtflag(void) const;
    void           fpmtflag(const int& fpmtflag);
    const bool&    use_d1(const int& id1) const;
    void           use_d1(const int& id1, const bool& use);
    const bool&    use_d2(const int& id2) const;
    void           use_d2(const int& id2, const bool& use);
    const GPhases& orbital_phases(void) const;
    void           orbital_phases(const GPhases& phases);
    void           orbital_period(const double& period, const GTime& time);
    double         orbital_phase(const GTime& time) const;
    const GPhases& pulsar_phases(void) const;
    void           pulsar_phases(const GPhases& phases);
    const GPulsar& pulsar(void) const;
    void           pulsar(const GPulsar& pulsar);
    bool           has_pulsar(void) const;
    void           read(const GFitsHDU& hdu);
    void           write(GFitsHDU& hdu) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMSelection& select);
    void free_members(void);

    // Protected members
    double                   m_e1_min;              //!< Minimum D1 energy deposit (MeV)
    double                   m_e1_max;              //!< Maximum D1 energy deposit (MeV)
    double                   m_e2_min;              //!< Minimum D2 energy deposit (MeV)
    double                   m_e2_max;              //!< Maximum D2 energy deposit (MeV)
    int                      m_tof_min;             //!< Minimum TOF window
    int                      m_tof_max;             //!< Maximum TOF window
    int                      m_psd_min;             //!< Minimum PSD window
    int                      m_psd_max;             //!< Maximum PSD window
    int                      m_reflag_min;          //!< Minimum rejection flag
    int                      m_reflag_max;          //!< Maximum rejection flag
    int                      m_vetoflag_min;        //!< Minimum veto flag
    int                      m_vetoflag_max;        //!< Maximum veto flag
    int                      m_fpmtflag;            //!< D2 PMT failures
                                                    //!< (0: exclude modules, 1: include modules,
                                                    //!<  2: use FPM information)
    bool                     m_use_d1[7];           //!< D1 module usage
    bool                     m_use_d2[14];          //!< D2 module usage
    GPhases                  m_orbital_phases;      //!< Phases for orbital phase selection
    GModelTemporalPhaseCurve m_orbital_phase_curve; //!< Orbital phase curve
    GPhases                  m_pulsar_phases;       //!< Phases for pulse phase selection
    GPulsar                  m_pulsar;              //!< Pulsar information

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
    mutable int m_num_reflag_min;      //!< Number of events below rejection flag threshold
    mutable int m_num_reflag_max;      //!< Number of events above rejection flag threshold
    mutable int m_num_vetoflag_min;    //!< Number of events below veto flag threshold
    mutable int m_num_vetoflag_max;    //!< Number of events above veto flag threshold
    mutable int m_num_no_scatter;      //!< Number of events without scatter angle
    mutable int m_num_invalid_modcom;  //!< Number of events with invalid minitelescopes
    mutable int m_num_d1module_off;    //!< Number of events excluded since D1 module off
    mutable int m_num_d2module_off;    //!< Number of events excluded since D2 module off
    mutable int m_num_fpmt;            //!< Number of events excluded due to failed PMT
    mutable int m_num_d1[7];           //!< Number of events per D1 module
    mutable int m_num_d2[14];          //!< Number of events per D2 module
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


/***********************************************************************//**
 * @brief Return minimum of D1 energy selection window
 *
 * @return Minimum of D1 energy selection window.
 ***************************************************************************/
inline
const double& GCOMSelection::e1_min(void) const
{
    return (m_e1_min);
}


/***********************************************************************//**
 * @brief Set minimum of D1 energy selection window
 *
 * @param[in] e1_min Minimum of D1 energy selection window.
 ***************************************************************************/
inline
void GCOMSelection::e1_min(const double& e1_min)
{
    m_e1_min = e1_min;
    return;
}


/***********************************************************************//**
 * @brief Return maximum of D1 energy selection window
 *
 * @return Maximum of D1 energy selection window.
 ***************************************************************************/
inline
const double& GCOMSelection::e1_max(void) const
{
    return (m_e1_max);
}


/***********************************************************************//**
 * @brief Set maximum of D1 energy selection window
 *
 * @param[in] e1_max Maximum of D1 energy selection window.
 ***************************************************************************/
inline
void GCOMSelection::e1_max(const double& e1_max)
{
    m_e1_max = e1_max;
    return;
}


/***********************************************************************//**
 * @brief Return minimum of D2 energy selection window
 *
 * @return Minimum of D2 energy selection window.
 ***************************************************************************/
inline
const double& GCOMSelection::e2_min(void) const
{
    return (m_e2_min);
}


/***********************************************************************//**
 * @brief Set minimum of D2 energy selection window
 *
 * @param[in] e2_min Minimum of D2 energy selection window.
 ***************************************************************************/
inline
void GCOMSelection::e2_min(const double& e2_min)
{
    m_e2_min = e2_min;
    return;
}


/***********************************************************************//**
 * @brief Return maximum of D2 energy selection window
 *
 * @return Maximum of D2 energy selection window.
 ***************************************************************************/
inline
const double& GCOMSelection::e2_max(void) const
{
    return (m_e2_max);
}


/***********************************************************************//**
 * @brief Set maximum of D2 energy selection window
 *
 * @param[in] e2_max Maximum of D2 energy selection window.
 ***************************************************************************/
inline
void GCOMSelection::e2_max(const double& e2_max)
{
    m_e2_max = e2_max;
    return;
}


/***********************************************************************//**
 * @brief Return minimum of ToF selection window
 *
 * @return Minimum of ToF selection window.
 ***************************************************************************/
inline
const int& GCOMSelection::tof_min(void) const
{
    return (m_tof_min);
}


/***********************************************************************//**
 * @brief Set minimum of ToF selection window
 *
 * @param[in] tof_min Minimum of ToF selection window.
 ***************************************************************************/
inline
void GCOMSelection::tof_min(const int& tof_min)
{
    m_tof_min = tof_min;
    return;
}


/***********************************************************************//**
 * @brief Return maximum of ToF selection window
 *
 * @return Maximum of ToF selection window.
 ***************************************************************************/
inline
const int& GCOMSelection::tof_max(void) const
{
    return (m_tof_max);
}


/***********************************************************************//**
 * @brief Set maximum of ToF selection window
 *
 * @param[in] tof_max Maximum of ToF selection window.
 ***************************************************************************/
inline
void GCOMSelection::tof_max(const int& tof_max)
{
    m_tof_max = tof_max;
    return;
}


/***********************************************************************//**
 * @brief Return minimum of PSD selection window
 *
 * @return Minimum of PSD selection window.
 ***************************************************************************/
inline
const int& GCOMSelection::psd_min(void) const
{
    return (m_psd_min);
}


/***********************************************************************//**
 * @brief Set minimum of PSD selection window
 *
 * @param[in] psd_min Minimum of PSD selection window.
 ***************************************************************************/
inline
void GCOMSelection::psd_min(const int& psd_min)
{
    m_psd_min = psd_min;
    return;
}


/***********************************************************************//**
 * @brief Return maximum of PSD selection window
 *
 * @return Maximum of PSD selection window.
 ***************************************************************************/
inline
const int& GCOMSelection::psd_max(void) const
{
    return (m_psd_max);
}


/***********************************************************************//**
 * @brief Set maximum of PSD selection window
 *
 * @param[in] psd_max Maximum of PSD selection window.
 ***************************************************************************/
inline
void GCOMSelection::psd_max(const int& psd_max)
{
    m_psd_max = psd_max;
    return;
}


/***********************************************************************//**
 * @brief Return minimum of Rejection Flag selection window
 *
 * @return Minimum of Rejection Flag selection window.
 ***************************************************************************/
inline
const int& GCOMSelection::reflag_min(void) const
{
    return (m_reflag_min);
}


/***********************************************************************//**
 * @brief Set minimum of Rejection Flag selection window
 *
 * @param[in] reflag_min Minimum of Rejection Flag selection window.
 ***************************************************************************/
inline
void GCOMSelection::reflag_min(const int& reflag_min)
{
    m_reflag_min = reflag_min;
    return;
}


/***********************************************************************//**
 * @brief Return maximum of Rejection Flag selection window
 *
 * @return Maximum of Rejection Flag selection window.
 ***************************************************************************/
inline
const int& GCOMSelection::reflag_max(void) const
{
    return (m_reflag_max);
}


/***********************************************************************//**
 * @brief Set maximum of Rejection Flag selection window
 *
 * @param[in] reflag_max Maximum of Rejection Flag selection window.
 ***************************************************************************/
inline
void GCOMSelection::reflag_max(const int& reflag_max)
{
    m_reflag_max = reflag_max;
    return;
}


/***********************************************************************//**
 * @brief Return minimum of Veto Flag selection window
 *
 * @return Minimum of Veto Flag selection window.
 ***************************************************************************/
inline
const int& GCOMSelection::vetoflag_min(void) const
{
    return (m_vetoflag_min);
}


/***********************************************************************//**
 * @brief Set minimum of Veto Flag selection window
 *
 * @param[in] vetoflag_min Minimum of Veto Flag selection window.
 ***************************************************************************/
inline
void GCOMSelection::vetoflag_min(const int& vetoflag_min)
{
    m_vetoflag_min = vetoflag_min;
    return;
}


/***********************************************************************//**
 * @brief Return maximum of Veto Flag selection window
 *
 * @return Maximum of Veto Flag selection window.
 ***************************************************************************/
inline
const int& GCOMSelection::vetoflag_max(void) const
{
    return (m_vetoflag_max);
}


/***********************************************************************//**
 * @brief Set maximum of Veto Flag selection window
 *
 * @param[in] vetoflag_max Maximum of Veto Flag selection window.
 ***************************************************************************/
inline
void GCOMSelection::vetoflag_max(const int& vetoflag_max)
{
    m_vetoflag_max = vetoflag_max;
    return;
}


/***********************************************************************//**
 * @brief Return failed PMT flag for D2 modules
 *
 * @return Failed PMT flag for D2 modules.
 *
 * Returns the failed PMT flag for D2 modules. The following values can be
 * returned:
 * - 0: excluded D2 modules with failed PMTs
 * - 1: include D2 modules with failed PMTs
 * - 2: include D2 modules with failed PMTs by excluding zone around failed
 *      PMTs
 ***************************************************************************/
inline
const int& GCOMSelection::fpmtflag(void) const
{
    return (m_fpmtflag);
}


/***********************************************************************//**
 * @brief Return orbital phases
 *
 * @return Orbital phases
 ***************************************************************************/
inline
const GPhases& GCOMSelection::orbital_phases(void) const
{
    return (m_orbital_phases);
}


/***********************************************************************//**
 * @brief Set orbital phases
 *
 * @param[in] phases Orbital phases.
 ***************************************************************************/
inline
void GCOMSelection::orbital_phases(const GPhases& phases)
{
    m_orbital_phases = phases;
    return;
}


/***********************************************************************//**
 * @brief Return orbital phase for a given time
 *
 * @param[in] time Time.
 * @return Orbital phase.
 ***************************************************************************/
inline
double GCOMSelection::orbital_phase(const GTime& time) const
{
    // Compute phase
    double phase = m_orbital_phase_curve.phase(time);

    // Return phase
    return phase;
}


/***********************************************************************//**
 * @brief Return pulsar phases
 *
 * @return Pulsar phases
 ***************************************************************************/
inline
const GPhases& GCOMSelection::pulsar_phases(void) const
{
    return (m_pulsar_phases);
}


/***********************************************************************//**
 * @brief Set pulsar phases
 *
 * @param[in] phases Pulsar phases.
 ***************************************************************************/
inline
void GCOMSelection::pulsar_phases(const GPhases& phases)
{
    m_pulsar_phases = phases;
    return;
}


/***********************************************************************//**
 * @brief Return pulsar
 *
 * @return Pulsar
 ***************************************************************************/
inline
const GPulsar& GCOMSelection::pulsar(void) const
{
    return (m_pulsar);
}


/***********************************************************************//**
 * @brief Set pulsar
 *
 * @param[in] pulsar Pulsar.
 ***************************************************************************/
inline
void GCOMSelection::pulsar(const GPulsar& pulsar)
{
    m_pulsar = pulsar;
    return;
}


/***********************************************************************//**
 * @brief Signals that pulsar selection should be performed
 *
 * @return True if pulsar selection should be performed, false otherwise
 *
 * This method returns true if the pulsar and the pulsar phase selection are
 * not empty.
 ***************************************************************************/
inline
bool GCOMSelection::has_pulsar(void) const
{
    return (!pulsar().is_empty() && !pulsar_phases().is_empty());
}

#endif /* GCOMSELECTION_HPP */
