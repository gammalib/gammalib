/***************************************************************************
 *               GCOMEventAtom.hpp - COMPTEL event atom class              *
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
 * @file GCOMEventAtom.hpp
 * @brief COMPTEL event atom class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMEVENTATOM_HPP
#define GCOMEVENTATOM_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GEventAtom.hpp"
#include "GCOMInstDir.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GCOMEventAtom
 *
 * @brief COMPTEL event atom class
 ***************************************************************************/
class GCOMEventAtom : public GEventAtom {

public:
    // Constructors and destructors
    GCOMEventAtom(void);
    GCOMEventAtom(const GCOMEventAtom& atom);
    virtual ~GCOMEventAtom(void);

    // Operators
    GCOMEventAtom& operator=(const GCOMEventAtom& atom);

    // Implemented pure virtual base class methods
    void               clear(void);
    GCOMEventAtom*     clone(void) const;
    std::string        classname(void) const;
    const GCOMInstDir& dir(void) const;
    const GEnergy&     energy(void) const;
    const GTime&       time(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void         dir(const GCOMInstDir& dir);
    void         energy(const GEnergy& energy);
    void         time(const GTime& time);
    void         time(const int& tjd, const int& tics);
    void         phibar(const float& phibar);
    const float& phibar(void) const;
    void         phi(const float& phi);
    const float& phi(void) const;
    void         theta(const float& theta);
    const float& theta(void) const;
    void         eha(const float& eha);
    const float& eha(void) const;
    void         e1(const float& e1);
    const float& e1(void) const;
    void         e2(const float& e2);
    const float& e2(void) const;
    void         psd(const float& psd);
    const float& psd(void) const;
    void         tof(const float& tof);
    const float& tof(void) const;
    void         x_d2(const float& x_d2);
    const float& x_d2(void) const;
    void         y_d2(const float& y_d2);
    const float& y_d2(void) const;
    void         modcom(const int& modcom);
    const int&   modcom(void) const;
    void         reflag(const int& reflag);
    const int&   reflag(void) const;
    void         veto(const int& veto);
    const int&   veto(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMEventAtom& atom);
    void free_members(void);

    // Protected members
    GCOMInstDir m_dir;    //!< Event direction
    GEnergy     m_energy; //!< Event energy
    GTime       m_time;   //!< Event time
    float       m_e1;     //!< D1 energy deposit (MeV)
    float       m_e2;     //!< D2 energy deposit (MeV)
    float       m_phibar; //!< Compton scatter angle (deg)
    float       m_theta;  //!< Zenith angle of scatter direction (deg)
    float       m_phi;    //!< Azimuth angle of scatter direction (deg)
    float       m_eha;    //!< Earth horizon angle (deg)
    float       m_psd;    //!< PSD value (channel)
    float       m_tof;    //!< Time of flight value (channel)
    float       m_x_d2;   //!< D2 model X position (mm)
    float       m_y_d2;   //!< D2 model X position (mm)
    int         m_modcom; //!< Mini telescope number
    int         m_reflag; //!< Rejection flag
    int         m_veto;   //!< Veto flag
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMEventAtom").
 ***************************************************************************/
inline
std::string GCOMEventAtom::classname(void) const
{
    return ("GCOMEventAtom");
}


/***********************************************************************//**
 * @brief Return event instrument direction
 *
 * @return Event instrument direction.
 *
 * Returns the direction of the event.
 ***************************************************************************/
inline
const GCOMInstDir& GCOMEventAtom::dir(void) const
{
    return m_dir;
}


/***********************************************************************//**
 * @brief Set event instrument direction
 *
 * @param[in] dir Event instrument direction.
 *
 * Sets the event instrument direction.
 ***************************************************************************/
inline
void GCOMEventAtom::dir(const GCOMInstDir& dir)
{
    m_dir = dir;
    return;
}


/***********************************************************************//**
 * @brief Return event energy
 *
 * @return Event energy.
 *
 * Returns the reconstructed energy of the photon on the sky.
 ***************************************************************************/
inline
const GEnergy& GCOMEventAtom::energy(void) const
{
    return m_energy;
}


/***********************************************************************//**
 * @brief Set event energy
 *
 * @param[in] energy Event energy.
 *
 * Sets the event energy.
 ***************************************************************************/
inline
void GCOMEventAtom::energy(const GEnergy& energy)
{
    m_energy = energy;
    return;
}


/***********************************************************************//**
 * @brief Return event time
 *
 * @return Event time.
 *
 * Returns the event triggering time.
 ***************************************************************************/
inline
const GTime& GCOMEventAtom::time(void) const
{
    return m_time;
}


/***********************************************************************//**
 * @brief Set event time
 *
 * @param[in] time Event time.
 *
 * Sets the event time.
 ***************************************************************************/
inline
void GCOMEventAtom::time(const GTime& time)
{
    m_time = time;
    return;
}


/***********************************************************************//**
 * @brief Set Compton scatter angle
 *
 * @param[in] phibar Compton scatter angle (deg).
 *
 * Sets Compton scatter angle.
 ***************************************************************************/
inline
void GCOMEventAtom::phibar(const float& phibar)
{
    m_phibar = phibar;
    return;
}


/***********************************************************************//**
 * @brief Return Compton scatter angle
 *
 * @return Compton scatter angle (deg).
 *
 * Returns Compton scatter angle.
 ***************************************************************************/
inline
const float& GCOMEventAtom::phibar(void) const
{
    return m_phibar;
}


/***********************************************************************//**
 * @brief Set scatter azimuth angle
 *
 * @param[in] phi Scatter azimuth angle (deg).
 *
 * Sets scatter azimuth angle.
 ***************************************************************************/
inline
void GCOMEventAtom::phi(const float& phi)
{
    m_phi = phi;
    return;
}


/***********************************************************************//**
 * @brief Return scatter azimuth angle
 *
 * @return Scatter azimuth angle (deg).
 *
 * Returns scatter azimuth angle.
 ***************************************************************************/
inline
const float& GCOMEventAtom::phi(void) const
{
    return m_phi;
}


/***********************************************************************//**
 * @brief Set scatter zenith angle
 *
 * @param[in] theta Scatter zenith angle (deg).
 *
 * Sets scatter zenith angle.
 ***************************************************************************/
inline
void GCOMEventAtom::theta(const float& theta)
{
    m_theta = theta;
    return;
}


/***********************************************************************//**
 * @brief Return scatter zenith angle
 *
 * @return Scatter zenith angle (deg).
 *
 * Returns scatter zenith angle.
 ***************************************************************************/
inline
const float& GCOMEventAtom::theta(void) const
{
    return m_theta;
}


/***********************************************************************//**
 * @brief Set Earth horizon angle
 *
 * @param[in] eha Earth horizon angle (deg).
 *
 * Sets Earth horizon angle.
 ***************************************************************************/
inline
void GCOMEventAtom::eha(const float& eha)
{
    m_eha = eha;
    return;
}


/***********************************************************************//**
 * @brief Return Earth horizon angle
 *
 * @return Earth horizon angle (deg).
 *
 * Returns Earth horizon angle.
 ***************************************************************************/
inline
const float& GCOMEventAtom::eha(void) const
{
    return m_eha;
}


/***********************************************************************//**
 * @brief Set D1 module energy deposit
 *
 * @param[in] e1 D1 module energy deposit (MeV).
 *
 * Sets D1 module energy deposit.
 ***************************************************************************/
inline
void GCOMEventAtom::e1(const float& e1)
{
    m_e1 = e1;
    return;
}


/***********************************************************************//**
 * @brief Return D1 module energy deposit
 *
 * @return D1 module energy deposit (MeV).
 *
 * Returns D1 module energy deposit.
 ***************************************************************************/
inline
const float& GCOMEventAtom::e1(void) const
{
    return m_e1;
}


/***********************************************************************//**
 * @brief Set D2 module energy deposit
 *
 * @param[in] e2 D2 module energy deposit (MeV).
 *
 * Sets D2 module energy deposit.
 ***************************************************************************/
inline
void GCOMEventAtom::e2(const float& e2)
{
    m_e2 = e2;
    return;
}


/***********************************************************************//**
 * @brief Return D2 module energy deposit
 *
 * @return D2 module energy deposit (MeV).
 *
 * Returns D2 module energy deposit.
 ***************************************************************************/
inline
const float& GCOMEventAtom::e2(void) const
{
    return m_e2;
}


/***********************************************************************//**
 * @brief Set PSD value
 *
 * @param[in] psd PSD value (channel).
 *
 * Sets PSD value.
 ***************************************************************************/
inline
void GCOMEventAtom::psd(const float& psd)
{
    m_psd = psd;
    return;
}

/***********************************************************************//**
 * @brief Return PSD value
 *
 * @return PSD value (channel).
 *
 * Returns the Pulse Shape Discriminator (PSD) channel value of the event.
 * The PSD value is used for the distinction between gammas and neutrons.
 ***************************************************************************/
inline
const float& GCOMEventAtom::psd(void) const
{
    return m_psd;
}

/***********************************************************************//**
 * @brief Set TOF value
 *
 * @param[in] tof TOF value (channel).
 *
 * Sets TOF value.
 ***************************************************************************/
inline
void GCOMEventAtom::tof(const float& tof)
{
    m_tof = tof;
    return;
}


/***********************************************************************//**
 * @brief Return TOF value
 *
 * @return TOF value (channel).
 *
 * Returns the Time Of Flight (TOF) channel value of the event. The TOF
 * value is used for the distinction between forward and backward scattering
 * events.
 ***************************************************************************/
inline
const float& GCOMEventAtom::tof(void) const
{
    return m_tof;
}


/***********************************************************************//**
 * @brief Set D2 module X value
 *
 * @param[in] x_d2 D2 module X value (mm).
 *
 * Sets D2 module X value.
 ***************************************************************************/
inline
void GCOMEventAtom::x_d2(const float& x_d2)
{
    m_x_d2 = x_d2;
    return;
}


/***********************************************************************//**
 * @brief Return D2 module X value
 *
 * @return D2 module X value (mm).
 *
 * Returns the D2 module X value of the event.
 ***************************************************************************/
inline
const float& GCOMEventAtom::x_d2(void) const
{
    return m_x_d2;
}


/***********************************************************************//**
 * @brief Set D2 module Y value
 *
 * @param[in] y_d2 D2 module Y value (mm).
 *
 * Sets D2 module Y value.
 ***************************************************************************/
inline
void GCOMEventAtom::y_d2(const float& y_d2)
{
    m_y_d2 = y_d2;
    return;
}


/***********************************************************************//**
 * @brief Return D2 module Y value
 *
 * @return D2 module Y value (mm).
 *
 * Returns the D2 module Y value of the event.
 ***************************************************************************/
inline
const float& GCOMEventAtom::y_d2(void) const
{
    return m_y_d2;
}


/***********************************************************************//**
 * @brief Set mini telescope
 *
 * @param[in] modcom Mini telescope.
 *
 * Sets the mini telescope of the event.
 ***************************************************************************/
inline
void GCOMEventAtom::modcom(const int& modcom)
{
    m_modcom = modcom;
    return;
}


/***********************************************************************//**
 * @brief Return mini telescope
 *
 * @return Mini telescope.
 *
 * Returns the mini telescope of the event.
 ***************************************************************************/
inline
const int& GCOMEventAtom::modcom(void) const
{
    return m_modcom;
}


/***********************************************************************//**
 * @brief Set rejection flag
 *
 * @param[in] reflag Rejection flag.
 *
 * Sets the rejection flag of the event.
 ***************************************************************************/
inline
void GCOMEventAtom::reflag(const int& reflag)
{
    m_reflag = reflag;
    return;
}


/***********************************************************************//**
 * @brief Return rejection flag
 *
 * @return Rejection flag.
 *
 * Returns the rejection flag of the event.
 ***************************************************************************/
inline
const int& GCOMEventAtom::reflag(void) const
{
    return m_reflag;
}


/***********************************************************************//**
 * @brief Set veto flag
 *
 * @param[in] veto Veto flag.
 *
 * Sets the veto flag of the event.
 ***************************************************************************/
inline
void GCOMEventAtom::veto(const int& veto)
{
    m_veto = veto;
    return;
}


/***********************************************************************//**
 * @brief Return veto flag
 *
 * @return Veto flag.
 *
 * Returns the veto flag of the event.
 ***************************************************************************/
inline
const int& GCOMEventAtom::veto(void) const
{
    return m_veto;
}

#endif /* GCOMEVENTATOM_HPP */
