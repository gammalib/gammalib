/***************************************************************************
 *          GMWLDatum.hpp - Multi-wavelength spectral point class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GMWLDatum.hpp
 * @brief Multi-wavelength spectral point class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMWLDATUM_HPP
#define GMWLDATUM_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventBin.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GMWLInstDir.hpp"


/***********************************************************************//**
 * @class GMWLDatum
 *
 * @brief Multi-wavelength spectral point class
 *
 * This class defines a spectral data point for the multi-wavelength
 * interface. It derives from the abstract GEventBin base class.
 ***************************************************************************/
class GMWLDatum : public GEventBin {

    // Friend classes
    friend class GMWLSpectrum;   //!< Needs access to load data

public:
    // Constructors and destructors
    GMWLDatum(void);
    GMWLDatum(const GMWLDatum& datum);
    virtual ~GMWLDatum(void);

    // Operators
    virtual GMWLDatum& operator= (const GMWLDatum& datum);

    // Implemented pure virtual base class methods
    virtual void               clear(void);
    virtual GMWLDatum*         clone(void) const;
    virtual std::string        classname(void) const;
    virtual double             size(void) const;
    virtual const GMWLInstDir& dir(void) const;
    virtual const GEnergy&     energy(void) const;
    virtual const GTime&       time(void) const;
    virtual double             counts(void) const;
    virtual double             error(void) const;
    virtual void               counts(const double& flux);
    virtual std::string        print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const GEnergy& energy_err(void) const;
    const double&  flux(void) const;
    const double&  flux_err(void) const;
    void           flux(const double& flux);
    void           flux_err(const double& error);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLDatum& datum);
    void free_members(void);

    // Protected members
    GMWLInstDir m_dir;       //!< Instrument direction of spectral point (not used)
    GTime       m_time;      //!< Time of spectral point (not used)
    GEnergy     m_eng;       //!< Energy of spectral point
    GEnergy     m_eng_err;   //!< Uncertainty in energy
    double      m_flux;      //!< Flux of spectral point (ph/cm2/s/MeV)
    double      m_flux_err;  //!< Uncertainty in flux (ph/cm2/s/MeV)

};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GMWLDatum").
 ***************************************************************************/
inline
std::string GMWLDatum::classname(void) const
{
    return ("GMWLDatum");
}


/***********************************************************************//**
 * @brief Return size of spectral bins
 *
 * @return Size of spectrl bin (always 1).
 ***************************************************************************/
inline
double GMWLDatum::size(void) const
{
    return 1.0;
}


/***********************************************************************//**
 * @brief Return instrument direction (dummy method)
 *
 * @return Instrument direction.
 ***************************************************************************/
inline
const GMWLInstDir& GMWLDatum::dir(void) const
{
    return m_dir;
}


/***********************************************************************//**
 * @brief Return energy of spectral bin
 *
 * @return Energy of spectral bin.
 ***************************************************************************/
inline
const GEnergy& GMWLDatum::energy(void) const
{
    return m_eng;
}


/***********************************************************************//**
 * @brief Return time of spectral bin
 *
 * @return Time of spectral bin.
 ***************************************************************************/
inline
const GTime& GMWLDatum::time(void) const
{
    return m_time;
}


/***********************************************************************//**
 * @brief Return flux of spectral bin
 *
 * @return Flux of spectral bin.
 ***************************************************************************/
inline
double GMWLDatum::counts(void) const
{
    return m_flux;
}


/***********************************************************************//**
 * @brief Set flux of spectral bin
 *
 * @param[in] flux Flux of spectral bin.
 ***************************************************************************/
inline
void GMWLDatum::counts(const double& flux)
{
    m_flux = flux;
    return;
}


/***********************************************************************//**
 * @brief Return flux of spectral bin
 *
 * @return Flux of spectral bin.
 ***************************************************************************/
inline
const double& GMWLDatum::flux(void) const
{
    return m_flux;
}


/***********************************************************************//**
 * @brief Return flux error of spectral bin
 *
 * @return Flux error of spectral bin.
 ***************************************************************************/
inline
const double& GMWLDatum::flux_err(void) const
{
    return m_flux_err;
}


/***********************************************************************//**
 * @brief Set flux of spectral bin
 *
 * @param[in] flux Flux of spectral bin.
 ***************************************************************************/
inline
void GMWLDatum::flux(const double& flux)
{
    m_flux = flux;
    return;
}


/***********************************************************************//**
 * @brief Set flux error of spectral bin
 *
 * @param[in] error Flux error of spectral bin.
 ***************************************************************************/
inline
void GMWLDatum::flux_err(const double& error)
{
    m_flux_err = error;
    return;
}


/***********************************************************************//**
 * @brief Return energy error of spectral bin
 *
 * @return Energy error of spectral bin.
 ***************************************************************************/
inline
const GEnergy& GMWLDatum::energy_err(void) const
{
    return m_eng_err;
}

#endif /* GMWLDATUM_HPP */
