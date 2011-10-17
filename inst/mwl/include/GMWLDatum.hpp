/***************************************************************************
 *         GMWLDatum.hpp  -  Multi-wavelength spectral point class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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
 * @author J. Knodlseder
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
    virtual double             size(void) const { return 1.0; }
    virtual const GInstDir&    dir(void) const { return m_dir; }
    virtual const GEnergy&     energy(void) const { return m_eng; }
    virtual const GTime&       time(void) const { return m_time; }
    virtual double             counts(void) const { return m_flux; }
    virtual double             error(void) const;
    virtual void               counts(const double& counts) { m_flux=counts; }
    virtual std::string        print(void) const;

    // Other methods
    GEnergy energy_err(void) const { return m_eng_err; }
    double  flux(void) const { return m_flux; }
    double  flux_err(void) const { return m_flux_err; }
    void    flux(const double& flux) { m_flux=flux; }
    void    flux_err(const double& flux_err) { m_flux_err=flux_err; }

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

#endif /* GMWLDATUM_HPP */
