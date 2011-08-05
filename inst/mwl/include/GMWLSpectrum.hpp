/***************************************************************************
 *          GMWLSpectrum.hpp  -  Multi-wavelength spectrum class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GMWLSpectrum.hpp
 * @brief Multi-wavelength spectrum class interface definition
 * @author J. Knodlseder
 */

#ifndef GMWLSPECTRUM_HPP
#define GMWLSPECTRUM_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEventCube.hpp"
#include "GMWLDatum.hpp"
#include "GFitsTable.hpp"
#include "GEnergy.hpp"


/***********************************************************************//**
 * @class GMWLSpectrum
 *
 * @brief Multi-wavelength spectrum class interface
 *
 * This class defines a multi-wavelength spectrum and is a container for
 * spectral points of type GMWLDatum. It derives from the abstract
 * GEventCube base class.
 ***************************************************************************/
class GMWLSpectrum : public GEventCube {

public:
    // Constructors and destructors
    GMWLSpectrum(void);
    explicit GMWLSpectrum(const std::string& filename);
    GMWLSpectrum(const GMWLSpectrum& spec);
    virtual ~GMWLSpectrum(void);

    // Operators
    virtual GMWLSpectrum&    operator=(const GMWLSpectrum& spec);
    virtual GMWLDatum*       operator[](const int& index);
    virtual const GMWLDatum* operator[](const int& index) const;

    // Implemented pure virtual methods
    virtual void          clear(void);
    virtual GMWLSpectrum* clone(void) const;
    virtual int           size(void) const { return m_data.size(); }
    virtual int           dim(void) const { return 1; }
    virtual int           naxis(int axis) const { return m_data.size(); }
    virtual void          load(const std::string& filename);
    virtual void          save(const std::string& filename, bool clobber = false) const;
    virtual void          read(const GFits& file);
    virtual void          write(GFits& file) const;
    virtual int           number(void) const;
    virtual std::string   print(void) const;

    // Other methods
    void                  load(const std::string& filename, const std::string& extname);
    void                  load(const std::string& filename, int extno);
    void                  read(const GFits& file, const std::string& extname);
    void                  read(const GFits& file, int extno);
    std::string           telescope(void) const { return m_telescope; }
    std::string           instrument(void) const { return m_instrument; }

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GMWLSpectrum& spec);
    void         free_members(void);
    virtual void set_energies(void) { return; }
    virtual void set_times(void) { return; }
    void         set_ebounds(void);
    void         read_fits(const GFitsTable* table);
    GEnergy      conv_energy(const double& energy, const std::string& unit);
    double       conv_flux(const GEnergy& energy, const double& flux, const std::string& unit);

    // Protected members
    std::string            m_telescope;   //!< Telescope name
    std::string            m_instrument;  //!< Instrument name
    std::vector<GMWLDatum> m_data;        //!< Spectral data
};

#endif /* GMWLSPECTRUM_HPP */
