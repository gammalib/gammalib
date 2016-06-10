/***************************************************************************
 *          GCTACubeExposure.hpp - CTA cube analysis exposure class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Chia-Chun Lu                                *
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
 * @file GCTACubeExposure.hpp
 * @brief CTA cube analysis exposure class definition
 * @author Chia-Chun Lu
 */

#ifndef GCTACUBEEXPOSURE_HPP
#define GCTACUBEEXPOSURE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFilename.hpp"
#include "GSkyMap.hpp"
#include "GEnergies.hpp"
#include "GNodeArray.hpp"
#include "GGti.hpp"

/* __ Forward declarations _______________________________________________ */
class GFits;
class GFitsHDU;
class GObservations;
class GCTAObservation;
class GCTAEventCube;


/***********************************************************************//**
 * @class GCTACubeExposure
 *
 * @brief CTA exposure cube class
 *
 * This class implements a CTA exposure cube which provides the average
 * exposure for binned analysis as function of sky position and energy.
 ***************************************************************************/
class GCTACubeExposure : public GBase {

public:
   
    // Constructors and destructors
    GCTACubeExposure(void);
    GCTACubeExposure(const GCTACubeExposure& cube);
    explicit GCTACubeExposure(const GFilename& filename);
    explicit GCTACubeExposure(const GCTAEventCube& cube);
    GCTACubeExposure(const std::string&   wcs,
                     const std::string&   coords,
                     const double&        x,
                     const double&        y,
                     const double&        dx,
                     const double&        dy,
                     const int&           nx,
                     const int&           ny,
                     const GEnergies&     energies);
    virtual ~GCTACubeExposure(void);

    // Operators
    GCTACubeExposure& operator=(const GCTACubeExposure& exp);
    double            operator()(const GSkyDir& dir, const GEnergy& energy) const;

    // Methods
    void               clear(void);
    GCTACubeExposure*  clone(void) const;
    std::string        classname(void) const;
    void               set(const GCTAObservation& obs);
    void               fill(const GObservations& obs, GLog* log = NULL);
    const GSkyMap&     cube(void) const;
    const GEnergies&   energies(void) const;
    const GGti&        gti(void) const;
    const GNodeArray&  elogmeans(void) const;
    const double&      livetime(void) const;
    const double&      ontime(void) const;
    double             deadc(void) const;
    void               read(const GFits& fits);
    void               write(GFits& file) const;
    void               load(const GFilename& filename);
    void               save(const GFilename& filename,
                            const bool&      clobber = false) const;
    const GFilename&   filename(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Methods
    void init_members(void);
    void copy_members(const GCTACubeExposure& exp);
    void free_members(void);
    void update(const double& logE) const;
    void set_eng_axis(void);
    void read_attributes(const GFitsHDU& hdu);
    void write_attributes(GFitsHDU& hdu) const;

    // Members
    mutable GFilename m_filename;  //!< Filename
    GSkyMap           m_cube;      //!< Average Exposure cube
    GEnergies         m_energies;  //!< Energy values for the Exposure cube
    GNodeArray        m_elogmeans; //!< Mean energy for the Exposure cube
    GGti              m_gti;       //!< Good time interval for the Exposure cube

    // Exposure attributes
    double            m_livetime;  //!< Livetime (sec)

private:
    // Response table computation cache for 1D access
    mutable int    m_inx_left;       //!< Index of left node
    mutable int    m_inx_right;      //!< Index of right node
    mutable double m_wgt_left;       //!< Weight of left node
    mutable double m_wgt_right;      //!< Weight of right node
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTACubeExposure").
 ***************************************************************************/
inline
std::string GCTACubeExposure::classname(void) const
{
    return ("GCTACubeExposure");
}


/***********************************************************************//**
 * @brief Return exposure cube
 *
 * @return Exposure cube.
 *
 * Returns the GSkyMap object that is used to store the exposure cube
 * information.
 ***************************************************************************/
inline
const GSkyMap& GCTACubeExposure::cube(void) const
{
    return (m_cube);
}


/***********************************************************************//**
 * @brief Return energies
 *
 * @return Energies
 ***************************************************************************/
inline
const GEnergies& GCTACubeExposure::energies(void) const
{
    return (m_energies);
}


/***********************************************************************//**
 * @brief Return arithmetic mean of log10 energies
 *
 * @return Arithmetic mean of log10 energies.
 ***************************************************************************/
inline
const GNodeArray& GCTACubeExposure::elogmeans(void) const
{
    return (m_elogmeans);
}


/***********************************************************************//**
 * @brief Return Good Time Intervals
 *
 * @return Good Time Intervals.
 ***************************************************************************/
inline
const GGti& GCTACubeExposure::gti(void) const
{
    return (m_gti);
}


/***********************************************************************//**
 * @brief Return livetime
 *
 * @return Livetime (seconds).
 ***************************************************************************/
inline
const double& GCTACubeExposure::livetime(void) const
{
    return (m_livetime);
}


/***********************************************************************//**
 * @brief Return ontime
 *
 * @return Ontime (seconds).
 ***************************************************************************/
inline
const double& GCTACubeExposure::ontime(void) const
{
    return (m_gti.ontime());
}


/***********************************************************************//**
 * @brief Return deadtime correction
 *
 * @return Deadtime correction factor.
 ***************************************************************************/
inline
double GCTACubeExposure::deadc(void) const
{
    return ((m_gti.ontime() > 0.0) ? m_livetime/m_gti.ontime() : 1.0);
}


/***********************************************************************//**
 * @brief Return exposure cube filename
 *
 * @return Exposure cube filename.
 *
 * Returns the filename from which the exposure cube was loaded or into which
 * the exposure cube has been saved.
 ***************************************************************************/
inline
const GFilename& GCTACubeExposure::filename(void) const
{
    return (m_filename);
}

#endif /* GCTACUBEEXPOSURE_HPP */
