/***************************************************************************
 *                GCTAExposure.hpp - CTA exposure cube class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Chia-Chun Lu                                     *
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
 * @file GCTAExposure.hpp
 * @brief CTA exposure cube class definition
 * @author Chia-Chun Lu
 */

#ifndef GCTAEXPOSURE_HPP
#define GCTAEXPOSURE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFits.hpp"
#include "GSkymap.hpp"
#include "GObservations.hpp"
#include "GCTAObservation.hpp"


/***********************************************************************//**
 * @class GCTAExposure
 *
 * @brief CTA exposure cube class
 *
 * This class implements a CTA exposure cube which provides the average
 * exposure for binned analysis as function of sky position and log10
 * energy.
 ***************************************************************************/
class GCTAExposure : public GBase {

public:
   
    // Constructors and destructors
    GCTAExposure(void);
    GCTAExposure(const GCTAExposure& cube);
    explicit GCTAExposure(const std::string& filename);
    GCTAExposure(const std::string&   wcs,
                 const std::string&   coords,
                 const double&        x,
                 const double&        y,
                 const double&        dx,
                 const double&        dy,
                 const int&           nx,
                 const int&           ny,
                 const GEbounds&      ebounds);
    virtual ~GCTAExposure(void);

    // Operators
    GCTAExposure& operator=(const GCTAExposure& exp);
    double        operator()(const GSkyDir& dir, const GEnergy& energy) const;

    // Methods
    void               clear(void);
    GCTAExposure*      clone(void) const;
    void               set(const GCTAObservation& obs);
    void               fill(const GObservations& obs);
    const GSkymap&     cube(void) const;
    const GEbounds&    ebounds(void) const;
    const GNodeArray&  elogmeans(void) const;
    void               read(const GFits& fits);
    void               write(GFits& file) const;
    void               load(const std::string& filename);
    void               save(const std::string& filename,
                            const bool& clobber = false) const;
    const std::string& filename(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Methods
    void init_members(void);
    void copy_members(const GCTAExposure& exp);
    void free_members(void);
    void clear_cube(void);
    void update(const double& logE) const;
    void set_eng_axis(void);

    // Data
    mutable std::string m_filename;  //!< Filename
    GSkymap             m_cube;      //!< Average Exposure cube
    GEbounds            m_ebounds;   //!< Energy bounds for the Exposure cube
    GNodeArray          m_elogmeans; //!< Mean energy for the Exposure cube

private:
    // Response table computation cache for 1D access
    mutable int    m_inx_left;       //!< Index of left node
    mutable int    m_inx_right;      //!< Index of right node
    mutable double m_wgt_left;       //!< Weight of left node
    mutable double m_wgt_right;      //!< Weight of right node
};


/***********************************************************************//**
 * @brief Return exposure cube
 *
 * @return Exposure cube.
 *
 * Returns the GSkymap object that is used to store the exposure cube
 * information.
 ***************************************************************************/
inline
const GSkymap& GCTAExposure::cube(void) const
{
    return (m_cube);
}


/***********************************************************************//**
 * @brief Return energy boundaries
 *
 * @return Energy boundaries
 ***************************************************************************/
inline
const GEbounds& GCTAExposure::ebounds(void) const
{
    return (m_ebounds);
}


/***********************************************************************//**
 * @brief Return arithmetic mean of log10 energies
 *
 * @return Arithmetic mean of log10 energies.
 ***************************************************************************/
inline
const GNodeArray& GCTAExposure::elogmeans(void) const
{
    return (m_elogmeans);
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
const std::string& GCTAExposure::filename(void) const
{
    return (m_filename);
}

#endif /* GCTAEXPOSURE_HPP */
