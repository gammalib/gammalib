/***************************************************************************
 *              GCTACubeBackground.hpp - CTA cube background class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Michael Mayer                     *
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
 * @file GCTACubeBackground.hpp
 * @brief CTA cube background class definition
 * @author Michael Mayer
 */

#ifndef GCTACUBEBACKGROUND_HPP
#define GCTACUBEBACKGROUND_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GSkymap.hpp"
#include "GModelSpectralNodes.hpp"
#include "GEbounds.hpp"

/* __ Forward declarations _______________________________________________ */
class GFits;
class GFitsBinTable;


/***********************************************************************//**
 * @class GCTACubeBackground
 *
 * @brief CTA cube background class
 ***************************************************************************/
class GCTACubeBackground : public GBase {

public:
    // Constructors and destructors
    GCTACubeBackground(void);
    explicit GCTACubeBackground(const std::string& filename);
    GCTACubeBackground(const GCTACubeBackground& bgd);
    virtual ~GCTACubeBackground(void);

    // Operators
    GCTACubeBackground& operator=(const GCTACubeBackground& bgd);
    double            operator()(const GSkyDir& dir, const GEnergy& energy) const;

    // Implemented pure virtual methods
    void                       clear(void);
    GCTACubeBackground*          clone(void) const;
    std::string                classname(void) const;
    void               set(const GSkymap& cube, const GEbounds& ebounds);
    void                       load(const std::string& filename);
    void                       read(const GFits& fits);
    void                       write(GFits& file) const;
    const std::string&                filename(void) const;
    GSkyDir               mc(const GEnergy& energy,
                                  const GTime& time,
                                  GRan& ran) const;
    const GModelSpectralNodes& spectrum(void) const;
    std::string                print(const GChatter& chatter = NORMAL) const;

    // Methods
    const GSkymap&     cube(void) const;
    const GEbounds&    ebounds(void) const;
    const GNodeArray&  elogmeans(void) const;
    void                       save(const std::string& filename,
                                    const bool& clobber = false) const;
    void                       set_mc_cone(const GSkyDir& centre, const double&  radius);
    double                  integral(const double& logE) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTACubeBackground& bgd);
    void free_members(void);
    void set_eng_axis(void);
    void update(const double& logE) const;
    void init_mc_cache();

    // Members
    mutable std::string       m_filename;    //!< Name of background response file
    GSkymap             m_cube;      //!< Background cube
    GEbounds            m_ebounds;   //!< Energy bounds for the Exposure cube
    GNodeArray          m_elogmeans; //!< Mean energy for the Exposure cube

    // Response table computation cache for 1D access
    mutable int    m_inx_left;       //!< Index of left node
    mutable int    m_inx_right;      //!< Index of right node
    mutable double m_wgt_left;       //!< Weight of left node
    mutable double m_wgt_right;      //!< Weight of right node

    // Monte Carlo cache
    mutable std::vector<double> m_mc_cache;    //!< Monte Carlo cache
    mutable GModelSpectralNodes m_mc_spectrum; //!< Response cube spectrum

};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTCubeABackground").
 ***************************************************************************/
inline
std::string GCTACubeBackground::classname(void) const
{
    return ("GCTACubeBackground");
}


/***********************************************************************//**
 * @brief Return exposure cube
 *
 * @return Exposure cube.
 *
 * Returns the GSkymap object that is used to store the exposure cube
 * information.
 ***************************************************************************/
inline
const GSkymap& GCTACubeBackground::cube(void) const
{
    return (m_cube);
}


/***********************************************************************//**
 * @brief Return energy boundaries
 *
 * @return Energy boundaries
 ***************************************************************************/
inline
const GEbounds& GCTACubeBackground::ebounds(void) const
{
    return (m_ebounds);
}


/***********************************************************************//**
 * @brief Return arithmetic mean of log10 energies
 *
 * @return Arithmetic mean of log10 energies.
 ***************************************************************************/
inline
const GNodeArray& GCTACubeBackground::elogmeans(void) const
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
const std::string& GCTACubeBackground::filename(void) const
{
    return (m_filename);
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
const GModelSpectralNodes& GCTACubeBackground::spectrum(void) const
{
    return (m_mc_spectrum);
}

#endif /* GCTACUBEBACKGROUND_HPP */
