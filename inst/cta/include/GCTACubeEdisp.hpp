/***************************************************************************
 *      GCTACubeEdisp.hpp - CTA cube analysis energy dispersion class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Michael Mayer                                    *
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
 * @file GCTACubeEdisp.hpp
 * @brief CTA cube analysis energy disperson class definition
 * @author Michael Mayer
 */

#ifndef GCTACUBEEDISP_HPP
#define GCTACUBEEDISP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GSkyMap.hpp"
#include "GEbounds.hpp"
#include "GEnergies.hpp"
#include "GNodeArray.hpp"

/* __ Forward declarations _______________________________________________ */
class GFits;
class GFilename;
class GObservations;
class GCTAEventCube;
class GCTAObservation;


/***********************************************************************//**
 * @class GCTACubeEdisp
 *
 * @brief CTA energy dispersion for stacked analysis
 *
 * This class implements a mean CTA energy dispersion which provides the
 * average energy dispersion probability density function for stacked
 * analysis as function of sky position, energy migration, and true log10
 * energy.
 ***************************************************************************/
class GCTACubeEdisp : public GBase {

public:
    // Constructors and destructors
	GCTACubeEdisp(void);
	GCTACubeEdisp(const GCTACubeEdisp& cube);
    explicit GCTACubeEdisp(const GFilename& filename);
    GCTACubeEdisp(const GCTAEventCube& cube,
                  const double&        mmax,
                  const int&           nmbins);
    GCTACubeEdisp(const std::string&   wcs,
                  const std::string&   coords,
                  const double&        x,
                  const double&        y,
                  const double&        dx,
                  const double&        dy,
                  const int&           nx,
                  const int&           ny,
                  const GEnergies&     energies,
                  const double&        mmax,
                  const int&           nmbins);
    virtual ~GCTACubeEdisp(void);

    // Operators
    GCTACubeEdisp& operator=(const GCTACubeEdisp& cube);
    double         operator()(const GSkyDir& dir,
                              const double&  migra,
                              const GEnergy& energy) const;

    // Methods
    void              clear(void);
    GCTACubeEdisp*    clone(void) const;
    std::string       classname(void) const;
    void              set(const GCTAObservation& obs);
    void              fill(const GObservations& obs, GLog* log = NULL);
    const GSkyMap&    map(void) const;
    const GEnergies&  energies(void) const;
    const GNodeArray& migras(void) const;
    double            migra_max(void) const;
    int               offset(const int& imigra, const int& iebin) const;
    GEbounds          ebounds(const GEnergy& obsEng) const;
    void              read(const GFits& fits);
    void              write(GFits& file) const;
    void              load(const GFilename& filename);
    void              save(const GFilename& filename,
                           const bool&      clobber = false) const;
    const GFilename&  filename(void) const;
    std::string       print(const GChatter& chatter = NORMAL) const;

protected:
    // Methods
    void init_members(void);
    void copy_members(const GCTACubeEdisp& cube);
    void free_members(void);
    void clear_cube(void);
    void update(const double& migra, const double& logEsrc) const;
    void set_eng_axis(void);
    void compute_ebounds(void) const;

    // Members
    mutable GFilename m_filename;  //!< Filename
    GSkyMap           m_cube;      //!< Energy dispersion cube
    GEnergies         m_energies;  //!< True energy values of cube
    GNodeArray        m_elogmeans; //!< Mean log10TeV energy for the Edisp cube
    GNodeArray        m_migras;    //!< Migra bins for the Edisp cube

private:
    // Response table computation cache for 2D access
    mutable int                   m_inx1;    //!< Index of upper left node
    mutable int                   m_inx2;    //!< Index of lower left node
    mutable int                   m_inx3;    //!< Index of upper right node
    mutable int                   m_inx4;    //!< Index of lower right node
    mutable double                m_wgt1;    //!< Weight of upper left node
    mutable double                m_wgt2;    //!< Weight of lower left node
    mutable double                m_wgt3;    //!< Weight of upper right node
    mutable double                m_wgt4;    //!< Weight of lower right node
    mutable std::vector<GEbounds> m_ebounds; //!< Energy boundaries
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTACubeEdisp").
 ***************************************************************************/
inline
std::string GCTACubeEdisp::classname(void) const
{
    return ("GCTACubeEdisp");
}


/***********************************************************************//**
 * @brief Return energy dispersion cube sky map
 *
 * @return Energy dispersion cube sky map.
 *
 * Returns the energy dispersion cube sky map.
 ***************************************************************************/
inline
const GSkyMap& GCTACubeEdisp::map(void) const
{
    return (m_cube);
}


/***********************************************************************//**
 * @brief Return energies for energy dispersion cub
 *
 * @return Energies for energy dispersion cube
 *
 * Returns the true energies of the energy dispersion cube.
 ***************************************************************************/
inline
const GEnergies& GCTACubeEdisp::energies(void) const
{
    return (m_energies);
}


/***********************************************************************//**
 * @brief Return migra fractions of true and measured photon energy
 *
 * @return Offset migra between fractions of true and measured photon energy
 ***************************************************************************/
inline
const GNodeArray& GCTACubeEdisp::migras(void) const
{
    return (m_migras);
}


/***********************************************************************//**
 * @brief Return maximum migra value
 *
 * @return Maximum migra value.
 ***************************************************************************/
inline
double GCTACubeEdisp::migra_max(void) const
{
    // Get maximum delta value
    double migra_max = (m_migras.size() > 0) ? m_migras[m_migras.size()-1] : 0.0;
    
    // Return
    return (migra_max);
}


/***********************************************************************//**
 * @brief Return edisp cube filename
 *
 * @return Energy dispersion cube filename.
 *
 * Returns the filename from which the energy dispersion cube was loaded or into which
 * the cube has been saved.
 ***************************************************************************/
inline
const GFilename& GCTACubeEdisp::filename(void) const
{
    return (m_filename);
}


/***********************************************************************//**
 * @brief Return map offset
 *
 * @param[in] imigra Migration bin index.
 * @param[in] iebin Energy bin index.
 * @return Map offset.
 *
 * Returns the offset of the energy dispersion in the sky map for a given
 * migration bin index and a given energy bin index.
 ***************************************************************************/
inline
int GCTACubeEdisp::offset(const int& imigra, const int& iebin) const
{
    return (imigra + iebin*m_migras.size());
}

#endif /* GCTACUBEEDISP_HPP */
