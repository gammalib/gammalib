/***************************************************************************
 *             GCTACubeBackground.hpp - CTA cube background class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2018 by Michael Mayer                               *
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
#include "GFilename.hpp"
#include "GSkyMap.hpp"
#include "GNodeArray.hpp"
#include "GEbounds.hpp"

/* __ Forward declarations _______________________________________________ */
class GFits;
class GFitsBinTable;
class GRan;
class GObservations;
class GCTAInstDir;


/***********************************************************************//**
 * @class GCTACubeBackground
 *
 * @brief CTA cube background class
 ***************************************************************************/
class GCTACubeBackground : public GBase {

public:
    // Constructors and destructors
    GCTACubeBackground(void);
    explicit GCTACubeBackground(const GFilename& filename);
    explicit GCTACubeBackground(const GCTAEventCube& cube);
    GCTACubeBackground(const GCTACubeBackground& bgd);
    GCTACubeBackground(const std::string&   wcs,
                       const std::string&   coords,
                       const double&        x,
                       const double&        y,
                       const double&        dx,
                       const double&        dy,
                       const int&           nx,
                       const int&           ny,
                       const GEbounds&      ebounds);
    virtual ~GCTACubeBackground(void);

    // Operators
    GCTACubeBackground& operator=(const GCTACubeBackground& bgd);
    double              operator()(const GCTAInstDir& dir,
                                   const GEnergy&     energy) const;

    // Methods
    void                clear(void);
    GCTACubeBackground* clone(void) const;
    std::string         classname(void) const;
    void                fill(const GObservations& obs, GLog* log = NULL);
    double              integral(const double& logE) const;
    void                read(const GFits& fits);
    void                write(GFits& file) const;
    void                load(const GFilename& filename);
    void                save(const GFilename& filename,
                             const bool&      clobber = false) const;
    const GSkyMap&      cube(void) const;
    const GEbounds&     ebounds(void) const;
    const GFilename&    filename(void) const;
    std::string         print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTACubeBackground& bgd);
    void free_members(void);
    void set_eng_axis(void);
    void update(const double& logE) const;

    // Members
    mutable GFilename m_filename;  //!< Name of background response file
    GSkyMap           m_cube;      //!< Background cube
    GEbounds          m_ebounds;   //!< Energy boundaries of background cube
    GNodeArray        m_elogmeans; //!< Mean energy for the background cube

    // Response table computation cache for 1D access
    mutable int    m_inx_left;       //!< Index of left node
    mutable int    m_inx_right;      //!< Index of right node
    mutable double m_wgt_left;       //!< Weight of left node
    mutable double m_wgt_right;      //!< Weight of right node

};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTACubeBackground").
 ***************************************************************************/
inline
std::string GCTACubeBackground::classname(void) const
{
    return ("GCTACubeBackground");
}


/***********************************************************************//**
 * @brief Return background cube
 *
 * @return Sky map holding the background cube.
 *
 * Returns the GSkyMap object that is used to store the background cube
 * information.
 ***************************************************************************/
inline
const GSkyMap& GCTACubeBackground::cube(void) const
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
 * @brief Return background cube filename
 *
 * @return Background cube filename.
 *
 * Returns the filename from which the background cube was loaded or into
 * which the background cube has been saved.
 ***************************************************************************/
inline
const GFilename& GCTACubeBackground::filename(void) const
{
    return (m_filename);
}

#endif /* GCTACUBEBACKGROUND_HPP */
