/***************************************************************************
 *                  GCTAPointing.hpp - CTA pointing class                  *
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
 * @file GCTAPointing.hpp
 * @brief CTA pointing class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAPOINTING_HPP
#define GCTAPOINTING_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GSkyDir.hpp"
#include "GTime.hpp"
#include "GTimeReference.hpp"
#include "GMatrix.hpp"
#include "GHorizDir.hpp"
#include "GNodeArray.hpp"
#include "GFitsTable.hpp"
#include "GCTAInstDir.hpp"


/***********************************************************************//**
 * @class GCTAPointing
 *
 * @brief CTA pointing class.
 *
 * This class implements a CTA pointing. For the time being it is assumed
 * that the pointing direction is time-independent.
 *
 * @todo No transformation from sky coordinates to geographic coordinates
 *       has so far been implemented. The azimuth and zenith angle are not
 *       meaningful.
 ***************************************************************************/
class GCTAPointing : public GBase {

public:
    // Constructors and destructors
    GCTAPointing(void);
    explicit GCTAPointing(const GSkyDir& dir);
    GCTAPointing(const std::string& filename,
                 const std::string& extname = "POINTING");
    GCTAPointing(const GCTAPointing& pnt);
    virtual ~GCTAPointing(void);

    // Operators
    virtual GCTAPointing& operator=(const GCTAPointing& pnt);

    // Implemented pure virtual methods
    virtual void          clear(void);
    virtual GCTAPointing* clone(void) const;
    virtual std::string   classname(void) const;
    virtual std::string   print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const GSkyDir& dir(void) const;
    void           dir(const GSkyDir& dir);
    GCTAInstDir    instdir(const GSkyDir& skydir) const;
    GSkyDir        skydir(const GCTAInstDir& instdir) const;
    const GMatrix& rot(void) const;
    const double&  zenith(void) const;
    const double&  azimuth(void) const;
    void           zenith(const double& zenith);  
    void           azimuth(const double& azimuth); 
    GHorizDir      dir_horiz(const GTime& time) const;
    void           load(const std::string& filename,
                        const std::string& extname = "POINTING");
    void           read(const GFitsTable& table);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAPointing& pnt);
    void free_members(void);
    void update(void) const;

    // Protected members
    GSkyDir             m_dir;         //!< Pointing direction in sky coordinates
    double              m_zenith;      //!< Pointing zenith angle
    double              m_azimuth;     //!< Pointing azimuth angle
    bool                m_has_table;   //!< Table is loaded
    GNodeArray          m_table_nodes; //!< Pointing nodes
    std::vector<double> m_table_az;    //!< Table of azimuths (rad)
    std::vector<double> m_table_alt;   //!< Table of altitudes (rad)
    GTime               m_table_tmin;  //!< Min time bound in table
    GTime               m_table_tmax;  //!< Max time bound in table
    GTimeReference      m_reference;   //!< Time reference

    // Cached members
    mutable bool    m_has_cache;  //!< Has transformation cache
    mutable GMatrix m_Rback;      //!< Rotation matrix
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAPointing").
 ***************************************************************************/
inline
std::string GCTAPointing::classname(void) const
{
    return ("GCTAPointing");
}


/***********************************************************************//**
 * @brief Return pointing sky direction
 *
 * @return Pointing sky direction.
 ***************************************************************************/
inline
const GSkyDir& GCTAPointing::dir(void) const
{
    return m_dir;
}


/***********************************************************************//**
 * @brief Return pointing zenith angle
 *
 * @return Pointing zenith angle.
 ***************************************************************************/
inline
const double& GCTAPointing::zenith(void) const
{
    return m_zenith;
}


/***********************************************************************//**
 * @brief Return pointing azimuth angle
 *
 * @return Pointing zenith angle.
 ***************************************************************************/
inline
const double& GCTAPointing::azimuth(void) const
{
    return m_azimuth;
}

/***********************************************************************//**
 * @brief assign zenith angle
 *
 * @param[in] zenith The zenith angle. (deg)
 ***************************************************************************/
inline
void GCTAPointing::zenith(const double& zenith)
{
    m_zenith = zenith;
}


/***********************************************************************//**
 * @brief assign azimuth angle
 *
 * @param[in] azimuth The azimuth angle. (deg)
 ***************************************************************************/
inline
void GCTAPointing::azimuth(const double& azimuth)
{
    m_azimuth = azimuth;
}
#endif /* GCTAPOINTING_HPP */
