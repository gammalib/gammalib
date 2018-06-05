/***************************************************************************
 *    GModelSpatialRadialProfile.hpp - Radial profile source model class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file GModelSpatialRadialProfile.hpp
 * @brief Radial profile model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALRADIALPROFILE_HPP
#define GMODELSPATIALRADIALPROFILE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialRadial.hpp"
#include "GNodeArray.hpp"
#include "GSkyRegionCircle.hpp"

/* __ Forward declaration ________________________________________________ */
class GXmlElement;
class GSkyDir;


/**************************************************************************
 * @class GModelSpatialRadialProfile
 *
 * @brief Radial profile source model class
 *
 * This class implements the spatial component of the factorised source
 * model for an arbitrary radial profile.
 ***************************************************************************/
class GModelSpatialRadialProfile : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelSpatialRadialProfile(void);
    explicit GModelSpatialRadialProfile(const GXmlElement& xml);
    GModelSpatialRadialProfile(const GModelSpatialRadialProfile& model);
    virtual ~GModelSpatialRadialProfile(void);

    // Operators
    virtual GModelSpatialRadialProfile& operator=(const GModelSpatialRadialProfile& model);

    // Pure virtual methods
    virtual void                        clear(void) = 0;
    virtual GModelSpatialRadialProfile* clone(void) const = 0;
    virtual std::string                 classname(void) const = 0;
    virtual std::string                 type(void) const = 0;
    virtual double                      theta_min(void) const = 0;
    virtual double                      theta_max(void) const = 0;
    virtual std::string                 print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented pure virtual base class methods
    virtual double      eval(const double&  theta,
                             const GEnergy& energy,
                             const GTime&   time,
                             const bool&    gradients = false) const;
    virtual GSkyDir     mc(const GEnergy& energy,
                           const GTime&   time,
                           GRan&          ran) const;
    virtual bool        contains(const GSkyDir& dir,
                                 const double&  margin = 0.0) const;
    virtual GSkyRegion* region(void) const;

    // Implement other methods
    int  num_nodes(void) const;
    void num_nodes(const int& number);

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GModelSpatialRadialProfile& model);
    void           free_members(void);
    int            cache_index(void) const;
    virtual double profile_value(const double& theta) const = 0;
    void           set_region(void) const;

    // Protected members
    bool                     m_coord_indep; //!< True if model independent
                                            //!< of sky coordinates
    int                      m_num_nodes;   //!< Number of profile nodes
    mutable GSkyRegionCircle m_region;      //!< Bounding circle

    // Pre-computed radial profile
    struct profile {
        std::vector<double> pars;           //!< Profile parameters
        GNodeArray          nodes;          //!< Profile nodes
        std::vector<double> values;         //!< Profile values
        std::vector<double> mc;             //!< Profile for MC
        double              mc_max;         //!< Maximum of profile for MC
    };

    // Pre-computation cache
    mutable std::vector<profile> m_profile;     //!< Pre-computation cache
};


/***********************************************************************//**
 * @brief Return number of nodes
 *
 * @return Number of nodes.
 *
 * Returns the number of nodes in the radial profile.
 ***************************************************************************/
inline
int GModelSpatialRadialProfile::num_nodes(void) const
{
    return (m_num_nodes);
}


/***********************************************************************//**
 * @brief Set number of nodes
 *
 * @param[in] number Number of nodes.
 *
 * Sets the number of nodes in the radial profile.
 ***************************************************************************/
inline
void GModelSpatialRadialProfile::num_nodes(const int& number)
{
    m_num_nodes = number;
    return;
}


/***********************************************************************//**
 * @brief Return boundary sky region
 *
 * @return Boundary sky region.
 *
 * Returns a sky region that fully encloses the spatial model component.
 ***************************************************************************/
inline
GSkyRegion* GModelSpatialRadialProfile::region(void) const
{
    set_region();
    return (&m_region);
}

#endif /* GMODELSPATIALRADIALPROFILE_HPP */
