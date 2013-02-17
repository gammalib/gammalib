/***************************************************************************
 *            GModelSpatialMap.hpp  -  Spatial map model class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GModelSpatialMap.hpp
 * @brief Spatial map model class interface definition
 * @author J. Knoedlseder
 */

#ifndef GMODELSPATIALMAP_HPP
#define GMODELSPATIALMAP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GSkymap.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpatialMap
 *
 * @brief Spatial map model
 *
 * This class implements the spatial component of the factorised source
 * model for a skymap.
 ***************************************************************************/
class GModelSpatialMap : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialMap(void);
    explicit GModelSpatialMap(const GXmlElement& xml);
    explicit GModelSpatialMap(const std::string& filename);
    GModelSpatialMap(const GModelSpatialMap& model);
    virtual ~GModelSpatialMap(void);

    // Operators
    virtual GModelSpatialMap& operator= (const GModelSpatialMap& model);

    // Implemented pure virtual methods
    virtual void              clear(void);
    virtual GModelSpatialMap* clone(void) const;
    virtual std::string       type(void) const { return "SpatialMap"; }
    virtual double            eval(const GSkyDir& srcDir) const;
    virtual double            eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir           mc(GRan& ran) const;
    virtual void              read(const GXmlElement& xml);
    virtual void              write(GXmlElement& xml) const;
    virtual std::string       print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialMap& model);
    void free_members(void);
    void load_map(const std::string& filename);

    // Protected members
    GModelPar           m_value;        //!< Value
    GSkymap             m_map;          //!< Skymap
    std::string         m_filename;     //!< Name of skymap
    std::vector<double> m_mc_cache;     //!< Monte Carlo cache
};

#endif /* GMODELSPATIALMAP_HPP */
