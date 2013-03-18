/***************************************************************************
 *     GModelSpatialPointSource.hpp - Spatial point source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialPointSource.hpp
 * @brief Point source spatial model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALPOINTSOURCE_HPP
#define GMODELSPATIALPOINTSOURCE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpatialPointSource
 *
 * @brief Point source spatial model
 *
 * This class implements a point source as the spatial component of the
 * factorised source model. The point source has two parameters: the Right
 * Ascension and Declination of the point source location.
 *
 * The model is of type "SkyDirFunction".
 ***************************************************************************/
class GModelSpatialPointSource : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialPointSource(void);
    explicit GModelSpatialPointSource(const GSkyDir& dir);
    explicit GModelSpatialPointSource(const GXmlElement& xml);
    GModelSpatialPointSource(const GModelSpatialPointSource& model);
    virtual ~GModelSpatialPointSource(void);

    // Operators
    virtual GModelSpatialPointSource& operator=(const GModelSpatialPointSource& model);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialPointSource* clone(void) const;
    virtual std::string               type(void) const { return "SkyDirFunction"; }
    virtual double                    eval(const GSkyDir& srcDir) const;
    virtual double                    eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir                   mc(GRan& ran) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;
    virtual std::string               print(void) const;

    // Other methods
    double  ra(void) const { return m_ra.value(); }
    double  dec(void) const { return m_dec.value(); }
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialPointSource& model);
    void free_members(void);

    // Protected members
    GModelPar m_ra;          //!< Right Ascension (deg)
    GModelPar m_dec;         //!< Declination (deg)
};

#endif /* GMODELSPATIALPOINTSOURCE_HPP */
