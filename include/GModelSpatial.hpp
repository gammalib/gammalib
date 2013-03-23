/***************************************************************************
 *          GModelSpatial.hpp - Spatial model abstract base class          *
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
 * @file GModelSpatial.hpp
 * @brief Abstract spatial model base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIAL_HPP
#define GMODELSPATIAL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"
#include "GRan.hpp"


/***********************************************************************//**
 * @class GModelSpatial
 *
 * @brief Abstract spatial model base class
 *
 * This class implements the spatial component of the factorized gamma-ray
 * source model.
 ***************************************************************************/
class GModelSpatial : public GBase {

public:
    // Constructors and destructors
    GModelSpatial(void);
    GModelSpatial(const GModelSpatial& model);
    virtual ~GModelSpatial(void);

    // Operators
    virtual GModelSpatial&   operator=(const GModelSpatial& model);
    virtual GModelPar&       operator[](const int& index);
    virtual const GModelPar& operator[](const int& index) const;
    virtual GModelPar&       operator[](const std::string& name);
    virtual const GModelPar& operator[](const std::string& name) const;

    // Pure virtual methods
    virtual void           clear(void) = 0;
    virtual GModelSpatial* clone(void) const = 0;
    virtual std::string    type(void) const = 0;
    virtual double         eval(const GSkyDir& srcDir) const = 0;
    virtual double         eval_gradients(const GSkyDir& srcDir) const = 0;
    virtual GSkyDir        mc(GRan& ran) const = 0;
    virtual void           read(const GXmlElement& xml) = 0;
    virtual void           write(GXmlElement& xml) const = 0;
    virtual std::string    print(void) const = 0;

    // Methods
    int  size(void) const;
    void autoscale(void);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatial& model);
    void free_members(void);

    // Proteced members
    std::vector<GModelPar*> m_pars;  //!< Parameter pointers
};

/***********************************************************************//**
 * @brief Return number of parameters
 *
 * @return Number of parameters in spatial model component.
 *
 * Returns the number of parameters in the spatial model component.
 ***************************************************************************/
inline
int GModelSpatial::size(void) const
{
    return (m_pars.size());
}

#endif /* GMODELSPATIAL_HPP */
