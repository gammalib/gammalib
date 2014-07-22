/***************************************************************************
 *          GModelSpatial.hpp - Spatial model abstract base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
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
#include "GPhoton.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GXmlElement.hpp"
#include "GRan.hpp"


/***********************************************************************//**
 * @class GModelSpatial
 *
 * @brief Abstract spatial model base class
 *
 * This class implements the spatial component of the factorized gamma-ray
 * source model. The spatial component is given by
 *
 * \f[
 *    S_{\rm p}(\vec{p} | E, t)
 * \f]
 *
 * where
 * - \f$\vec{p}\f$ is the true photon arrival direction,
 * - \f$E\f$ is the true photon energy, and
 * - \f$t\f$ is the true photon arrival time.
 *
 * The spatial component describes the energy and time dependent morphology
 * of the source. It satisfies
 * \f[
 *    \int_{\Omega} S_{\rm p}(\vec{p} | E, t) d\Omega = 1
 * \f]
 * for all \f$E\f$ and \f$t\f$, hence the spatial component does not
 * impact the spatially integrated spectral and temporal properties of the
 * source.
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
    virtual GClassCode     code(void) const = 0;
    virtual double         eval(const GPhoton& photon) const = 0;
    virtual double         eval_gradients(const GPhoton& photon) const = 0;
    virtual GSkyDir        mc(const GEnergy& energy, const GTime& time,
                              GRan& ran) const = 0;
    virtual double         norm(const GSkyDir& dir,
                                const double&  radius) const = 0;
    virtual void           read(const GXmlElement& xml) = 0;
    virtual void           write(GXmlElement& xml) const = 0;
    virtual std::string    print(const GChatter& chatter = NORMAL) const = 0;

    // Methods
    GModelPar&       at(const int& index);
    const GModelPar& at(const int& index) const;
    int              size(void) const;
    void             autoscale(void);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatial& model);
    void free_members(void);

    // Proteced members
    std::vector<GModelPar*> m_pars;  //!< Parameter pointers
};


/***********************************************************************//**
 * @brief Returns model parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Model parameter.
 *
 * Returns model parameter without @p index range checking.
 ***************************************************************************/
inline
GModelPar& GModelSpatial::operator[](const int& index)
{
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns model parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Model parameter.
 *
 * Returns model parameter without @p index range checking.
 ***************************************************************************/
inline
const GModelPar& GModelSpatial::operator[](const int& index) const
{
    return *(m_pars[index]);
}


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
