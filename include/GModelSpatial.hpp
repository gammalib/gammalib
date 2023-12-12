/***************************************************************************
 *          GModelSpatial.hpp - Spatial model abstract base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2023 by Juergen Knoedlseder                         *
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
#include "GSkyRegionCircle.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GXmlElement.hpp"
#include "GRan.hpp"
#include "GFunction.hpp"
#include "GIntegral.hpp"

/* __ Forward declarations _______________________________________________ */
class GSkyRegion;


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
    virtual std::string    classname(void) const = 0;
    virtual GClassCode     code(void) const = 0;
    virtual double         eval(const GPhoton& photon,
                                const bool& gradients = false) const = 0;
    virtual GSkyDir        mc(const GEnergy& energy, const GTime& time,
                              GRan& ran) const = 0;
    virtual double         mc_norm(const GSkyDir& dir,
                                   const double&  radius) const = 0;
    virtual bool           contains(const GSkyDir& dir,
                                    const double&  margin = 0.0) const = 0;
    virtual void           read(const GXmlElement& xml) = 0;
    virtual void           write(GXmlElement& xml) const = 0;
    virtual std::string    print(const GChatter& chatter = NORMAL) const = 0;

    // Virtual methods
    virtual double flux(const GSkyRegion& region,
                        const GEnergy&    srcEng  = GEnergy(),
                        const GTime&      srcTime = GTime()) const;

    // Methods
    std::string       type(void) const;
    void              type(const std::string& type);
    GModelPar&        at(const int& index);
    const GModelPar&  at(const int& index) const;
    bool              has_par(const std::string& name) const;
    bool              has_free_pars(void) const;
    int               size(void) const;
    void              autoscale(void);
    const GSkyRegion* region(void) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelSpatial& model);
    void         free_members(void);
    virtual void set_region(void) const = 0;

    // Kernel for spatial model radial integration
    class circle_int_kern_rho : public GFunction {
    public:
        circle_int_kern_rho(const GModelSpatial*    model,
                            const GSkyRegion&       region,
                            const GSkyDir&          centre,
                            const GEnergy&          srcEng,
                            const GTime&            srcTime) :
                            m_model(model),
                            m_region(region),
                            m_centre(centre),
                            m_srcEng(srcEng),
                            m_srcTime(srcTime) { }
        double eval(const double& rho);
    public:
        const GModelSpatial* m_model;     //!< Spatial model
        const GSkyRegion&    m_region;    //!< Sky region
        const GSkyDir&       m_centre;    //!< Model centre
        GEnergy              m_srcEng;    //!< Photon energy
        GTime                m_srcTime;   //!< Photon time
    };

    // Kernel for spatial model azimuth angle integration
    class circle_int_kern_omega : public GFunction {
    public:
        circle_int_kern_omega(const GModelSpatial*    model,
                              const GSkyRegion&       region,
                              const GSkyDir&          centre,
                              const double&           rho,
                              const GEnergy&          srcEng,
                              const GTime&            srcTime) :
                              m_model(model),
                              m_region(region),
                              m_centre(centre),
                              m_rho(rho),
                              m_srcEng(srcEng),
                              m_srcTime(srcTime) { }
        double eval(const double& omega);
    public:
        const GModelSpatial* m_model;   //!< Spatial model
        const GSkyRegion&    m_region;  //!< Sky region
        const GSkyDir&       m_centre;  //!< Model centre
        double               m_rho;     //!< Offset from center of the region
        GEnergy              m_srcEng;  //!< Photon energy
        GTime                m_srcTime; //!< Photon time
    };

    // Proteced members
    std::string             m_type;   //!< Spatial model type
    GSkyRegionCircle        m_region; //!< Bounding circle
    std::vector<GModelPar*> m_pars;   //!< Parameter pointers
};


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the spatial model.
 ***************************************************************************/
inline
std::string GModelSpatial::type(void) const
{
    return (m_type);
}


/***********************************************************************//**
 * @brief Set model type
 *
 * @param[in] type Model type.
 *
 * Set the type of the spatial model.
 ***************************************************************************/
inline
void GModelSpatial::type(const std::string& type)
{
    m_type = type;
    return;
}


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
    return (int)m_pars.size();
}


/***********************************************************************//**
 * @brief Return boundary sky region
 *
 * @return Boundary sky region.
 *
 * Returns a sky region that fully encloses the spatial model component.
 ***************************************************************************/
inline
const GSkyRegion* GModelSpatial::region(void) const
{
    const_cast<GModelSpatial*>(this)->set_region();
    return (&m_region);
}

#endif /* GMODELSPATIAL_HPP */
