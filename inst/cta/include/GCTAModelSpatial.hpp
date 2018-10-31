/***************************************************************************
 *         GCTAModelSpatial.hpp - Spatial model abstract base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GCTAModelSpatial.hpp
 * @brief Abstract spatial model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELSPATIAL_HPP
#define GCTAMODELSPATIAL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFunction.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GModelPar.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GObservation;
class GXmlElement;
class GCTAInstDir;
class GCTAObservation;


/***********************************************************************//**
 * @class GCTAModelSpatial
 *
 * @brief Abstract spatial model class
 *
 * This class implements the spatial component of the CTA background model.
 ***************************************************************************/
class GCTAModelSpatial : public GBase {

public:
    // Constructors and destructors
    GCTAModelSpatial(void);
    GCTAModelSpatial(const GCTAModelSpatial& model);
    virtual ~GCTAModelSpatial(void);

    // Operators
    virtual GCTAModelSpatial& operator=(const GCTAModelSpatial& model);
    virtual GModelPar&        operator[](const int& index);
    virtual const GModelPar&  operator[](const int& index) const;
    virtual GModelPar&        operator[](const std::string& name);
    virtual const GModelPar&  operator[](const std::string& name) const;

    // Pure virtual methods
    virtual void              clear(void) = 0;
    virtual GCTAModelSpatial* clone(void) const = 0;
    virtual std::string       classname(void) const = 0;
    virtual std::string       type(void) const = 0;
    virtual double            eval(const GCTAInstDir& dir,
                                   const GEnergy&     energy,
                                   const GTime&       time,
                                   const bool&        gradients = false) const = 0;
    virtual double            mc_max_value(const GCTAObservation& obs) const = 0;
    virtual void              read(const GXmlElement& xml) = 0;
    virtual void              write(GXmlElement& xml) const = 0;
    virtual std::string       print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented virtual methods
    virtual GCTAInstDir       mc(const GEnergy&         energy,
                                 const GTime&           time,
                                 const GCTAObservation& obs,
                                 GRan&                  ran) const;

    // Methods
    int                       size(void) const;
    virtual double            npred(const GEnergy&      energy,
                                    const GTime&        time,
                                    const GObservation& obs) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelSpatial& model);
    void free_members(void);

    // RoI integration kernel over theta
    class npred_roi_kern_theta : public GFunction {
    public:
        npred_roi_kern_theta(const GCTAModelSpatial* spatial,
                             const GEnergy&          energy,
                             const GTime&            time,
                             const int&              min_iter,
                             const int&              max_iter) :
                             m_spatial(spatial),
                             m_energy(energy),
                             m_time(time),
                             m_min_iter(min_iter),
                             m_max_iter(max_iter) { }
        double eval(const double& theta);
    protected:
        const GCTAModelSpatial* m_spatial;  //!< Pointer to spatial component
        GEnergy                 m_energy;   //!< Energy
        GTime                   m_time;     //!< Time
        int                     m_min_iter; //!< Minimum number of Romberg iterations
        int                     m_max_iter; //!< Maximum number of Romberg iterations
    };

    // RoI integration kernel over phi
    class npred_roi_kern_phi : public GFunction {
    public:
        npred_roi_kern_phi(const GCTAModelSpatial* spatial,
                           const GEnergy&          energy,
                           const GTime&            time,
                           const double&           theta) :
                           m_spatial(spatial),
                           m_energy(energy),
                           m_time(time),
                           m_theta(theta) { }
        double eval(const double& phi);
    protected:
        const GCTAModelSpatial* m_spatial;  //!< Pointer to spatial component
        GEnergy                 m_energy;   //!< Energy
        GTime                   m_time;     //!< Time
        double                  m_theta;    //!< Offset angle (radians)
    };

    // Proteced members
    std::vector<GModelPar*> m_pars;  //!< Parameter pointers
};


/***********************************************************************//**
 * @brief Return number of model parameters
 *
 * @return Number of model parameters.
 ***************************************************************************/
inline
int GCTAModelSpatial::size(void) const
{
    return ((int)m_pars.size());
}

#endif /* GCTAMODELSPATIAL_HPP */
