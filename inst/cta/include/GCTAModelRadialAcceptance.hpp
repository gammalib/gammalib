/***************************************************************************
 *     GCTAModelRadialAcceptance.hpp  -  Radial acceptance model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2012 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadialAcceptance.hpp
 * @brief Radial acceptance model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELRADIALACCEPTANCE_HPP
#define GCTAMODELRADIALACCEPTANCE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <cmath>
#include "GModelData.hpp"
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GEvent.hpp"
#include "GObservation.hpp"
#include "GXmlElement.hpp"
#include "GFunction.hpp"
#include "GCTAEventList.hpp"
#include "GCTAModelRadial.hpp"


/***********************************************************************//**
 * @class GCTAModelRadialAcceptance
 *
 * @brief Radial acceptance model class
 *
 * This class implements a radial acceptance model for CTA.
 ***************************************************************************/
class GCTAModelRadialAcceptance : public GModelData {

public:
    // Constructors and destructors
    GCTAModelRadialAcceptance(void);
    explicit GCTAModelRadialAcceptance(const GXmlElement& xml);
    explicit GCTAModelRadialAcceptance(const GCTAModelRadial& radial,
                                       const GModelSpectral& spectral);
    GCTAModelRadialAcceptance(const GCTAModelRadialAcceptance& model);
    virtual ~GCTAModelRadialAcceptance(void);

    // Operators
    virtual GCTAModelRadialAcceptance& operator=(const GCTAModelRadialAcceptance& model);

    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GCTAModelRadialAcceptance* clone(void) const;
    virtual std::string                type(void) const { return "RadialAcceptance"; }
    virtual double                     eval(const GEvent& event,
                                            const GObservation& obs) const;
    virtual double                     eval_gradients(const GEvent& event,
                                                      const GObservation& obs) const;
    virtual double                     npred(const GEnergy& obsEng, const GTime& obsTime,
                                             const GObservation& obs) const;
    virtual GCTAEventList*             mc(const GObservation& obs, GRan& ran) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;
    virtual std::string                print(void) const;

    // Other methods
    GCTAModelRadial* radial(void)   const { return m_radial; }
    GModelSpectral*  spectral(void) const { return m_spectral; }
    GModelTemporal*  temporal(void) const { return m_temporal; }

protected:
    // Protected methods
    void             init_members(void);
    void             copy_members(const GCTAModelRadialAcceptance& model);
    void             free_members(void);
    void             set_pointers(void);
    bool             valid_model(void) const;
    GCTAModelRadial* xml_radial(const GXmlElement& radial) const;
    GModelSpectral*  xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal*  xml_temporal(const GXmlElement& temporal) const;

    // ROI integration kernel
    class roi_kern : public GFunction {
    public:
        roi_kern(const GCTAModelRadial* parent, const double& roi, const double& dist) :
                 m_parent(parent),
                 m_roi(roi),
                 m_cosroi(std::cos(roi)),
                 m_dist(dist),
                 m_cosdist(std::cos(dist)),
                 m_sindist(std::sin(dist)) { }
        double eval(double r);
    protected:
        const GCTAModelRadial* m_parent;   //!< Pointer to radial model
        double                 m_roi;      //!< ROI radius in radians
        double                 m_cosroi;   //!< Cosine of ROI radius
        double                 m_dist;     //!< Distance between pointing and ROI centre in radians
        double                 m_cosdist;  //!< Cosine of distance
        double                 m_sindist;  //!< Sinus of distance
    };

    // Proteced data members
    GCTAModelRadial* m_radial;       //!< Radial model
    GModelSpectral*  m_spectral;     //!< Spectral model
    GModelTemporal*  m_temporal;     //!< Temporal model
};

#endif /* GCTAMODELRADIALACCEPTANCE_HPP */
