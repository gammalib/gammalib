/***************************************************************************
 *      GCTAModelBackground.hpp - generic background model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @author Michael Mayer
 */

#ifndef GCTAMODELBACKGROUND_HPP
#define GCTAMODELBACKGROUND_HPP

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
#include "GModelSpatial.hpp"


/***********************************************************************//**
 * @class GCTAModelBackground
 *
 * @brief CTA background model class
 *
 * This class implements a general background model for CTA.
 ***************************************************************************/
class GCTAModelBackground : public GModelData {

public:
    // Constructors and destructors
	GCTAModelBackground(void);
    explicit GCTAModelBackground(const GXmlElement& xml);
    explicit GCTAModelBackground(const GModelSpatial& spatial,
                                       const GModelSpectral& spectral);
    GCTAModelBackground(const GCTAModelBackground& model);
    virtual ~GCTAModelBackground(void);

    // Operators
    virtual GCTAModelBackground& operator=(const GCTAModelBackground& model);

    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GCTAModelBackground* clone(void) const;
    virtual std::string                type(void) const { return "CTABackground"; }
    virtual double                     eval(const GEvent& event,
                                            const GObservation& obs) const;
    virtual double                     eval_gradients(const GEvent& event,
                                                      const GObservation& obs) const;
    virtual double                     npred(const GEnergy& obsEng, const GTime& obsTime,
                                             const GObservation& obs) const;
    virtual GCTAEventList*             mc(const GObservation& obs, GRan& ran) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;
    virtual std::string                print(const GChatter& chatter = NORMAL) const;

    // Other methods
    GModelSpatial* spatial(void)   const { return m_spatial; }
    GModelSpectral*  spectral(void) const { return m_spectral; }
    GModelTemporal*  temporal(void) const { return m_temporal; }

protected:
    // Protected methods
    void             init_members(void);
    void             copy_members(const GCTAModelBackground& model);
    void             free_members(void);
    void             set_pointers(void);
    bool             valid_model(void) const;
    GModelSpatial* xml_spatial(const GXmlElement& spatial) const;
    GModelSpectral*  xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal*  xml_temporal(const GXmlElement& temporal) const;

    class npred_roi_kern_theta : public GFunction {
    public:
    	npred_roi_kern_theta(const GModelSpatial*   model,
                                     const GEnergy&         obsEng,
                                     const GTime&           obsTime,
                                     const GMatrix&         rot,
                                     double roi,
                                     double dist,
                                     double omega0) :
                                     m_model(model),
                                     m_obsEng(obsEng),
                                     m_obsTime(obsTime),
                                     m_rot(rot),
                                     m_roi(roi),
                                     m_cosroi(std::cos(roi)),
                                     m_dist(dist),
                                     m_cosdist(std::cos(dist)),
                                     m_sindist(std::sin(dist)),
                                     m_omega0(omega0) { }
        double eval(double theta);
    protected:
        const GModelSpatial*   m_model;      //!< Spatial model
        const GEnergy&         m_obsEng;     //!< True photon energy
        const GTime&           m_obsTime;    //!< True photon arrival time
        const GMatrix&         m_rot;        //!< Rotation matrix
        double                 m_roi;      //!< ROI radius in radians
		double                 m_cosroi;   //!< Cosine of ROI radius
		double                 m_dist;     //!< Distance between pointing and ROI centre in radians
		double                 m_cosdist;  //!< Cosine of distance
		double                 m_sindist;  //!< Sinus of distance
	    const double&                  m_omega0;     //!< Position angle of ROI
    };

    class npred_roi_kern_phi : public GFunction {
    public:
    	npred_roi_kern_phi(const GModelSpatial*   model,
                                   const GEnergy&         obsEng,
                                   const GTime&           obsTime,
                                   const GMatrix&         rot,
                                   const double&          theta,
                                   const double&          sin_theta) :
                                   m_model(model),
                                   m_obsEng(obsEng),
                                   m_obsTime(obsTime),
                                   m_rot(rot),
                                   m_theta(theta),
                                   m_cos_theta(std::cos(theta)),
                                   m_sin_theta(sin_theta) { }
        double eval(double phi);
    protected:
        const GModelSpatial*   m_model;      //!< Spatial model
        const GEnergy&         m_obsEng;     //!< True photon energy
        const GTime&           m_obsTime;    //!< True photon arrival time
        const GMatrix&         m_rot;        //!< Rotation matrix
        const double&          m_theta;      //!< Offset angle (radians)
        double                 m_cos_theta;  //!< Cosine of offset angle
        const double&          m_sin_theta;  //!< Sine of offset angle
    };

    // Proteced data members
    GModelSpatial* m_spatial;       //!< spatial model
    GModelSpectral*  m_spectral;     //!< Spectral model
    GModelTemporal*  m_temporal;     //!< Temporal model
    GMatrix m_rot;             //!!< Rotation matrix from model system to skydir
};

#endif /* GMODELSPATIAL_HPP */
