/***************************************************************************
 *      GCTAModelCubeBackground.hpp - CTA cube background model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Michael Mayer                               *
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
 * @file GCTAModelCubeBackground.hpp
 * @brief CTA cube background model class interface definition
 * @author Michael Mayer
 */

#ifndef GCTAMODELCUBEBACKGROUND_HPP
#define GCTAMODELCUBEBACKGROUND_HPP

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
#include "GCTAObservation.hpp"
#include "GCTAEventList.hpp"
#include "GModelSpatial.hpp"


/***********************************************************************//**
 * @class GCTAModelCubeBackground
 *
 * @brief CTA cube background model class
 *
 * This class implements a cube background model for CTA.
 ***************************************************************************/
class GCTAModelCubeBackground : public GModelData {

public:
    // Constructors and destructors
    GCTAModelCubeBackground(void);
    explicit GCTAModelCubeBackground(const GXmlElement& xml);
    GCTAModelCubeBackground(const GModelSpatial& spatial,
                            const GModelSpectral& spectral);
    GCTAModelCubeBackground(const GCTAModelCubeBackground& model);
    GCTAModelCubeBackground(const GCTAObservation& obs,
                            const std::string&     filename,
                            const GModelSpectral&  spectral,
                            const int&             nx_sky   = 0,
                            const int&             ny_sky   = 0,
                            const int&             n_energy = 0);
    virtual ~GCTAModelCubeBackground(void);

    // Operators
    virtual GCTAModelCubeBackground& operator=(const GCTAModelCubeBackground& model);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GCTAModelCubeBackground* clone(void) const;
    virtual std::string              type(void) const;
    virtual bool                     is_constant(void) const;
    virtual double                   eval(const GEvent& event,
                                          const GObservation& obs) const;
    virtual double                   eval_gradients(const GEvent& event,
                                                    const GObservation& obs) const;
    virtual double                   npred(const GEnergy& obsEng,
                                           const GTime& obsTime,
                                           const GObservation& obs) const;
    virtual GCTAEventList*           mc(const GObservation& obs, GRan& ran) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;
    virtual std::string              print(const GChatter& chatter = NORMAL) const;

    // Other methods
    GModelSpatial*  spatial(void) const;
    GModelSpectral* spectral(void) const;
    GModelTemporal* temporal(void) const;

protected:
    // Protected methods
    void            set_spatial(const GCTAObservation& obs,
                                const std::string& filename,
                                const int& nx_sky   = 0,
                                const int& ny_sky   = 0,
                                const int& n_energy = 0);
    void            init_members(void);
    void            copy_members(const GCTAModelCubeBackground& model);
    void            free_members(void);
    void            set_pointers(void);
    bool            valid_model(void) const;
    GModelSpatial*  xml_spatial(const GXmlElement& spatial) const;
    GModelSpectral* xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal* xml_temporal(const GXmlElement& temporal) const;

    class npred_roi_kern_theta : public GFunction {
    public:
    	npred_roi_kern_theta(const GModelSpatial* model,
                             const GEnergy&       obsEng,
                             const GTime&         obsTime,
                             const GMatrix&       rot,
                             const double&        roi,
                             const double&        dist,
                             const double&        omega0) :
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
        double eval(const double& theta);
    protected:
        const GModelSpatial* m_model;    //!< Spatial model
        const GEnergy&       m_obsEng;   //!< True photon energy
        const GTime&         m_obsTime;  //!< True photon arrival time
        const GMatrix&       m_rot;      //!< Rotation matrix
        const double&        m_roi;      //!< ROI radius in radians
		double               m_cosroi;   //!< Cosine of ROI radius
		const double&        m_dist;     //!< Distance between pointing and ROI centre in radians
		double               m_cosdist;  //!< Cosine of distance
		double               m_sindist;  //!< Sinus of distance
	    const double&        m_omega0;   //!< Position angle of ROI
    };

    class npred_roi_kern_phi : public GFunction {
    public:
    	npred_roi_kern_phi(const GModelSpatial* model,
                           const GEnergy&       obsEng,
                           const GTime&         obsTime,
                           const GMatrix&       rot,
                           const double&        theta,
                           const double&        sin_theta) :
                           m_model(model),
                           m_obsEng(obsEng),
                           m_obsTime(obsTime),
                           m_rot(rot),
                           m_theta(theta),
                           m_cos_theta(std::cos(theta)),
                           m_sin_theta(sin_theta) { }
        double eval(const double& phi);
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
    GModelSpatial*  m_spatial;   //!< Spatial model
    GModelSpectral* m_spectral;  //!< Spectral model
    GModelTemporal* m_temporal;  //!< Temporal model
    GMatrix         m_rot;       //!< Rotation matrix from model system to skydir

    // Npred cache
    mutable std::vector<std::string> m_npred_names;    //!< Model names
    mutable std::vector<GEnergy>     m_npred_energies; //!< Model energy
    mutable std::vector<GTime>       m_npred_times;    //!< Model time
    mutable std::vector<double>      m_npred_values;   //!< Model values
};


/***********************************************************************//**
 * @brief Return data model type
 *
 * @return Data model type.
 *
 * Returns the type of the data model.
 ***************************************************************************/
inline
std::string GCTAModelCubeBackground::type(void) const
{
    return ("CTACubeBackground");
}


/***********************************************************************//**
 * @brief Signals if sky model is temporally constant
 *
 * @return True if sky model is temporally constant, false otherwise.
 *
 * Signals if the sky model is temporally constant. A temporally constant
 * model is a model that has a temporal component of type "Constant".
 ***************************************************************************/
inline
bool GCTAModelCubeBackground::is_constant(void) const
{
    return (m_temporal != NULL && m_temporal->type() == "Constant");
}


/***********************************************************************//**
 * @brief Return spatial model component
 *
 * @return Pointer to spatial model component.
 *
 * Returns a pointer to the spatial model component of the model. The pointer
 * is of type GModelSpatial. Note that a NULL pointer may be returned if the
 * sky model has no spatial model component.
 ***************************************************************************/
inline
GModelSpatial* GCTAModelCubeBackground::spatial(void) const
{
    return (m_spatial);
}


/***********************************************************************//**
 * @brief Return spectral model component
 *
 * @return Pointer to spectral model component.
 *
 * Returns a pointer to the spectral model component of the model. The
 * pointer is of type GModelSpectral. Note that a NULL pointer may be
 * returned if the sky model has no spectral model component.
 ***************************************************************************/
inline
GModelSpectral* GCTAModelCubeBackground::spectral(void) const
{
    return (m_spectral);
}


/***********************************************************************//**
 * @brief Return temporal model component
 *
 * @return Pointer to temporal model component.
 *
 * Returns a pointer to the temporal model component of the model. The
 * pointer is of type GModelTemporal. Note that a NULL pointer may be
 * returned if the sky model has no temporal model component.
 ***************************************************************************/
inline
GModelTemporal* GCTAModelCubeBackground::temporal(void) const
{
    return (m_temporal);
}

#endif /* GCTAMODELCUBEBACKGROUND_HPP */
