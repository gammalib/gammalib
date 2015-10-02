/***************************************************************************
 *       GCTAModelAeffBackground.hpp - CTA Aeff background model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Michael Mayer                                    *
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
 * @file GCTAModelAeffBackground.hpp
 * @brief CTA Aeff background model class definition
 * @author Michael Mayer
 */

#ifndef GCTAMODELAEFFBACKGROUND_HPP
#define GCTAMODELAEFFBACKGROUND_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelData.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GFunction.hpp"
#include "GCTAEventList.hpp"
#include "GCTAAeff.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GCTAModelIrfBackground
 *
 * @brief CTA IRF background model class
 ***************************************************************************/
class GCTAModelAeffBackground : public GModelData {

public:
    // Constructors and destructors
    GCTAModelAeffBackground(void);
    explicit GCTAModelAeffBackground(const GXmlElement& xml);
    explicit GCTAModelAeffBackground(const GModelSpectral& spectral);
    GCTAModelAeffBackground(const GCTAModelAeffBackground& bgd);
    virtual ~GCTAModelAeffBackground(void);

    // Operators
    GCTAModelAeffBackground& operator=(const GCTAModelAeffBackground& bgd);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GCTAModelAeffBackground* clone(void) const;
    virtual std::string              classname(void) const;
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
    GModelSpectral* spectral(void) const;
    GModelTemporal* temporal(void) const;

private:
    // Methods
    void            init_members(void);
    void            copy_members(const GCTAModelAeffBackground& bgd);
    void            free_members(void);
    void            set_pointers(void);
    bool            valid_model(void) const;
    GModelSpectral* xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal* xml_temporal(const GXmlElement& temporal) const;
    double          aeff_integral(const GObservation& obs, const double& logE) const;

    // ROI integration kernel over theta
    class npred_roi_kern_theta : public GFunction {
    public:
        npred_roi_kern_theta(const GCTAAeff* aeff,
                             const double&   logE,
                             const int&      iter) :
                             m_aeff(aeff),
                             m_logE(logE),
                             m_iter(iter) { }
        double eval(const double& theta);
    protected:
        const GCTAAeff* m_aeff;  //!< Pointer to effectve area
        const double&   m_logE;  //!< Log10 of energy
        const int&      m_iter;  //!< Romberg iterations
    };

    // ROI integration kernel over phi
    class npred_roi_kern_phi : public GFunction {
    public:
        npred_roi_kern_phi(const GCTAAeff* aeff,
                           const double&   logE,
                           const double&   theta) :
                           m_aeff(aeff),
                           m_logE(logE),
                           m_theta(theta) { }
        double eval(const double& phi);
    protected:
        const GCTAAeff* m_aeff;   //!< Pointer to effective area
        const double&   m_logE;   //!< Log10 of energy
        const double&   m_theta;  //!< Offset angle (radians)
    };

    // Members
    GModelSpectral* m_spectral;      //!< Spectral model
    GModelTemporal* m_temporal;      //!< Temporal model
    int             m_n_mc_energies; //!< Energy sampling for MC spectrum

    // Npred cache
    mutable std::vector<std::string> m_npred_names;    //!< Model names
    mutable std::vector<GEnergy>     m_npred_energies; //!< Model energy
    mutable std::vector<GTime>       m_npred_times;    //!< Model time
    mutable std::vector<double>      m_npred_values;   //!< Model values
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAModelAeffBackground").
 ***************************************************************************/
inline
std::string GCTAModelAeffBackground::classname(void) const
{
    return ("GCTAModelAeffBackground");
}


/***********************************************************************//**
 * @brief Return data model type
 *
 * @return Data model type.
 *
 * Returns the type of the data model.
 ***************************************************************************/
inline
std::string GCTAModelAeffBackground::type(void) const
{
    return ("CTAAeffBackground");
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
bool GCTAModelAeffBackground::is_constant(void) const
{
    return (m_temporal != NULL && m_temporal->type() == "Constant");
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
GModelSpectral* GCTAModelAeffBackground::spectral(void) const
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
GModelTemporal* GCTAModelAeffBackground::temporal(void) const
{
    return (m_temporal);
}

#endif /* GCTAMODELAEFFBACKGROUND_HPP */
