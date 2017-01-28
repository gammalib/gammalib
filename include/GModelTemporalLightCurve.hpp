/***************************************************************************
 *     GModelTemporalLightCurve.hpp - Temporal light curve model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GModelTemporalLightCurve.hpp
 * @brief Light curve model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELTEMPORALLIGHTCURVE_HPP
#define GMODELTEMPORALLIGHTCURVE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModelPar.hpp"
#include "GModelTemporal.hpp"
#include "GNodeArray.hpp"
#include "GFilename.hpp"
#include "GTimeReference.hpp"
#include "GTime.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GXmlElement;


/***********************************************************************//**
 * @class GModelTemporalLightCurve
 *
 * @brief Light curve model class
 *
 * This class implements a light curve defined by nodes given specified in
 * a FITS file. The model is defined by
 *
 * \f[
 *    S_{\rm t}(t) = r(t) \times {\tt m\_norm}
 * \f]
 *
 * where
 * \f$r(t)$\f is the rate defined by linear interpolation between the nodes
 * in a FITS file, and
 * \f${\tt m\_norm}$\f is a normalisation constant.
 ***************************************************************************/
class GModelTemporalLightCurve : public GModelTemporal {

public:
    // Constructors and destructors
    GModelTemporalLightCurve(void);
    explicit GModelTemporalLightCurve(const GXmlElement& xml);
    GModelTemporalLightCurve(const GFilename& filename,
                             const double&    norm = 1.0);
    GModelTemporalLightCurve(const GModelTemporalLightCurve& model);
    virtual ~GModelTemporalLightCurve(void);

    // Operators
    virtual GModelTemporalLightCurve& operator=(const GModelTemporalLightCurve& model);

    // Implemented virtual base class methods
    virtual void                      clear(void);
    virtual GModelTemporalLightCurve* clone(void) const;
    virtual std::string               classname(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const GTime& srcTime,
                                           const bool& gradients = false) const;
    virtual GTimes                    mc(const double& rate, const GTime& tmin,
                                         const GTime& tmax, GRan& ran) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;
    virtual std::string               print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const GFilename& filename(void) const;
    void             filename(const GFilename& filename);
    double           norm(void) const;
    void             norm(const double& norm);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelTemporalLightCurve& model);
    void free_members(void);
    void load_nodes(const GFilename& filename);
    void mc_update(const GTime& tmin, const GTime& tmax) const;

    // Protected members
    GModelPar           m_norm;       //!< Normalization factor
    mutable GNodeArray  m_nodes;      //!< Time nodes of function
    std::vector<double> m_values;     //!< Function values at nodes
    GFilename           m_filename;   //!< Name of file function
    GTimeReference      m_timeref;    //!< Time reference
    GTime               m_tmin;       //!< Minimum time of model
    GTime               m_tmax;       //!< Maximum time of model

    // Cached members for MC
    mutable GTime               m_mc_tmin;         //!< Minimum time
    mutable GTime               m_mc_tmax;         //!< Maximum time
    mutable double              m_mc_eff_duration; //!< Effective duration
    mutable std::vector<double> m_mc_cum;          //!< Cumulative distribution
    mutable std::vector<double> m_mc_slope;        //!< Slope of interval
    mutable std::vector<double> m_mc_offset;       //!< Offset of interval
    mutable std::vector<double> m_mc_time;         //!< Start time of interval
    mutable std::vector<double> m_mc_dt;           //!< Length of interval
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelTemporalLightCurve").
 ***************************************************************************/
inline
std::string GModelTemporalLightCurve::classname(void) const
{
    return ("GModelTemporalLightCurve");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "LightCurve".
 *
 * Returns the type of the temporal model.
 ***************************************************************************/
inline
std::string GModelTemporalLightCurve::type(void) const
{
    return "LightCurve";
}


/***********************************************************************//**
 * @brief Return normalization factor
 *
 * @return Normalization factor.
 *
 * Returns the normalization factor.
 ***************************************************************************/
inline
double GModelTemporalLightCurve::norm(void) const
{
    return (m_norm.value());
}


/***********************************************************************//**
 * @brief Set normalization factor 
 *
 * @param[in] norm Normalization factor.
 *
 * Sets the normalization factor.
 ***************************************************************************/
inline
void GModelTemporalLightCurve::norm(const double& norm)
{
    m_norm.value(norm);
    return;
}


/***********************************************************************//**
 * @brief Return file name
 *
 * @return Name of node file.
 *
 * Returns the name of the file function node file.
 ***************************************************************************/
inline
const GFilename& GModelTemporalLightCurve::filename(void) const
{
    return (m_filename);
}


/***********************************************************************//**
 * @brief Loads nodes from node file and set filename
 *
 * @param[in] filename Node file name.
 *
 * Loads the nodes from a file function node file and sets the filename.
 ***************************************************************************/
inline
void GModelTemporalLightCurve::filename(const GFilename& filename)
{
    load_nodes(filename);
    return;
}

#endif /* GMODELTEMPORALLIGHTCURVE_HPP */
