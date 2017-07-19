/***************************************************************************
 *     GModelTemporalPhaseCurve.hpp - Temporal phase curve model class     *
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
 * @file GModelTemporalPhaseCurve.hpp
 * @brief Temporal phase curve model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELTEMPORALPHASECURVE_HPP
#define GMODELTEMPORALPHASECURVE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModelPar.hpp"
#include "GModelTemporal.hpp"
#include "GNodeArray.hpp"
#include "GFilename.hpp"
#include "GTime.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GXmlElement;


/***********************************************************************//**
 * @class GModelTemporalPhaseCurve
 *
 * @brief Temporal phase curve model class
 *
 * This class implements a temporal variation based on the phase
 * \f$\Phi(t)\f$ that is computed using
 *
 * \f[
 *    \Phi(t) = \Phi_0 + f(t-t_0) + \frac{1}{2}\dot{f} (t-t_0)^2 +
 *                                  \frac{1}{6}\ddot{f} (t-t_0)^3
 * \f]
 *
 * where
 * \f$t_0\f$ is a reference time,
 * \f$\Phi_0\f$ is the phase at the reference time,
 * \f$f\f$ is the variation frequency at the reference time,
 * \f$\dot{f}\f$ is the first derivative of the variation frequency at the
 * reference time, and
 * \f$\ddot{f}\f$ is the second derivative of the variation frequency at the
 * reference time.
 *
 * The temporal variation is then computed using
 *
 * \f[
 *    S_{\rm t}(t) = r(\Phi(t)) \times {\tt m\_norm}
 * \f]
 *
 * where \f$r(\Phi(t))\f$ is the rate defined by linear interpolation between
 * the nodes in a FITS file, and \f${\tt m\_norm}\f$ is a normalisation
 * constant.
 ***************************************************************************/
class GModelTemporalPhaseCurve : public GModelTemporal {

public:
    // Constructors and destructors
    GModelTemporalPhaseCurve(void);
    explicit GModelTemporalPhaseCurve(const GXmlElement& xml);
    GModelTemporalPhaseCurve(const GFilename& filename,
                             const GTime&     mjd,
                             const double&    phase,
                             const double&    f0,
                             const double&    f1,
                             const double&    f2,
                             const double&    norm = 1.0,
                             const bool&      normalize = true);
    GModelTemporalPhaseCurve(const GModelTemporalPhaseCurve& model);
    virtual ~GModelTemporalPhaseCurve(void);

    // Operators
    virtual GModelTemporalPhaseCurve& operator=(const GModelTemporalPhaseCurve& model);

    // Implemented virtual base class methods
    virtual void                      clear(void);
    virtual GModelTemporalPhaseCurve* clone(void) const;
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
    GTime            mjd(void) const;
    void             mjd(const GTime& time);
    double           phase(void) const;
    double           phase(const GTime& time) const;
    void             phase(const double& phase);
    double           value(const double& phase) const;
    double           f0(void) const;
    void             f0(const double& f0);
    double           f1(void) const;
    void             f1(const double& f1);
    double           f2(void) const;
    void             f2(const double& f2);
    bool             normalize(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelTemporalPhaseCurve& model);
    void free_members(void);
    void load_nodes(const GFilename& filename);
    void normalize_nodes(void);

    // Protected parameter members
    GModelPar           m_norm;          //!< Normalization factor
    GModelPar           m_mjd;           //!< Reference MJD
    GModelPar           m_phase;         //!< Phase at reference MJD
    GModelPar           m_f0;            //!< Frequency at reference MJD
    GModelPar           m_f1;            //!< First freq. derivative at reference MJD
    GModelPar           m_f2;            //!< Second freq. derivative at reference MJD
    bool                m_normalize;     //!< Normalize phase curve (default: true)
    bool                m_has_normalize; //!< XML has normalize attribute

    // Protected members
    mutable GNodeArray  m_nodes;    //!< Phase values of nodes
    std::vector<double> m_values;   //!< Function values at nodes
    GFilename           m_filename; //!< Name of phase file function
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelTemporalPhaseCurve").
 ***************************************************************************/
inline
std::string GModelTemporalPhaseCurve::classname(void) const
{
    return ("GModelTemporalPhaseCurve");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "PhaseCurve".
 *
 * Returns the type of the temporal model.
 ***************************************************************************/
inline
std::string GModelTemporalPhaseCurve::type(void) const
{
    return ("PhaseCurve");
}


/***********************************************************************//**
 * @brief Return file name
 *
 * @return Name of node file.
 *
 * Returns the name of the file function node file.
 ***************************************************************************/
inline
const GFilename& GModelTemporalPhaseCurve::filename(void) const
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
void GModelTemporalPhaseCurve::filename(const GFilename& filename)
{
    load_nodes(filename);
    return;
}


/***********************************************************************//**
 * @brief Return normalization factor
 *
 * @return Normalization factor.
 *
 * Returns the normalization factor.
 ***************************************************************************/
inline
double GModelTemporalPhaseCurve::norm(void) const
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
void GModelTemporalPhaseCurve::norm(const double& norm)
{
    m_norm.value(norm);
    return;
}


/***********************************************************************//**
 * @brief Return reference Modified Julian Day
 *
 * @return Reference Modified Julian Day.
 *
 * Returns the reference Modified Julian Day.
 ***************************************************************************/
inline
GTime GModelTemporalPhaseCurve::mjd(void) const
{
    GTime time;
    time.mjd(m_mjd.value());
    return (time);
}


/***********************************************************************//**
 * @brief Set reference Modified Julian Day
 *
 * @param[in] mjd Reference Modified Julian Day.
 *
 * Sets the reference Modified Julian Day.
 ***************************************************************************/
inline
void GModelTemporalPhaseCurve::mjd(const GTime& mjd)
{
    m_mjd.value(mjd.mjd());
    return;
}


/***********************************************************************//**
 * @brief Return phase at reference Modified Julian Day
 *
 * @return Phase at reference Modified Julian Day.
 *
 * Returns the phase at reference Modified Julian Day.
 ***************************************************************************/
inline
double GModelTemporalPhaseCurve::phase(void) const
{
    return (m_phase.value());
}


/***********************************************************************//**
 * @brief Set phase at reference Modified Julian Day
 *
 * @param[in] phase Phase at reference Modified Julian Day.
 *
 * Sets the phase at reference Modified Julian Day.
 ***************************************************************************/
inline
void GModelTemporalPhaseCurve::phase(const double& phase)
{
    m_phase.value(phase);
    return;
}


/***********************************************************************//**
 * @brief Return frequency at reference Modified Julian Day
 *
 * @return Frequency at reference Modified Julian Day.
 *
 * Returns the frequency at reference Modified Julian Day.
 ***************************************************************************/
inline
double GModelTemporalPhaseCurve::f0(void) const
{
    return (m_f0.value());
}


/***********************************************************************//**
 * @brief Set frequency at reference Modified Julian Day
 *
 * @param[in] f0 Frequency at reference Modified Julian Day.
 *
 * Sets the frequency at reference Modified Julian Day.
 ***************************************************************************/
inline
void GModelTemporalPhaseCurve::f0(const double& f0)
{
    m_f0.value(f0);
    return;
}


/***********************************************************************//**
 * @brief Return first frequency derivative at reference Modified Julian Day
 *
 * @return First frequency derivative at reference Modified Julian Day.
 *
 * Returns the first frequency derivative at reference Modified Julian Day.
 ***************************************************************************/
inline
double GModelTemporalPhaseCurve::f1(void) const
{
    return (m_f1.value());
}


/***********************************************************************//**
 * @brief Set first frequency derivative at reference Modified Julian Day
 *
 * @param[in] f1 First frequency derivative at reference Modified Julian Day.
 *
 * Sets the first frequency derivative at reference Modified Julian Day.
 ***************************************************************************/
inline
void GModelTemporalPhaseCurve::f1(const double& f1)
{
    m_f1.value(f1);
    return;
}


/***********************************************************************//**
 * @brief Return second frequency derivative at reference Modified Julian Day
 *
 * @return Second frequency derivative at reference Modified Julian Day.
 *
 * Returns the second frequency derivative at reference Modified Julian Day.
 ***************************************************************************/
inline
double GModelTemporalPhaseCurve::f2(void) const
{
    return (m_f2.value());
}


/***********************************************************************//**
 * @brief Set second frequency derivative at reference Modified Julian Day
 *
 * @param[in] f2 Second frequency derivative at reference Modified Julian Day.
 *
 * Sets the second frequency derivative at reference Modified Julian Day.
 ***************************************************************************/
inline
void GModelTemporalPhaseCurve::f2(const double& f2)
{
    m_f2.value(f2);
    return;
}


/***********************************************************************//**
 * @brief Return normalization flag
 *
 * @return True if the phase curve has been normalized, false otherwise.
 *
 * Signals whether a phase curve has been normalized or not.
 ***************************************************************************/
inline
bool GModelTemporalPhaseCurve::normalize(void) const
{
    return (m_normalize);
}

#endif /* GMODELTEMPORALPHASECURVE_HPP */
