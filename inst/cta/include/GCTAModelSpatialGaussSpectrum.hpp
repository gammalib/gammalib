/***************************************************************************
 *  GCTAModelSpatialGaussSpectrum.hpp - Spatial energy dependent Gaussian  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019 by Juergen Knoedlseder                              *
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
 * @file GCTAModelSpatialGaussSpectrum.hpp
 * @brief Spatial energy dependent Gaussian interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELSPATIALGAUSSSPECTRUM_HPP
#define GCTAMODELSPATIALGAUSSSPECTRUM_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpectral.hpp"
#include "GCTAModelSpatial.hpp"

/* __ Forward declarations _______________________________________________ */
class GXmlElement;
class GCTAInstDir;
class GEnergy;
class GTime;
class GCTAObservation;


/***********************************************************************//**
 * @class GCTAModelSpatialGaussSpectrum
 *
 * @brief Spatial energy dependent Gaussian model class
 ***************************************************************************/
class GCTAModelSpatialGaussSpectrum  : public GCTAModelSpatial {

public:
    // Constructors and destructors
    GCTAModelSpatialGaussSpectrum(void);
    explicit GCTAModelSpatialGaussSpectrum(const double& sigma);
    explicit GCTAModelSpatialGaussSpectrum(const GModelSpectral& sigma);
    explicit GCTAModelSpatialGaussSpectrum(const GXmlElement& xml);
    GCTAModelSpatialGaussSpectrum(const GCTAModelSpatialGaussSpectrum& model);
    virtual ~GCTAModelSpatialGaussSpectrum(void);

    // Operators
    virtual GCTAModelSpatialGaussSpectrum& operator=(const GCTAModelSpatialGaussSpectrum& model);

    // Implemented pure virtual methods
    virtual void                           clear(void);
    virtual GCTAModelSpatialGaussSpectrum* clone(void) const;
    virtual std::string                    classname(void) const;
    virtual std::string                    type(void) const;
    virtual double                         eval(const GCTAInstDir& dir,
                                                const GEnergy&     energy,
                                                const GTime&       time,
                                                const bool&        gradients = false) const;
    virtual double                         mc_max_value(const GCTAObservation& obs) const;
    virtual void                           read(const GXmlElement& xml);
    virtual void                           write(GXmlElement& xml) const;
    virtual std::string                    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const GModelSpectral* sigma(void) const;
    void                  sigma(const double& sigma);
    void                  sigma(const GModelSpectral& sigma);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelSpatialGaussSpectrum& model);
    void free_members(void);
    void set_pointers(void);


    // Protected members
    GModelSpectral* m_sigma;  // Pointer to sigma spectrum
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAModelSpatialGaussSpectrum").
 ***************************************************************************/
inline
std::string GCTAModelSpatialGaussSpectrum::classname(void) const
{
    return ("GCTAModelSpatialGaussSpectrum");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type "EnergyDependentGaussian".
 ***************************************************************************/
inline
std::string GCTAModelSpatialGaussSpectrum::type(void) const
{
    return ("EnergyDependentGaussian");
}


/***********************************************************************//**
 * @brief Return maximum function value for Monte Carlo simulations
 *
 * @param[in] obs CTA Observation.
 * @return Maximum function value for Monte Carlo simulations.
 *
 * This method always returns 1.
 ***************************************************************************/
inline
double GCTAModelSpatialGaussSpectrum::mc_max_value(const GCTAObservation& obs) const
{
    return 1.0;
}


/***********************************************************************//**
 * @brief Return pointer to sigma spectrum
 *
 * @return Pointer to sigma spectrum.
 ***************************************************************************/
inline
const GModelSpectral* GCTAModelSpatialGaussSpectrum::sigma(void) const
{
    return (m_sigma);
}

#endif /* GCTAMODELSPATIALGAUSSSPECTRUM_HPP */
