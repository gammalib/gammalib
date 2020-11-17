/***************************************************************************
 * GCTAModelSpatialMultiplicative.hpp - Multiplicative spatial model class *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018-2020 by Juergen Knoedlseder                         *
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
 * @file GCTAModelSpatialMultiplicative.hpp
 * @brief Multiplicative spatial model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELSPATIALMULTIPLICATIVE_HPP
#define GCTAMODELSPATIALMULTIPLICATIVE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GCTAModelSpatial.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GXmlElement;
class GCTAInstDir;
class GEnergy;
class GTime;
class GCTAInstDir;
class GCTAObservation;


/***********************************************************************//**
 * @class GCTAModelSpatialMultiplicative
 *
 * @brief Multiplicative spatial model class
 *
 * This class implements a multiplicative spatial component of the CTA
 * background model.
 ***************************************************************************/
class GCTAModelSpatialMultiplicative : public GCTAModelSpatial {

public:
    // Constructors and destructors
    GCTAModelSpatialMultiplicative(void);
    explicit GCTAModelSpatialMultiplicative(const GXmlElement& xml);
    GCTAModelSpatialMultiplicative(const GCTAModelSpatialMultiplicative& model);
    virtual ~GCTAModelSpatialMultiplicative(void);

    // Operators
    virtual GCTAModelSpatialMultiplicative& operator=(const GCTAModelSpatialMultiplicative& model);

    // Pure virtual methods
    virtual void                            clear(void);
    virtual GCTAModelSpatialMultiplicative* clone(void) const;
    virtual std::string                     classname(void) const;
    virtual std::string                     type(void) const;
    virtual double                          eval(const GCTAInstDir& dir,
                                                 const GEnergy&     energy,
                                                 const GTime&       time,
                                                 const bool&        gradients = false) const;
    virtual double                          mc_max_value(const GCTAObservation& obs) const;
    virtual void                            read(const GXmlElement& xml);
    virtual void                            write(GXmlElement& xml) const;
    virtual std::string                     print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void                    append(const GCTAModelSpatial& spatial,
                                   const std::string&      name="");
    int                     components(void) const;
    const GCTAModelSpatial* component(const int& index) const;
    const GCTAModelSpatial* component(const std::string& name) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelSpatialMultiplicative& model);
    void free_members(void);

    // Protected members
    std::string                    m_type;       //!< Model type
    std::vector<GCTAModelSpatial*> m_spatial;    //!< Container of spatial models
    std::vector<std::string>       m_components; //!< Names of components
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAModelSpatialMultiplicative").
 ***************************************************************************/
inline
std::string GCTAModelSpatialMultiplicative::classname(void) const
{
    return ("GCTAModelSpatialMultiplicative");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the spatial multiplicative model.
 ***************************************************************************/
inline
std::string GCTAModelSpatialMultiplicative::type(void) const
{
    return (m_type);
}


/***********************************************************************//**
 * @brief Return number of spatial components
 *
 * @return Number of model components.
 *
 * Returns the number of spatial components.
 ***************************************************************************/
inline
int GCTAModelSpatialMultiplicative::components(void) const
{
    return ((int)m_spatial.size());
}

#endif /* GCTAMODELSPATIALMULTIPLICATIVE_HPP */
