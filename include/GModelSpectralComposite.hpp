/***************************************************************************
 *         GModelSpectralComposite.hpp - Spectral composite model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Michael Mayer                                    *
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
 * @file GModelSpectralComposite.hpp
 * @brief Composite spectral model class interface definition
 * @author Michael Mayer
 */

#ifndef GMODELSPECTRALCOMPOSITE_HPP
#define GMODELSPECTRALCOMPOSITE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpectral.hpp"
#include "GModelPar.hpp"
#include "GEnergy.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GTime;
class GXmlElement;


/***********************************************************************//**
 * @class GModelSpectralComposite
 *
 * @brief Composite spectral model class
 *
 * This class implements a composite spectrum. The model is defined by the sum
 * of several individual spectral models. Each spectral model can be added via
 * the XML interface or using the append() method.
 *
 ***************************************************************************/
class GModelSpectralComposite : public GModelSpectral {

public:
    // Constructors and destructors
	GModelSpectralComposite(void);
    explicit GModelSpectralComposite(const GXmlElement& xml);
    GModelSpectralComposite(const GModelSpectralComposite& model);
    virtual ~GModelSpectralComposite(void);

    // Operators
    virtual GModelSpectralComposite& operator=(const GModelSpectralComposite& model);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GModelSpectralComposite* clone(void) const;
    virtual std::string              classname(void) const;
    virtual std::string              type(void) const;
    virtual double                   eval(const GEnergy& srcEng,
                                         const GTime&   srcTime = GTime(),
                                         const bool&    gradients = false) const;
    virtual double                   flux(const GEnergy& emin,
                                         const GEnergy& emax) const;
    virtual double                   eflux(const GEnergy& emin,
                                         const GEnergy& emax) const;
    virtual GEnergy                  mc(const GEnergy& emin,
                                        const GEnergy& emax,
                                        const GTime&   time,
                                        GRan&          ran) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;
    virtual std::string              print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void                  append(const GModelSpectral& spec,
                                 const std::string&    name="");
    int                   components(void) const;
    const GModelSpectral* component(const int& index) const;
    const GModelSpectral* component(const std::string& name) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralComposite& model);
    void free_members(void);
    void update_mc_cache(const GEnergy& emin, const GEnergy& emax) const;

    // Protected members
    std::string                  m_type;        //!< Model type
    std::vector<GModelSpectral*> m_spectral;    //!< Container of spectral models
    std::vector<std::string>     m_components;  //!< Names of components

    // MC cache
    mutable double               m_mc_flux;     //!< Flux cache
    mutable std::vector<double>  m_mc_probs;    //!< Probailities of individual components
    mutable GEnergy              m_mc_emin;     //!< Last minimum energy
    mutable GEnergy              m_mc_emax;     //!< Last maximum energy
    mutable std::vector<double>  m_mc_values;   //!< Parameter values
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralComposite").
 ***************************************************************************/
inline
std::string GModelSpectralComposite::classname(void) const
{
    return ("GModelSpectralComposite");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the spectral composite model.
 ***************************************************************************/
inline
std::string GModelSpectralComposite::type(void) const
{
    return (m_type);
}


/***********************************************************************//**
 * @brief Return number of spectral components
 *
 * @return Number of model components.
 *
 * Returns the number of spectral components.
 ***************************************************************************/
inline
int GModelSpectralComposite::components(void) const
{
    return (m_spectral.size());
}

#endif /* GMODELSPECTRALCOMPOSITE_HPP */
