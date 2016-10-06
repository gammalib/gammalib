/***************************************************************************
 *    GModelSpectralMultiplicative.hpp - Spectral power law model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Michael Mayer                               *
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
 * @file GModelSpectralMultiplicative.hpp
 * @brief Multiplicative spectral model class interface definition
 * @author Michael Mayer
 */

#ifndef GMODELSPECTRALMULTIPLICATIVE_HPP
#define GMODELSPECTRALMULTIPLICATIVE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpectral.hpp"
#include "GModelPar.hpp"
#include "GEnergy.hpp"
#include "GModelSpectralNodes.hpp"
#include "GFunction.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GTime;
class GXmlElement;



/***********************************************************************//**
 * @class GModelSpectralMultiplicative
 *
 * @brief Multiplicative spectral model class
 *
 * This class implements a Multiplicative spectrum. The model is defined by the sum
 * of several individual spectral models. Each spectral model can be added via
 * the XML interface or using the append() method.
 *
 ***************************************************************************/
class GModelSpectralMultiplicative : public GModelSpectral {

public:
    // Constructors and destructors
	GModelSpectralMultiplicative(void);
    explicit GModelSpectralMultiplicative(const GXmlElement& xml);
    GModelSpectralMultiplicative(const GModelSpectralMultiplicative& model);
    virtual ~GModelSpectralMultiplicative(void);

    // Operators
    virtual GModelSpectralMultiplicative& operator=(const GModelSpectralMultiplicative& model);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GModelSpectralMultiplicative* clone(void) const;
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
    void                  append(const GModelSpectral& spec, const std::string& name="");
    int                   components(void) const;
    const GModelSpectral* component(const int& index) const;
    const GModelSpectral* component(const std::string& name) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralMultiplicative& model);
    void free_members(void);
    void add_component(const GXmlElement& spec);
    void update_mc_cache(const GEnergy& emin, const GEnergy& emax) const;

    // Class to determine to the integral photon flux
	class flux_kern : public GFunction {
	public:
		// Constructor
		flux_kern(std::vector<GModelSpectral*> spec) :
				  m_spec(spec){};

		// Method
		double eval(const double& x) {
			GEnergy energy(x, "MeV");
			double value = m_spec[0]->eval(energy);
			for(int i=1;i<m_spec.size();++i){
				value*=m_spec[i]->eval(energy);
			}
			return value;
		}

		// Members
	protected:
		std::vector<GModelSpectral*> m_spec; //!< Spectral models
	};

	// Class to determine the integral energy flux, derviation of flux_kern
	class eflux_kern : public flux_kern {
	public:
		// Constructor
		eflux_kern(std::vector<GModelSpectral*> spec) :
				   flux_kern(spec) {};

		// Method
		double eval(const double& x) {
			return x * flux_kern::eval(x);
		}
	};

    // Protected members
    std::string                  m_type;        //!< Model type
    std::vector<GModelSpectral*> m_spectral;    //!< Container of spectral models
    std::vector<std::string>     m_components;  //!< Names of components

    // MC cache
    mutable GModelSpectralNodes  m_mc_spectrum; //!< MC spectrum cache
    mutable int                  m_n_mc_nodes;  //!< Probailities of individual components
    mutable GEnergy              m_mc_emin;     //!< Last minimum energy
    mutable GEnergy              m_mc_emax;     //!< Last maximum energy

};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpectralMultiplicative").
 ***************************************************************************/
inline
std::string GModelSpectralMultiplicative::classname(void) const
{
    return ("GModelSpectralMultiplicative");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the spectral Multiplicative model.
 ***************************************************************************/
inline
std::string GModelSpectralMultiplicative::type(void) const
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
int GModelSpectralMultiplicative::components(void) const
{
    return (m_spectral.size());
}


#endif /* GMODELSPECTRALMULTIPLICATIVE_HPP */
