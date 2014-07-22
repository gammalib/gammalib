/***************************************************************************
 *         GModelSpectral.hpp - Abstract spectral model base class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
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
 * @file GModelSpectral.hpp
 * @brief Abstract spectral model base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPECTRAL_HPP
#define GMODELSPECTRAL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GModelPar.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GRan.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpectral
 *
 * @brief Abstract spectral model base class
 *
 * This class implements the spectral component of the factorized source
 * model
 *
 * \f[
 *    S_{\rm E}(E | t)
 * \f]
 *
 * where
 * - \f$E\f$ is the true photon energy, and
 * - \f$t\f$ is the true photon arrival time.
 *
 * The spectral component describes the spatially integrated time dependent
 * spectral distribution of the source. It satisfies
 * \f[
 *    \int_{E} S_{\rm E}(E | t) dE = \Phi
 * \f]
 * for all \f$t\f$, where \f$\Phi\f$ is the spatially and spectrally
 * integrated total source flux. The spectral component does not impact
 * the temporal properties of the integrated flux \f$\Phi\f$.
 ***************************************************************************/
class GModelSpectral : public GBase {

public:
    // Constructors and destructors
    GModelSpectral(void);
    GModelSpectral(const GModelSpectral& model);
    virtual ~GModelSpectral(void);

    // Operators
    virtual GModelSpectral&  operator=(const GModelSpectral& model);
    virtual GModelPar&       operator[](const int& index);
    virtual const GModelPar& operator[](const int& index) const;
    virtual GModelPar&       operator[](const std::string& name);
    virtual const GModelPar& operator[](const std::string& name) const;

    // Pure virtual methods
    virtual void            clear(void) = 0;
    virtual GModelSpectral* clone(void) const = 0;
    virtual std::string     type(void) const = 0;
    virtual double          eval(const GEnergy& srcEng,
                                 const GTime& srcTime) const = 0;
    virtual double          eval_gradients(const GEnergy& srcEng,
                                           const GTime& srcTime) = 0;
    virtual double          flux(const GEnergy& emin,
                                 const GEnergy& emax) const = 0;
    virtual double          eflux(const GEnergy& emin,
                                  const GEnergy& emax) const = 0;
    virtual GEnergy         mc(const GEnergy& emin, const GEnergy& emax,
                               const GTime& time, GRan& ran) const = 0;
    virtual void            read(const GXmlElement& xml) = 0;
    virtual void            write(GXmlElement& xml) const = 0;
    virtual std::string     print(const GChatter& chatter = NORMAL) const = 0;

    // Methods
    GModelPar&       at(const int& index);
    const GModelPar& at(const int& index) const;
    int              size(void) const;
    void             autoscale(void);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectral& model);
    void free_members(void);

    // Proteced members
    std::vector<GModelPar*> m_pars;  //!< Parameter pointers
};


/***********************************************************************//**
 * @brief Returns model parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Model parameter.
 *
 * Returns model parameter without @p index range checking.
 ***************************************************************************/
inline
GModelPar& GModelSpectral::operator[](const int& index)
{
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns model parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Model parameter.
 *
 * Returns model parameter without @p index range checking.
 ***************************************************************************/
inline
const GModelPar& GModelSpectral::operator[](const int& index) const
{
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Return number of parameters
 *
 * @return Number of parameters in spectral model component.
 *
 * Returns the number of parameters in the spectral model component.
 ***************************************************************************/
inline
int GModelSpectral::size(void) const
{
    return (m_pars.size());
}

#endif /* GMODELSPECTRAL_HPP */
