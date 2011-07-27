/***************************************************************************
 *        GModelSpectralFunc.hpp  -  Spectral function model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GModelSpectralFunc.hpp
 * @brief Spectral function model class definition
 * @author J. Knodlseder
 */

#ifndef GMODELSPECTRALFUNC_HPP
#define GMODELSPECTRALFUNC_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"
#include "GNodeArray.hpp"


/***********************************************************************//**
 * @class GModelSpectralFunc
 *
 * @brief Spectral function model class
 *
 * This class implements an arbitrary function as the spectral component of
 * the model. The function is defined as
 * \f[I(E)=norm f(E)\f]
 * where
 * \f$norm\f$ is the normalization of the function.
 ***************************************************************************/
class GModelSpectralFunc : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralFunc(void);
    explicit GModelSpectralFunc(const std::string& filename);
    explicit GModelSpectralFunc(const GXmlElement& xml);
    GModelSpectralFunc(const GModelSpectralFunc& model);
    virtual ~GModelSpectralFunc(void);

    // Operators
    virtual GModelSpectralFunc& operator=(const GModelSpectralFunc& model);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpectralFunc* clone(void) const;
    virtual std::string         type(void) const { return "FileFunction"; }
    virtual double              eval(const GEnergy& srcEng) const;
    virtual double              eval_gradients(const GEnergy& srcEng) const;
    virtual double              flux(const GEnergy& emin, const GEnergy& emax) const;
    virtual double              eflux(const GEnergy& emin, const GEnergy& emax) const;
    virtual GEnergy             mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(void) const;

    // Other methods
    double norm(void) const { return m_norm.real_value(); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralFunc& model);
    void free_members(void);
    void load_nodes(const std::string& filename);

    // Protected members
    GModelPar           m_norm;      //!< Normalization factor
    GNodeArray          m_nodes;     //!< Nodes of function
    std::vector<double> m_values;    //!< Function values at nodes
    std::string         m_filename;  //!< Name of file function
};

#endif /* GMODELSPECTRALFUNC_HPP */
