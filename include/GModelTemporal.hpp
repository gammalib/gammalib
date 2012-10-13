/***************************************************************************
 *        GModelTemporal.hpp  -  Abstract temporal model base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2012 by Juergen Knoedlseder                         *
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
 * @file GModelTemporal.hpp
 * @brief Abstract temporal model base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELTEMPORAL_HPP
#define GMODELTEMPORAL_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include <iostream>
#include "GBase.hpp"
#include "GLog.hpp"
#include "GModelPar.hpp"
#include "GTime.hpp"
#include "GRan.hpp"


/***********************************************************************//**
 * @class GModelTemporal
 *
 * @brief Abstract temporal model base class
 *
 * This class implements the temporal component of the factorised model.
 * The temporal component of the factorised model is supposed to describe the
 * relative variation of the source flux with respect to the mean value
 * that is given by the spectral component. Normally, this model will have
 * a mean value of 1.
 ***************************************************************************/
class GModelTemporal : public GBase {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GModelTemporal& model);
    friend GLog&         operator<<(GLog& log,        const GModelTemporal& model);

public:
    // Constructors and destructors
    GModelTemporal(void);
    GModelTemporal(const GModelTemporal& model);
    virtual ~GModelTemporal(void);

    // Operators
    virtual GModelTemporal&  operator=(const GModelTemporal& model);
    virtual GModelPar&       operator[](const int& index);
    virtual const GModelPar& operator[](const int& index) const;
    virtual GModelPar&       operator[](const std::string& name);
    virtual const GModelPar& operator[](const std::string& name) const;

    // Virtual methods
    virtual void            clear(void) = 0;
    virtual GModelTemporal* clone(void) const = 0;
    virtual std::string     type(void) const = 0;
    virtual double          eval(const GTime& srcTime) const = 0;
    virtual double          eval_gradients(const GTime& srcTime) const = 0;
    virtual GTimes          mc(const double& rate, const GTime& tmin,
                               const GTime& tmax, GRan& ran) const = 0;
    virtual void            read(const GXmlElement& xml) = 0;
    virtual void            write(GXmlElement& xml) const = 0;
    virtual std::string     print(void) const = 0;

    // Methods
    int size(void) const { return m_pars.size(); }
    
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelTemporal& model);
    void free_members(void);

    // Proteced members
    std::vector<GModelPar*> m_pars;  //!< Parameter pointers
};

#endif /* GMODELTEMPORAL_HPP */
