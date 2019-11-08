/***************************************************************************
 *              GCTAResponseCache.i - CTA response cache class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018-2019 by Juergen Knoedlseder                         *
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
 * @file GCTAResponseCache.i
 * @brief CTA response cache class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAResponseCache.hpp"
%}

/* __ Typemaps ___________________________________________________________ */
/*
 * Typemap so that the contains() method returns two values, first the
 * flag indicating whether the value exists or not, and second the
 * value, i.e.
 *
 * flag, value = cache.contains('Crab', GEnergy(1.0,'TeV'), GEnergy(1.0,'TeV'))
 *
 */
%include "typemaps.i"
%apply double *OUTPUT { double* value };


/***********************************************************************//**
 * @class GCTAResponseCache
 *
 * @brief CTA response cache class
 ***************************************************************************/
class GCTAResponseCache : public GBase {

public:
    // Constructors and destructors
    GCTAResponseCache(void);
    GCTAResponseCache(const GCTAResponseCache& cache);
    virtual ~GCTAResponseCache(void);

    // Methods
    void               clear(void);
    GCTAResponseCache* clone(void) const;
    std::string        classname(void) const;
    int                size(void) const;
    bool               is_empty(void) const;
    void               set(const std::string& name,
                           const GEnergy&     ereco,
                           const GEnergy&     etrue,
                           const double&      value);
    void               set(const std::string& name,
                           const GCTAInstDir& dir,
                           const GEnergy&     ereco,
                           const GEnergy&     etrue,
                           const double&      value);
    void               remove(const std::string& name);
    bool               contains(const std::string& name,
                                const GEnergy&     ereco,
                                const GEnergy&     etrue,
                                double*            value) const;
    bool               contains(const std::string& name,
                                const GCTAInstDir& dir,
                                const GEnergy&     ereco,
                                const GEnergy&     etrue,
                                double*            value) const;
};


/***********************************************************************//**
 * @brief GCTAResponseCache class extension
 ***************************************************************************/
%extend GCTAResponseCache {
    GCTAResponseCache copy() {
        return (*self);
    }
};
