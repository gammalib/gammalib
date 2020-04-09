/***************************************************************************
 *                 GResponseCache.i - Response cache class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GResponseCache.i
 * @brief Response cache class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GResponseCache.hpp"
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
 * @class GResponseCache
 *
 * @brief Response cache class
 ***************************************************************************/
class GResponseCache : public GBase {

public:
    // Constructors and destructors
    GResponseCache(void);
    GResponseCache(const GResponseCache& cache);
    virtual ~GResponseCache(void);

    // Methods
    void            clear(void);
    GResponseCache* clone(void) const;
    std::string     classname(void) const;
    int             size(void) const;
    bool            is_empty(void) const;
    void            set(const std::string& name,
                        const GEnergy&     ereco,
                        const GEnergy&     etrue,
                        const double&      value);
    void            set(const std::string& name,
                        const GInstDir&    dir,
                        const GEnergy&     ereco,
                        const GEnergy&     etrue,
                        const double&      value);
    void            remove(const std::string& name);
    bool            contains(const std::string& name,
                             const GEnergy&     ereco,
                             const GEnergy&     etrue,
                             double*            value) const;
    bool            contains(const std::string& name,
                             const GInstDir&    dir,
                             const GEnergy&     ereco,
                             const GEnergy&     etrue,
                             double*            value) const;
};


/***********************************************************************//**
 * @brief GResponseCache class extension
 ***************************************************************************/
%extend GResponseCache {
    GResponseCache copy() {
        return (*self);
    }
};
