/***************************************************************************
 *           GResponseVectorCache.i - Response vector cache class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020-2022 by Juergen Knoedlseder                         *
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
 * @file GResponseVectorCache.i
 * @brief Response vector cache class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GResponseVectorCache.hpp"
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
%typemap(argout) GVector *OUTPUT {
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GVector, 0 |  0);
}

/***********************************************************************//**
 * @class GResponseVectorCache
 *
 * @brief Response cache class
 ***************************************************************************/
class GResponseVectorCache : public GBase {

public:
    // Constructors and destructors
    GResponseVectorCache(void);
    GResponseVectorCache(const GResponseVectorCache& cache);
    virtual ~GResponseVectorCache(void);

    // Methods
    void                  clear(void);
    GResponseVectorCache* clone(void) const;
    std::string           classname(void) const;
    bool                  is_empty(void) const;
    int                   size(void) const;
    void                  set(const std::string& cache_id,
                              const GVector&     vector);
    void                  remove(const std::string& cache_id);
    bool                  contains(const std::string& cache_id,
                                   GVector*           irfs = NULL) const;
    void                  load(const GFilename& filename);
    void                  save(const GFilename& filename,
                               const bool& clobber = false) const;
    void                  read(const GFitsTable& table);
};


/***********************************************************************//**
 * @brief GResponseVectorCache class extension
 ***************************************************************************/
%extend GResponseVectorCache {
    GResponseVectorCache copy() {
        return (*self);
    }
};
