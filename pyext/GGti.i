/***************************************************************************
 *                    GGti.i - Good time interval class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GGti.i
 * @brief Good time interval class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GGti.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GGti
 *
 * @brief Good Time Interval class
 ***************************************************************************/
class GGti : public GBase {
public:
    // Constructors and destructors
    GGti(void);
    GGti(const GGti& gti);
    explicit GGti(const GTimeReference& ref);
    virtual ~GGti(void);

    // Methods
    void                  clear(void);
    GGti*                 clone(void) const;
    int                   size(void) const;
    bool                  isempty(void) const;
    void                  append(const GTime& tstart, const GTime& tstop);
    void                  insert(const GTime& tstart, const GTime& tstop);
    void                  merge(void);
    void                  merge(const GTime& tstart, const GTime& tstop);
    void                  reduce(const GTime& tstart, const GTime& tstop);
    void                  remove(const int& index);
    void                  reserve(const int& num);
    void                  extend(const GGti& gti);
    void                  load(const std::string& filename,
                               const std::string& extname = "GTI");
    void                  save(const std::string& filename, bool clobber,
                               const std::string& extname = "GTI") const;
    void                  read(const GFitsTable* hdu);
    void                  write(GFits* file,
                                const std::string& extname = "GTI") const;
    const GTime&          tstart(void) const;
    const GTime&          tstop(void) const;
    const GTime&          tstart(const int& index) const;
    const GTime&          tstop(const int& index) const;
    const double&         telapse(void) const;
    const double&         ontime(void) const;
    void                  reference(const GTimeReference& ref);
    const GTimeReference& reference(void) const;
    bool                  contains(const GTime& time) const;
};


/***********************************************************************//**
 * @brief GGti class extension
 ***************************************************************************/
%extend GGti {
    char *__str__() {
        return tochar(self->print());
    }
    GGti copy() {
        return (*self);
    }
};
