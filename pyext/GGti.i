/***************************************************************************
 *                   GGti.i  -  Good time interval class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @brief Good time interval class Python interface definition
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
 * @brief Interface for the GTI class.
 *
 * This class holds a list of Good Time Intervals, i.e. time intervals that
 * are valid for science analysis. Times are stored using the GTime class.
 ***************************************************************************/
class GGti : public GBase {
public:
    // Constructors and destructors
    GGti(void);
    GGti(const GGti& gti);
    virtual ~GGti(void);

    // Methods
    void   clear(void);
    GGti*  clone(void) const;
    int    size(void) const;
    void   add(const GTime& tstart, const GTime& tstop);
    void   append(const GTime& tstart, const GTime& tstop);
    void   insert(const GTime& tstart, const GTime& tstop);
    void   reduce(const GTime& tstart, const GTime& tstop);
	void   load(const std::string& filename,
                const std::string& extname = "GTI");
	void   save(const std::string& filename, bool clobber,
                const std::string& extname = "GTI") const;
    void   read(GFitsTable* hdu);
    void   write(GFits* file, const std::string& extname = "GTI") const;
	GTime  tstart(void) const;
	GTime  tstop(void) const;
	GTime  tstart(int inx) const;
	GTime  tstop(int inx) const;
	double telapse(void) const;
	double ontime(void) const;
    double mjdref(void) const;
    bool   isin(const GTime& time) const;
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
