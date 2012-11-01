/***************************************************************************
 *             GEbounds.i  -  Energy boundary class python I/F             *
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
 * @file GEbounds.i
 * @brief GEbounds class python interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEbounds.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GEbounds
 *
 * @brief Interface for the GEbounds class.
 *
 * This class holds a list of energy intervals.
 ***************************************************************************/
class GEbounds : public GBase {
public:
    // Constructors and destructors
    GEbounds(void);
    GEbounds(const GEbounds& ebds);
    explicit GEbounds(const std::string& filename, const std::string& extname);
    virtual ~GEbounds(void);

    // Methods
    void      clear(void);
    GEbounds* clone(void) const;
    void      append(const GEnergy& emin, const GEnergy& emax);
    void      insert(const GEnergy& emin, const GEnergy& emax);
    void      setlin(const GEnergy& emin, const GEnergy& emax, const int& num);
    void      setlog(const GEnergy& emin, const GEnergy& emax, const int& num);
	void      load(const std::string& filename,
                   const std::string& extname);
	void      save(const std::string& filename, bool clobber,
                   const std::string& extname) const;
    void      read(GFitsTable* hdu);
    void      write(GFits* file, const std::string& extname) const;
    int       index(const GEnergy& eng) const;
    int       size(void) const;
    GEnergy   emin(void) const;
    GEnergy   emax(void) const;
    GEnergy   emin(int inx) const;
    GEnergy   emax(int inx) const;
    GEnergy   emean(int inx) const;
    GEnergy   elogmean(int inx) const;
    bool      contains(const GEnergy& eng) const;
};


/***********************************************************************//**
 * @brief GEbounds class extension
 ***************************************************************************/
%extend GEbounds {
    char *__str__() {
        return tochar(self->print());
    }
    GEbounds copy() {
        return (*self);
    }
};
