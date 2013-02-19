/***************************************************************************
 *                  GGti.hpp - Good time interval class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GGti.hpp
 * @brief Good time interval class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GGTI_HPP
#define GGTI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GContainer.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GTime.hpp"
#include "GTimeReference.hpp"


/***********************************************************************//**
 * @class GGti
 *
 * @brief Good Time Interval class
 *
 * This class holds a list of Good Time Intervals, i.e. time intervals that
 * are valid for science analysis. Times are stored using the GTime class.
 * The class also holds information about the time reference, which will
 * be retained when reading and used when writing, so that GTIs are always
 * written in the specified time reference.
 ***************************************************************************/
class GGti : public GContainer {

public:
    // Constructors and destructors
    GGti(void);
    GGti(const GGti& gti);
    explicit GGti(const GTimeReference& ref);
    virtual ~GGti(void);

    // Operators
    GGti& operator=(const GGti& gti);

    // Methods
    void                  clear(void);
    GGti*                 clone(void) const;
    int                   size(void) const;
    bool                  isempty(void) const;
    void                  add(const GTime& tstart, const GTime& tstop);
    void                  append(const GTime& tstart, const GTime& tstop);
    void                  insert(const GTime& tstart, const GTime& tstop);
    void                  reduce(const GTime& tstart, const GTime& tstop);
    void                  pop(const int& index);
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
    std::string           print(void) const;

protected:
    // Protected methods
    void  init_members(void);
    void  copy_members(const GGti& gti);
    void  free_members(void);
    void  set_attributes(void);
    void  insert_gti(const int& index, const GTime& tstart, const GTime& tstop);
    void  merge_gtis(void);

    // Protected data area
    int             m_num;       //!< Number of intervals
    GTime           m_tstart;    //!< Start of observation
    GTime           m_tstop;     //!< Stop of observation
    double          m_ontime;    //!< Sum of GTI durations (in seconds)
    double          m_telapse;   //!< Time between start of first GTI and stop of last GTI (in seconds)
    GTime          *m_start;     //!< Array of start times
    GTime          *m_stop;      //!< Array of stop times
    GTimeReference  m_reference; //!< Time reference
};

#endif /* GGTI_HPP */
