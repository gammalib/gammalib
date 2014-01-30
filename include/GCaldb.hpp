/***************************************************************************
 *                 GCaldb.hpp - Calibration database class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file GCaldb.hpp
 * @brief Calibration database class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCALDB_HPP
#define GCALDB_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GCaldb
 *
 * @brief Interface for the Calibration database class
 *
 * This class holds the definition of a calibration database. If the void
 * constructor is invoked it tries setting the calibration database access
 * path using the GAMMALIB_CALDB or CALDB environment variables (the former
 * takes precedence over the latter). Alternatively, the calibration database
 * path can be passed to the constructor to override the environment variable
 * or set using the dir() method.
 *
 * It is assumed that the calibration data are found under
 *
 *     $CALDB/data/<mission>[/<instrument>]
 *
 * or
 *
 *     $GAMMALIB_CALDB/data/<mission>[/<instrument>]
 *
 * and that the Calibration Index File (CIF) is located at
 *
 *     $CALDB/data/<mission>[/<instrument>]/caldb.indx
 *
 * or
 *
 *     $GAMMALIB_CALDB/data/<mission>[/<instrument>]/caldb.indx
 *
 * where \<mission\> is the name of the mission and \<instrument\> is the
 * optional instrument name (all lower case).
 *
 * The calibration database for a given mission and instrument is opened
 * using the open() method. Once opened, database information can be
 * accessed. After usage, the database is closed using the close() method.
 ***************************************************************************/
class GCaldb : public GBase {

public:
    // Constructors and destructors
    GCaldb(void);
    GCaldb(const GCaldb& caldb);
    explicit GCaldb(const std::string& pathname);
    GCaldb(const std::string& mission, const std::string& instrument);
    virtual ~GCaldb(void);

    // Operators
    GCaldb& operator= (const GCaldb& caldb);

    // Methods
    void               clear(void);
    GCaldb*            clone(void) const;
    int                size(void) const;
    const std::string& dir(void) const;
    void               dir(const std::string& pathname);
    void               open(const std::string& mission,
                            const std::string& instrument = "");
    void               close(void);
    std::string        filename(const std::string& detector,
                                const std::string& filter,
                                const std::string& codename,
                                const std::string& date,
                                const std::string& time,
                                const std::string& expr);
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GCaldb& caldb);
    void        free_members(void);
    std::string rootdir(void) const;
    void        set_database(const std::string& pathname);
    std::string path(const std::string& mission, const std::string& instrument = "");
    std::string cifname(const std::string& mission, const std::string& instrument = "");

    // Protected data area
    std::string m_caldb;        //!< CALDB root directory
    std::string m_mission;      //!< Mission of opened database
    std::string m_instrument;   //!< Instrument of opened database
    std::string m_cifname;      //!< CIF filename of opened database
    GFits       m_fits;         //!< CIF FITS file
    GFitsTable* m_cif;          //!< Pointer to CIF table
};


/***********************************************************************//**
 * @brief Return calibration directory
 *
 * @return Calibration directory.
 ***************************************************************************/
inline
const std::string& GCaldb::dir(void) const
{
    return m_caldb;
}

#endif /* GCALDB_HPP */
