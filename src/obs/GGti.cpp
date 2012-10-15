/***************************************************************************
 *                 GGti.cpp  -  Good time interval class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @file GGti.cpp
 * @brief Good time interval class implementation.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GGti.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableDoubleCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_TSTART                                          "GGti::tstart(int)"
#define G_TSTOP                                            "GGti::tstop(int)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GGti::GGti(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] gti Good Time Intervals.
 ***************************************************************************/
GGti::GGti(const GGti& gti)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(gti);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GGti::~GGti(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] gti Good Time Intervals.
 ***************************************************************************/
GGti& GGti::operator= (const GGti& gti)
{
    // Execute only if object is not identical
    if (this != &gti) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(gti);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear Good Time Intervals
 ***************************************************************************/
void GGti::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Add one Good Time Interval
 *
 * @param[in] tstart Start time of interval.
 * @param[in] tstop Stop time of interval.
 *
 * This method adds a new GTI to the object. Adding means that if the
 * specified time interval overlaps with an existing time interval, the
 * existing time interval will be extended. If overlap with more than a
 * single interval exists, the intervals will be merged. Otherwise a new
 * interval will be inserted, respecting the time ordering of the intervals.
 *
 * If the time interval is not valid (tstop <= tstart), nothing is done.
 *
 * @todo Verify that both times have the same MJD reference, and that this
 *       reference corrsponds to the reference of the GTI.
 ***************************************************************************/
void GGti::add(const GTime& tstart, const GTime& tstop)
{
    // Continue only if time interval is valid
    if (tstop > tstart) {

        // Determine index at which GTI should be inserted
        int inx = 0;
        for (int i = 0; i < m_num; ++i) {
            if (tstart < m_start[i])
                break;
        }

        // Insert GTI
        insert_gti(inx, tstart, tstop);

        // Merge any overlapping GTIs
        merge_gtis();

    } // endif: Time interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append one Good Time Interval
 *
 * @param[in] tstart Start time of interval.
 * @param[in] tstop Stop time of interval.
 *
 * This method append a new GTI at the end of the object. In contrast to the
 * add() method, no checking of time ordering or interval overlapping will be
 * done.
 *
 * If the time interval is not valid (tstop <= tstart), nothing is done.
 *
 * @todo Verify that both times have the same MJD reference, and that this
 *       reference corrsponds to the reference of the GTI.
 ***************************************************************************/
void GGti::append(const GTime& tstart, const GTime& tstop)
{
    // Continue only if time interval is valid
    if (tstop > tstart) {

        // Insert GTI at end of list
        insert_gti(m_num, tstart, tstop);

    } // endif: Time interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert one Good Time Interval
 *
 * @param[in] tstart Start time of interval.
 * @param[in] tstop Stop time of interval.
 *
 * This method inserts a new GTI at the appropriate time ordered position in
 * the object. Time ordering will be done based on the start time of the GTI.
 * No check for interval overlap will be done.
 *
 * If the time interval is not valid (tstop <= tstart), nothing is done.
 *
 * @todo Verify that both times have the same MJD reference, and that this
 *       reference corrsponds to the reference of the GTI.
 ***************************************************************************/
void GGti::insert(const GTime& tstart, const GTime& tstop)
{
    // Continue only if time interval is valid
    if (tstop > tstart) {

        // Determine index at which GTI should be inserted
        int inx = 0;
        for (int i = 0; i < m_num; ++i) {
            if (tstart < m_start[i])
                break;
        }

        // Insert GTI
        insert_gti(inx, tstart, tstop);

    } // endif: Time interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Reduce Good Time Intervals to specified interval
 *
 * @param[in] tstart Start time of interval.
 * @param[in] tstop Stop time of interval.
 *
 * This method reduces the Good Time Intervals to the specified interval.
 * Reducing means that all GTIs are dropped that fall outside the interval,
 * and GTIs will be limited to within the interval if they overlap.
 *
 * If the time interval is not valid (tstop <= tstart), nothing is done.
 *
 * @todo Verify that both times have the same MJD reference, and that this
 *       reference corrsponds to the reference of the GTI.
 ***************************************************************************/
void GGti::reduce(const GTime& tstart, const GTime& tstop)
{
    // Continue only if time interval is valid
    if (tstop > tstart) {

        // Adjust existing GTIs. This will limit all GTIs to [tstart,tstop].
        // All GTIs outside [tstart,tstop] will have start > stop. The number
        // of valid GTIs will also be determined.
        int num = 0;
        for (int i = 0; i < m_num; ++i) {
            if (m_start[i] < tstart) {
                m_start[i] = tstart;
            }
            if (m_stop[i] > tstop) {
                m_stop[i] = tstop;
            }
            if (m_start[i] <= m_stop[i]) {
                num++;
            }
        }

        // If we still have valid GTIs then allocate memory for them, copy
        // over the start and stop times, and update the attributes
        if (num > 0) {

            // Allocate new intervals
            GTime* start = new GTime[num];
            GTime* stop  = new GTime[num];

            // Copy valid intervals
            for (int i = 0; i < m_num; ++i) {
                if (m_start[i] <= m_stop[i]) {
                    start[i] = m_start[i];
                    stop[i]  = m_stop[i];
                }
            }

            // Free old memory
            if (m_start != NULL) delete [] m_start;
            if (m_stop  != NULL) delete [] m_stop;

            // Set new memory
            m_start = start;
            m_stop  = stop;

            // Set attributes
            m_num = num;
            set_attributes();

        } // endif: there were still valid GTIs

        // ... otherwise we remove all GTIs
        else {
            clear();
        }

    } // endif: Time interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load GTIs from FITS file.
 *
 * @param[in] filename FITS filename.
 * @param[in] extname GTI extension name (default is "GTI")
 *
 * This method loads the Good Time Intervals from a FITS file.
 ***************************************************************************/
void GGti::load(const std::string& filename, const std::string& extname)
{
    // Allocate FITS file
    GFits file;

    // Open FITS file
    file.open(filename);

    // Get GTI HDU
    GFitsTable* hdu = file.table(extname);

    // Read GTI from HDU
    read(hdu);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save GTI intervals to FITS file.
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite any existing GTI extension?
 * @param[in] extname GTI extension name (default is "GTI")
 *
 * The method saves GTI into a FITS file. If the file does not exist it is
 * created. If the file exists the GTI is appended as extension. If another
 * GTI exists already it is overwritten if clobber=true.
 ***************************************************************************/
void GGti::save(const std::string& filename, bool clobber,
                const std::string& extname) const
{
    // Allocate FITS file
    GFits file;

    // Write GTI to FITS file
    write(&file, extname);

    // Save to file
    file.saveto(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read GTI intervals from file.
 *
 * @param[in] hdu Pointer to FITS table.
 *
 * @todo Method assumes that times are in MET.
 * @todo Read header keywords.
 ***************************************************************************/
void GGti::read(GFitsTable* hdu)
{
    // Free members
    free_members();

    // Initialise attributes
    init_members();

    // Extract GTI information from FITS file
    m_num = hdu->integer("NAXIS2");
    if (m_num > 0) {

        // Set GTIs
        m_start = new GTime[m_num];
        m_stop  = new GTime[m_num];
        for (int i = 0; i < m_num; ++i) {
            m_start[i].met((*hdu)["START"].real(i));
            m_stop[i].met((*hdu)["STOP"].real(i));
        }

        // Set attributes
        set_attributes();

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write GTI intervals into FITS file.
 *
 * @param[in] file Pointer to FITS file.
 * @param[in] extname GTI extension name (default is "GTI")
 *
 * The method saves GTI into a FITS file. If the file does not exist it is
 * created. If the file exists the GTI is appended as extension. If another
 * GTI exists already it is overwritten if clobber=true.
 *
 * @todo Method assumes that times are in MET.
 * @todo Write header keywords.
 * @todo Implement clobber method for overwriting of existing GTIs.
 ***************************************************************************/
void GGti::write(GFits* file, const std::string& extname) const
{
    // Create GTI columns
    GFitsTableDoubleCol cstart = GFitsTableDoubleCol("START", m_num);
    GFitsTableDoubleCol cstop  = GFitsTableDoubleCol("STOP", m_num);

    // Fill GTI columns
    for (int i = 0; i < m_num; ++i) {
        cstart(i) = m_start[i].met();
        cstop(i)  = m_stop[i].met();
    }

    // Create GTI table
    GFitsBinTable* table = new GFitsBinTable(m_num);
    table->append_column(cstart);
    table->append_column(cstop);
    table->extname(extname);

    // Write to FITS file
    file->append(*table);

    // Free table
    delete table;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns start time of specific GTI
 *
 * @param[in] inx Index of GTI (0 ... size()-1).
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
GTime GGti::tstart(int inx) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (inx < 0 || inx >= m_num)
        throw GException::out_of_range(G_TSTART, inx, 0, m_num-1);
    #endif

    // Return
    return (m_start[inx]);
}


/***********************************************************************//**
 * @brief Returns stop time of specific GTI
 *
 * @param[in] inx Index of GTI (0 ... size()-1).
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
GTime GGti::tstop(int inx) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (inx < 0 || inx >= m_num)
        throw GException::out_of_range(G_TSTOP, inx, 0, m_num-1);
    #endif

    // Return
    return (m_stop[inx]);
}


/***********************************************************************//**
 * @brief Tests whether time is within GTI intervals
 *
 * @param[in] time Time to be tested.
 ***************************************************************************/
bool GGti::isin(const GTime& time) const
{
    // Initialise test
    bool found = false;

    // Test all GTIs
    for (int i = 0; i < m_num; ++i) {
        if (time >= m_start[i] && time <= m_stop[i]) {
            found = true;
            break;
        }
    }

    // Return result
    return found;
}


/***********************************************************************//**
 * @brief Print Good Time Intervals
 ***************************************************************************/
std::string GGti::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GGti ===");
    result.append("\n"+parformat("Number of intervals")+str(size()));
    result.append("\n"+parformat("Ontime")+str(ontime())+" sec");
    result.append("\n"+parformat("Elapsed time")+str(telapse())+" sec");
    result.append("\n"+parformat("MJD reference date")+str(mjdref())+" days");
    result.append("\n"+parformat("Time range"));
    result.append(tstart().print()+" - "+tstop().print());

    // Return result
    return result;
}
    

/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GGti::init_members(void)
{
    // Initialise members
    m_num     = 0;
    m_tstart.clear();
    m_tstop.clear();
    m_ontime  = 0.0;
    m_telapse = 0.0;
    m_mjdref  = 0.0;
    m_start   = NULL;
    m_stop    = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] gti GGti members which should be copied.
 ***************************************************************************/
void GGti::copy_members(const GGti& gti)
{
    // Copy attributes
    m_num     = gti.m_num;
    m_tstart  = gti.m_tstart;
    m_tstop   = gti.m_tstop;
    m_ontime  = gti.m_ontime;
    m_telapse = gti.m_telapse;
    m_mjdref  = gti.m_mjdref;

    // Copy start/stop times
    if (m_num > 0) {
        m_start = new GTime[m_num];
        m_stop  = new GTime[m_num];
        for (int i = 0; i < m_num; ++i) {
            m_start[i] = gti.m_start[i];
            m_stop[i]  = gti.m_stop[i];
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GGti::free_members(void)
{
    // Free memory
    if (m_start != NULL) delete [] m_start;
    if (m_stop  != NULL) delete [] m_stop;

    // Signal free pointers
    m_start = NULL;
    m_stop  = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set class attributes
 ***************************************************************************/
void GGti::set_attributes(void)
{
    // Continue only if there are GTIs
    if (m_num > 0) {

        // Set attributes
        m_tstart  = m_start[0];
        m_tstop   = m_stop[m_num-1];
        m_telapse = m_tstop.met() - m_tstart.met();
        m_ontime  = 0.0;
        for (int i = 0; i < m_num; ++i) {
            m_ontime += (m_stop[i].met() - m_start[i].met());
        }

    } // endif: there were GTIs

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone GTI
 *
 * Cloning provides a copy of the GTIs. Cloning is used to allocate
 * derived classes into a base class pointer.
 ***************************************************************************/
GGti* GGti::clone(void) const
{
    // Clone this image
    return new GGti(*this);
}


/***********************************************************************//**
 * @brief Insert Good Time Interval
 *
 * @param[in] inx Index at which GTI is to be inserted.
 * @param[in] tstart Start time of interval to be inserted.
 * @param[in] tstop Stop time of interval to be inserted.
 *
 * Inserts a GTI at a given position in the GTI array. This method does not
 * assure the time ordering of the GTIs, this has to be done by the client
 * that determines the appropriate value of inx.
 * Invalid parameters do not produce any exception, but are handled
 * transparently. If the GTI is invalid (tstop <= tstart), no interval
 * will be inserted. If the index is out of the valid range, the index
 * will be adjusted to either the first or the last element.
 ***************************************************************************/
void GGti::insert_gti(int inx, const GTime& tstart, const GTime& tstop)
{
    // Continue only if time interval is valid
    if (tstop > tstart) {

        // If inx is out of range then adjust it
        if (inx < 0)     inx = 0;
        if (inx > m_num) inx = m_num;

        // Allocate new intervals
        int    num   = m_num+1;
        GTime* start = new GTime[num];
        GTime* stop  = new GTime[num];

        // Copy intervals before GTI to be inserted
        for (int i = 0; i < inx; ++i) {
            start[i] = m_start[i];
            stop[i]  = m_stop[i];
        }

        // Insert GTI
        start[inx] = tstart;
        stop[inx]  = tstop;

        // Copy intervals after GTI to be inserted
        for (int i = inx+1; i < num; ++i) {
            start[i] = m_start[i-1];
            stop[i]  = m_stop[i-1];
        }

        // Free memory
        if (m_start != NULL) delete [] m_start;
        if (m_stop  != NULL) delete [] m_stop;

        // Set new memory
        m_start = start;
        m_stop  = stop;

        // Set attributes
        m_num = num;
        set_attributes();

    } // endif: Time interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Merge all overlapping Good Time Intervals
 *
 * Merges all overlapping GTIs into a single GTI. This method assumes that
 * the GTI are ordered correctly by start time. Note that the method does
 * not actually reduce the memory size but just updates the information on
 * the number of elements in the array.
 ***************************************************************************/
void GGti::merge_gtis(void)
{
    // Find overlaps
    int i   = 0;
    int num = m_num;
    while (i < num-1) {

        // If GTI overlaps with following one then merge both GTIs, move
        // all remaining GTIs one position up, and reduce number of elements
        if (m_start[i+1] <= m_stop[i]) {
            m_start[i] = (m_start[i] < m_start[i+1]) ? m_start[i] : m_start[i+1];
            m_stop[i]  = (m_stop[i]  > m_stop[i+1])  ? m_stop[i]  : m_stop[i+1];
            for (int k = i+2; k < num; ++k) {
                m_start[k-1] = m_start[k];
                m_stop[k-1]  = m_stop[k];
            }
            num--;
        }

        // Otherwise increment GTI index
        else
            i++;

    } // endwhile: there were still GTIs to check

    // Update number of elements in GTI
    m_num = num;

    // Return
    return;
}
