/***************************************************************************
 *                  GGti.cpp - Good time interval class                    *
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
#define G_POP                                               "GGti::pop(int&)"
#define G_TSTART                                         "GGti::tstart(int&)"
#define G_TSTOP                                           "GGti::tstop(int&)"

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
    // Initialise class members
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
    // Initialise class members
    init_members();

    // Copy members
    copy_members(gti);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Time reference constructor
 *
 * @param[in] ref Time reference.
 *
 * Constructs Good Time Intervals using a specific time reference. The time
 * reference will be used when writing the Good Time Intervals into a FITS
 * file.
 ***************************************************************************/
GGti::GGti(const GTimeReference& ref)
{
    // Initialise class members
    init_members();

    // Set time reference
    this->reference(ref);

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
 * @return Good Time Intervals.
 ***************************************************************************/
GGti& GGti::operator= (const GGti& gti)
{
    // Execute only if object is not identical
    if (this != &gti) {

        // Free members
        free_members();

        // Initialise private members
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
 * @brief Clone Good Time Intervals
 *
 * @return Pointer to deep copy of Good Time Intervals.
 ***************************************************************************/
GGti* GGti::clone(void) const
{
    return new GGti(*this);
}


/***********************************************************************//**
 * @brief Append Good Time Interval
 *
 * @param[in] tstart Start time of interval.
 * @param[in] tstop Stop time of interval.
 *
 * Appends a Good Time Interval at the end of the container.
 *
 * The method does nothing if the Good Time Interval is not valid
 * (i.e. @p tstart >= @p tstop).
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
 * @brief Insert Good Time Interval
 *
 * @param[in] tstart Start time of interval.
 * @param[in] tstop Stop time of interval.
 *
 * Inserts a Good Time Interval into the container after the first interval
 * that has a start time smaller than @p tstart. The method implicitely
 * assumes that the Good Time Intervals are ordered by increasing start time.
 *
 * The method does nothing if the Good Time Interval is not valid
 * (i.e. @p tstart >= @p tstop).
 ***************************************************************************/
void GGti::insert(const GTime& tstart, const GTime& tstop)
{
    // Continue only if time interval is valid
    if (tstop > tstart) {

        // Determine index at which GTI should be inserted
        int inx = 0;
        for (int i = 0; i < m_num; ++i) {
            if (tstart < m_start[i]) {
                break;
            }
        }

        // Insert interval
        insert_gti(inx, tstart, tstop);

    } // endif: Time interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Merge all overlapping Good Time Intervals
 *
 * Merges all overlapping or connecting successive Good Time Intervals. The
 * method implicitely assumes that the intervals are ordered by increasing
 * start time.
 *
 * Note that the method does not actually reduce the memory size but just
 * updates the information on the number of elements in the array.
 ***************************************************************************/
void GGti::merge(void)
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
        else {
            i++;
        }

    } // endwhile: there were still GTIs to check

    // Update number of elements in GTI
    m_num = num;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Merge Good Time Interval into container
 *
 * @param[in] tstart Start time of interval.
 * @param[in] tstop Stop time of interval.
 *
 * Inserts a Good Time Interval into the container after the first interval
 * that has a start time smaller than @p tstart and then merges any
 * overlapping or connecting Good Time Intervals. The method implicitely
 * assumes that the intervals are ordered by increasing start time.
 *
 * The method does nothing if the Good Time Interval is not valid
 * (i.e. @p tstart >= @p tstop).
 ***************************************************************************/
void GGti::merge(const GTime& tstart, const GTime& tstop)
{
    // Continue only if time interval is valid
    if (tstop > tstart) {

        // Determine index at which GTI should be inserted
        int inx = 0;
        for (int i = 0; i < m_num; ++i) {
            if (tstart < m_start[i]) {
                break;
            }
        }

        // Insert GTI
        insert_gti(inx, tstart, tstop);

        // Merge any overlapping GTIs
        merge();

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
 * Reduces the Good Time Intervals to the specified interval. Reducing means
 * that all Good Time Intervals are dropped that fall outside the specified
 * interval [@p tstart, @p tstop], and Good Time Intervals will be limited
 * to >@p tstart and <=@p tstop in case that their boundaries are outside
 * [@p tstart, @p tstop].
 *
 * The method does nothing if the Good Time Interval is not valid
 * (i.e. @p tstart >= @p tstop).
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
 * @brief Remove Good Time Interval
 *
 * @param[in] index Good Time Interval index (0 ... size()-1).
 *
 * Removes Good Time Interval at @p index from the container. All intervals
 * after the specified @p index are moved forward by one position.
 *
 * Note that the method does not actually reduce the memory size but just
 * updates the information on the number of elements in the array.
 ***************************************************************************/
void GGti::pop(const int& index)
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_POP, index, 0, m_num-1);
    }
    #endif

    // Move all elements located after index forward
    for (int i = index+1; i < m_num; ++i) {
        m_start[i-1] = m_start[i];
        m_stop[i-1]  = m_stop[i];
    }

    // Reduce number of elements by one
    m_num--;

    // Update attributes
    set_attributes();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Reserve space for Good Time Intervals
 *
 * @param[in] num Number of elements.
 *
 * This method does nothing (it is needed to satisfy the GContainer
 * interface).
 ***************************************************************************/
void GGti::reserve(const int& num)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Append Good Time Intervals
 *
 * @param[in] gti Good Time Intervals.
 *
 * Append Good Time Intervals to the container. The method performs automatic
 * time reference conversion in case that the specified Good Time Intervals
 * @p gti have a time reference that differs from that of the current
 * instance.
 ***************************************************************************/
void GGti::extend(const GGti& gti)
{
    // Do nothing if Good Time Intervals are empty
    if (!gti.isempty()) {

        // Allocate new intervals
        int    num   = m_num+gti.size();
        GTime* start = new GTime[num];
        GTime* stop  = new GTime[num];

        // Initialise index
        int inx = 0;
        
        // Copy existing intervals
        for (; inx < m_num; ++inx) {
            start[inx] = m_start[inx];
            stop[inx]  = m_stop[inx];
        }

        // Append intervals. Convert to GTI reference on the fly.
        for (int i = 0; i < gti.size(); ++i, ++inx) {
            double tstart = gti.m_start[i].convert(gti.reference());
            double tstop  = gti.m_stop[i].convert(gti.reference());
            start[inx].set(tstart, this->reference());
            stop[inx].set(tstop,   this->reference());
        }

        // Free memory
        if (m_start != NULL) delete [] m_start;
        if (m_stop  != NULL) delete [] m_stop;

        // Set new memory
        m_start = start;
        m_stop  = stop;

        // Set number of elements
        m_num = num;

        // Set attributes
        set_attributes();

    } // endif: Good Time Intervals were not empty

    // Return
    return;
}



/***********************************************************************//**
 * @brief Load Good Time Intervals from FITS file
 *
 * @param[in] filename FITS filename.
 * @param[in] extname GTI extension name (defaults to "GTI")
 *
 * Loads the Good Time Intervals from FITS file.
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
 * @brief Save Good Time Intervals intervals to FITS file
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite any existing GTI extension?
 * @param[in] extname GTI extension name (defaults to "GTI")
 *
 * Saves Good Time Intervals into extension @p extname of a FITS file. If the
 * file does not exist it is created. If the file exists the GTI is appended
 * as extension. If another GTI exists already it is overwritten if
 * @p clobber=true.
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
 * @brief Read Good Time Intervals and time reference from FITS table
 *
 * @param[in] hdu Pointer to FITS table.
 *
 * Reads the Good Time Intervals and time reference from a FITS table.
 ***************************************************************************/
void GGti::read(const GFitsTable* hdu)
{
    // Free members
    free_members();

    // Initialise attributes
    init_members();

    // Read time reference
    m_reference.read(hdu);

    // Extract GTI information from FITS file
    m_num = hdu->integer("NAXIS2");
    if (m_num > 0) {

        // Set GTIs
        m_start = new GTime[m_num];
        m_stop  = new GTime[m_num];
        for (int i = 0; i < m_num; ++i) {
            m_start[i].set((*hdu)["START"].real(i), m_reference);
            m_stop[i].set((*hdu)["STOP"].real(i), m_reference);
        }

        // Set attributes
        set_attributes();

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write Good Time Intervals and time reference into FITS object
 *
 * @param[in] file Pointer to FITS file.
 * @param[in] extname GTI extension name (defaults to "GTI")
 *
 * Saves Good Time Intervals and time reference into a FITS object. If the
 * file does not exist it is created. If the file exists the GTI is appended
 * as extension. If another GTI exists already it is overwritten if
 * @p clobber=true.
 *
 * @todo Implement clobber method for overwriting of existing GTIs.
 ***************************************************************************/
void GGti::write(GFits* file, const std::string& extname) const
{
    // Create GTI columns
    GFitsTableDoubleCol cstart = GFitsTableDoubleCol("START", m_num);
    GFitsTableDoubleCol cstop  = GFitsTableDoubleCol("STOP", m_num);

    // Fill GTI columns in specified time reference
    for (int i = 0; i < m_num; ++i) {
        cstart(i) = m_start[i].convert(m_reference);
        cstop(i)  = m_stop[i].convert(m_reference);
    }

    // Create GTI table
    GFitsBinTable* table = new GFitsBinTable(m_num);
    table->append_column(cstart);
    table->append_column(cstop);
    table->extname(extname);

    // Write time reference
    m_reference.write(table);

    // Write to FITS file
    file->append(*table);

    // Free table
    delete table;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns earliest start time in Good Time Intervals
 *
 * @return Earliest start time in Good Time Intervals.
 ***************************************************************************/
const GTime& GGti::tstart(void) const
{
    // Return
    return m_tstart;
}


/***********************************************************************//**
 * @brief Returns latest stop time in Good Time Intervals
 *
 * @return Latest stop time in Good Time Intervals.
 ***************************************************************************/
const GTime& GGti::tstop(void) const
{
    // Return
    return m_tstop;
}


/***********************************************************************//**
 * @brief Returns start time for a given Good Time Interval
 *
 * @param[in] index Good Time Interval index (0 ... size()-1).
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
const GTime& GGti::tstart(const int& index) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_TSTART, index, 0, m_num-1);
    }
    #endif

    // Return
    return (m_start[index]);
}


/***********************************************************************//**
 * @brief Returns stop time for a given Good Time Interval
 *
 * @param[in] index Good Time Interval index (0 ... size()-1).
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
const GTime& GGti::tstop(const int& index) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_TSTOP, index, 0, m_num-1);
    }
    #endif

    // Return
    return (m_stop[index]);
}


/***********************************************************************//**
 * @brief Returns elapsed time
 *
 * @return Elapsed time [seconds].
 *
 * Returns the elapsed time in seconds. The elapsed time is defined as the
 * time difference between the latest stop time and the earliest start time
 * in the Good Time Intervals.
 ***************************************************************************/
const double& GGti::telapse(void) const
{
    // Return
    return m_telapse;
}


/***********************************************************************//**
 * @brief Returns ontime
 *
 * @return Ontime [seconds].
 *
 * Returns the ontime in seconds. The ontime is defined as the sum of all
 * Good Time Intervals.
 ***************************************************************************/
const double& GGti::ontime(void) const
{
    // Return
    return m_ontime;
}


/***********************************************************************//**
 * @brief Set time reference for Good Time Intervals
 *
 * @param[in] ref Time reference.
 *
 * Sets the time reference of the Good Time Intervals. This defines the
 * reference time which will be writted into the FITS file upon saving of
 * the Good Time Intervals.
 ***************************************************************************/
void GGti::reference(const GTimeReference& ref)
{
    // Set time reference
    m_reference = ref;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return time reference for Good Time Intervals
 *
 * @return Time reference.
 *
 * Returns the time reference of the Good Time Intervals.
 ***************************************************************************/
const GTimeReference& GGti::reference(void) const
{
    // Return time reference
    return m_reference;
}


/***********************************************************************//**
 * @brief Checks whether Good Time Intervals contain time
 *
 * @param[in] time Time to be checked.
 *
 * Checks if a given @p time falls in at least one of the Good Time
 * Intervals. The method exits when the first matching interval has been
 * found.
 ***************************************************************************/
bool GGti::contains(const GTime& time) const
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
 *
 * @return String containing Good Time Interval information
 ***************************************************************************/
std::string GGti::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GGti ===");

    // Append GTI information
    result.append("\n"+parformat("Number of intervals")+str(size()));
    result.append("\n"+parformat("Ontime")+str(ontime())+" sec");
    result.append("\n"+parformat("Elapsed time")+str(telapse())+" sec");
    result.append("\n"+parformat("Time range"));
    result.append(str(tstart().convert(m_reference)));
    result.append(" - ");
    result.append(str(tstop().convert(m_reference)));
    result.append(" "+reference().timeunit());
    result.append(" ("+reference().timesys()+")");
    result.append("\n"+parformat("Reference MDJ"));
    result.append(str(reference().mjdref()));

    // Append time reference information
    //result.append("\n"+reference().print());

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
    m_start   = NULL;
    m_stop    = NULL;

    // Initialise time reference with native reference
    GTime time;
    m_reference = time.reference();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] gti Good Time Intervals.
 ***************************************************************************/
void GGti::copy_members(const GGti& gti)
{
    // Copy attributes
    m_num       = gti.m_num;
    m_tstart    = gti.m_tstart;
    m_tstop     = gti.m_tstop;
    m_ontime    = gti.m_ontime;
    m_telapse   = gti.m_telapse;
    m_reference = gti.m_reference;

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
 *
 * Compute the following class attributes:
 *
 *     m_tstart  - Earliest start time of GTIs
 *     m_stop    - Latest stop time of GTIs
 *     m_telapse - Latest stop time minus earliest start time of GTIs [sec]
 *     m_ontime  - Sum of all intervals [sec]
 ***************************************************************************/
void GGti::set_attributes(void)
{
    // If there are intervals then determine the start and stop time
    // from these intervals ...
    if (m_num > 0) {
        m_tstart = m_start[0];
        m_tstop  = m_stop[0];
        for (int i = 1; i < m_num; ++i) {
            if (m_start[i] < m_tstart) m_tstart = m_start[i];
            if (m_stop[i]  > m_tstop)  m_tstop  = m_stop[i];
        }
    }

    // ... otherwise clear the start and stop time
    else {
        m_tstart.clear();
        m_tstop.clear();
    }

    // Set attributes
    m_telapse = m_tstop.secs() - m_tstart.secs();
    m_ontime  = 0.0;
    for (int i = 0; i < m_num; ++i) {
        m_ontime += (m_stop[i].secs() - m_start[i].secs());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert Good Time Interval
 *
 * @param[in] index Index after which interval is inserted.
 * @param[in] tstart Start time of interval.
 * @param[in] tstop Stop time of interval.
 *
 * Inserts a Good Time Interval after the specified @p index in the Good
 * Time Intervals. The method does not reorder the intervals by time,
 * instead the client needs to determine the approriate @p index.
 *
 * Invalid parameters do not produce any exception, but are handled
 * transparently. If the interval is invalid (i.e. @p tstart >= @p tstop)
 * then nothing is done. If the @p index is out of the valid range, the
 * index will be adjusted to either the first or the last element.
 ***************************************************************************/
void GGti::insert_gti(const int& index, const GTime& tstart, const GTime& tstop)
{
    // Continue only if time interval is valid
    if (tstop > tstart) {

        // Set index
        int inx = index;

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

        // Set number of elements
        m_num = num;
        
        // Set attributes
        set_attributes();

    } // endif: Time interval was valid

    // Return
    return;
}
