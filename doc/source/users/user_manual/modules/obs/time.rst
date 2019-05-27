.. _um_obs_time:

Times in GammaLib
=================

Times in GammaLib are implemented by the :doxy:`GTime` class that provides
transparent handling of times independent of their time reference, time
system and time unit.
Time is stored in :doxy:`GTime` in seconds in a GammaLib native
reference system, which has zero time at January 1, 2010, 00:00:00
Terrestrial Time (TT).
With respect to Coordinated Universal Time (UTC), TT time is greater than
UTC time by 66.184 sec at January 1, 2010, 00:00:00.
The difference is composed of leap seconds that synchronize the Earth
rotation variations with a uniform clock, and of a fixed offset of 32.184
seconds between TT and International Atomic Time (TAI).

The value of a :doxy:`GTime` instance can be set in native seconds
using :doxy:`GTime::secs` or in native days using :doxy:`GTime::days`.
It can furthermore also be set in Julian Days (TT) using :doxy:`GTime::jd`,
in Modified Julian Days (TT) using :doxy:`GTime::mjd` and as a string
in the ISO 8601 standard YYYY-MM-DDThh:mm:ss.s (UTC) using 
:doxy:`GTime::utc`.
Equivalent methods exist for retrieving the time in various formats,
allowing thus conversion from one format to the other.
In addition, the :doxy:`GTime::now` method sets the time to current time of the
computer, the method :doxy:`GTime::lmst` returns the local mean sidereal time
and the method :doxy:`GTime::last` the local apparent sidereal time for a given
geographic longitude.

As instrument times are generally given in a local reference, conversion
between different time reference systems is also supported.
Time references are specified by the :doxy:`GTimeReference` class that
define the Modified Julian Days (MJD) reference in TT, the time unit
(seconds or days), the time system (TT or UTC) and the time reference
(LOCAL).
A time can be converted into a reference using the :doxy:`GTime::convert`
method and set to a value specified in a given reference using the
:doxy:`GTime::set` method.
The native GammaLib reference can be retrieved using the
:doxy:`GTime::reference` method.

Instances of :doxy:`GTime` can be added, subtracted and compared using the
usual mathematical operators. Subtracting two instances of time gives the
differences in seconds.

Instances of :doxy:`GTime` can be collected in the :doxy:`GTimes` container
class.
This class is for the moment implemented as a minimal container class
without support for reading from and writing to files.

:doxy:`GTime` objects are also used in the definition of Good Time
Intervals (GTIs), which are intervals of contiguous time during which
data are valid for science analysis. GTIs are implemented by the
:doxy:`GGti` class which is a container of time intervals, formed
by a start and a stop value. These values can be accessed using the
:doxy:`GGti::tstart(int)` and :doxy:`GGti::tstop(int)` methods, both 
returning a :doxy:`GTime` object.
The summed length of all intervals is known
as the ontime which is returned by the :doxy:`GGti::ontime` method in
units of seconds. The elapsed time, returned by :doxy:`GGti::telapse`
in seconds is the difference between the last stop time and the first
start time. Time intervals are appended or inserted using the
:doxy:`GGti::append` and :doxy:`GGti::insert` methods. These methods
do not check whether intervals overlap in time, which may lead to an
errornous ontime value. To remove overlaps, the :doxy:`GGti::merge`
method can be used that will merge all overlapping intervals.
GTIs can be limited to a specific interval by applying the 
:doxy:`GGti::reduce` method.
GTIs can be written to or read from FITS files. The file format is the
standard OGIP format used for many high energy missions. The time
reference of the stored values will be defined by a :doxy:`GTimeReference`
object that can be set and retrieved using the :doxy:`GGti::reference`
method. Time reference information is written to the FITS file in
OGIP compliant header keywords.
