Constants and utility functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :doxy:`GTools.hpp` header defines a number of constants and utility 
functions that are widely used in GammaLib.

Constants
^^^^^^^^^

The following constants are available:

========================== ================================= =======
Constant                   Value                             Purpose
========================== ================================= =======
``gammalib::MeV2erg``      :math:`1.6021765 \times 10^{-6}`  Converts MeV to erg
``gammalib::erg2MeV``      :math:`624150.96`                 Converts erg to MeV
``gammalib::MeV2Angstrom`` :math:`1.239841875\times 10^{-2}` Converts MeV to Angstrom
``gammalib::pc2cm``        :math:`3.08568025 \times 10^{18}` Converts pc to cm
``gammalib::sec_in_day``   :math:`86400.0`                   Number of seconds in one day
``gammalib::sec2day``      :math:`1/86400`                   Converts seconds to days
``gammalib::tai2tt``       :math:`32.184`                    Converts TAI to TT system
``gammalib::mec2``         :math:`0.5109989461`              Electron rest mass in MeV
========================== ================================= =======


Functions
^^^^^^^^^

The following functions are available:

===================================== ===========
Function                              Description
===================================== ===========
``gammalib::strip_whitespace``        Strips all leading and trailing whitespace from string.
``gammalib::strip_chars``             Strips all leading and trailing characters from string.
``gammalib::rstrip_chars``            Strips all trailing characters from string.
``gammalib::replace_segment``         Replace string segment by another segment.
``gammalib::expand_env``              Replace any environment variables in string by its value.
``gammalib::filepath``                Build path from filename and path.
``gammalib::str``                     Conversion of C-types to strings.
``gammalib::tochar``                  Conversion of string to ``char``.
``gammalib::toshort``                 Conversion of string to ``short``.
``gammalib::toushort``                Conversion of string to ``unsigned short``.
``gammalib::toint``                   Conversion of string to ``int``.
``gammalib::touint``                  Conversion of string to ``unsigned int``.
``gammalib::tolong``                  Conversion of string to ``long``.
``gammalib::toulong``                 Conversion of string to ``unsigned long``.
``gammalib::tolonglong``              Conversion of string to ``long long``.
``gammalib::toulonglong``             Conversion of string to ``unsigned long long``.
``gammalib::tofloat``                 Conversion of string to ``float``.
``gammalib::todouble``                Conversion of string to ``double``.
``gammalib::toupper``                 Conversion of string to upper case letters.
``gammalib::tolower``                 Conversion of string to lower case letters.
``gammalib::split``                   Split string in vector of strings.
``gammalib::fill``                    Fill string with a number of replications of a string.
``gammalib::left``                    Left justify string to achieve a given length of characters.
``gammalib::right``                   Right justify string to achieve a given length of characters.
``gammalib::centre``                  Centre string to achieve a given length of characters.
``gammalib::parformat``               Format string for parameter value display.
``gammalib::number``                  Append a ``s`` to a noun if the number is larger than one.
``gammalib::plaw_photon_flux``        Compute photon flux under a power law.
``gammalib::plaw_energy_flux``        Compute energy flux under a power law.
``gammalib::elogmean``                Computes geometric mean of energy.
``gammalib::dir_exists``              Check whether a directory exists.
``gammalib::is_infinite``             Check whether a double precision value is infinite.
``gammalib::is_notanumber``           Check whether a double precision value is not a number.
``gammalib::contains``                Check whether a string contains a sub-string.
``gammalib::warning``                 Dump warning in console.
``gammalib::xml2str``                 Converts XML to string.
``gammalib::str2xml``                 Converts string to XML.
``gammalib::xml_has_par``             Checks is XML file has parameter.
``gammalib::xml_need_par``            Require specific parameter in XML file.
``gammalib::xml_get_par``             Get parameter from XML file.
``gammalib::xml_get_attr``            Get attribute from XML file.
``gammalib::xml_check_par``           Check parameter in XML file.
``gammalib::xml_file_expand``         Expand file name in XML file.
``gammalib::xml_file_reduce``         Reduce file name in XML file.
``gammalib::xml_get_name_value_pair`` Get name/value pair from XML node.
``gammalib::recv``                    Receive on socket with timeout.
``gammalib::roi_arclength``           Compute arc length of intersection with Region of Interest.
===================================== ===========
