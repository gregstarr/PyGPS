# PyGPS

This project is merging with a more versatile ionosphere data processing project, found here:
[GeoScience Ionospheric Tool](https://github.com/aldebaran1/gsit)

PyGPS is a python codebase for working with GPS data.
It can be used for:

* reading RINEX files from observation and navigation data
* processing TEC
* finding satellite positions in various coordinate systems
* calculating pierce point locations
* removing satellite and receiver bias.

All this code was written under the larger goal of integrating GPS data processing into [GeoData](https://www.github.com/jswoboda/GeoDataPython) which is John Swoboda’s open source ionosphere data processing project.
The code culminates in the function `GDfromRinex()` which runs through
the entire data processing algorithm from RINEX file to GeoData object for one receiver.
The basic procedure is highlighted below:

    GDfromRinex()

1. Read RINEX obs, nav and satellite bias files
2. calculate TEC for each satellite at each time, remove satellite bias
3. Mark hardware detected loss of lock events and software detected cycle slips
4. Calculate satellite position at each time
5. Produce mapping valued (function of elevation)
6. Put everything in giant Panel4D
7. Calculate receiver bias
8. Transfer data into a form that GeoData understands

Along the way the rinex data, which takes a long time to read in its natural format, can be turned into an HDF5 file which reads in a few seconds.
The RINEX reading was adapted from [PyRinex by Michael Hirsch](https://www.github.com/scienceopen/pyrinex) and most of the rest of the code was adapted from code by Bill Rideout.
