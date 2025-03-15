from __future__ import print_function, division

import numpy as np

from jplephem.spk import SPK

from astropy.constants import c
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

KERNEL = SPK.open('de430.bsp')

def bjd_tdb(jd_utc, sky_position):
    """
    Return BJD_TDB from JD_UTC at geocenter for object at desired RA/Dec.
    Parameters
    ----------
    jd_utc : flt or ``astropy.Time`` or ``astropy.Quatity``
        One or more times, as Julian Date in the UTC time scale. UTC is what
        any networked computer's timestamp records.
    sky_position : ``astropy.coordinates.SkyCoord``
        Position on the sky o Fourier transform object of interest. Must be
        a single coordinate.
    Returns
    -------
    bjd_tdb : `astropy.Quantity`
        BJD in the TDB time scale for each of the input JD_UTC.
    """
    # Calculate jd_tdb from jd_utc for each of the input jd_utc
    jd_utc = Time(jd_utc, format='jd', scale='utc')
    jd_tdb = jd_utc.tdb
    # Calculate the position of the sun at each jd_tdb
    sun_earth_displacement = KERNEL[0, 3].compute(jd_tdb.jd) * u.km

    sum_axis = sun_earth_displacement.ndim - 1
    # Calculate the time offset to the barycenter using the plane wave approximation
    delta_t_utc = (sun_earth_displacement.T *
                   sky_position.cartesian.xyz).sum(axis=sum_axis)/c

    # Add that offset to jd_tdb to obtain bjd_tdb
    bjd_tdb = jd_tdb + delta_t_utc

    return bjd_tdb