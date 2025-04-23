import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

#ABle to read the observations
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

#Convert the observation into a format which enables the rebinning method
from pysynphot import observation
from pysynphot import spectrum

#Fitting packages for the fit to normalise in the end
from scipy import optimize
from scipy.optimize import curve_fit

#Packages to convert time to BJD
from astropy import time, coordinates as coord, units as u
from astropy.coordinates import FK5
# from astropy.constants import c

#package to extract names from a directory
import glob
import pathlib
import os

pd.set_option('float_format', '{:f}'.format)



# hier nog een optie inzetten om apo te gebruiken
def BJD_calculator(FITStime, source_coordinates):
    "Function which calculates the BJD based on a timestamp and coordinates"
    "Only works for la palma observations"

    # In: '2018-10-13T06:16:52.584086' , out: 2458404.764557102

    # Import of the time, makes a Time format.
    time_variable = Time(FITStime, format='fits')

    # Creating a skycoordinate,
    sc = coord.SkyCoord(source_coordinates, frame=FK5)

    # Importing the location of lapalma observatorium.
    palma = coord.EarthLocation.of_site('lapalma')

    # a time and location, scale=utc: universal time code
    times = time.Time(time_variable, format='jd', location=palma, scale='utc')

    # Calculating the correction for the time in jd, output= TimeDelta_object
    ltt_bary = times.light_travel_time(sc)

    # barycorr = sc.radial_velocity_correction(obstime=Time('2016-6-4'), location=keck)
    barycorr = sc.radial_velocity_correction(kind='barycentric', obstime=times)
    barycorr_kms = barycorr.to(u.km / u.s)

    # Adding both times for barycentre time
    barycentric_time = times + ltt_bary.tcb

    return (round(barycentric_time.value, 6), round(barycorr_kms.value, 6))