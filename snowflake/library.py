""" Useful functions for working with icecube software
"""
from I3Tray import I3Tray
from icecube import icetray
from icecube.icetray import I3Units
from icecube import astro, millipede
from icecube.hdfwriter import I3HDFWriter
from collections import namedtuple, defaultdict
import numpy as np

# namedtuple class for storing charges and times for each dom
DomInfo = namedtuple('DomInfo', 'dom times charges')

def get_rde_map(resource):
    """ Builds a map of the rde based on ice-models/resources/
    """
    DomEff = namedtuple('DomEff', 'rde grp')
    data = np.loadtxt(resource)
    data_map = {}
    for row in data:
        data_map[icetray.OMKey(int(row[0]), int(row[1]))] = DomEff(row[2], row[3])

    return data_map


def update_dom_cal(frame, rde_map):
    """ Updates the relative dom efficiencies in dom calibration of I3Calibration
    """
    cal = frame['I3Calibration']
    dom_cal = cal.dom_cal
    for dom in rde_map:
        # dom_cal[dom].relative_dom_eff = rde*(1+0.35*grp)
        if dom_cal.has_key(dom):
            dom_cal[dom].relative_dom_eff *= rde_map[dom].rde


def excluded_doms(frame, exclude_list, keep_partial=True):
    """Returns a list of excluded doms for the current frame based on the
    exclude_list. Partially excluded DOMs are kept by default but if
    'keep_partial' is flagged as False then those DOMs will be
    excluded as well.

    return: [OMKey(...), ...]

    """
    excluded = []
    for category in exclude_list:
        if frame.Has(category):
            for k in frame[category]:
                if isinstance(k, icetray.OMKey):
                    excluded.append(k)
                elif not keep_partial and isinstance(k[0], icetray.OMKey):
                    excluded.append(k[0])

    return excluded


def missing(key):
    """function to check if frame has key. Can be passed to IceTray
    modules as If parameter
    """
    return lambda frame: not frame.Has(key)


def print_event(particle, header):
    """ Prints out conversions of the particle info
    """
    mjd = header.start_time.mod_julian_day_double
    eqtr = astro.I3GetEquatorialFromDirection(particle.dir,
                                              header.start_time)
    print header
    print particle
    print '------'
    print 'mjd = {:.2f}'.format(mjd)
    print 'dec, ra = {:.2f}, {:.2f}'.format(eqtr.dec/I3Units.deg,
                                            eqtr.ra/I3Units.deg)


def stringify(all_pulses, min_q=0):
    """sorts pulses into a dictionary with string numbers as keys and the
    doms on that string with their pulses as values

    all_pulses: an I3RecoPulseSeriesMap
    min_q: minimum total charge requirement

    returns: {i3string:[Pulse1,...]}
    """
    i3strings = defaultdict(list)
    # loop over all pulses on doms across all strings
    for dom in all_pulses.iterkeys():
        pulses = all_pulses[dom]
        times, charges = [], []
        for pulse in pulses:
            times.append(pulse.time)
            charges.append(pulse.charge)
        times = np.array(times)
        charges = np.array(charges)
        # append each dom's pulse onto string
        if np.sum(charges) > min_q:
            i3strings[dom.string].append(DomInfo(dom,
                                                 times,
                                                 charges))

    return i3strings


def hdfwriter(inp, out, subeventstreams=None, keys=None, types=None):
    """ Tabulates inp data into an HDF5 file
    """
    tray = I3Tray()
    tray.Add('I3Reader', filename=inp)
    tray.Add(I3HDFWriter,
             output=out,
             keys=keys,
             SubEventStreams=subeventstreams,
             types=types)
    tray.Execute()
    tray.Finish()
