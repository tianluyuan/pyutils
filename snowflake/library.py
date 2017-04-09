"""Useful functions for working with icecube software

N.B. don't delete millipede import!!! It is needed to book the
millipede fit params
"""
import os
import re
from collections import namedtuple, defaultdict
from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses, astro, millipede
from icecube.hdfwriter import I3HDFWriter
import numpy as np

# namedtuple class for storing charges and times for each dom
DomInfo = namedtuple('DomInfo', 'dom times charges')

def get_rde_map(resource=os.path.join(os.path.expandvars('$I3_SRC'),
                                      'ice-models',
                                      'resources', 'models',
                                      'spice-latest-full',
                                      'eff-f2k')):
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
    if isinstance(inp, str):
        inp = [inp]

    tray.Add('I3Reader', Filenamelist=inp)
    tray.Add(I3HDFWriter,
             output=out,
             keys=keys,
             SubEventStreams=subeventstreams,
             types=types)
    tray.Execute()
    tray.Finish()


def hese_names(run_id, event_id=None):
    """ stolen from ckopper's file. this is probably incomplete.

    *run_id* can be either run_id or an I3EventHeader
    *event_id* needs to be specified if header is a run_id
    """
    if isinstance(run_id, dataclasses.I3EventHeader):
        event_id = run_id.event_id
        run_id = run_id.run_id

    names = {
        (118178,66452255):'Camilla the Chicken',	#3
        (118283,9445773):'Rowlf',	#5
        (118381,19162840):'Beaker',	#7
        (118435,58198553):'Mr. Snuffleupagus',	#9
        (118545,63733662):'Bert',	#11
        (118549,11722208):'Sweetums',	#13
        (118602,23096391):'Dr. Bunsen Honeydew',	#15
        (118607,40435683):'Lew Zealand',	#17
        (119196,37713300):'Zoot',	#19
        (119214,8606380):'Link Hogthrob',	#21
        (119316,36556705):'Ernie',	#23
        (119352,56498321):'Guy Smiley',	#25
        (119404,80750561):'Miss Piggy',	#27
        (119470,48424887):'Gonzo the Great',	#29
        (119474,33152537):'Rizzo',	#31
        (119595,30769232):'Count von Count',	#33
        (119674,8449256):'Kermit',	#35
        (119842,82622124):'Sam Eagle',	#37
        (120045,22615214):'Dr. Teeth',	#39
        (115994,2538090):'Fozzy',	#41
        (115994,29874216):'Scooter',	#43
        (116528,52433389):'Animal',	#45
        (116698,10198436):'Swedish Chef',	#47
        (116876,63208734):'Dr. Strangepork',	#49
        (116878,34919596):'Statler and Waldorf',	#51
        (117322,7422546):'Pepe the King Prawn',	#53
        (117371,31623515):'Crazy Harry',	#55
        (117782,49441871):'Floyd',	#57
        (118145,5142726):'Oscar'	#59
    }
    try:
        return names[(run_id, event_id)]
    except KeyError:
        return '{} {}'.format(run_id, event_id)

def parse_input(inputfile, eventnum=None):
    """ return the run number parsed from inputfile
    """
    filename = os.path.basename(inputfile)
    m = re.match(r'.*([0-9]{6,8}).*', filename)
    runnum = m.group(1).lstrip('0')
    if eventnum is not None and eventnum > 0:
        runnum += '_{}'.format(eventnum)
    return runnum
