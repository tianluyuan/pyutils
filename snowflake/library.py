"""Useful functions for working with icecube software
"""
import os
import re
import copy
from collections import namedtuple, defaultdict
from . import standalone
from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, astro, simclasses, millipede, dataio
from icecube.hdfwriter import I3HDFWriter, I3SimHDFWriter
from icecube.frame_object_diff import segments
from icecube.BadDomList.BadDomListTraySegment import BadDomList
import numpy as np

# Constants
THC = dataclasses.I3Constants.theta_cherenkov
CCC = dataclasses.I3Constants.c

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


def update_dom_bad(frame, baddoms):
    """ Appends extra bad doms to BadDomsList
    """
    for bad in baddoms:
        frame['BadDomsList'].append(bad)


def update_dom_eff(frame, rde_map):
    """ Updates the relative dom efficiencies in dom calibration of I3Calibration
    """
    cal = frame['I3Calibration']
    dom_cal = cal.dom_cal
    for dom in rde_map:
        # dom_cal[dom].relative_dom_eff = rde*(1+0.35*grp)
        if dom in dom_cal:
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
    print(header)
    print(particle)
    print('------')
    print('mjd = {:.2f}'.format(mjd))
    print('dec, ra = {:.2f}, {:.2f}'.format(eqtr.dec/I3Units.deg,
                                            eqtr.ra/I3Units.deg))


def stringify(all_pulses, min_q=0):
    """sorts pulses into a dictionary with string numbers as keys and the
    doms on that string with their pulses as values

    all_pulses: an I3RecoPulseSeriesMap
    min_q: minimum total charge requirement

    returns: {i3string:[Pulse1,...]}
    """
    i3strings = defaultdict(list)
    # loop over all pulses on doms across all strings
    for dom in all_pulses.keys():
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


def qtot(all_pulses):
    """ returns simple qtot summing over charges in pulse series
    """
    if isinstance(all_pulses, simclasses.I3MCPESeriesMap):
        return sum([pulse.npe for pulses in all_pulses.values() for pulse in pulses])
    else:
        return sum([pulse.charge for pulses in all_pulses.values() for pulse in pulses])


def hdfwriter(inp, out, subeventstreams=None, keys=None, types=None, fn=None):
    """ Tabulates inp data into an HDF5 file
    """
    tray = I3Tray()
    if isinstance(inp, str):
        inp = [inp]

    tray.Add('I3Reader', Filenamelist=inp)
    if fn is not None:
        tray.Add(fn)
    tray.Add(I3HDFWriter,
             output=out,
             keys=keys,
             SubEventStreams=subeventstreams,
             types=types)
    tray.Execute()
    tray.Finish()


def simhdfwriter(inp, out, runnumber=0, subeventstreams=None, keys=None, types=None):
    """ Tabulates Q-frames without event headers into an HDF5 file
    """
    tray = I3Tray()
    if isinstance(inp, str):
        inp = [inp]

    tray.Add('I3Reader', Filenamelist=inp)
    tray.Add(I3SimHDFWriter,
             output=out,
             keys=keys,
             RunNumber=runnumber,
             types=types)
    tray.Execute()
    tray.Finish()


def filter_event(inp, out, run, event, subevent=None):
    def is_event(frame):
        return frame['I3EventHeader'].run_id == run and frame['I3EventHeader'].event_id == event and (subevent is None or frame['I3EventHeader'].sub_event_id == subevent)

    filter_func(inp, out, is_event)


def filter_func(inp, out, fn=lambda frame:True,
                filter_streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ]):
    tray = I3Tray()
    if isinstance(inp, str):
        inp = [inp]
    
    tray.Add('I3Reader', Filenamelist=inp)
    tray.Add(fn,
             streams=filter_streams)
    tray.Add('I3Writer', 'writer', filename=out,
             streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.Physics, icetray.I3Frame.DAQ])
    tray.Execute()

    
def hese_names(run_id, event_id=None):
    """ stolen from ckopper's file. this is probably incomplete.

    *run_id* can be either run_id or an I3EventHeader
    *event_id* needs to be specified if header is a run_id
    """
    if isinstance(run_id, dataclasses.I3EventHeader):
        event_id = run_id.event_id
        run_id = run_id.run_id
    return standalone.simple_hese_names(run_id, event_id)


def parse_input(inputfile, eventnum=None):
    """ return the run number parsed from inputfile
    """
    filename = os.path.basename(inputfile)
    m = re.match(r'.*([0-9]{6,8}).*', filename)
    runnum = m.group(1).lstrip('0')
    if eventnum is not None and eventnum > 0:
        runnum += '_{}'.format(eventnum)
    return runnum


def ikeep(hcenter, frame, dom):
    """ hcenter is a histogram in time

    return indices of hcenter that fall outside of the saturation window and calibration errata
    """
    iskips = [np.array([False]*len(hcenter))]
    if frame.Has('CalibrationErrata') and dom in frame['CalibrationErrata']:
        for cerrata in frame['CalibrationErrata'][dom]:
            iskips.append(
                abs(hcenter-(cerrata.start+cerrata.stop)/(2*I3Units.microsecond)) < (cerrata.stop-cerrata.start)/(2*I3Units.microsecond))

    if frame.Has('SaturationWindows') and dom in frame['SaturationWindows']:
        for swindow in frame['SaturationWindows'][dom]:
            iskips.append(
                abs(hcenter-(swindow.start+swindow.stop)/(2*I3Units.microsecond)) < (swindow.stop-swindow.start)/(2*I3Units.microsecond))

    ikeep = ~np.any(iskips, axis=0)
    return ikeep


# #############################
# Here broken events are removed
# #############################
class RemoveBrokenEvents(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox("OutBox")
        self.AddParameter("Mask", "Name of the bitmask you want to test")

    def Configure(self):
        self.mask = self.GetParameter("Mask")

    def Physics(self, fr):
        try:
            fr[self.mask].apply(fr)
        except Exception:
            pass
        else:
            self.PushFrame(fr)

    def DAQ(self, fr):
        self.PushFrame(fr)


def get_deposit_energy(mctree):
    losses = 0
    for p in mctree:
        if not p.is_cascade: continue
        if not p.location_type == dataclasses.I3Particle.InIce: continue
        if p.shape == p.Dark: continue
        if p.type in [p.Hadrons, p.PiPlus, p.PiMinus, p.NuclInt]:
            #hadlosses += p.energy
            if p.energy < 1*I3Units.GeV:
                losses += 0.8*p.energy
            else:
                energyScalingFactor = 1.0 + ((p.energy/I3Units.GeV/0.399)**-0.130)*(0.467 - 1)
                losses += energyScalingFactor*p.energy
        else:
            #emlosses += p.energy
            losses += p.energy

    return losses


def rebuild_gcd(gcddiff, gcdout='GCD.i3.zst', runid=0, issim=False, writeqp=False):
    tray = I3Tray()
    tray.Add('I3Reader', Filename=gcddiff)

    tray.Add(segments.uncompress,
             base_path='/data/exp/IceCube/2016/internal-system/PoleBaseGCDs/')
    tray.Add(BadDomList,
             Simulation=issim,
             RunId=runid)
    streams = [icetray.I3Frame.Geometry,
               icetray.I3Frame.Calibration,
               icetray.I3Frame.DetectorStatus]
    if writeqp:
        streams.extend([icetray.I3Frame.DAQ, icetray.I3Frame.Physics])
    tray.Add('I3Writer', Filename=gcdout,
             Streams=streams)

    tray.Execute()
    tray.Finish()


def splitQ(infile, outdir):
    """ split Qframes into individual i3 files
    """
    with dataio.I3File(infile) as f:
        while f.more():
            qfr = f.pop_daq()
            eve = qfr['I3EventHeader']
            outfile = dataio.I3File(
                os.path.join(outdir, f'{eve.run_id}.{eve.event_id}.i3.zst'), 'w')
            outfile.push(qfr)
            outfile.close()


def splitP(infile, outdir, subeventstreams=None, append_fname=False):
    """ split Pframes into individual i3 files
    """
    with dataio.I3File(infile) as f:
        while f.more():
            pfr = f.pop_physics()
            eve = pfr['I3EventHeader']
            if (subeventstreams is not None and
                (eve.sub_event_stream != subeventstreams or
                 eve.sub_event_stream not in subeventstreams)):
                continue
            parents = f.get_mixed_frames()
            outfile = dataio.I3File(
                os.path.join(outdir, f'{eve.run_id}.{eve.event_id}.{os.path.basename(infile) if append_fname else "i3.zst"}'), 'w')
            for parent_frame in parents:
                outfile.push(parent_frame)
            outfile.push(pfr)
            outfile.close()


def refine_vertex_time(vertex, time, direction, pulses, omgeo):
    min_d = np.inf
    min_t = time
    adj_d = 0
    for om in pulses.keys():
        rvec = omgeo[om].position-vertex
        _l = -rvec*direction
        _d = np.sqrt(rvec.mag2-_l**2) # closest approach distance
        if _d < min_d: # closest om
            min_d = _d
            min_t = pulses[om][0].time
            adj_d = _l+_d/np.tan(THC)-_d/(np.cos(THC)*np.sin(THC)) # translation distance
    if np.isinf(min_d):
        return time
    return min_t + adj_d/CCC


################## pulse cleaning
def _weighted_quantile_arg(values, weights, q=0.5):
    indices = np.argsort(values)
    sorted_indices = np.arange(len(values))[indices]
    medianidx = (weights[indices].cumsum()/weights[indices].sum()).searchsorted(q)
    if (0 <= medianidx) and (medianidx < len(values)):
        return sorted_indices[medianidx]
    else:
        return np.nan

def weighted_quantile(values, weights, q=0.5):
    if len(values) != len(weights):
        raise ValueError("shape of `values` and `weights` don't match!")
    index = _weighted_quantile_arg(values, weights, q=q)
    if not np.isnan(index):
        return values[index]
    else:
        return np.nan

def weighted_median(values, weights):
    return weighted_quantile(values, weights, q=0.5)

def late_pulse_cleaning(frame, Pulses, Residual=1.5e3*I3Units.ns):
    pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, Pulses)
    mask = dataclasses.I3RecoPulseSeriesMapMask(frame, Pulses)
    counter, charge = 0, 0
    qtot = 0
    times = dataclasses.I3TimeWindowSeriesMap()
    for omkey, ps in pulses.items():
        if len(ps) < 2:
            if len(ps) == 1:
                qtot += ps[0].charge
            continue
        ts = np.asarray([p.time for p in ps])
        cs = np.asarray([p.charge for p in ps])
        median = weighted_median(ts, cs)
        qtot += cs.sum()
        ### DEBUG
        # if cs.sum()>200:
        #     from matplotlib import pyplot as plt
        #     plt.figure()
        #     plt.hist(ts, bins=np.arange(median-0.5*Residual, median+3*Residual, 50), weights=cs, histtype='step')
        #     [plt.vlines(_, 0, 10) for _ in [median-Residual, median, median+Residual]]
        #     plt.title(omkey)
        #     plt.yscale('log')
        #     plt.savefig(f'out/misc/pulses/{omkey.string}_{omkey.om}.png')
        for p in ps:
            if p.time >= (median+Residual):
                if not times.has_key(omkey):
                    ts = dataclasses.I3TimeWindowSeries()
                    ts.append(dataclasses.I3TimeWindow(median+Residual, np.inf)) # this defines the **excluded** time window
                    times[omkey] = ts
                mask.set(omkey, p, False)
                counter += 1
                charge += p.charge
    frame[Pulses+"LatePulseCleaned"] = mask
    frame[Pulses+"LatePulseCleanedTimeWindows"] = times
    try:
        frame[Pulses+"LatePulseCleanedTimeRange"] = copy.deepcopy(frame[Pulses+"TimeRange"])
    except KeyError:
        frame[Pulses+"LatePulseCleanedTimeRange"] = copy.deepcopy(frame["CalibratedWaveformRange"])
