"""Useful functions for working with icecube software
"""
import os
import re
from collections import namedtuple, defaultdict
from . import standalone
from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, simclasses, dataio
from icecube.hdfwriter import I3HDFWriter, I3SimHDFWriter
from icecube.frame_object_diff import segments
from icecube.BadDomList.BadDomListTraySegment import BadDomList
from icecube.sim_services import ShowerParameters
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
    """
    Builds a map of the relative DOM efficiencies (RDE) based on ice-model resources.

    Parameters
    ----------
    resource : str, optional
        Path to the resource file containing the RDE data.

    Returns
    -------
    dict
        A mapping of OMKeys to their respective relative DOM efficiencies and group values.
    """
    DomEff = namedtuple('DomEff', 'rde grp')
    data = np.loadtxt(resource)
    data_map = {}
    for row in data:
        data_map[icetray.OMKey(int(row[0]), int(row[1]))] = DomEff(row[2], row[3])

    return data_map


def update_dom_bad(frame, baddoms):
    """
    Appends additional bad DOMs to the `BadDomsList` in the frame.

    Parameters
    ----------
    frame : I3Frame
        The frame to update.
    baddoms : list
        List of additional bad DOMs to append.
    """
    for bad in baddoms:
        frame['BadDomsList'].append(bad)


def update_dom_eff(frame, rde_map):
    """
    Updates the relative DOM efficiencies in the DOM calibration of `I3Calibration`.

    Parameters
    ----------
    frame : I3Frame
        The frame containing the calibration data.
    rde_map : dict
        Mapping of OMKeys to their respective relative DOM efficiencies (RDE) values.
    """
    cal = frame['I3Calibration']
    dom_cal = cal.dom_cal
    for dom in rde_map:
        if dom in dom_cal:
            dom_cal[dom].relative_dom_eff = rde_map[dom].rde * (1 + 0.35 * rde_map[dom].grp)


def excluded_doms(frame, exclude_list, keep_partial=True):
    """
    Returns a list of excluded DOMs for the current frame based on the exclude list.

    Returns a list of excluded doms for the current frame based on the exclude_list.
    Partially excluded DOMs are kept by default but if 'keep_partial' is flagged as
    False then those DOMs will be excluded as well.

    Parameters
    ----------
    frame : I3Frame
        The frame to check for excluded DOMs.
    exclude_list : list
        Categories of DOMs to exclude.
    keep_partial : bool, optional
        Whether to keep partially excluded DOMs. Default is True.

    Returns
    -------
    list
        List of excluded OMKeys.
    """
    excluded = []
    for category in exclude_list:
        if frame.Has(category):
            if keep_partial and isinstance(frame[category], dataclasses.I3TimeWindowSeriesMap):
                continue
            for k in frame[category]:
                if isinstance(k, icetray.OMKey):
                    excluded.append(k)

    return excluded


def missing(key):
    """
    Checks if a frame is missing a specific key.

    Parameters
    ----------
    key : str
        Key to check in the frame.

    Returns
    -------
    function
        A lambda function to use as an `If` parameter in IceTray modules.
    """
    return lambda frame: not frame.Has(key)


def print_event(particle, header):
    """
    Prints out conversions of the particle information, including MJD and equatorial coordinates.

    Parameters
    ----------
    particle : I3Particle
        The particle to print information for.
    header : I3EventHeader
        The event header containing the start time.
    """
    from icecube import astro
    mjd = header.start_time.mod_julian_day_double
    eqtr = astro.I3GetEquatorialFromDirection(particle.dir, header.start_time)
    print(header)
    print(particle)
    print('------')
    print(f'mjd = {mjd:.2f}')
    print(f'dec, ra = {eqtr.dec / I3Units.deg:.2f}, {eqtr.ra / I3Units.deg:.2f}')


def stringify(all_pulses, min_q=0):
    """
    Sorts pulses into a dictionary with string numbers as keys and DOMs on that string with their pulses as values.

    Parameters
    ----------
    all_pulses : I3RecoPulseSeriesMap
        The pulse series map to process.
    min_q : float, optional
        Minimum total charge requirement. Default is 0.

    Returns
    -------
    dict
        A dictionary mapping string numbers to lists of `DomInfo` objects.
    """
    i3strings = defaultdict(list)
    for dom in all_pulses.keys():
        pulses = all_pulses[dom]
        times, charges = [], []
        for pulse in pulses:
            times.append(pulse.time)
            charges.append(pulse.charge)
        times = np.array(times)
        charges = np.array(charges)
        if np.sum(charges) > min_q:
            i3strings[dom.string].append(DomInfo(dom, times, charges))

    return i3strings


def qtot(all_pulses):
    """
    Computes the total charge (qtot) by summing over charges in a pulse series.

    Parameters
    ----------
    all_pulses : I3RecoPulseSeriesMap or I3MCPESeriesMap
        The pulse series map to process.

    Returns
    -------
    float
        Total charge in the pulse series.
    """
    if isinstance(all_pulses, simclasses.I3MCPESeriesMap):
        return sum([pulse.npe for pulses in all_pulses.values() for pulse in pulses])
    else:
        return sum([pulse.charge for pulses in all_pulses.values() for pulse in pulses])


def hdfwriter(inp, out, subeventstreams=None, keys=None, types=None, fn=None):
    """
    Writes input data into an HDF5 file.

    Parameters
    ----------
    inp : str or list
        Input file or list of input files.
    out : str
        Output HDF5 file path.
    subeventstreams : list, optional
        Sub-event streams to include in the HDF5 file.
    keys : list, optional
        Specific keys to include in the HDF5 file.
    types : list, optional
        Specific types to include in the HDF5 file.
    fn : callable, optional
        Optional function to apply during data processing.

    Returns
    -------
    None
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
    """
    Writes simulation Q-frames without event headers into an HDF5 file.

    Parameters
    ----------
    inp : str or list
        Input file or list of input files.
    out : str
        Output HDF5 file path.
    runnumber : int, optional
        Run number to include in the output. Default is 0.
    subeventstreams : list, optional
        Sub-event streams to include in the HDF5 file.
    keys : list, optional
        Specific keys to include in the HDF5 file.
    types : list, optional
        Specific types to include in the HDF5 file.

    Returns
    -------
    None
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
    """
    Filters and extracts a specific event from input files and writes it to an output file.

    Parameters
    ----------
    inp : str or list
        Input file or list of input files.
    out : str
        Output file path.
    run : int
        Run ID of the event to filter.
    event : int
        Event ID of the event to filter.
    subevent : int, optional
        Sub-event ID of the event to filter. Default is None.

    Returns
    -------
    None
    """
    def is_event(frame):
        return (frame.Has('I3EventHeader') and
                frame['I3EventHeader'].run_id == run and
                frame['I3EventHeader'].event_id == event and
                (subevent is None or frame['I3EventHeader'].sub_event_id == subevent))

    filter_func(inp, out, is_event)


def filter_func(inp, out, fn=lambda frame: True,
                filter_streams=[icetray.I3Frame.Physics, icetray.I3Frame.DAQ]):
    """
    Filters frames from input files based on a provided filtering function and writes to an output file.

    Parameters
    ----------
    inp : str or list
        Input file or list of input files.
    out : str
        Output file path.
    fn : callable, optional
        Filtering function that determines which frames to include. Default is a function that includes all frames.
    filter_streams : list, optional
        Streams to filter. Default is Physics and DAQ streams.

    Returns
    -------
    None
    """
    tray = I3Tray()
    if isinstance(inp, str):
        inp = [inp]

    tray.Add('I3Reader', Filenamelist=inp)
    tray.Add(fn, streams=filter_streams)
    tray.Add('I3Writer', 'writer', filename=out,
             streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.Physics, icetray.I3Frame.DAQ])
    tray.Execute()


def hese_names(run_id, event_id=None):
    """
    Generates a human-readable name for a high-energy starting event (HESE).

    Parameters
    ----------
    run_id : int or I3EventHeader
        Run ID or I3EventHeader object.
    event_id : int, optional
        Event ID. Only required if `run_id` is not an I3EventHeader object.

    Returns
    -------
    str
        Human-readable HESE name.

    Reference
    ---------
    Based on a file from C. Kopper
    """
    if isinstance(run_id, dataclasses.I3EventHeader):
        event_id = run_id.event_id
        run_id = run_id.run_id
    return standalone.simple_hese_names(run_id, event_id)


def parse_input(inputfile, eventnum=None):
    """
    Parses the run number from an input file name.

    Parameters
    ----------
    inputfile : str
        Path to the input file.
    eventnum : int, optional
        Event number to include in the output. Default is None.

    Returns
    -------
    str
        Parsed run number, optionally appended with the event number.
    """
    filename = os.path.basename(inputfile)
    m = re.match(r'.*([0-9]{6,8}).*', filename)
    runnum = m.group(1).lstrip('0')
    if eventnum is not None and eventnum > 0:
        runnum += f'_{eventnum}'
    return runnum


def ikeep(hcenter, frame, dom):
    """
    Identifies indices of a histogram that fall outside of saturation windows and calibration errata.

    Parameters
    ----------
    hcenter : ndarray
        Histogram center values in time.
    frame : I3Frame
        Current frame containing calibration errata and saturation windows.
    dom : OMKey
        OMKey of the DOM to check.

    Returns
    -------
    ndarray
        Boolean array indicating which indices to keep.
    """
    iskips = [np.array([False] * len(hcenter))]
    if frame.Has('CalibrationErrata') and dom in frame['CalibrationErrata']:
        for cerrata in frame['CalibrationErrata'][dom]:
            iskips.append(
                abs(hcenter - (cerrata.start + cerrata.stop) / (2 * I3Units.microsecond)) <
                (cerrata.stop - cerrata.start) / (2 * I3Units.microsecond))

    if frame.Has('SaturationWindows') and dom in frame['SaturationWindows']:
        for swindow in frame['SaturationWindows'][dom]:
            iskips.append(
                abs(hcenter - (swindow.start + swindow.stop) / (2 * I3Units.microsecond)) <
                (swindow.stop - swindow.start) / (2 * I3Units.microsecond))

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


def isem(particle):
    """
    Checks if a particle is an electromagnetic (EM) particle.

    Parameters
    ----------
    particle : I3Particle
        The particle to check.

    Returns
    -------
    bool
        True if the particle is an EM particle, otherwise False.
    """
    return particle.type in [
        dataclasses.I3Particle.ParticleType.EMinus,
        dataclasses.I3Particle.ParticleType.EPlus,
        dataclasses.I3Particle.ParticleType.Brems,
        dataclasses.I3Particle.ParticleType.DeltaE,
        dataclasses.I3Particle.ParticleType.PairProd,
        dataclasses.I3Particle.ParticleType.Gamma,
        dataclasses.I3Particle.ParticleType.Pi0,  # Pi0 decays to 2 gammas and produce EM showers
    ]


def ishadron(particle):
    """
    Checks if a particle is a hadron.

    Parameters
    ----------
    particle : I3Particle
        The particle to check.

    Returns
    -------
    bool
        True if the particle is a hadron, otherwise False.
    """
    return particle.type in [
        dataclasses.I3Particle.ParticleType.Hadrons,
        dataclasses.I3Particle.ParticleType.Neutron,
        dataclasses.I3Particle.ParticleType.PiPlus,
        dataclasses.I3Particle.ParticleType.PiMinus,
        dataclasses.I3Particle.ParticleType.K0_Long,
        dataclasses.I3Particle.ParticleType.KPlus,
        dataclasses.I3Particle.ParticleType.KMinus,
        dataclasses.I3Particle.ParticleType.PPlus,
        dataclasses.I3Particle.ParticleType.PMinus,
        dataclasses.I3Particle.ParticleType.K0_Short,
        dataclasses.I3Particle.ParticleType.Eta,
        dataclasses.I3Particle.ParticleType.Lambda,
        dataclasses.I3Particle.ParticleType.SigmaPlus,
        dataclasses.I3Particle.ParticleType.Sigma0,
        dataclasses.I3Particle.ParticleType.SigmaMinus,
        dataclasses.I3Particle.ParticleType.Xi0,
        dataclasses.I3Particle.ParticleType.XiMinus,
        dataclasses.I3Particle.ParticleType.OmegaMinus,
        dataclasses.I3Particle.ParticleType.NeutronBar,
        dataclasses.I3Particle.ParticleType.LambdaBar,
        dataclasses.I3Particle.ParticleType.SigmaMinusBar,
        dataclasses.I3Particle.ParticleType.Sigma0Bar,
        dataclasses.I3Particle.ParticleType.SigmaPlusBar,
        dataclasses.I3Particle.ParticleType.Xi0Bar,
        dataclasses.I3Particle.ParticleType.XiPlusBar,
        dataclasses.I3Particle.ParticleType.OmegaPlusBar,
        dataclasses.I3Particle.ParticleType.DPlus,
        dataclasses.I3Particle.ParticleType.DMinus,
        dataclasses.I3Particle.ParticleType.D0,
        dataclasses.I3Particle.ParticleType.D0Bar,
        dataclasses.I3Particle.ParticleType.DsPlus,
        dataclasses.I3Particle.ParticleType.DsMinusBar,
        dataclasses.I3Particle.ParticleType.LambdacPlus,
        dataclasses.I3Particle.ParticleType.WPlus,
        dataclasses.I3Particle.ParticleType.WMinus,
        dataclasses.I3Particle.ParticleType.Z0,
        dataclasses.I3Particle.ParticleType.NuclInt,
        dataclasses.I3Particle.ParticleType.WeakInt,
    ]


def get_deposit_energy(mctree, tstop=np.inf):
    r""" computes the true deposited EM-equivalent energy from the MCTree

    Parameters
    ----------
    mctree : I3MCTree object
        Contains the true particles produced during signal and background
        generation.
    tstop : float, (default np.inf) optional
        Cut off time for inclusion in the calculation. Can be used to compute
        the deposited energy up to e.g. tau decay.

    Returns
    -------
    float
        Sum over all EM-equivalent, InIce, non-Dark energy losses
    """
    losses = 0
    for p in mctree:
        if not p.is_cascade: continue
        if not p.location_type == dataclasses.I3Particle.InIce: continue
        if p.shape == p.Dark: continue
        if p.time > tstop: continue
        losses += ShowerParameters(p.type, p.energy).emScale * p.energy

    return losses


def rebuild_gcd(gcddiff, gcdout='GCD.i3.zst', runid=0, issim=False, writeqp=False):
    """
    Rebuilds GCD (Geometry, Calibration, and DetectorStatus) file from GCDDiff file.

    Parameters
    ----------
    gcddiff : str
        Path to the GCD difference file.
    gcdout : str, optional
        Path to the output GCD file. Default is 'GCD.i3.zst'.
    runid : int, optional
        Run ID for the output file. Default is 0.
    issim : bool, optional
        Whether the file is for simulation data. Default is False.
    writeqp : bool, optional
        Whether to include Q and P frames in the output. Default is False.

    Returns
    -------
    None
    """
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
    """
    Splits Q-frames into individual I3 files.

    Parameters
    ----------
    infile : str
        Path to the input file.
    outdir : str
        Directory to store the output Q-frame files.

    Returns
    -------
    None
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
    """
    Splits P-frames into individual I3 files.

    Parameters
    ----------
    infile : str
        Path to the input file.
    outdir : str
        Directory to store the output P-frame files.
    subeventstreams : list, optional
        Specific sub-event streams to include. Default is None.
    append_fname : bool, optional
        Whether to append the input file name to the output file names. Default is False.

    Returns
    -------
    None
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
    """
    Refines the vertex time based on the closest approach distance.

    Parameters
    ----------
    vertex : ndarray
        The vertex position as a numpy array.
    time : float
        Initial vertex time.
    direction : ndarray
        Direction vector of the event
    """
    min_d = np.inf
    min_t = time
    adj_d = 0
    for om in pulses.keys():
        rvec = omgeo[om].position-vertex
        _l = -rvec*direction
        if _l > 0.:  # only look ahead from vertex
            continue
        _d = np.sqrt(rvec.mag2-_l**2) # closest approach distance
        if _d < min_d: # closest om
            min_d = _d
            min_t = pulses[om][0].time
            adj_d = _l+_d/np.tan(THC)-_d/(np.cos(THC)*np.sin(THC)) # translation distance
    if np.isinf(min_d):
        return time
    return min_t + adj_d/CCC
