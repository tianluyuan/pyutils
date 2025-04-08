from icecube import icetray, dataclasses, millipede
from icecube import phys_services
import numpy
from snowflake import library
import logging
logger = logging.getLogger(__name__)


class Unfold(icetray.I3Module):
    """ A class originally gifted by K. Jero (from JvS?)
    """

    def __init__(self, ctx):
        super(Unfold, self).__init__(ctx)
        self.AddParameter('Loss_Vector_Name',
                          'Name of the loss vector or I3MCTree',
                          'I3MCTree')
        self.AddParameter('Pulses', 'Name of the Pulses', 'InIcePulses')
        self.AddParameter('FitName',
                          'Name of the fit to use for closest approach distance calculations',
                          'SeedTrack')
        self.AddParameter('CascadePhotonicsService',
                          'Photonics service for cascades',
                          None)
        self.AddParameter('ExcludedDOMs',
                          'DOMs to exclude',
                          ['BadDomsList', 'CalibrationErrata', 'SaturationWindows'])
        self.AddParameter('PhotonsPerBin',
                          'Number of photoelectrons to include in each timeslice',
                          0)
        self.AddParameter('BinSigma',
                          'Bayesian blocking sigma',
                          0)

        self.AddOutBox('OutBox')

    def Configure(self):
        self.input_loss_vect_name = self.GetParameter('Loss_Vector_Name')
        self.pulses = self.GetParameter('Pulses')
        self.fitname = self.GetParameter('FitName')
        self.cscd_service = self.GetParameter('CascadePhotonicsService')
        self.exclude_doms = self.GetParameter('ExcludedDOMs')
        self.ppb = self.GetParameter('PhotonsPerBin')
        self.bs = self.GetParameter('BinSigma')

    def Physics(self, frame):
        if not frame.Has(self.input_loss_vect_name):
            logger.warning(f'no unfolding for {self.input_loss_vect_name}')
            self.PushFrame(frame)
            return True

        self.millipede = millipede.PyPyMillipede(self.context)
        caddict = {}
        ObQtotdict = {}
        ObQdict = {}
        ExQdict = {}
        TBdict = {}

        if self.input_loss_vect_name == 'I3MCTree':
            I3MCTree = frame['I3MCTree']
            sources = [
                loss for loss in I3MCTree.get_daughters(
                    I3MCTree.first_child(
                        I3MCTree[0]))]
        elif isinstance(frame[self.input_loss_vect_name], dataclasses.I3Particle) :
            sources = [frame[self.input_loss_vect_name]]
        else:
            sources = frame[self.input_loss_vect_name]
        for s in sources:
            logger.debug(f'time: {s.time} energy: {s.energy}')

        # This line needs to call get_photonics not the service itself
        self.millipede.SetParameter('CascadePhotonicsService', self.cscd_service)
        self.millipede.SetParameter('ExcludedDOMs', self.exclude_doms)
        self.millipede.SetParameter('Pulses', self.pulses)
        self.millipede.SetParameter('PhotonsPerBin', self.ppb)
        self.millipede.SetParameter('BinSigma', self.bs)
        self.millipede.DatamapFromFrame(frame)
        response = self.millipede.GetResponseMatrix(sources)
        logger.info(f'Fit Statistics For Losses: {self.millipede.FitStatistics(sources, response, params=None)}')
        edeps = [p.energy for p in sources]
        responsemat = response.to_I3Matrix()
        expectations = numpy.inner(responsemat, edeps)
        thisreco = frame[self.fitname]
        I3OMGeo = frame['I3Geometry'].omgeo
        I3EH = frame['I3EventHeader']
        try:
            ps = frame[self.pulses].apply(frame)
        except:
            ps = frame[self.pulses]

        itera = -1
        observedqtmp = []
        expectedqtmp = []
        qtimebinstmp = []
        hitinfotmp = []
        ExcludedDOMlist = library.excluded_doms(frame, self.exclude_doms)
        for k, dc in self.millipede.domCache.items():
            valid = dc.valid
            if k in ExcludedDOMlist:
                for v in valid:
                    if v:
                        itera += 1
                continue
            cad = phys_services.I3Calculator.closest_approach_distance(
                thisreco, I3OMGeo[k].position)

            ObQ = 0
            if k in ps:
                for p in ps[k]:
                    ObQ += p.charge

            caddict[k] = cad
            ObQtotdict[k] = ObQ
            thisobdomq = []
            thisexdomq = []
            edgeiter = 0
            timebins = []
            for v in valid:
                if v:
                    itera += 1
                    thischarge = 0
                    tl = dc.time_bin_edges[edgeiter]
                    th = dc.time_bin_edges[edgeiter + 1]
                    if k in ps:
                        for p in ps[k]:
                            if p.time < th and p.time >= tl:
                                thischarge += p.charge
                    timebins.append(tl)  # [tl,th]
                    try:
                        thisexdomq.append(expectations[itera])
                    except:
                        logger.warning('Problem extracting expected Q')
                        pass
                    thisobdomq.append(thischarge)
                    cad = phys_services.I3Calculator.closest_approach_distance(
                        thisreco, I3OMGeo[k].position)
                    cap = phys_services.I3Calculator.closest_approach_position(
                        thisreco, I3OMGeo[k].position)
                    # print 'FitName=', self.fitname
                    # print k, itera
                    # print 'Time Edges', dc.time_bin_edges[edgeiter], dc.time_bin_edges[edgeiter + 1]
                    # print 'Expected charge', expectations[itera], 'Observed charge', thischarge
                    # print 'Closest Approach=', cad
                    # print 'Closest ApproachPosition=', cap
                    # print
                    edgeiter += 1

            timebins.append(th) # last bin edge
            ObQdict[k] = thisobdomq
            ExQdict[k] = thisexdomq
            TBdict[k] = timebins
            if len(thisexdomq) == 0:
                continue
            hitinfotmp.append([I3OMGeo[k].position.x,
                               I3OMGeo[k].position.y,
                               I3OMGeo[k].position.z,
                               cad, I3EH.run_id,
                               I3EH.event_id, k[0], k[1]])
            expectedqtmp.append(thisexdomq)
            observedqtmp.append(thisobdomq)
            qtimebinstmp.append(timebins)
        frame[self.fitname + '_' + self.input_loss_vect_name +
              '_cads'] = dataclasses.I3MapKeyDouble(caddict)
        frame[self.fitname + '_' + self.input_loss_vect_name +
              '_ObQtots'] = dataclasses.I3MapKeyDouble(ObQtotdict)
        frame[self.fitname + '_' + self.input_loss_vect_name +
              '_ObQ'] = dataclasses.I3MapKeyVectorDouble(ObQdict)
        frame[self.fitname + '_' + self.input_loss_vect_name +
              '_ExQ'] = dataclasses.I3MapKeyVectorDouble(ExQdict)
        frame[self.fitname + '_' + self.input_loss_vect_name +
              '_TB'] = dataclasses.I3MapKeyVectorDouble(TBdict)
        self.PushFrame(frame)
        logger.info(f'Unfold done for run {I3EH.run_id} event {I3EH.event_id}')
