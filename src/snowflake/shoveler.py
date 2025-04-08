""" Useful functions for working with steamshovel
"""
from I3Tray import I3Units
from icecube import dataclasses

def cascade_reshape(cascade):
    """ A cascade will by default be rendered as a sphere. This gives it a direction and speed.
    """
    cascade.shape = dataclasses.I3Particle.InfiniteTrack
    cascade.type = dataclasses.I3Particle.Nu
    cascade.speed = 3.*10**8 * I3Units.m/I3Units.s
    cascade.length = 0.1 * I3Units.m


def add_track(frame, name, params):
    """ Adds a cascade with params to the frame
    """
    particle = dataclasses.I3Particle(dataclasses.I3Particle.InfiniteTrack, dataclasses.I3Particle.MuMinus)
    particle.speed = 0.3
    particle.pos = dataclasses.I3Position(*params[:3])
    particle.dir = dataclasses.I3Direction(params[3]*I3Units.deg, params[4]*I3Units.deg)
    particle.energy = params[5]
    particle.time = params[6]

    frame[name] = particle
