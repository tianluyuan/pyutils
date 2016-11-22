""" Useful functions for working with steamshovel
"""
from I3Tray import I3Units
from icecube import dataclasses

def cascade_reshape(cascade):
    """ A cascade will by default be rendered as a sphere. This gives it a direction and speed.
    """
    cascade.shape = dataclasses.I3Particle.InfiniteTrack
    cascade.type = dataclasses.I3Particle.unknown
    cascade.speed = 3.*10**8 * I3Units.m/I3Units.s
    cascade.length = 0.1 * I3Units.m
