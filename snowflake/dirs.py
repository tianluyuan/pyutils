""" Module containing settings for icecube dirs
"""
import os

FIG = os.path.join(os.getenv('HOME'), 'public_html', 'share', 'fig')
DATA = os.path.join(os.getenv('HOME'), 'projects', 'icecube', 'data')
CONDOR_JOB = os.path.join('/scratch', os.getenv('USER'), 'condor')
CONDOR_OUT = os.path.join('/data', 'user', os.getenv('USER'), 'condor')
