""" Module containing settings for icecube dirs
"""
import os

FIG = os.path.join('/home', os.getenv('USER'), 'public_html', 'share', 'fig')
STORE = os.path.join('/data', 'user', os.getenv('USER'), 'store')
DATA = os.path.join('/home', os.getenv('USER'), 'projects', 'icecube', 'data')
