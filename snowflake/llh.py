import os
from collections import namedtuple
from common.calculator import vmf_stats, mean_ang, med_ang_res
import pandas as pd
import numpy as np


def qtots(dat):
    """ returns dom-by-dom charge on ppc simulated hit files
    """
    if os.path.splitext(dat)[-1] == '.csv':
        dat = pd.read_csv(dat, index_col=0)
        qtots = [len(dat.loc[(dat['str']==s) & (dat['om'] == om)]) for s
                 in dat['str'].unique() for om in
                 dat.loc[dat['str']==s]['om'].unique()]
    else:
        dat = pd.read_csv(dat, sep=' ', header=None, names=('str', 'om', 'time', 'charge'))
        qtots = [dat.loc[(dat['str']==s) & (dat['om'] == om)]['charge'].sum() for s
                 in dat['str'].unique() for om in
                 dat.loc[dat['str']==s]['om'].unique()]

    return qtots
    

def read(llhout, llhcut=np.inf):
    """ Reads output file from llh and parses the useful bits into a pandas dataframe

    *llhout* is the output file from llh 10 or llh 16
    """
    llhdat = pd.read_csv(llhout, delim_whitespace=True, header=None,
                         names='l rlogl x y z zenith azimuth e t a b'.split(),
                         error_bad_lines=False)
    llhdat = llhdat.loc[(llhdat['l'].str.isdigit()) | (llhdat['l'] == '+')].apply(pd.to_numeric, errors='ignore')
    llhsteps = llhdat.loc[(llhdat['rlogl'] < llhcut)]
    return llhsteps


def llh_stats(finput, llhchoice='minlast', llhcut=1):
    """ returns stats based on llh output

    *llhchoice* can be
    'min' returns the absolute best fit... in this case there are no errors
    'cut' uses fits with rlogl<1
    'last' uses last 5% of steps
    'all' uses all steps
    'minlast' uses min rlogl half of last 10% of steps. used in llh.cxx
    """
    centerz = namedtuple('centerz', 'rlogl x y z zenith azimuth e t')
    errorz = namedtuple('errorz', 'dl dx dy dz dr dA de dt N')
    if llhchoice == 'min':
        llhsteps = read(finput, llhcut)
        rlogl, x, y, z, zenith, azimuth, e, t = llhsteps.loc[llhsteps['rlogl'].idxmin()].values[1:-2]
        dl = dx = dy = dz = dr = de = dt = 0
        kappa = np.inf
    else:
        if llhchoice == 'cut':
            llhsteps = read(finput, llhcut)
        elif llhchoice == 'last':
            llhfull = read(finput, np.inf)
            # keep only the last 5%
            keep = int(0.05*len(llhfull.index))
            llhsteps = llhfull.iloc[-keep:]
        elif llhchoice == 'minlast':
            llhfull = read(finput, np.inf)
            keep = int(0.1*len(llhfull.index))
            llhsteps = llhfull.iloc[-keep:].sort_values('rlogl').iloc[:keep/2]
        elif llhchoice == 'all':
            # avg over all steps
            llhsteps = read(finput, np.inf)

        rlogl, x, y, z, t = llhsteps[['rlogl', 'x', 'y', 'z', 't']].mean().values
        e = 10**np.log10(llhsteps['e']).mean()
        dl, dx, dy, dz, dt = llhsteps[['rlogl', 'x', 'y', 'z', 't']].std().values
        de = e*np.log10(llhsteps['e']).std()*np.log(10)
        dr = np.sqrt(llhsteps['x']**2+llhsteps['y']**2+llhsteps['z']**2).std()
        # zenith, azimuth, R, kappa = vmf_stats(np.radians(llhsteps['zenith']), np.radians(llhsteps['azimuth']))
        norm, zenith, azimuth = mean_ang(np.radians(llhsteps['zenith']), np.radians(llhsteps['azimuth']))
        dA = med_ang_res(zenith, azimuth, np.radians(llhsteps['zenith']), np.radians(llhsteps['azimuth']))

    centers = centerz(rlogl, x, y, z, np.degrees(zenith), np.degrees(azimuth), e, t)
    errors = errorz(dl, dx, dy, dz, dr, np.degrees(dA), de, dt, len(llhsteps))
    return centers, errors
