import numpy as np


def centers(x):
    return (x[:-1]+x[1:])*0.5


def edges(x):
    """ returns bin edges with centers that match x. approx "inverse" of center(x).
    """
    c = centers(x)
    return np.concatenate(([2*x[0]-c[0]], c, [2*x[-1]-c[-1]]))


def calc_nbins(x):
    n =  (np.max(x) - np.min(x)) / (2 * len(x)**(-1/3) * (np.percentile(x, 75) - np.percentile(x, 25)))
    return np.floor(n)


def calc_bins(x):
    nbins = calc_nbins(x)
    return np.linspace(np.min(x), np.max(x)+2, num=nbins+1)


def cart_to_sphe(x, y, z):
    r = np.sqrt(x*x+y*y+z*z)
    theta = np.arccos(z/r)
    phi = np.arctan(y/x)

    return r, theta, phi


def vmf_stats(thetas, phis, p=3):
    """Return the average of directional coordinates assuming they are
    drawn from a Von-Mises Fisher distribution (gaussian on a sphere)
    """
    xsum = np.sum(np.sin(thetas)*np.cos(phis))
    ysum = np.sum(np.sin(thetas)*np.sin(phis))
    zsum = np.sum(np.cos(thetas))

    norm = np.sqrt(xsum**2+ysum**2+zsum**2)
    r, theta, phi = cart_to_sphe(xsum/norm, ysum/norm, zsum/norm)
    R = norm/thetas.size
    # below are approximations
    kappa = R*(p-R**2)/(1-R**2)
    return theta, phi, R, kappa
