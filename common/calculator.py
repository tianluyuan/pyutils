import numpy as np
import healpy as hp


def centers(x):
    return (x[:-1]+x[1:])*0.5


def edges(x):
    """returns bin edges with centers that approximately match x. approx
    "inverse" of center(x). Note that it is impossible to ensure that
    centers(edges(x)) == x for all x using the functions defined in
    this module as the returned values of centers(e) are subject to
    constraints that do not necessarily exist for arbitrary center
    points.
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
    phi = np.arctan2(y, x)

    return r, theta, phi


def sphe_to_cart(r, theta, phi):
    x = r*np.sin(theta)*np.cos(phi)    
    y = r*np.sin(theta)*np.sin(phi)
    z  = r*np.cos(theta)

    return x, y, z


def center_angle(theta0, phi0, theta1, phi1):
    """given two directions in spherical coordinates calculate the angle
    between them

    angles should be in radians
    """
    pos1 = sphe_to_cart(1, theta0, phi0)
    pos2 = sphe_to_cart(1, theta1, phi1)
    dotp = np.sum(np.asarray(pos1).T*np.asarray(pos2).T, axis=-1)
    return np.arccos(dotp)

    
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
    return theta, phi, R, kappa, np.median(center_angle(theta, phi, thetas, phis))


def med_ang_res(thetas, phis, nside=64):
    """ Returns mode direction and median deviation from that as the angular res
    """
    if len(thetas) == 0:
        return np.nan, np.nan, np.nan
    pixs = hp.ang2pix(nside, thetas, phis)
    mode_pix = np.bincount(pixs).argmax()
    mode_th, mode_phi = hp.pix2ang(nside, mode_pix)
    return mode_th, mode_phi, np.median(center_angle(mode_th, mode_phi, thetas, phis))
