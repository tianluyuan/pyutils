import numpy as np
from scipy.stats import gamma


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
    n =  (np.max(x) - np.min(x)) / (2 * len(x)**(-1./3) * (np.percentile(x, 75) - np.percentile(x, 25)))
    return 10 if np.isnan(n) else int(n)


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


def sphe_to_kent(theta, phi):
    """ shift x,y,z->x2,x3,x1 which is the coordinate system used in Kent 1982

    return np.array((x1,x2,x3), ...) which is the format read-in by kent_distribution
    """
    xyz = sphe_to_cart(1, theta, phi)
    xs = np.asarray((xyz[2], xyz[0], xyz[1])).T
    return xs


def center_angle(theta0, phi0, theta1, phi1):
    """given two directions in spherical coordinates calculate the angle
    between them

    angles should be in radians
    """
    pos1 = sphe_to_cart(1, theta0, phi0)
    pos2 = sphe_to_cart(1, theta1, phi1)
    dotp = np.sum(np.asarray(pos1).T*np.asarray(pos2).T, axis=-1)
    return np.asfarray(np.arccos(np.clip(dotp,-1,1)))


def vmf_stats(thetas, phis, p=3):
    """Return the average of directional coordinates assuming they are
    drawn from a Von-Mises Fisher distribution (gaussian on a sphere)
    """
    norm, theta, phi = mean_ang(thetas, phis)
    
    R = norm/thetas.size
    # below are approximations
    kappa = R*(p-R**2)/(1-R**2)
    d = 1-np.sum(np.dot(np.array(sphe_to_cart(1, theta, phi)), np.array(sphe_to_cart(1, thetas, phis)))**2)/thetas.size
    sigma = np.sqrt(d/(thetas.size*R**2))
    return theta, phi, R, kappa, sigma


def mean_ang(thetas, phis):
    """ calculates mean angle given thetas and phis
    """
    xsum = np.sum(np.sin(thetas)*np.cos(phis))
    ysum = np.sum(np.sin(thetas)*np.sin(phis))
    zsum = np.sum(np.cos(thetas))

    norm, mean_th, mean_phi = cart_to_sphe(xsum, ysum, zsum)
    return norm, mean_th, mean_phi


def mode_ang(thetas, phis):
    """ Returns mode direction
    """
    import healpy as hp
    if len(thetas) == 0:
        return np.nan, np.nan, np.nan
    nside = 128
    pixs = hp.ang2pix(nside, thetas, phis)
    counts = np.bincount(pixs)
    mode_pixs = np.flatnonzero(counts==counts.max())
    while nside > 16 and len(mode_pixs) > 1:
        nside /= 2
        pixs = hp.ang2pix(nside, thetas, phis)
        counts = np.bincount(pixs)
        mode_pixs = np.flatnonzero(counts==counts.max())

    norm, mode_th, mode_phi = mean_ang(*hp.pix2ang(nside, mode_pixs))
    return norm, mode_th, mode_phi, 


def med_ang_res(theta, phi, thetas, phis):
    """ calculate median angular deviation from theta, phi
    """
    return np.median(center_angle(theta, phi, thetas, phis))


def poisson_llh(hdata, hexp):
    """ returns the poisson llh evaluated from hexp for the hdata

    *hdata* and *hexp* must have the same length
    """
    pllh = np.ma.masked_where(hexp==0,
                              gamma.logpdf(hexp, hdata+1))
    return -np.sum(pllh), pllh.count()


def dchi2(hdata, hexp):
    """ returns the unnormalized chi2 evaluated from hexp for hdata
    """
    return np.sum((hdata-hexp)**2), len(hdata)


def dima_llh(hdata, hexp, sigma=0.1, ns=1, nd=1):
    """ returns dima's llh
    """
    s = hexp*ns
    d = hdata*nd
    ls = np.ma.masked_where(s+d==0, (s+d)/(ns+nd))
    ld = ls
    tol = 1e-8

    if sigma > 0:
        while True:
            f = ns*ls-s-np.log(ld/ls)/sigma**2
            df = 1+(1/(nd*ld)+1/(ns*ls))/sigma**2
            dls = -f/(ns*df)
            ls += dls
            ld = (s+d-ns*ls)/nd
            if np.all(abs(dls) < tol*ls):
                break

    dllh = s*np.log(hexp/ls)+d*np.log(hdata/ld)+np.log(ld/ls)**2/(2*sigma**2)
    return np.sum(dllh), dllh.count()


def most_likely(arr, weights=None):
    """ return the densest region given a 1D array of data
    """
    binning = calc_bins(arr)
    harr = np.histogram(arr, binning, weights=weights)[0] 
    return centers(binning)[np.argmax(harr)]


def interval(arr, percentile=68.3):
    """returns the *percentile* shortest interval around the
    mode
    """
    if len(arr) == 0:
        return np.nan, np.nan, np.nan
    center = most_likely(arr)
    sarr = sorted(arr)
    delta = np.abs(sarr - center)
    curr_low = np.argmin(delta)
    curr_up = curr_low
    npoints = len(sarr)
    while curr_up - curr_low < percentile/100.*npoints:
        if curr_low == 0:
            curr_up += 1
        elif curr_up == npoints-1:
            curr_low -= 1
        elif sarr[curr_up]-sarr[curr_low-1] < sarr[curr_up+1]-sarr[curr_low]:
            curr_low -= 1
        elif sarr[curr_up]-sarr[curr_low-1] > sarr[curr_up+1]-sarr[curr_low]:
            curr_up += 1
        elif (curr_up - curr_low) % 2:
            # they are equal so step half of the time up and down
            curr_low -= 1
        else:
            curr_up += 1

    return sarr[curr_low], center, sarr[curr_up]
